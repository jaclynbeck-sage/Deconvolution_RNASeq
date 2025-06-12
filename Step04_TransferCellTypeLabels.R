# This script transfers cell type labels from the Seattle Reference Atlas to all
# single cell data sets, in order to get consistent broad and fine cell-type
# assignments across all the data sets. Some cell type labels are renamed and/or
# re-grouped based on expression similarity to reduce the number of sub-classes,
# excluded genes (mitochondrial and non-coding) are removed, and data sets are
# subset to only cells above a certain level of annotation confidence as output
# by MapMyCells.
#
# Labels are re-assigned as follows:
#   Broad class:
#     "Neuronal: GABAergic" => "Inhibitory"
#     "Neuronal: Glutamatergic" => "Excitatory"
#     "Non-neuronal and Non-neural" => the cell's sub class
#       Exception: sub classes "Endothelial" and "VLMC" => "Vascular"
#
#   Sub class:
#     "Endothelial", "VLMC" => "Vascular"
#     "Lamp5", "Lamp5 Lhx6" => "Lamp5 / Lamp5 Lhx6"
#     "L2/3 IT", "L6 IT" => "L2/3/6 IT"
#     "L4 IT", "L5 IT" => "L4/5 IT"
#     "L5 ET" => removed
#     "Pax6", "Sncg" => "Pax6 / Sncg"
#     "Sst", "Sst Chodl" => "Sst / Sst Chodl"
#     All others: as-is
#
# Group merges were determined by examining UMAPs and dendrograms of all data
# sets after label transfer. Notes:
#   * Although Endothelial and VLMC cells form semi-distinct clusters in all
#     data sets, the cell populations were so small that it made sense to merge
#     them.
#   * Lamp5 and Lamp5 Lhx6 are distinct populations in cain and seaRef but on
#     the dendrograms and the other datasets' UMAPS they intermix and are
#     smaller populations, so they were merged.
#   * Although L2/3 IT and L6 IT are distinct populations on the UMAPs of most
#     data sets, on the dendrograms L6 IT nests inside L2/3 IT supertypes.
#     Initial marker finding with L6 IT as a separate cluster resulted in
#     too few markers being found for L6 IT in many cases, due to the high level
#     of similarity between L6 IT and L2/3 IT. Therefore these two clusters were
#     merged.
#   * L4 IT and L5 IT do not separate by sub-class on UMAPs or dendrograms.
#     Rather, they form two distinct but closely-related clusters that contain a
#     mix of certain L4 IT and L5 IT supertypes together. The two subclasses
#     were merged rather than splitting by supertype in order to avoid having
#     too many highly-similar IT clusters from L2, 3, 4, 5, and 6.
#   * In the case of Pax6, Sncg, and Vip, these are somewhat distinct
#     populations that form a gradient on the UMAP with no clear borders
#     between the cell types. On the dendrograms of most data sets, Pax6 and
#     some Sncg supertypes intermix, and some Sngc and some Vip supertypes
#     intermix, making it difficult to separate them into clear clusters. As
#     both Pax6 and Sncg are small populations, they were merged together and
#     Vip was left separated to avoid Vip cells dominating the expression
#     profile of a combined cluster.
#   * Sst Chodl is an extremely small population in all data sets and
#     consistently clusters with Sst cells, so these were merged.
#   * L5 ET cells are a distinct population in all data sets but make up a
#     fraction of a percent of each sample on average, which makes it difficult
#     to confidently find markers for these cells or accurately predict in
#     deconvolution. The population does not cluster with any other populations,
#     so there was no clear subclass to combine it with and the L5 ET subclass
#     as a whole was removed.
#
# NOTE: Annotation files were obtained by uploading all pre-processed h5ad files
# from Step02 to MapMyCells (https://knowledge.brain-map.org/mapmycells/process/),
# aligning to "10x Human MTG SEA-AD (CCN20230505)" with algorithm "Deep
# Generative Mapping", and manually putting the resulting annotation files on
# Synapse.

library(synapser)
library(MatrixGenerics)
library(dplyr)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- all_singlecell_datasets()

synLogin()

# Annotation files from MapMyCells
anno_ids <- list(cain = "syn68239068.1",
                 lau = "syn68239074.1",
                 lengEC = "syn68252071.1",
                 lengSFG = "syn68252122.1",
                 mathys = "syn68239087.1")

# For threaded processing of seaRef
n_cores <- max(parallel::detectCores() / 2, 1)
DelayedArray::setAutoBPPARAM(BPPARAM = BiocParallel::MulticoreParam(n_cores))

# Process each data set
for (dataset in datasets) {
  sce <- Load_PreprocessedData(dataset, remove_excluded = TRUE)

  if (dataset == "seaRef") {
    # seaRef didn't need re-mapping so we don't have an annotation file for it.
    # We rename a few columns to match the annotation files for the other data
    # sets.
    new_metadata <- colData(sce) |>
      as.data.frame() |>
      dplyr::rename(supertype_name = supertype)
    colnames(new_metadata) <- str_replace(colnames(new_metadata),
                                          "confidence",
                                          "softmax_probability")
    new_metadata$class_name <- new_metadata$broad_class
    new_metadata$subclass_name <- new_metadata$sub_class
  } else {
    # All other data sets have an annotation file from MapMyCells
    anno_f <- synGet(anno_ids[[dataset]],
                     downloadLocation = dir_preprocessed,
                     ifcollision = "overwrite.local")

    anno <- read.csv(anno_f$path, comment.char = "#") |>
      select(-ends_with("label"))

    new_metadata <- merge(colData(sce), anno)
  }

  # All datasets -- adjust broad and sub-classes and rename columns
  new_metadata <- new_metadata |>
    as.data.frame() |>
    mutate(
      # Preserve any previous annotations for debugging
      original_broad_class = broad_class,
      original_sub_class = sub_class,

      # Shorten broad class names, assign actual glial cell names to non-neuronal cells
      broad_class = case_match(
        class_name,
        "Neuronal: GABAergic" ~ "Inhibitory",
        "Neuronal: Glutamatergic" ~ "Excitatory",
        .default = subclass_name
      ),

      # Reclassify Endothelial and VLMC as "Vascular"
      broad_class = case_match(
        broad_class,
        c("Endothelial", "VLMC") ~ "Vascular",
        .default = broad_class),

      # Glia are not split into supertypes, some neuronal subtypes are merged
      # together
      sub_class = case_match(
        subclass_name,
        c("Endothelial", "VLMC") ~ "Vascular",
        c("L2/3 IT", "L6 IT") ~ "L2/3/6 IT",
        c("L4 IT", "L5 IT") ~ "L4/5 IT",
        c("Lamp5", "Lamp5 Lhx6") ~ "Lamp5 / Lamp5 Lhx6",
        c("Pax6", "Sncg") ~ "Pax6 / Sncg",
        c("Sst", "Sst Chodl") ~ "Sst / Sst Chodl",
        .default = subclass_name),

      # Factor
      broad_class = factor(broad_class),
      sub_class = factor(sub_class)
    ) |>

    # Column renames
    dplyr::rename(
      broad_class_softmax_probability = class_softmax_probability,
      sub_class_softmax_probability = subclass_softmax_probability,
      supertype = supertype_name
    ) |>
    select(-class_name, -subclass_name) |>
    as("DataFrame")

  # Add new metadata back to the sce object
  rownames(new_metadata) <- new_metadata$cell_id
  colData(sce) <- new_metadata[colnames(sce), ]

  # Remove L5 ET cells
  sce <- sce[, sce$sub_class != "L5 ET"]
  sce$sub_class <- factor(sce$sub_class)

  # Only keep cells where the probability for the annotation is >= 0.95 for
  # broad and sub class, and > 0.5 for supertype. The min broad class confidence
  # for seaRef is 0.988, and 99.9% of seaRef cells have a sub class confidence
  # >= 0.95, so these cutoffs seem reasonable.
  sce$pass_QC <- sce$broad_class_softmax_probability >= 0.95 &
    sce$sub_class_softmax_probability >= 0.95

  sce <- sce[, sce$pass_QC]

  sce$tmm_factors <- edgeR::calcNormFactors(counts(sce),
                                            lib.size = sce$lib_size,
                                            method = "TMMwsp")

  Save_SingleCell(dataset, sce)

  # For determining cell type merges, create a dendrogram
  counts_log <- sageRNAUtils::simple_log2norm(counts(sce), pseudocount = 1)

  ct_mean <- sapply(sort(unique(sce$supertype)), function(ct) {
    counts_sub <- counts_log[, sce$supertype == ct]

    if (is.null(dim(counts_sub))) {
      counts_sub
    } else {
      MatrixGenerics::rowMeans(counts_sub)
    }
  })

  colnames(ct_mean) <- sort(unique(sce$supertype))

  # Remove genes where no cell types have > 1 log2-CPM expression
  ct_mean <- ct_mean[rowMaxs(ct_mean) >= 1, ]

  # Subset to most variable genes after low expression removal
  var_genes <- Seurat::FindVariableFeatures(ct_mean, nfeatures = 4000) |>
    subset(variable == TRUE)

  ct_dist <- as.dist(1 - cor(ct_mean[rownames(var_genes), ]))
  ct_clust <- hclust(ct_dist, method="average")
  plot(ct_clust, cex = 0.5)
  saveRDS(ct_clust, file.path(dir_tmp, str_glue("{dataset}_supertype_dend.rds")))
}
