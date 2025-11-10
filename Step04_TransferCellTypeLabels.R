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
#     "L4 IT_1" supertype, all "L5 IT" => "L5 IT"
#     "L4 IT" except "L4 IT_1" => "L4 IT"
#     "Pax6", "Sncg" => "Pax6 / Sncg"
#     "Sst", "Sst Chodl" => "Sst / Sst Chodl"
#     All others: as-is
#
# Group merges were determined by examining UMAPs, clusters, and dendrograms of
# all data sets after label transfer. Notes:
#   * The Allen Institute whole human brain nomenclature groups cells roughly as
#     follows, based on mapping of the cain dataset to both nomenclatures:
#       (Chandelier, Lamp5-Lhx6) => LAMP5-LHX6 and Chandelier
#       (L2/3 IT, most L4 IT, some L5 IT) => Upper-layer intratelencephalic
#       (Some L4 IT, most L5 IT, L6 IT, L6 IT Car3) => Deep-layer intratelencephalic
#       L5/6 NP => Deep-layer near-projecting
#       (L6 CT, L6b) => Deep-layer corticothalamic and 6b
#       (Lamp5, Pax6, Sncg, Vip) => CGE interneuron
#       (Pvalb, Sst, Sst Chodl) => MGE interneuron
#   * Although Endothelial and VLMC cells form semi-distinct clusters in all
#     data sets, the cell populations were so small and closely related that it
#     made sense to merge them.
#   * L2/3 IT and L6 IT form distinct populations on the UMAPs of most data
#     sets, but on the dendrograms, L6 IT clusters with L2/3 IT supertypes as
#     a slightly separate mini-cluster. Based on the fact that the whole brain
#     nomenclature separates upper-layer and deep-layer, L6 was left separated
#     but could merge with 2/3 if analysis shows issues with marker finding for
#     this population.
#   * L4 IT and L5 IT form two distinct but closely-related clusters that
#     contain a mix of certain L4 IT and L5 IT supertypes together. L4 IT types
#     2, 3, and 4 cluster together, while L5 IT cells cluster together with L4
#     IT_1. Two L5 IT supertypes form a semi-distinct cluster from the rest of
#     L5/L4_1. To avoid splitting groups too finely but still keep cells with
#     similar expression together, L4 IT 2, 3, and 4 were assigned as group "L4
#     IT" and L4 IT_1 plus all L5 IT were assigned to group "L5 IT".
#   * Lamp5 and Lamp5 Lhx6 appear as two mostly-distinct but connected clusters
#     in all UMAPs, but in dendrograms Lhx6_1 clusters inside of the Lamp5. As
#     These are closely-related and smaller populations,two subclasses were
#     merged.
#   * In the case of Pax6, Sncg, and Vip, these are somewhat distinct
#     populations that form a gradient on the UMAP with no clear borders between
#     the cell types. On the dendrograms of most data sets, Pax6 and 2-3 Sncg
#     supertypes intermix, and most Sngc and some Vip supertypes intermix, while
#     most Vip make their own cluster. It is difficult to separate them into
#     clear clusters. One option based on all the dendrograms is Pax6 by itself
#     and Sncg + Vip together, but both Pax6 and Sncg are small populations and
#     cluster nearly on top of each other in all UMAPs, so Pax6 and Sncg were
#     merged together into one group and Vip was left alone.
#   * Sst Chodl is an extremely small population in all data sets and
#     consistently clusters with Sst cells, so these were merged.
#   * L5 ET cells are a distinct population in all data sets but make up a
#     fraction of a percent of each sample on average, which makes it difficult
#     to confidently find markers for these cells or accurately predict in
#     deconvolution. The population does not cluster with any other populations,
#     so there was no clear subclass to combine it with. Therefore the L5 ET
#     subclass as a whole was removed.
#
# NOTE: Annotation files were obtained by uploading all pre-processed h5ad files
# from Step02 to MapMyCells (https://knowledge.brain-map.org/mapmycells/process/),
# aligning to "10x Human MTG SEA-AD (CCN20230505)" with algorithm "Deep
# Generative Mapping", and manually putting the resulting annotation files on
# Synapse.

library(synapser)
library(MatrixGenerics)
library(dplyr)
library(Seurat)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- all_singlecell_datasets()

synLogin()

# Annotation files from MapMyCells
anno_ids <- config::get("step04_label_transfer")

# Process each data set
for (dataset in datasets) {
  sce <- Load_PreprocessedData(dataset, remove_excluded = TRUE)

  sink(file.path(dir_tmp, str_glue("{dataset}_label_transfer.log")))

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
    anno_f <- synGet(anno_ids[[paste0(dataset, "_annotations")]],
                     downloadLocation = dir_tmp,
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
      sub_class_unmerged = subclass_name,

      sub_class = case_match(
        subclass_name,
        c("Endothelial", "VLMC") ~ "Vascular",
        "L4 IT" ~ ifelse(supertype_name == "L4 IT_1", "L5 IT", "L4 IT"),
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
  # broad and sub class. The min broad class confidence for seaRef is 0.988, and
  # 99.9% of seaRef cells have a sub class confidence >= 0.95, so these cutoffs
  # seem reasonable.
  sce$pass_QC <- sce$broad_class_softmax_probability >= 0.95 &
    sce$sub_class_softmax_probability >= 0.95

  cat("Pass QC:", "\n")
  print(table(sce$sub_class, sce$pass_QC))

  sce <- sce[, sce$pass_QC]

  # Remove "singleton" cells that cluster with a different subclass
  so <- CreateSeuratObject(counts(sce), meta.data = as.data.frame(colData(sce)))
  so <- so |>
    NormalizeData() |>
    FindVariableFeatures(nfeatures = 4000) |>
    ScaleData() |>
    RunPCA() |>
    FindNeighbors(dims = 1:20) |>
    FindClusters(resolution = 2)

  so$seurat_clusters <- paste0("C", so$seurat_clusters)

  so <- DietSeurat(so, layers = "data",
                   dimreducs = "pca",
                   features = VariableFeatures(so))

  min_cells <- 5 # At least 5 cells from a cell type must be in a cluster
  min_pct <- 1.0 # At least 1% of a cell type's total population must be in a cluster

  # Percentage of each cell type's total population that appears in each cluster
  clusters <- table(so$seurat_clusters, so$sub_class)
  clusters <- sweep(clusters, 2, colSums(clusters), "/") * 100

  cat("\nCluster counts:\n")
  print(table(so$seurat_clusters, so$sub_class))
  cat("\nCluster percents:\n")
  print(round(clusters, digits = 1))

  removals <- so@meta.data |>
    group_by(seurat_clusters) |>
    mutate(
      remove = as.logical(clusters[cur_group()$seurat_clusters, sub_class] < min_pct),
      remove = remove | table(sub_class)[sub_class] < min_cells
    ) |>
    subset(remove == TRUE)

  cat(str_glue("Removing {nrow(removals)} mis-clustered cells."), "\n")

  so$remove <- so$cell_id %in% removals$cell_id

  cat("\nCell removals:\n")
  print(table(so$sub_class, so$remove))
  sink()

  saveRDS(so, file.path(dir_tmp, str_glue("{dataset}_seurat.rds")))

  remove_ids <- so$cell_id[so$remove]
  sce <- sce[, !(colnames(sce) %in% remove_ids)]

  sce$tmm_factors <- edgeR::calcNormFactors(counts(sce),
                                            lib.size = sce$lib_size,
                                            method = "TMMwsp")

  Save_SingleCell(dataset, sce)

  rm(so)
  gc()

  # For determining cell type merges, create a dendrogram
  aggr_sce <- scuttle::aggregateAcrossCells(
    sce,
    ids = paste(sce$sample, "::", sce$supertype),
    statistics = c("sum")
  ) |>
    scuttle::logNormCounts(transform = "log", assay.type = "counts")

  aggr_sce <- scuttle::aggregateAcrossCells(
    aggr_sce,
    ids = aggr_sce$supertype,
    statistics = c("mean"),
    use.assay.type = "logcounts"
  )

  aggr <- assay(aggr_sce, "logcounts")

  pdf(file.path(dir_figures, str_glue("{dataset}_supertype_dend.pdf")),
      width = 16, height = 8)

  par(cex = 1) # Can't set cex in plot() directly for dendrograms

  exc_sub_class <- unique(sce$sub_class_unmerged[sce$broad_class == "Excitatory"])
  colors_exc <- viridis::turbo(n = length(exc_sub_class))
  names(colors_exc) <- as.character(exc_sub_class)

  inh_sub_class <- unique(sce$sub_class_unmerged[sce$broad_class == "Inhibitory"])
  colors_inh <- viridis::turbo(n = length(inh_sub_class))
  names(colors_inh) <- as.character(inh_sub_class)

  colors_all <- viridis::turbo(n = length(unique(sce$sub_class_unmerged)))
  names(colors_all) <- as.character(unique(sce$sub_class_unmerged))

  supertype_map <- as.character(aggr_sce$sub_class_unmerged)
  names(supertype_map) <- as.character(aggr_sce$supertype)

  color_leaf <- function(node, colors) {
    if (is.leaf(node)) {
      a <- attributes(node)
      labCol <- colors[supertype_map[a$label]] |> as.character()
      attr(node, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    node
  }

  exc_supertypes <- aggr_sce$supertype[aggr_sce$broad_class == "Excitatory"]
  inh_supertypes <- aggr_sce$supertype[aggr_sce$broad_class == "Inhibitory"]

  var_genes_exc <- rowVars(aggr[, colnames(aggr) %in% exc_supertypes]) |>
    sort(decreasing = TRUE)
  var_genes_exc <- names(var_genes_exc)[1:4000]

  ct_dist <- as.dist(1 - cor(aggr[var_genes_exc, colnames(aggr) %in% exc_supertypes]))
  ct_clust <- hclust(ct_dist, method = "average") |>
    as.dendrogram() |>
    dendrapply(color_leaf, colors = colors_exc)
  plot(ct_clust)

  saveRDS(ct_clust, file.path(dir_tmp, str_glue("{dataset}_supertype_dend_exc.rds")))

  var_genes_inh <- rowVars(aggr[, colnames(aggr) %in% inh_supertypes]) |>
    sort(decreasing = TRUE)
  var_genes_inh <- names(var_genes_inh)[1:4000]

  ct_dist <- as.dist(1 - cor(aggr[var_genes_inh, colnames(aggr) %in% inh_supertypes]))
  ct_clust <- hclust(ct_dist, method = "average") |>
    as.dendrogram() |>
    dendrapply(color_leaf, colors = colors_inh)
  plot(ct_clust)

  saveRDS(ct_clust, file.path(dir_tmp, str_glue("{dataset}_supertype_dend_inh.rds")))

  var_genes_all <- rowVars(aggr) |> sort(decreasing = TRUE)
  var_genes_all <- names(var_genes_all)[1:4000]

  ct_dist <- as.dist(1 - cor(aggr[var_genes_all, ]))
  ct_clust <- hclust(ct_dist, method = "average") |>
    as.dendrogram() |>
    dendrapply(color_leaf, colors = colors_all)

  par(cex = 0.5) # Can't set cex in plot() directly for dendrograms
  plot(ct_clust)

  saveRDS(ct_clust, file.path(dir_tmp, str_glue("{dataset}_supertype_dend_all.rds")))

  # Dendrogram of sub classes instead of supertypes
  aggr_sce <- scuttle::aggregateAcrossCells(
    sce,
    ids = paste(sce$sample, "::", sce$sub_class_unmerged),
    statistics = c("sum")
  ) |>
    scuttle::logNormCounts(transform = "log", assay.type = "counts")

  aggr_sce <- scuttle::aggregateAcrossCells(
    aggr_sce,
    ids = aggr_sce$sub_class_unmerged,
    statistics = c("mean"),
    use.assay.type = "logcounts"
  )

  aggr <- assay(aggr_sce, "logcounts")

  colors_sub <- viridis::turbo(n = length(aggr_sce$sub_class_unmerged))
  names(colors_sub) <- as.character(aggr_sce$sub_class_unmerged)

  # Re-using this variable name so color_leaf works
  supertype_map <- as.character(aggr_sce$sub_class_unmerged)
  names(supertype_map) <- as.character(aggr_sce$sub_class_unmerged)

  var_genes_all <- rowVars(aggr) |> sort(decreasing = TRUE)
  var_genes_all <- names(var_genes_all)[1:4000]

  ct_dist <- as.dist(1 - cor(aggr[var_genes_all, ]))
  ct_clust <- hclust(ct_dist, method = "average") |>
    as.dendrogram() |>
    dendrapply(color_leaf, colors = colors_sub)
  par(cex = 1)
  plot(ct_clust)

  dev.off()
}
