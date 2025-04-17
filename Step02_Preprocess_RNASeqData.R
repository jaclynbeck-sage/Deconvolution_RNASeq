# Load metadata and counts matrices from Synapse and create either a
# SingleCellExperiment (single cell data) or a Summarized experiment (bulk data)
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol. Gene names and cell type labels are also
# standardized across data sets.
library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(edgeR)

source(file.path("functions", "Step02_Preprocess_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

# Setup ------------------------------------------------------------------------

datasets <- c("cain", "lau", "leng", "mathys", "seaRef", # Single cell
              "Mayo", "MSBB", "ROSMAP") # Bulk

# Helper functions
is_bulk <- function(dataset) {
  return(dataset %in% c("Mayo", "MSBB", "ROSMAP"))
}

is_singlecell <- function(dataset) {
  return(!is_bulk(dataset))
}

# Download files and process data ----------------------------------------------

for (dataset in datasets) {
  message(str_glue("Creating data set for {dataset}..."))
  files <- DownloadData(dataset)

  # Allow threading of operations on DelayedMatrix for seaRef. This formula
  # assumes 8 GB of RAM per CPU.
  if (dataset == "seaRef") {
    n_cores <- max(parallel::detectCores() / 2, 1)
    setAutoBPPARAM(BPPARAM = BiocParallel::SnowParam(n_cores))
  }

  ## Read in metadata file -----------------------------------------------------

  metadata_list <- ReadMetadata(dataset, files)
  metadata <- metadata_list$metadata

  # For data privacy reasons, each data set's covariates are saved in a separate
  # file from the main SingleCellExperiment object and are not uploaded anywhere
  covariates <- metadata_list$covariates
  Save_Covariates(dataset, covariates)

  if (is_singlecell(dataset)) {
    colnames(metadata) <- c("cell_id", "sample", "diagnosis",
                            "broad_class", "sub_class")
    rownames(metadata) <- metadata$cell_id
  } else { # bulk
    colnames(metadata) <- c("sample", "diagnosis", "tissue")
    rownames(metadata) <- metadata$sample
  }


  ## Read in matrix of counts --------------------------------------------------

  counts <- ReadCounts(dataset, files)

  samples <- intersect(rownames(metadata), colnames(counts))
  metadata <- metadata[samples, ]
  counts <- counts[, samples]


  ## Remove sample outliers (bulk only) ----------------------------------------

  if (is_bulk(dataset)) {
    outliers <- FindOutliers_BulkData(dataset, covariates, counts,
                                      do_plot = FALSE, sd_threshold = 4)

    print(str_glue("{length(outliers)} outlier samples will be removed from {dataset}."))
    counts <- counts[, setdiff(colnames(counts), outliers)]
    metadata <- metadata[colnames(counts),]
  }


  ## Convert gene names --------------------------------------------------------

  # Bulk only -- Convert bulk data Ensembl IDs to gene symbols.
  if (is_bulk(dataset) || dataset == "mathys") {
    genes <- EnsemblIdToGeneSymbol(rownames(counts))
  } else {
    # Single cell only -- update potentially outdated gene symbols to the most
    # current version possible, and get the matching Ensembl IDs too.
    genes <- UpdateGeneSymbols(dataset, rownames(counts))
  }

  # Some symbols are not unique so we use the first entry in the list for
  # duplicate symbols, which is what Seurat does and is probably consistent with
  # most of the single cell data sets.
  dupes <- duplicated(genes$hgnc_symbol)
  genes <- genes[!dupes, ]

  # Puts 'genes' in the same gene order as counts so the counts matrix doesn't
  # get rearranged unnecessarily
  genes <- genes[intersect(rownames(counts), rownames(genes)), ]

  # Assign rownames to be the canonical hgnc symbol, applies to both bulk and
  # single cell
  counts <- counts[rownames(genes), ]
  rownames(counts) <- genes$hgnc_symbol
  rownames(genes) <- genes$hgnc_symbol


  ## Adjust for and remove mitochondrial/non-coding genes ----------------------

  # Keep track of original library size before excluding mito/nc genes
  metadata$lib_size <- colSums(counts)

  mt_genes <- grepl("^MT-", rownames(counts))
  pct_mt <- colSums(counts[mt_genes, ]) / metadata$lib_size
  metadata$percent_mito <- pct_mt

  # Remove non-coding genes -- grep pattern from Green et al 2023.
  nc_genes <- grepl("^(AC\\d+{3}|AL\\d+{3}|AP\\d+{3}|LINC\\d+{3})", rownames(counts))
  pct_nc <- colSums(counts[nc_genes, ]) / metadata$lib_size
  metadata$percent_noncoding <- pct_nc

  genes$exclude <- mt_genes | nc_genes

  if (is_singlecell(dataset)) {
    # Since all of these are single nucleus datasets, we expect the percent of
    # mitochondrial expression to be near-zero. We set the threshold to 0.05 to
    # allow for a small amount of noise or minor non-nuclear RNA contamination.
    # Ideally this threshold would be closer to 0.02 or 0.03, but some of the
    # datasets are too noisy for this to be realistic. The seaRef dataset is
    # pre-capped at 0.05%, so this is the threshold we chose.

    # TODO move to top of script as settable variable
    # Assumes a ratio of 8 GB per 1 CPU. This formula should be safe for all
    # data sets.
    n_cores <- max(parallel::detectCores() - 2, 1)
    counts <- QC_SingleCell(metadata, counts,
                            mt_threshold = 0.05,
                            dataset_name = dataset,
                            n_cores = n_cores)
  } else {
    # Remove bulk samples with > 35% mito genes by count. The percentages in
    # Mayo and ROSMAP especially can get pretty high and these samples should be
    # removed as low-quality. The 35% threshold was chosen by looking at the
    # value of median(pct_mt) + 3*mad(pct_mt), which is near 0.35 for both Mayo
    # and ROSMAP (MSBB is much lower), and visual inspection of the distribution
    # of pct_mt for each data set.
    mt_threshold <- 0.35

    if (any(pct_mt > mt_threshold)) {
      print(str_glue(paste("Removing {sum(pct_mt > mt_threshold)} samples from",
                           "{dataset} due to high mitochondrial gene expression.")))
    }
    counts <- counts[, pct_mt <= mt_threshold] # Exclude samples with high mitochondrial genes
  }

  # Remove genes that are expressed in less than 3 cells (or samples) after
  # filtering for outliers and high mitochondrial percentages
  ok <- rowSums(counts > 0) >= 3
  genes$exclude <- genes$exclude | !ok


  ## Final modifications to metadata -------------------------------------------

  # Make sure metadata has the same samples and is in the same order as counts
  metadata <- metadata[colnames(counts), ]

  #if (is_singlecell(dataset)) {
  #  metadata$broad_class <- RemapCelltypeNames(metadata$broad_class)
  #}

  # Ensure "sample" is a factor and not numeric (fixes mathys samples)
  metadata$sample <- factor(metadata$sample)

  for (col in colnames(metadata)) {
    if (!is.numeric(metadata[, col]) && col != "cell_id") {
      metadata[, col] <- factor(metadata[, col])
    }
  }

  # TMM normalization factors -- bulk only. Single cell factors will be
  # calculated after cell type mapping.
  if (is_bulk(dataset)) {
    tmm <- edgeR::calcNormFactors(counts[!genes$exclude, ],
                                  lib.size = metadata$lib_size,
                                  method = "TMM")
    metadata$tmm_factors <- tmm
  }

  gc()


  ## Bulk data -- create SummarizedExperiment and save -------------------------

  if (is_bulk(dataset)) {
    se <- SummarizedExperiment(assays = list(counts = counts),
                               colData = metadata,
                               rowData = genes)

    Save_PreprocessedData(dataset, se)
  }

  ## Single cell data -- create SingleCellExperiment object and save -----------

  else {
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = metadata,
                                rowData = genes)

    # If the counts matrix is a DelayedArray (i.e. from the SEA-AD files), the
    # sce file will contain a pointer to the original data file rather than
    # writing the full data to disk again
    Save_PreprocessedData(dataset, sce)

    # Save a metadata-free copy as an h5ad file for mapping. seaRef is already
    # mapped and doesn't need this step
    if (dataset != "seaRef") {
      ad <- import("anndata")

      obs <- data.frame(cell_id = metadata$cell_id,
                        row.names = metadata$cell_id)
      var <- data.frame(ensembl_gene_id = genes$ensembl_gene_id,
                        row.names = genes$ensembl_gene_id)

      adata <- ad$AnnData(X = t(counts),
                          obs = obs,
                          var = var)
      adata$write_h5ad(filename = file.path(dir_preprocessed,
                                            str_glue("{dataset}_preprocessed.h5ad")),
                       compression = "gzip")
    }
  }
}
