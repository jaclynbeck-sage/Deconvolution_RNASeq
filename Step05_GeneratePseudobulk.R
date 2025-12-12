# This script generates pseudobulk data sets from single cell data sets. The
# following pseudobulk data sets are created for both broad and fine cell types:
#   1. Pseudobulk of each sample
#   2. Pseudobulk of each cell type + sample combination, to create "pure"
#      pseudobulk samples of each cell type
#   3. Simulated bulk data that is pseudobulked from a random subset of cells,
#      which is used for Scaden only
#
# In the first two cases, the raw counts are summed together for each group, TMM
# factors are computed for each new sample, and the percentage of RNA and
# percentage of cells that contributed to each pseudobulked sample are
# calculated. The data is saved as a SummarizedExperiment object.
#
# For the last case, we generate a set of random cell type proportions for 1000
# simulated samples, randomly select cells from each cell type to meet those
# proportions, and sum the CPM or TMM values of those cells together to create
# each simulated sample. The data is saved as an anndata/h5ad file in the format
# Scaden expects.

library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(stringr)
library(dplyr)
library(edgeR)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- all_singlecell_datasets()

for (dataset in datasets) {
  sce <- Load_SingleCell(dataset, "broad_class", normalization = "counts")

  ## Pseudobulk by sample ------------------------------------------------------

  message(str_glue("{dataset}: creating pseudobulk by sample..."))
  pb <- scuttle::aggregateAcrossCells(sce, ids = sce$sample, statistics = "sum")

  colData(pb) <- colData(pb)[, c("sample", "diagnosis", "ncells")]
  pb$tmm_factors <- edgeR::normLibSizes(counts(pb), method = "TMMwsp")

  # The counts for broud and sub type pseudobulk sets are the same, only the
  # metadata changes. But we create each as a separate file to make looping
  # easier further down the pipeline.
  for (granularity in c("broad_class", "sub_class")) {
    sce$celltype <- colData(sce)[, granularity]
    propCells <- table(sce$sample, sce$celltype)
    propCells <- sweep(propCells, 1, rowSums(propCells), "/")

    pctRNA <- CalculatePercentRNA(sce, granularity)
    metadata(pb) <- list("propCells" = propCells,
                         "pctRNA" = pctRNA)

    Save_Pseudobulk(pb, dataset, "sc_samples", granularity)
  }


  ## Pseudobulk by sample x celltype -------------------------------------------

  message(str_glue("{dataset}: creating pseudobulk pure samples..."))

  # Here, the counts are different between broad and sub class
  for (granularity in c("broad_class", "sub_class")) {
    sce$celltype <- colData(sce)[, granularity]

    pb <- scuttle::aggregateAcrossCells(sce,
                                        ids = paste(sce$sample, sce$celltype),
                                        statistics = "sum")

    # Reassign "sample" to be the new ID, but preserve original sample values
    pb$sample_orig <- pb$sample
    pb$sample <- pb$ids

    colData(pb) <- colData(pb)[, c("sample", "sample_orig", "diagnosis", "celltype", "ncells")]

    pb$tmm_factors <- edgeR::normLibSizes(counts(pb), method = "TMMwsp")

    # There will be a single "1" value per row, no need to divide for percentages
    propCells <- table(pb$sample, pb$celltype)

    # Since these are pure samples, pctRNA and propCells = 1 where the cell type
    # matches the pure sample. pctRNA is therefore identical to propCells.
    metadata(pb) <- list("propCells" = propCells,
                         "pctRNA" = propCells)

    Save_PseudobulkPureSamples(pb, dataset, granularity)
  }


  ## Simulated pseudobulk for Scaden -------------------------------------------

  # This is functionally equivalent to `scaden_simulate`, but much faster and
  # allows for a seed to be set for reproducibility.
  message(str_glue("{dataset}: creating simulated pseudobulks for Scaden..."))

  for (granularity in c("broad_class", "sub_class")) {
    set.seed(sageRNAUtils::string_to_seed(paste(dataset, granularity, "scaden")))

    sce$celltype <- colData(sce)[, granularity]

    celltypes <- levels(sce$celltype)
    n_celltypes <- length(celltypes) # number of unique cell types
    n_samples <- 1000 # number of samples to generate
    n_cells <- 1000 # number of cells to sum for each sample

    # Generate 1000 pseudobulk samples -- 500 with all cell types and 500 that
    # are missing at least one cell type
    fractions <- runif(n_samples * n_celltypes) |>
      matrix(nrow = n_samples, ncol = n_celltypes,
             dimnames = list(paste0("s", 1:n_samples), celltypes))

    # How many cell types to drop for the dropout samples, always drop at least
    # 1 and keep at least one
    n_drop <- sample(1:(n_celltypes-1), size = n_samples / 2, replace = TRUE)

    # Keep sparse samples at the end of the dataset like Scaden does, even
    # though it might not affect anything
    sparse_start <- n_samples / 2
    for (ind in 1:length(n_drop)) {
      cts_drop <- sample(celltypes, size = n_drop[ind], replace = FALSE)
      fractions[sparse_start + ind, cts_drop] <- 0
    }

    # Force sample props to sum to 1
    fractions <- sweep(fractions, 1, rowSums(fractions), "/")

    # For each simulated sample, pick random cells from each cell type at the
    # fraction specified, and pseudobulk the result. Although the above
    # pseudobulk code adds counts together, the code below follows the method
    # used in the Scaden paper's pre-processing code, which normalizes every
    # cell by library size and scales the values to the median library size,
    # before adding the values together.
    med_lib <- median(colSums(counts(sce)))
    cpm_mat <- scuttle::calculateCPM(sce) * med_lib / 1e6
    tmm_mat <- scuttle::normalizeCounts(cpm_mat, sce$tmm_factors,
                                        log = FALSE,
                                        center.size.factors = FALSE)

    fractions_real <- fractions

    # Use identical cells for both cpm and tmm pseudobulking
    pb_cpm <- matrix(0, nrow = nrow(sce), ncol = n_samples,
                     dimnames = list(rownames(sce), rownames(fractions)))

    pb_tmm <- matrix(0, nrow = nrow(sce), ncol = n_samples,
                     dimnames = list(rownames(sce), rownames(fractions)))

    for (samp in rownames(fractions)) {
      # How many cells from each cell type to sample
      sample_num <- round(fractions[samp, ] * n_cells)

      cells <- lapply(celltypes, function(ct) {
        sample(colnames(sce)[sce$celltype == ct], size = sample_num[ct], replace = TRUE)
      }) |>
        unlist()

      pb_cpm[, samp] <- rowSums(cpm_mat[, cells])
      pb_tmm[, samp] <- rowSums(tmm_mat[, cells])

      # Adjust the fractions to reflect the actual number of cells that went into
      # the sample
      fractions_real[samp, ] <- sample_num / sum(sample_num)
    }

    fractions_real <- as.data.frame(fractions_real)

    # Turn this into a SummarizedExperiment. Preserve the original cell type
    # names as a metadata variable, because they get modified with make.names()
    # when fractions_real is converted to a DataFrame.
    pb_se <- SummarizedExperiment(list(counts = pb_cpm),
                                  colData = DataFrame(fractions_real),
                                  metadata = list(celltype_names = colnames(fractions)))

    save_SimulatedScadenData(pb_se, dataset, granularity, "cpm")

    # Save a separate file for TMM data
    assay(pb_se, "counts") <- pb_tmm
    save_SimulatedScadenData(pb_se, dataset, granularity, "tmm")
  }
}
