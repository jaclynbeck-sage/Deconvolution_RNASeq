# This script generates pseudobulk data sets from single cell data sets. The
# following pseudobulk data sets are created for both broad and fine cell types:
#   1. Pseudobulk of each sample
#   2. Pseudobulk of each cell type + sample combination, to create "pure"
#      pseudobulk samples of each cell type
#
# In each case, the raw counts are summed together for each group, TMM factors
# are computed for each new sample, and the percentage of RNA and percentage of
# cells that contributed to each pseudobulked sample are calculated.

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
}
