# Runs AutogeneS using reticulate to find marker genes for each cell type.
#
# This script will save 3 different versions of markers from one run of
# AutogeneS:
#   1) The set of markers with the lowest expression correlation between cell types
#   2) The set of markers with the highest expression distance between cell types
#   3) The set of markers with the best score when lowest correlation and highest
#      distance are equally weighted
#
# This script assumes there is already a python virtual environment called
# "r-omnideconv" that has been set up with all needed packages
# (see Step00_InitialSetupInstall.R).

library(reticulate)
library(SingleCellExperiment)
library(anndata)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

use_virtualenv("r-omnideconv")

FindMarkers_AutogeneS <- function(datasets, granularities) {
  for (dataset in datasets) {
    for (granularity in granularities) {
      message(str_glue("Finding markers for {dataset} / {granularity}..."))

      sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

      sc <- import("scanpy")
      ag <- import("autogenes")

      adata <- AnnData(X = t(counts(sce)),
                       obs = data.frame(colData(sce)),
                       var = data.frame(rowData(sce)))

      rm(sce)

      sc$pp$filter_genes(adata, min_cells = 10)
      sc$pp$highly_variable_genes(adata, flavor = "seurat_v3", n_top_genes = 6000)
      sc$pp$normalize_total(adata, target_sum = 1e6) # cpm

      seed <- sageRNAUtils::string_to_seed(paste(dataset, granularity, "AutogeneS")) |>
        as.integer()

      ag$init(adata, use_highly_variable = TRUE, celltype_key = "celltype")
      ag$optimize(ngen = 5000L, seed = seed, mode = "standard",
                  offspring_size = 100L, verbose = FALSE)

      wts <- list("correlation" = c(-1, 0),
                  "combined" = c(-1, 1),
                  "distance" = c(0, 1))

      # Get the marker genes for each type of weight, rank them by fold-change
      # between the target cell type and other cell types, and put them in ranked
      # order in a named list, separated by cell type.
      for (key in names(wts)) {
        wt <- wts[[key]]
        inds <- ag$select(weights = wts[[key]])
        markers <- ag$adata()$var_names[inds]

        avgs <- Load_SignatureMatrix(dataset, granularity, "log_cpm")
        marker_results <- Get_QualityMarkers(avgs, markers, granularity)

        print(str_glue("Markers for {dataset} / {granularity} cell types ({key}):"))
        print(lengths(marker_results$filtered))

        list_final <- list("all" = marker_results$all,
                           "filtered" = marker_results$filtered)

        Save_Markers(list_final, dataset, granularity,
                     marker_type = "autogenes",
                     marker_subtype = key)
      }
    }
  }
}
