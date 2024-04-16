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
# "autogenes_env" that has been set up with all needed packages
# (see Step00_InitialSetupInstall.R).

library(reticulate)
library(SingleCellExperiment)
library(anndata)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

use_virtualenv("autogenes_env")

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

      sc$pp$filter_cells(adata, min_genes = 200)
      sc$pp$filter_genes(adata, min_cells = 10)
      sc$pp$highly_variable_genes(adata, flavor = "seurat_v3", n_top_genes = 6000)
      sc$pp$normalize_total(adata, target_sum = 1e6) # cpm

      ag$init(adata, use_highly_variable = TRUE, celltype_key = "celltype")
      ag$optimize(ngen = 5000L, seed = 0L, mode = "standard",
                  offspring_size = 100L, verbose = FALSE)

      wts <- list("correlation" = c(-1, 0),
                  "combined" = c(-1, 1),
                  "distance" = c(0, 1))

      # Get the marker genes for each type of weight, rank them by fold-change
      # between the target cell type and other cell types, and put them in ranked
      # order in a named list, separated by cell type.
      for (key in names(wts)) {
        wt <- wts[[key]]
        inds <- ag$select(weights = wt)
        sc_means <- ag$adata()$transpose()
        sc_means <- sc_means[inds]$copy()

        X <- sc_means$X
        markers <- sc_means$obs_names

        # Which cell type has the highest expression for each gene
        maxs <- apply(X, 1, which.max)
        log2FC <- apply(log2(X + 1), 1, function(row) {
          sorted <- sort(row, decreasing = TRUE)
          return(sorted[1] - sorted[2]) # Log space is subtraction
        })

        sorted_logfc <- data.frame(celltype = sc_means$var_names[maxs],
                                   log2FC = log2FC,
                                   gene = markers) %>%
          dplyr::arrange(desc(log2FC))

        markers_filt <- subset(sorted_logfc, log2FC >= 1)

        print(str_glue("Markers for {dataset} / {granularity} cell types ({key}):"))
        print(table(markers_filt$celltype))

        markers_list <- sapply(sort(unique(sorted_logfc$celltype)), function(ct) {
          return(sorted_logfc$gene[sorted_logfc$celltype == ct])
        })

        markers_filt_list <- sapply(markers_list, function(M) {
          return(intersect(M, markers_filt$gene))
        })

        # Save all the markers just for reference, but we also want a filtered
        # list with genes where the log2FC between the highest and
        # second-highest expressing cell type is >= 1
        list_final <- list("all" = markers_list,
                           "filtered" = markers_filt_list)

        Save_Markers(list_final, dataset, granularity,
                     marker_type = "autogenes",
                     marker_subtype = key)
      }
    }
  }
}
