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
library(zellkonverter)

source(file.path("functions", "FileIO_HelperFunctions.R"))

use_virtualenv("autogenes_env")

FindMarkers_AutogeneS <- function(datasets, granularities) {
  for (dataset in datasets) {
    for (granularity in granularities) {
      print(str_glue("Finding markers for {dataset} / {granularity}..."))

      sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

      sc <- import("scanpy")
      ag <- import("autogenes")

      adata <- SCE2AnnData(sce)
      rm(sce)

      sc$pp$filter_cells(adata, min_genes=200)
      sc$pp$filter_genes(adata, min_cells=10)
      sc$pp$highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=4000)
      sc$pp$normalize_total(adata, target_sum=1e4)

      ag$init(adata, use_highly_variable=TRUE, celltype_key='celltype')
      ag$optimize(ngen=5000L, seed=0L, mode='standard', offspring_size=100L,
                  verbose=FALSE)

      wts <- list("correlation" = c(-1,0),
                  "combined" = c(-1,1),
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
        markers <- sc_means$obs_names$to_list()

        # Which cell type has the highest expression for each gene
        maxs <- apply(X, 1, which.max)
        cts <- sc_means$var_names[maxs-1]$to_list() # 0-index for python

        print(str_glue("Markers for {dataset} / {granularity} cell types ({key}):"))
        print(table(cts))

        # This isn't exactly FC because it uses the mean of means, but it's good
        # enough for figuring out a relative ordering within each cell type
        fc <- lapply(1:length(maxs), function(M) {
          ct <- maxs[M]
          val1 <- X[M, ct]
          non_ct <- setdiff(1:ncol(X), ct)
          val2 <- mean(X[M, non_ct])
          return(val1 / max(val2, 0.001))
        })

        # Data frame ordered by fold-change -> collapsed to one row per cell type
        # with a list of genes for each cell type -> dict
        markers_df <- data.frame("gene" = markers,
                                 "celltype" = cts,
                                 "fc" = unlist(fc))
        markers_df <- markers_df %>% arrange(desc(fc))
        markers_list <- sapply(sort(unique(markers_df$celltype)), function(ct) {
          return(markers_df$gene[markers_df$celltype == ct])
        })

        saveRDS(markers_list, file.path(dir_markers,
                                        str_glue("autogenes_markers_{dataset}_{granularity}_{key}.rds")))
      }
    }
  }
}
