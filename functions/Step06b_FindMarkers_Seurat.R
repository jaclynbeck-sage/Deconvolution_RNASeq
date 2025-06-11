# Finds marker genes for each cell type using Seurat. A gene is considered a
# marker gene if:
#   1) It is significantly up-regulated in the cell type at p <= 0.05
#   2) The log fold-change is > 0.5
#   3) The gene is not significantly up-regulated in any other cell type (broad
#      class) or is up-regulated in at most one other cell type (sub class)
#   4) The log fold-change between the highest-expressing cell type and the
#      second-highest expressing cell type is >= 1 (broad class), or the log
#      fold-change between the highest-expressing cell type and the highest
#      expressing cell type for which the gene is not a marker is >= 0.25 (sub
#      class)

library(SingleCellExperiment)
library(Seurat)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

FindMarkers_Seurat <- function(datasets, granularities) {
  options(future.globals.maxSize = 5 * 1024^3) # 5 GB

  for (dataset in datasets) {
    for (granularity in granularities) {
      message(str_glue("Finding markers for {dataset} / {granularity}..."))

      sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

      seurat <- as.Seurat(sce, data = NULL)
      rm(sce)

      seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

      Idents(seurat) <- seurat$celltype
      markers <- FindAllMarkers(seurat,
                                logfc.threshold = 0.5,
                                min.pct = 0.1,
                                only.pos = TRUE)

      markers <- subset(markers, p_val_adj <= 0.05) |>
        dplyr::rename(celltype = cluster)

      avgs <- AggregateExpression(seurat,
                                  features = unique(markers$gene),
                                  return.seurat = TRUE) |>
        GetAssayData(layer = "data")

      # Filter to markers that have a large enough gap between the target cell
      # type and non-target cell types
      res <- Get_QualityMarkers(avgs, markers, granularity)
      sorted_logfc <- res$sorted_logfc
      markers <- res$markers

      # Create one full list containing all genes that were identified as
      # markers, and one list filtered to (log2FC between highest and
      # second-highest expression) > 1
      markers_all <- sapply(levels(markers$celltype), function(ct) {
        return(subset(markers, celltype == ct)$gene)
      })

      markers_filt <- sapply(levels(markers$celltype), function(ct) {
        return(sorted_logfc$gene[sorted_logfc$celltype == ct])
      })

      print(str_glue("Markers for {dataset} / {granularity} cell types:"))
      print(lengths(markers_filt))

      markers_list <- list("seurat_results" = markers,
                           "all" = markers_all,
                           "filtered" = markers_filt)

      Save_Markers(markers_list, dataset, granularity, marker_type = "seurat")

      rm(seurat)
      gc()
    }
  }
}
