# Finds marker genes for each cell type using Seurat / MAST. A gene is
# considered a marker gene if:
#   1) It is significantly up-regulated in the cell type at p <= 0.05
#   2) The log fold-change is > 0.5
#   3) The gene is not significantly up-regulated in any other cell type
#   4) The log fold-change between the highest-expressing cell type and the
#      second-highest expressing cell type is >= 1.
#
# We use Seurat's FindAllMarkers instead of calling MAST directly, to take
# advantage of Seurat's built-in pre-filtering of genes that don't meet cluster
# expression percentage or logFC thresholds prior to calling MAST.
#
# We don't run multiple datasets in parallel for this script because one call to
# MAST needs and will use all available CPUs.

library(SingleCellExperiment)
library(Seurat)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

FindMarkers_Seurat <- function(datasets, granularities) {
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
                                only.pos = TRUE, test.use = "MAST",
                                latent.vars = c("nCount_originalexp"))

      markers <- subset(markers, p_val_adj <= 0.05)
      dupes <- markers$gene[duplicated(markers$gene)]

      markers <- subset(markers, !(gene %in% dupes))

      avgs <- AggregateExpression(seurat,
                                  features = unique(markers$gene),
                                  group.by = "celltype",
                                  return.seurat = TRUE)
      avgs <- GetAssayData(avgs, layer = "data") |>
        as.data.frame()

      getLog2FC <- function(cols) {
        sorted <- sort(cols, decreasing = TRUE)
        return(sorted[1] - sorted[2]) # Log space is subtraction
      }

      # Mark which cell type has the highest expression for each gene, get the
      # log2-fold-change between the highest expression and the 2nd-highest
      # expression for each gene, sort the dataframe with highest log2FC first.
      sorted_logfc <- avgs %>%
        dplyr::mutate(highest = colnames(avgs)[apply(., 1, which.max)],
                      log2FC = apply(., 1, getLog2FC),
                      gene = rownames(.)) %>%
        subset(log2FC >= 1) %>%
        dplyr::arrange(desc(log2FC)) %>%
        select(gene, highest)

      # Create one full list containing all genes that are unique to one cell
      # type, and one list filtered to (log2FC between highest and
      # second-highest expression) > 1
      markers_all <- sapply(levels(markers$cluster), function(ct) {
        return(subset(markers, cluster == ct)$gene)
      })
      names(markers_all) <- levels(markers$cluster)

      markers_filt <- sapply(colnames(avgs), function(ct) {
        return(sorted_logfc$gene[sorted_logfc$highest == ct])
      })

      print(str_glue("Markers for {dataset} / {granularity} cell types:"))
      print(lengths(markers_filt))

      markers_list <- list("mast_results" = markers,
                           "all" = markers_all,
                           "filtered" = markers_filt)

      Save_Markers(markers_list, dataset, granularity, marker_type = "seurat")

      rm(seurat)
      gc()
    }
  }
}
