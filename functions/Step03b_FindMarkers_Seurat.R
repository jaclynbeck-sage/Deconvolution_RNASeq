# Finds marker genes for each cell type using Seurat / MAST. A gene is
# considered a marker gene if:
#   1) It is significantly up-regulated in the cell type at p <= 0.01
#   2) The log fold-change is > 0.5
#   3) The gene is not significantly up-regulated in any other cell type
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

FindMarkers_Seurat <- function(datasets, granluarities) {
  for (dataset in datasets) {
    for (granularity in granularities) {
      print(str_glue("Finding markers for {dataset} / {granularity}..."))

      sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

      seurat <- as.Seurat(sce, data = NULL)
      rm(sce)

      seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

      Idents(seurat) <- seurat$celltype
      markers <- FindAllMarkers(seurat, logfc.threshold = 0.5,
                                only.pos = TRUE, test.use = "MAST",
                                latent.vars = c("nCount_originalexp"))
      markers <- subset(markers, p_val_adj <= 0.01)

      # Remove marker genes that are present in more than one cell type, so that
      # each cell type has a set of unique marker genes
      dupes <- markers$gene[duplicated(markers$gene)]
      markers <- subset(markers, !(gene %in% dupes))

      print(str_glue("Markers for {dataset} / {granularity} cell types:"))
      print(table(markers$cluster))

      # Format in a way dtangle/hspe understand
      markers <- markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC))

      markers_list <- sapply(levels(markers$cluster), function(ct) {
        m <- subset(markers, cluster == ct)
        return(m$gene)
      })

      saveRDS(markers_list, file.path(dir_markers,
                                      str_glue("seurat_markers_{dataset}_{granularity}.rds")))

      rm(seurat)
      gc()
    }
  }
}
