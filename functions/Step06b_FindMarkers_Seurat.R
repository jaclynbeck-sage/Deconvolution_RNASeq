# Finds marker genes for each cell type using Seurat / MAST. A gene is
# considered a marker gene if:
#   1) It is significantly up-regulated in the cell type at p <= 0.01
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
      print(str_glue("Finding markers for {dataset} / {granularity}..."))

      sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

      seurat <- as.Seurat(sce, data = NULL)
      rm(sce)

      seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

      Idents(seurat) <- seurat$celltype
      markers <- FindAllMarkers(seurat, logfc.threshold = 0.5,
                                only.pos = TRUE, test.use = "MAST",
                                latent.vars = c("nCount_originalexp"))
      markers <- subset(markers, p_val_adj <= 0.05)
      dupes <- markers$gene[duplicated(markers$gene)]

      markers2 <- subset(markers, !(gene %in% dupes))

      avgs <- AverageExpression(seurat, features = markers2$gene, slot = "data")
      avgs <- as.data.frame(avgs[[1]])

      getLog2FC <- function(cols) {
        sorted <- sort(cols, decreasing = TRUE)
        return(sorted[1] - sorted[2]) # Log space is subtraction
      }

      # Mark which cell type has the highest expression for each gene, get the
      # log2-fold-change between the highest expression and the 2nd-highest
      # expression for each gene, sort the dataframe with highest log2FC first.
      avgs2 <- avgs %>% mutate(highest = colnames(avgs)[apply(., 1, which.max)],
                               log2FC = apply(., 1, getLog2FC),
                               gene = rownames(.)) %>%
                subset(log2FC >= 1) %>% arrange(desc(log2FC)) %>%
                select(gene, highest)

      # Format in a way dtangle/HSPE understand -- one full list containing all
      # genes that are unique to one cell type, one list filtered to (log2FC
      # between highest and second-highest expression) > 1
      markers_all <- sapply(levels(markers2$cluster), function(ct) {
        return(subset(markers2, cluster == ct)$gene)
      })
      names(markers_all) <- levels(markers2$cluster)

      markers_filt <- sapply(colnames(avgs), function(ct) {
        return(avgs2$gene[avgs2$highest == ct])
      })

      print(str_glue("Markers for {dataset} / {granularity} cell types:"))
      print(lengths(markers_filt))

      markers_list <- list("mast_results" = markers,
                           "all" = markers_all,
                           "filtered" = markers_filt)

      saveRDS(markers_list, file.path(dir_markers,
                                      str_glue("seurat_markers_{dataset}_{granularity}.rds")))

      rm(seurat)
      gc()
    }
  }
}
