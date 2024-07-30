# Finds genes which are significantly affected by diagnosis, making them bad
# candidates for cell marker genes. This script is very similar to Step06b,
# except that it compares differential expression between diagnoses within
# each celltype. A gene is marked as significantly affected if for any cell
# type, for any diagnosis vs the other diagnoses:
#   1) It is significantly up- or down-regulated in the cell type + diagnosis at p <= 0.05
#   2) The log fold-change is > 0.5 or < -0.5
#
# We don't run multiple datasets in parallel for this script because one call to
# FindAllMarkers() needs and will use all available CPUs.

library(SummarizedExperiment)
library(Seurat)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

FindMarkerExclusions <- function(datasets, granularities) {
  for (dataset in datasets) {
    if (dataset == "seaRef") { # seaRef only has control samples
      next
    }

    gene_list <- lapply(granularities, function(granularity) {
      message(str_glue("Finding markers affected by diagnosis for {dataset} / {granularity}..."))

      sce <- Load_SingleCell(dataset, granularity, "counts")
      seurat <- as.Seurat(sce, data = NULL)
      rm(sce)
      gc()

      seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

      ident_2 <- "Control"

      # Lau is the only data set that doesn't use "Control" as the name for
      # control samples
      if (dataset == "lau") {
        ident_2 <- "healthy control"
      }

      Idents(seurat) <- seurat$diagnosis

      s_obj <- SplitObject(seurat, split.by = "celltype")
      markers_all <- lapply(s_obj, FindMarkers,
                            ident.1 = "AD",
                            ident.2 = ident_2,
                            logfc.threshold = 0.5,
                            min.pct = 0.1,
                            test.use = "MAST",
                            latent.vars = c("nCount_originalexp"))

      markers_all <- lapply(markers_all, function(X) {
        X$gene <- rownames(X)
        return(X)
      })

      markers <- do.call(rbind, markers_all)
      markers <- subset(markers, p_val_adj <= 0.05)
      genes <- unique(markers$gene)

      print(str_glue("{length(genes)} affected genes found for {dataset} / {granularity}"))

      rm(seurat)
      gc()

      return(list("genes" = genes,
                  "mast_results" = markers))
    })

    names(gene_list) <- granularities

    final_list <- list("genes" = unique(unlist(sapply(gene_list, "[[", "genes"))),
                       "mast_results" = lapply(gene_list, "[[", "mast_results"))

    print(str_glue("{length(final_list$genes)} total affected genes found for {dataset}"))

    saveRDS(final_list, file.path(dir_metadata,
                                  str_glue("{dataset}_excluded_genes.rds")))
  }
}
