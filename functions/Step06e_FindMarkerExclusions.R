# Finds genes which are significantly affected by diagnosis, making them bad
# candidates for cell marker genes. This script is very similar to Step06d,
# except that it compares differential expression between diagnoses within
# each celltype. A gene is marked as significantly affected if for any cell
# type, for any diagnosis vs the other diagnoses:
#   1) It is significantly up- or down-regulated in the cell type + diagnosis at p <= 0.05
#   2) The log fold-change is > 0.5 or < -0.5
#
# We don't run multiple datasets in parallel for this script because one call to
# DESeq() needs and will use all available CPUs.

library(SummarizedExperiment)
library(DESeq2)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

FindMarkerExclusions <- function(datasets, granularities) {
  for (dataset in datasets) {
    if (dataset == "seaRef") { # seaRef only has control samples
      next
    }

    for (granularity in granularities) {
      message(str_glue("Finding markers affected by diagnosis for {dataset} / {granularity}..."))

      pb <- Load_PseudobulkPureSamples(dataset, granularity,
                                       output_type = "counts")

      dds <- DESeqDataSet(pb, design = ~ diagnosis * celltype)
      dds <- DESeq(dds, test = "Wald",
                   fitType = "parametric",
                   sfType = "poscounts")

      mod <- model.matrix(design(dds), data = colData(pb))

      dds_res <- lapply(levels(pb$celltype), function(ct) {
        res_ct <- lapply(levels(pb$diagnosis), function(dg) {
          # This provides a contrast for (ct and dg) vs (ct and not dg)
          cont <- colMeans(mod[dds$celltype == ct & dds$diagnosis == dg, ]) -
            colMeans(mod[dds$celltype == ct & dds$diagnosis != dg, ])

          res <- results(dds, contrast = cont, alpha = 0.05)
          res <- subset(res, padj <= 0.05 & abs(log2FoldChange) > 0.5)

          if (nrow(res) == 0) {
            return(NULL)
          }

          res <- data.frame(res) %>% dplyr::arrange(desc(log2FoldChange))
          res$celltype <- ct
          res$diagnosis <- dg
          res$gene <- rownames(res)
          return(res)
        })

        return(do.call(rbind, res_ct))
      })

      dds_res <- do.call(rbind, dds_res)
      genes <- unique(dds_res$gene)

      print(str_glue("{length(genes)} affected genes found for {dataset} / {granularity}"))

      genes_list <- list("genes" = genes,
                         "deseq_results" = dds_res)

      saveRDS(genes_list, file.path(dir_metadata,
                                    str_glue("{dataset}_{granularity}_excluded_genes.rds")))

      rm(dds)
      gc()
    }
  }
}
