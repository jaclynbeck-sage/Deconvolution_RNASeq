# Finds marker genes for each cell type using DESeq2 on pseudobulk samples of
# individual cell types. A gene is considered a marker gene if:
#   1) It is significantly up-regulated in the cell type at p <= 0.05
#   2) The log fold-change is > 0.5
#   3) The gene is not significantly up-regulated in any other cell type
#   4) The log fold-change between the highest-expressing cell type and the
#      second-highest expressing cell type is >= 1.
#
# Two methods are implemented, "DESeq2" which uses the default Wald test with
# a parametric fit, and "glmGamPoi", which uses a likelihood ratio test fitted
# with a gamma poisson distribution.
#
# NOTE: For broad_class celltypes, the results from both fit types are identical
# except that more genes are filtered out as duplicates with "glmGamPoi". Since
# this method takes a much longer time to run than the Wald test, we do not use
# "glmGamPoi" for broad_class marker finding. Verification of results on
# sub_class marker finding is still pending.
#
# We don't run multiple datasets in parallel for this script because one call to
# DESeq() needs and will use all available CPUs.

library(SummarizedExperiment)
library(DESeq2)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

FindMarkers_DESeq2 <- function(datasets, granularities) {
  for (dataset in datasets) {
    for (granularity in granularities) {
      for (fit_type in c("DESeq2")) { #, "glmGamPoi")) {
        print(str_glue("Finding markers for {dataset} / {granularity} / {fit_type}..."))

        pb <- Load_PseudobulkPureSamples(dataset, granularity,
                                         output_type = "counts")

        # For most datasets, we model the interaction between diagnosis and
        # celltype as explanatory for gene expression. For the seaRef dataset,
        # there is no diagnosis, so the model has to be just '~ celltype'.
        if (length(levels(pb$diagnosis)) > 1) { # data with diagnosis info
          dds <- DESeqDataSet(pb, design = ~ diagnosis * celltype)
        }
        else { # seaRef
          dds <- DESeqDataSet(pb, design = ~ celltype)
        }

        if (fit_type == "DESeq2") {
          dds <- DESeq(dds, test = "Wald", fitType = "parametric",
                       sfType = "poscounts")
        }
        else { # glmGamPoi
          dds <- DESeq(dds, test = "LRT", fitType = "glmGamPoi",
                       sfType = "poscounts", reduced = ~1)
        }

        mod <- model.matrix(design(dds), data = colData(pb))

        dds_res <- lapply(levels(pb$celltype), function(ct) {
          # This provides a contrast for (ct) vs (all other ct), controlling for
          # the effect of diagnosis
          cont <- colMeans(mod[dds$celltype == ct,]) -
                    colMeans(mod[dds$celltype != ct,])

          res <- results(dds, contrast = cont, alpha = 0.05)
          res <- subset(res, padj <= 0.05 & log2FoldChange > 0.5)
          res <- data.frame(res) %>% arrange(desc(log2FoldChange))
          res$celltype <- ct
          res$gene <- rownames(res)
          return(res)
        })
        dds_res <- do.call(rbind, dds_res)

        dupes <- unique(dds_res$gene[duplicated(dds_res$gene)])
        unique_markers <- setdiff(dds_res$gene, dupes)

        # Create a list of markers per cell type, with genes that are only
        # differentially expressed in one cell type
        markers_all <- sapply(sort(unique(dds_res$celltype)), function(ct) {
          res <- subset(dds_res, celltype == ct & gene %in% unique_markers)
          return(res$gene)
        })

        # Find genes where the logFC of the highest-expressing cell type vs
        # the second-highest expressing cell type is > 1
        data_norm <- varianceStabilizingTransformation(dds, blind = TRUE)
        data_means <- sapply(levels(data_norm$celltype), function(ct) {
          rowMeans(assay(data_norm)[unique_markers, data_norm$celltype == ct])
        })

        data_means <- as.data.frame(data_means)

        getLog2FC <- function(cols) {
          sorted <- sort(cols, decreasing = TRUE)
          return(sorted[1] - sorted[2]) # Log space is subtraction
        }

        avgs <- data_means %>%
                  mutate(highest = colnames(data_means)[apply(., 1, which.max)],
                         log2FC = apply(., 1, getLog2FC),
                         gene = rownames(.)) %>%
                  subset(log2FC >= 1) %>% arrange(desc(log2FC)) %>%
                  select(gene, highest)

        # Filter genes to those where (log2FC between highest and second-highest
        # expression) > 1
        markers_filt <- sapply(colnames(data_means), function(ct) {
          return(avgs$gene[avgs$highest == ct])
        })

        print(str_glue("Markers for {dataset} / {granularity} / {fit_type} cell types:"))
        print(lengths(markers_filt))

        markers_list <- list("deseq_results" = dds_res,
                             "all" = markers_all,
                             "filtered" = markers_filt)

        saveRDS(markers_list,
                file.path(dir_markers,
                          str_glue("deseq2_markers_{dataset}_{granularity}_{fit_type}.rds")))

        rm(dds)
        gc()
      }
    }
  }
}
