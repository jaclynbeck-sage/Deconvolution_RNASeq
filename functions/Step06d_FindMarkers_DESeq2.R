# Finds marker genes for each cell type using DESeq2 on pseudobulk samples of
# individual cell types. A gene is considered a marker gene if:
#   1) It is significantly up-regulated in the cell type at p <= 0.05
#   2) The log fold-change is > 0.5
#   3) The gene is not significantly up-regulated in any other cell type
#   4) The log fold-change between the highest-expressing cell type and the
#      second-highest expressing cell type is >= 1.

library(SummarizedExperiment)
library(DESeq2)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

FindMarkers_DESeq2 <- function(datasets, granularities, n_cores) {
  for (dataset in datasets) {
    for (granularity in granularities) {
      message(str_glue("Finding markers for {dataset} / {granularity}..."))

      pb <- Load_PseudobulkPureSamples(dataset, granularity,
                                       output_type = "counts")

      # Use only genes that express > 10 counts in > 10% of at least one cell
      # type's cells (or > 5 cells if the cell type has < 50 cells)
      over_10 <- assay(pb, "counts") > 10
      expr_genes <- sapply(levels(pb$celltype), function(ct) {
        rowSums(over_10[, pb$celltype == ct]) > max(0.1 * sum(pb$celltype == ct), 5)
      })

      expr_genes <- rownames(expr_genes)[rowSums(expr_genes) > 0]

      dds <- DESeqDataSet(pb[expr_genes, ], design = ~celltype)

      dds <- DESeq(dds, test = "Wald",
                   fitType = "parametric",
                   sfType = "poscounts",
                   parallel = (n_cores > 1),
                   BPPARAM = BiocParallel::SnowParam(n_cores))

      mod <- model.matrix(design(dds), data = colData(dds))

      # This can take several minutes per cell type even when running in parallel
      dds_res <- lapply(levels(pb$celltype), function(ct) {
        # This provides a contrast for (ct) vs (all other ct), controlling for
        # the effect of diagnosis
        cont <- colMeans(mod[dds$celltype == ct, ]) -
          colMeans(mod[dds$celltype != ct, ])

        res <- results(dds, contrast = cont, alpha = 0.05,
                       parallel = (n_cores > 1),
                       BPPARAM = BiocParallel::SnowParam(n_cores))

        res <- subset(res, padj <= 0.05 & log2FoldChange > 0.5)
        res <- data.frame(res) |> dplyr::arrange(desc(log2FoldChange))
        res$celltype <- ct
        res$gene <- rownames(res)
        return(res)
      })
      dds_res <- do.call(rbind, dds_res)

      # Find genes where the logFC of the highest-expressing cell type vs
      # the second-highest expressing cell type is > 1
      data_norm <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
      data_means <- sapply(levels(data_norm$celltype), function(ct) {
        tmp <- assay(data_norm)[unique(dds_res$gene), ]
        return(rowMeans(tmp[, data_norm$celltype == ct]))
      })

      res <- Get_QualityMarkers(data_means, dds_res, granularity)
      sorted_logfc <- res$sorted_logfc
      dds_res <- res$markers

      # Create a list of markers per cell type, with genes that are only
      # differentially expressed in one (or two) cell type
      markers_all <- sapply(sort(unique(dds_res$celltype)), function(ct) {
        res <- subset(dds_res, celltype == ct)
        return(res$gene)
      })

      # Filter genes to those where (log2FC between highest and second-highest
      # expression) > 1
      markers_filt <- sapply(colnames(data_means), function(ct) {
        return(sorted_logfc$gene[sorted_logfc$celltype == ct])
      })

      print(str_glue("Markers for {dataset} / {granularity} cell types:"))
      print(lengths(markers_filt))

      markers_list <- list("deseq_results" = dds_res,
                           "all" = markers_all,
                           "filtered" = markers_filt)

      Save_Markers(markers_list, dataset, granularity, marker_type = "deseq2")

      rm(dds)
      gc()
    }
  }
}
