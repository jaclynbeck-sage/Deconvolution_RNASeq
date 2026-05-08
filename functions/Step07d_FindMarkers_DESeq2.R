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
                                       normalization = "counts")

      # Use only genes that express > 10 counts in at least 3 samples
      over_10 <- rowSums(assay(pb, "counts") > 10)
      expr_genes <- names(over_10)[over_10 >= 3]

      # For most datasets, we model the interaction between diagnosis and
      # celltype as explanatory for gene expression. For the seaRef dataset,
      # there is no diagnosis, so the model has to be just '~ group'.
      if (dataset == "seaRef") {
        formula <- as.formula("~ group")
      } else {
        formula <- as.formula("~ group * diagnosis")
      }

      dds_base <- DESeqDataSet(pb[expr_genes, ], design = ~ 1)

      # To prevent warnings printing out
      dds_base$diagnosis <- factor(make.names(as.character(dds_base$diagnosis)))

      dds_res <- lapply(levels(pb$celltype), function(ct) {
        print(str_glue("Getting results for {dataset} / {ct}..."))

        dds <- dds_base
        dds$group <- factor(paste0("celltype_", dds$celltype == ct))
        design(dds) <- formula

        dds <- DESeq(dds, test = "Wald",
                     fitType = "parametric",
                     sfType = "poscounts",
                     quiet = TRUE,
                     parallel = (n_cores > 1),
                     BPPARAM = BiocParallel::MulticoreParam(n_cores))

        # Differences between this cell type and all other cells, controlling
        # for diagnosis (if applicable)
        res <- results(dds,
                       contrast = c("group", "celltype_TRUE", "celltype_FALSE"),
                       alpha = 0.05,
                       parallel = (n_cores > 1),
                       BPPARAM = BiocParallel::MulticoreParam(n_cores))

        res <- subset(res, padj <= 0.05 & log2FoldChange >= 0.5)
        res <- data.frame(res) |> dplyr::arrange(desc(log2FoldChange))
        res$celltype <- ct
        res$gene <- rownames(res)
        return(res)
      })

      dds_res <- do.call(rbind, dds_res)

      avgs <- Load_SignatureMatrix(dataset, granularity, "log_cpm")

      marker_results <- Get_QualityMarkers(avgs, dds_res$gene, granularity)

      print(str_glue("Markers for {dataset} / {granularity} cell types:"))
      print(lengths(marker_results$filtered))

      all_genes <- as.character(unlist(marker_results$all))

      markers_list <- list(
        "genes" = all_genes,
        "deseq_results" = dds_res,
        "all" = lapply(marker_results$all, match, all_genes), # index into all_genes
        "filtered" = lapply(marker_results$filtered, match, all_genes) # index into all_genes
      )

      Save_Markers(markers_list, dataset, granularity, marker_type = "deseq2")
    }
  }
}
