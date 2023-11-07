library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

batch_vars <- list("Mayo" = "flowcell",
                   "MSBB" = "sequencingBatch",
                   "ROSMAP" = "final_batch")

for (dataset in datasets) {
  bulk_raw_counts <- Load_BulkData(dataset, output_type = "log_cpm",
                                   regression_method = "none")

  bulk_deseq2 <- Load_BulkData(dataset, output_type = "log_cpm",
                               regression_method = "deseq2")

  bulk_dream <- Load_BulkData(dataset, output_type = "log_cpm",
                              regression_method = "dream")

  bulk_edger <- Load_BulkData(dataset, output_type = "log_cpm",
                              regression_method = "edger")

  pca_raw <- prcomp(t(assay(bulk_raw_counts, "counts")))
  pca_deseq2 <- prcomp(t(assay(bulk_deseq2, "counts")))
  pca_dream <- prcomp(t(assay(bulk_dream, "counts")))
  pca_edger <- prcomp(t(assay(bulk_edger, "counts")))

  metadata <- colData(bulk_raw_counts)
  covariates <- Load_Covariates(dataset)
  metadata <- Clean_BulkCovariates(dataset, metadata, covariates)

  metadata$batch <- metadata[,batch_vars[[dataset]]]

  df_raw <- merge(pca_raw$x[,1:2], metadata, by = "row.names", sort = FALSE)
  df_deseq2 <- merge(pca_deseq2$x[,1:2], metadata, by = "row.names", sort = FALSE)
  df_dream <- merge(pca_dream$x[,1:2], metadata, by = "row.names", sort = FALSE)
  df_edger <- merge(pca_edger$x[,1:2], metadata, by = "row.names", sort = FALSE)

  plt1 = ggplot(df_raw, aes(x = PC1, y = PC2, color = batch)) +
            geom_point() + ggtitle("Raw counts")
  plt2 = ggplot(df_deseq2, aes(x = PC1, y = PC2, color = batch)) +
            geom_point() + ggtitle("DESeq2")
  plt3 = ggplot(df_dream, aes(x = PC1, y = PC2, color = batch)) +
            geom_point() + ggtitle("Dream")
  plt4 = ggplot(df_edger, aes(x = PC1, y = PC2, color = batch)) +
            geom_point() + ggtitle("Edger")
  plt5 = ggplot(df_raw, aes(x = PC1, y = PC2, color = tissue)) +
            geom_point() + ggtitle("Raw counts")
  plt6 = ggplot(df_deseq2, aes(x = PC1, y = PC2, color = tissue)) +
            geom_point() + ggtitle("DESeq2")
  plt7 = ggplot(df_dream, aes(x = PC1, y = PC2, color = tissue)) +
            geom_point() + ggtitle("Dream")
  plt8 = ggplot(df_edger, aes(x = PC1, y = PC2, color = tissue)) +
            geom_point() + ggtitle("Edger")

  print(wrap_plots(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, nrow = 2))

  plt1 = ggplot(df_raw, aes(x = PC1, y = PC2, color = diagnosis)) +
    geom_point() + ggtitle("Raw counts")
  plt2 = ggplot(df_deseq2, aes(x = PC1, y = PC2, color = diagnosis)) +
    geom_point() + ggtitle("DESeq2")
  plt3 = ggplot(df_dream, aes(x = PC1, y = PC2, color = diagnosis)) +
    geom_point() + ggtitle("Dream")
  plt4 = ggplot(df_edger, aes(x = PC1, y = PC2, color = diagnosis)) +
    geom_point() + ggtitle("Edger")

  plt5 = ggplot(df_raw, aes(x = PC1, y = PC2, color = sex)) +
    geom_point() + ggtitle("Raw counts")
  plt6 = ggplot(df_deseq2, aes(x = PC1, y = PC2, color = sex)) +
    geom_point() + ggtitle("DESeq2")
  plt7 = ggplot(df_dream, aes(x = PC1, y = PC2, color = sex)) +
    geom_point() + ggtitle("Dream")
  plt8 = ggplot(df_edger, aes(x = PC1, y = PC2, color = sex)) +
    geom_point() + ggtitle("Edger")
  print(wrap_plots(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, nrow = 2))
}
