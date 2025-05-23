# This script looks at correlations between the bulk data and its covariates.
# This is not part of the main pipeline. Rather, it is used to manually examine
# what the covariates look like.
library(SummarizedExperiment)
library(scuttle)
library(corrplot)
library(dplyr)
library(sageseqr)

# TODO this is really out of date

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

batch_vars <- list("Mayo" = "flowcell",
                   "MSBB" = "sequencingBatch",
                   "ROSMAP" = "final_batch")

for (dataset in datasets) {
  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)

  covariates <- Load_Covariates(dataset)
  metadata <- Clean_BulkCovariates(dataset, colData(bulk), covariates)

  metadata <- data.frame(metadata[colnames(bulk),]) %>%
                select(-sample, -percent_noncoding, -tmm_factors)

  expr_norm <- calculateCPM(bulk)
  expr_norm <- normalizeCounts(expr_norm, bulk$tmm_factors,
                                 log = TRUE, center.size.factors = FALSE)

  # Association/correlation of covariates with each other
  res <- get_association_statistics(metadata)
  print(corrplot(res$estimate, p.mat = res$pval,
                 sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
                 pch.cex = 0.2, tl.cex = 0.5))

  # res$plot has nonsig correlations set to 0 for easier viewing
  print(corrplot(res$plot, p.mat = res$pval,
                 sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
                 pch.cex = 0.2, tl.cex = 0.5))

  res <- run_pca_and_plot_correlations(expr_norm, metadata, scaled = FALSE)
  print(res$pc_results)
}
