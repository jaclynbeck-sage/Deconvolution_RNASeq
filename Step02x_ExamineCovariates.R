# This script looks at correlations between the bulk data and its covariates.
# All or nearly all of the covariates are significantly correlated with gene
# expression (as determined by correlation with the top PCAs of log-normalized
# data), but several covariates are highly correlated with each other. So this
# script also runs a stepwise regression to determine the combination of
# covariates that best fits the data.
#
# This is not part of the main pipeline. Rather, it is used to manually examine
# what the covariates look like and determine what models to use for the
# normalization step further down the pipeline.
library(SummarizedExperiment)
library(scuttle)
library(corrplot)
library(dplyr)
library(sageseqr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

batch_vars <- list("Mayo" = "flowcell",
                   "MSBB" = "sequencingBatch",
                   "ROSMAP" = "final_batch")

for (dataset in datasets) {
  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)

  covariates <- Load_Covariates(dataset)
  covariates <- subset(covariates, specimenID %in% colnames(bulk))

  # Fix some column types
  covariates[, batch_vars[[dataset]]] <- factor(covariates[, batch_vars[[dataset]]])

  for (col in c("sex", "race", "spanish", "ethnicity", "individualID", "apoe4_allele")) {
    if (col %in% colnames(covariates)) {
      covariates[,col] <- factor(covariates[,col])
    }
  }
  covariates$age_death[covariates$age_death == "90+"] = 90
  covariates$age_death <- as.numeric(covariates$age_death)

  for (colname in colnames(covariates)) {
    if (is.numeric(covariates[,colname]) & colname != "apoe4_allele") {
      covariates[,colname] <- scale(covariates[,colname])
    }
  }

  # Remove duplicate columns that already exist in colData
  covariates <- covariates %>% select(-diagnosis, -tissue)

  # Merge covariates into the metadata
  metadata <- merge(colData(bulk), covariates, by.x = "sample",
                    by.y = "specimenID", sort = FALSE)
  rownames(metadata) <- metadata$sample

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

  sig_covars <- res$significant_covariates

  cqn_counts <- list(E = expr_norm) # Fake voom structure needed for stepwise_regression

  # evaluate best model using "~diagnosis + tissue" as the base model and
  # iteratively adding variables. I force tissue to be in the base model because
  # it is nearly 100% correlated with batch, and batch always gets added but
  # tissue gets left out for Mayo. Since edger and deseq2 don't use random
  # effects, we need to ensure that the tissue fixed effect is in the model.
  res <- stepwise_regression(metadata, primary_variable = "diagnosis",
                             cqn_counts = cqn_counts,
                             model_variables = sig_covars,
                             random_effect = c("individualID", batch_vars[[dataset]]),
                             add_model = c("tissue"))

  saveRDS(res, file.path(dir_tmp,
                         str_glue("stepwise_regression_{dataset}.rds")))
}
