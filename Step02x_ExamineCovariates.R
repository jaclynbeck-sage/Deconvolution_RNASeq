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

  # Save the whole output for debugging purposes
  saveRDS(res, file.path(dir_tmp,
                         str_glue("stepwise_regression_{dataset}.rds")))

  # Save the formula as both a fixed and a mixed effect formula
  vars <- all.vars(res$formula)
  fixed <- setdiff(vars, c("individualID", batch_vars[[dataset]]))
  mixed <- paste0("(1|", setdiff(vars, fixed), ")")

  form_fixed <- paste("~", paste(fixed, collapse = " + "))
  form_mixed <- paste(form_fixed, "+", paste(mixed, collapse = " + "))

  forms <- list("formula_fixed" = form_fixed,
                "formula_mixed" = form_mixed)

  Save_ModelFormulas(dataset, forms)
}
