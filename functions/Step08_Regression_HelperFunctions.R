# Helper functions for regressing out technical effects in bulk data.
#
# For finding the best model, we first remove covariates that are highly
# correlated with other covariates (R^2 > 0.5) to avoid over-fitting. To
# determine which covariates to remove, I converted all factors to numeric
# values and calculated the correlations between all pairs of covariates. For
# pairs of covariates with R^2 > 0.5 (correlation ~ +/-0.7):
#   1. RIN2 was removed in every data set due to R^2 > 0.99 with RIN
#   2. Any covariate that paired with diagnosis or tissue, which have to be in
#      the model, was removed
#   3. For remaining pairs, the covariate with the higher mean R^2 with all
#      other covariates was removed.

# Runs stepwise regression, iteratively adding variables to the model until it
# can't be improved any more. This is done to create both a fixed-effects-only
# model and a mixed effect model.
#
# Arguments:
#   dataset - the name of the data set
#   bulk_se - a SummarizedExperiment object containing an assay called "counts"
#   covariates - a samples x variables data frame with covariates to consider.
#                Must contain only factors or numerical values.
#   batch_vars - a named list of batch variables, where the name is a data set
#                and the list item is the name of the column in 'covariates'
#                that corresponds to batch
Find_BestModel <- function(dataset, bulk_se, covariates, batch_vars) {
  expr_norm <- scuttle::calculateCPM(assay(bulk_se, "counts"))
  expr_norm <- scuttle::normalizeCounts(expr_norm,
                                        size.factors = bulk_se$tmm_factors,
                                        transform = "log",
                                        pseudo.count = 0.5,
                                        center.size.factors = FALSE)

  # Remove sample and tmm_factors from covariates as they should not be in the
  # model. RIN2 is also removed due to extremely high correlation with RIN.
  covariates <- data.frame(covariates) %>%
    dplyr::select(-sample, -tmm_factors, -RIN2)

  # ROSMAP needs projid removed as well as it's redundant with individualID
  if (dataset == "ROSMAP") {
    covariates <- dplyr::select(covariates, -projid)
  }

  valid_covariates <- Filter_Covariates(covariates)
  print(paste("Removing covariates:",
              paste(setdiff(colnames(covariates), valid_covariates),
                collapse = ", ")))

  # Get significant covariates -- all of the variables in Mayo, MSBB, and ROSMAP
  # end up being significantly correlated with the first couple of PCs so this
  # part is not strictly necessary. It's mostly here in case we add another data
  # set in which some variables aren't significantly correlated to PCs, or if we
  # want to add more covariates to existing data sets.
  res <- run_pca_and_plot_correlations(expr_norm,
                                       covariates[, valid_covariates],
                                       scaled = TRUE)
  sig_covars <- res$significant_covariates
  print(paste("Significant covariates:", paste(sig_covars, collapse = ", ")))

  cqn_counts <- list(E = expr_norm) # Fake voom structure needed for stepwise_regression

  rand_vars <- intersect(sig_covars, c("individualID", batch_vars[[dataset]]))

  # IndividualID and batch can't be fixed effects
  fixed_vars <- setdiff(sig_covars, rand_vars)

  # Fixed effect model: evaluate best model using "~diagnosis + tissue" as the
  # base model and iteratively adding fixed effect variables.
  res_fixed <- stepwise_regression(covariates,
                                   primary_variable = "diagnosis",
                                   cqn_counts = cqn_counts,
                                   model_variables = fixed_vars,
                                   random_effect = NULL,
                                   add_model = "tissue")

  # rand_vars <- paste0("(1|", rand_vars, ")")
  # model <- mvIC::mvForwardStepwise(exprObj = expr_norm,
  #                                 baseFormula = "~ diagnosis + tissue",
  #                                 data = covariates,
  #                                 variables = c(fixed_vars, rand_vars))

  # Mixed effect model includes random variables
  res_mixed <- stepwise_regression(covariates,
                                   primary_variable = "diagnosis",
                                   cqn_counts = cqn_counts,
                                   model_variables = sig_covars,
                                   random_effect = rand_vars,
                                   add_model = c("tissue"))

  # Save the formulas as strings instead of formula objects, so they can be
  # manipulated easier elsewhere.
  vars_fixed <- all.vars(res_fixed$formula)
  formula_fixed <- paste("~", paste(vars_fixed, collapse = " + "))

  vars_mixed <- all.vars(res_mixed$formula)
  rand_vars <- intersect(vars_mixed, c("individualID", batch_vars[[dataset]]))
  fixed_vars <- setdiff(vars_mixed, rand_vars)

  if (length(rand_vars) > 0) {
    rand_vars <- paste0("(1 | ", rand_vars, ")")
  }

  formula_mixed <- paste("~", paste(c(fixed_vars, rand_vars), collapse = " + "))

  formulas <- list("formula_fixed" = formula_fixed,
                   "formula_mixed" = formula_mixed)
  print(formulas)

  # Save in case we need to re-run the regression without having to run this
  # step again
  Save_ModelFormulas(dataset, formulas)

  return(formulas)
}


# Removes covariates with a high level of correlation with other covariates, as
# described at the top of this file.
#
# Arguments:
#   covariates - a samples x variable data frame containing covariates, which
#                may be factors or numeric values. Plain strings (non-factored)
#                are not allowed.
#
# Returns:
#   a vector of variable names that should be removed from consideration when
#   trying to create the linear model
Filter_Covariates <- function(covariates) {
  # Convert all factors to numeric values
  cov_numeric <- covariates %>%
    dplyr::mutate(across(where(is.factor), as.numeric)) %>%
    as.matrix()

  corr <- cor(cov_numeric)
  r2 <- corr^2
  diag(r2) <- 0
  mean_r2 <- rowMeans(r2)

  r2_melt <- r2
  r2_melt[upper.tri(r2_melt)] <- 0  # Avoids picking up both (a vs b) and (b vs a)
  r2_melt <- melt(r2_melt) %>%  # Data frame with columns Var1, Var2, Value (R^2)
    subset(value >= 0.5) %>%    # R^2 > 0.5 only
    dplyr::arrange(desc(value))

  to_remove <- c()
  for (R in 1:nrow(r2_melt)) {
    vars <- as.character(c(r2_melt$Var1[R], r2_melt$Var2[R]))

    if (any(vars %in% to_remove)) {
      # No need to re-check if we're already removing one of these variables
      next
    } else if (any(vars %in% c("diagnosis", "tissue"))) {
      # Remove the variable that isn't diagnosis or tissue
      to_remove <- c(to_remove,
                     vars[!(vars %in% c("diagnosis", "tissue"))])
    } else {
      # Remove the variable with the largest mean R^2 with other variables
      to_remove <- c(to_remove,
                     vars[which.max(mean_r2[vars])])
    }
  }

  to_remove <- unique(to_remove)
  return(setdiff(colnames(cov_numeric), to_remove))
}
