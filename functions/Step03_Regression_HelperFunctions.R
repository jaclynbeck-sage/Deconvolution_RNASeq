# Helper functions for regressing out technical effects in bulk data.
#
# For finding the best model, we first remove covariates that are highly
# correlated with other covariates (R^2 > 0.5) to avoid over-fitting. To
# determine which covariates to remove, I calculated the correlations between
# all pairs of numeric covariates. For pairs of covariates with R^2 > 0.5
# (correlation ~ +/-0.7), the covariate with the higher mean R^2 with all other
# remaining covariates was removed.

# Runs stepwise regression, iteratively adding variables to the model until it
# can't be improved any more. This is done to create both a fixed-effects-only
# model and a mixed effect model.
#
# Arguments:
#   dataset - the name of the data set
#   bulk_se - a SummarizedExperiment object containing an assay called "counts"
#   covariates - a samples x variables data frame with covariates to consider.
#                Must contain only factors or numerical values.
Find_BestModel <- function(dataset, tissue, bulk_se, covariates) {
  expr_norm <- sageRNAUtils::simple_log2norm(as.matrix(assay(bulk_se, "counts")),
                                             size_factors = bulk_se$tmm_factors,
                                             pseudocount = 0.5)

  # Remove variables that should not be in the model. RIN2 is also removed due
  # to extremely high correlation with RIN. We explicitly keep some biological
  # variables and batch.
  covariates <- data.frame(covariates) |>
    dplyr::select(
      where(is.numeric),
      any_of(c("diagnosis", "sex", "race", "ethnicity", "ageDeath")), # keep
      batch, # keep
      -any_of(
        c("lib_size", "tmm_factors", "flowcell", "final_batch", "RIN2",
          "sequencingBatch", "ad_cerad", "high_braak", "projid", "educ",
          "yearsEducation")
      ) # numeric variables to remove
    )

  sink(file.path(dir_tmp, str_glue("{dataset}_{tissue}_formulas.log")))

  to_remove <- Filter_Covariates(covariates)
  print(paste("Removing covariates:", paste(to_remove, collapse = ", ")))

  valid_covariates <- setdiff(colnames(covariates), to_remove)

  rand_vars <- intersect(valid_covariates, "batch")

  # batch can't be fixed effects, and diagnosis is already in the model
  fixed_vars <- setdiff(valid_covariates, c(rand_vars, "diagnosis"))

  # Fixed effect model: evaluate best model using "~ diagnosis" as the base
  # model and iteratively adding fixed effect variables.
  res_fixed <- mvIC::mvForwardStepwise(exprObj = expr_norm,
                                       baseFormula = "~ diagnosis",
                                       data = covariates,
                                       variables = fixed_vars)

  rand_vars <- paste0("(1|", rand_vars, ")")
  res_mixed <- mvIC::mvForwardStepwise(exprObj = expr_norm,
                                       baseFormula = "~ diagnosis",
                                       data = covariates,
                                       variables = c(fixed_vars, rand_vars))

  # Save the formulas as strings instead of formula objects, so they can be
  # manipulated easier elsewhere.
  vars_fixed <- all.vars(res_fixed$formula)
  formula_fixed <- paste("~", paste(vars_fixed, collapse = " + "))

  vars_mixed <- all.vars(res_mixed$formula)
  rand_vars <- intersect(vars_mixed, "batch")
  fixed_vars <- setdiff(vars_mixed, rand_vars)

  if (length(rand_vars) > 0) {
    rand_vars <- paste0("(1 | ", rand_vars, ")")
  }

  formula_mixed <- paste("~", paste(c(fixed_vars, rand_vars), collapse = " + "))

  formulas <- list("formula_fixed" = formula_fixed,
                   "formula_mixed" = formula_mixed)
  print(formulas)

  sink()

  # Save in case we need to re-run the regression without having to run this
  # step again
  Save_ModelFormulas(str_glue("{dataset}_{tissue}"), formulas)

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
  # Subset to numeric values only
  cov_numeric <- covariates |>
    select(where(is.numeric)) |>
    as.matrix()

  corr <- cor(cov_numeric)
  na_vars <- colSums(is.na(cov_numeric))

  # Remove columns with NA values
  corr <- corr[na_vars == 0, na_vars == 0]

  to_remove <- remove_correlated_variables(corr, na_vars, 0.5)
  to_remove <- c(to_remove, names(na_vars)[na_vars > 0])

  return(to_remove)
}


remove_correlated_variables <- function(cor_mat, na_vars, R2_threshold = 0.5) {
  r2 <- cor_mat^2
  diag(r2) <- NA

  r2_melt <- r2
  r2_melt[upper.tri(r2_melt, diag = TRUE)] <- 0  # Avoids picking up both (a vs b) and (b vs a)
  r2_melt <- r2_melt |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(cols = -var1,
                        names_to = "var2", values_to = "value") |>
    subset(value >= R2_threshold) |>  # pairs with R^2 > 0.5 (cor ~ 0.7) only
    dplyr::arrange(desc(value))

  to_remove <- c()

  for (R in 1:nrow(r2_melt)) {
    vars <- as.character(c(r2_melt$var1[R], r2_melt$var2[R]))

    if (any(vars %in% to_remove)) {
      # No need to re-check if we're already removing one of these variables
      next
    } else {
      # If either variable has any NA values, and they don't have the same number
      # of NAs, remove the one with the most NAs
      if (any(na_vars[vars] > 0) && (na_vars[vars[1]] != na_vars[vars[2]])) {
        to_remove <- c(to_remove,
                       vars[which.max(na_vars[vars])])
      } else {
        cur_vars <- setdiff(rownames(r2), to_remove)
        mean_r2 <- rowMeans(r2[cur_vars, cur_vars], na.rm = TRUE)

        # Remove the variable with the largest mean R^2 with the other remaining variables
        to_remove <- c(to_remove,
                       vars[which.max(mean_r2[vars])])
      }
    }
  }

  return(unique(to_remove))
}
