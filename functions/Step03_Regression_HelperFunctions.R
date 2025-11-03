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
Find_BestModel <- function(dataset, tissue, bulk_se, covariates, plot_var_explained = FALSE) {
  expr_norm <- sageRNAUtils::simple_log2norm(
    as.matrix(assay(bulk_se, "counts")),
    size_factors = bulk_se$tmm_factors,
    pseudocount = 0.5
  )

  sink(file.path(dir_tmp, str_glue("{dataset}_{tissue}_formulas.log")))

  # Remove variables that should not be in the model, and variables that have
  # NA values or are all the same value.
  covariates_clean <- data.frame(covariates) |>
    sageRNAUtils::remove_unusable_covariates() |>
    dplyr::select(
      -any_of(c("lib_size", "tmm_factors")),
      -where(~any(is.na(.x))) # Remove any columns with at least 1 NA value
    )

  plot_title <- paste(dataset, tissue)
  vars_keep <- remove_correlated_vars(covariates_clean, expr_norm, n_cores,
                                      plot_var_explained, plot_title)
  print(paste(tissue, "- recommended variables:",
              paste(vars_keep, collapse = ", ")))

  rand_vars <- intersect(vars_keep, "batch")

  # batch can't be in fixed effects, and diagnosis is already in the model
  fixed_vars <- setdiff(vars_keep, c(rand_vars, "diagnosis"))

  base_form_fixed <- "~ diagnosis"
  base_form_mixed <- base_form_fixed

  # Special case: MSBB IFG needs "batch" added to the fixed formula in order for
  # edgeR output to have minimal batch effects, but mvIC doesn't add it.
  if (tissue == "IFG") {
    base_form_fixed <- paste(base_form_fixed, "+ batch")
  }

  # Fixed effect model: evaluate best model using "~ diagnosis" as the base
  # model and iteratively adding fixed effect variables.
  res_fixed <- mvIC::mvForwardStepwise(exprObj = expr_norm,
                                       baseFormula = base_form_fixed,
                                       data = covariates_clean,
                                       variables = c(fixed_vars, rand_vars))

  rand_vars <- paste0("(1|", rand_vars, ")")
  res_mixed <- mvIC::mvForwardStepwise(exprObj = expr_norm,
                                       baseFormula = base_form_mixed,
                                       data = covariates_clean,
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


vp_summary <- function(var_part) {
  var_part |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    select(-Residuals) |>
    tidyr::pivot_longer(-gene, names_to = "covariate", values_to = "value") |>
    group_by(covariate) |>
    summarize(mean_pct = mean(value),
              median_pct = median(value)) |>
    mutate(rank_mean = rank(-mean_pct),
           rank_median = rank(-median_pct),
           total_rank = (rank_mean + rank_median) / 2) |>
    arrange(desc(mean_pct))
}


remove_correlated_vars <- function(covariates_clean, expr_norm, n_cores,
                                   plot_var_explained, plot_title) {
  form_cc <- paste("~", paste0(colnames(covariates_clean), collapse = " + "))
  cc <- canCorPairs(form_cc, covariates_clean, showWarnings = FALSE)
  cc[upper.tri(cc, diag = TRUE)] <- NA

  print(cc)

  # Use most variable genes only for extractVarPart
  var_genes <- rowVars(expr_norm) |> sort(decreasing = TRUE) |> names()
  var_genes <- var_genes[1:5000]

  varpart <- fitExtractVarPartModel(
    expr_norm[var_genes, ], form_cc, covariates_clean,
    BPPARAM = MulticoreParam(n_cores)
  ) |>
    sortCols()

  if (plot_var_explained) {
    print(plotVarPart(varpart, main = plot_title))
  }

  print(colMeans(varpart))

  vp <- vp_summary(varpart)
  cor_pairs <- cc |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "cor_value") |>
    subset(cor_value > 0.7) |>
    arrange(desc(cor_value)) |>
    rowwise() |>
    mutate(
      var1_rank = vp$rank_mean[vp$covariate == var1],
      var2_rank = vp$rank_mean[vp$covariate == var2],
      lowest_pct = which.max(c(var1_rank, var2_rank)),
      var_remove = c(var1, var2)[lowest_pct]
    ) |>
    ungroup()

  print(cor_pairs)

  to_remove <- c()
  for (N in 1:nrow(cor_pairs)) {
    pair <- c(cor_pairs$var1[N], cor_pairs$var2[N])

    # Only remove one of the two variables if neither has been removed yet
    if (!any(pair %in% to_remove)) {
      to_remove <- c(to_remove, cor_pairs$var_remove[N])
    }
  }

  print(paste("Removing", length(to_remove), "variables:",
              paste(to_remove, collapse = ", ")))

  vars_keep <- vp |>
    arrange(desc(mean_pct)) |>
    subset(!(covariate %in% to_remove)) |> # mean_pct > 0.01 &
    pull(covariate)

  if (plot_var_explained & !all(colnames(covariates_clean) %in% vars_keep)) {
    form_rec <- paste("~", paste(vars_keep, collapse = " + "))
    varpart2 <- fitExtractVarPartModel(
      expr_norm[var_genes, ], form_rec, covariates_clean,
      BPPARAM = MulticoreParam(n_cores)
    ) |>
      sortCols()

    print(plotVarPart(varpart2, main = paste(plot_title, "(Recommended variables)")))
  }

  return(vars_keep)
}
