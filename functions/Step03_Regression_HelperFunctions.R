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
Find_BestModel <- function(dataset, tissue, bulk_se, covariates, n_cores, plot_var_explained = FALSE) {
  expr_norm <- sageRNAUtils::simple_log2norm(
    as.matrix(assay(bulk_se, "counts")),
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
  cat(paste(tissue, "- recommended variables:",
            paste(vars_keep, collapse = ", ")), "\n")
  cat("\n", "\n")

  rand_vars <- intersect(vars_keep, "batch")

  # batch can't be in fixed effects, and diagnosis is already in the model
  fixed_vars <- setdiff(vars_keep, c(rand_vars, "diagnosis"))

  # Fixed effect model: evaluate best model using "~ diagnosis" as the base
  # model and iteratively adding fixed effect variables.
  res_fixed <- mvIC::mvForwardStepwise(exprObj = expr_norm,
                                       baseFormula = "~ diagnosis",
                                       data = covariates_clean,
                                       variables = c(fixed_vars, rand_vars))

  rand_vars <- paste0("(1|", rand_vars, ")")
  res_mixed <- mvIC::mvForwardStepwise(exprObj = expr_norm,
                                       baseFormula = "~ diagnosis",
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
  cat("\n")

  sink()

  # Save in case we need to re-run the regression without having to run this
  # step again
  Save_ModelFormulas(str_glue("{dataset}_{tissue}"), formulas)

  return(formulas)
}


# Clean_BulkCovariates - takes a covariates dataframe for one of the bulk
# datasets, makes sure that categorical variables are factors, scales numerical
# variables, and merges the cleaned covariates dataframe with the metadata
# dataframe.
#
# Arguments:
#   metadata - the metadata dataframe (colData()) from a SummarizedExperiment
#   covariates - a dataframe of covariates, where rows are samples and columns
#                are the covariates
#   scale_numerical - TRUE or FALSE, whether to scale numeric columns
#
# Returns:
#   a dataframe with merged metadata and cleaned covariates
Clean_BulkCovariates <- function(metadata, covariates, scale_numerical = TRUE) {
  covariates <- subset(covariates, specimenID %in% rownames(metadata))

  # Only keep these variables of interest and drop all other categorical
  # variables. "diagnosis" and "tissue" are excluded from this list because they
  # already exist in `metadata` and don't need to be copied over from
  # `covariates`
  cat_cols_keep <- c("specimenID", "sex", "ageDeath", "batch")

  num_cols_keep <- c("RIN", "PMI", "picard_PERCENT_DUPLICATION",
                     "rsem_uniquely_aligned_percent")

  # Bin ages in 5-year intervals
  covariates <- covariates |>
    mutate(
      ageDeath_num = suppressWarnings(as.numeric(ageDeath)),
      ageDeath = case_when(
        !is.na(ageDeath_num) ~ cut(
          ageDeath_num,
          breaks = c(0, 65, 70, 75, 80, 85, 90),
          labels = c("Under 65", "65 - 69", "70 - 74", "75 - 79", "80 - 84", "85 - 89"),
          right = FALSE,
          include.lowest = TRUE
        ),
        ageDeath == "90_or_over" ~ "90+", # Fix for Mayo
        .default = ageDeath # 90+ or NA
      )
    ) |>
    select(any_of(cat_cols_keep), any_of(num_cols_keep))

  # Re-order the age factors so "Under 65" is first
  covariates$ageDeath <- factor(covariates$ageDeath,
                                levels = c("Under 65", "65 - 69", "70 - 74",
                                           "75 - 79", "80 - 84", "85 - 89", "90+"))

  for (col in cat_cols_keep) {
    if (col %in% colnames(covariates)) {
      covariates[, col] <- factor(as.character(covariates[, col]))
    }
  }

  # Merge covariates into the metadata
  sample_order <- rownames(metadata)
  metadata <- merge(metadata, covariates,
                    by.x = "sample", by.y = "specimenID",
                    sort = FALSE)
  rownames(metadata) <- metadata$sample

  # Put the data frame back in the original order, as merge might change it
  metadata <- data.frame(metadata[sample_order, ])

  # Scale numerical covariates if scale == TRUE
  if (scale_numerical) {
    for (colname in colnames(metadata)) {
      if (is.numeric(metadata[, colname])) {
        metadata[, colname] <- as.numeric(scale(metadata[, colname]))
      }
    }
  }

  return(metadata)
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

  cat("Correlation:", "\n")
  print(cc)
  cat("\n", "\n")

  # Use most variable genes only for extractVarPart
  var_genes <- rowVars(expr_norm) |> sort(decreasing = TRUE) |> names()
  var_genes <- var_genes[1:5000]

  n_pairs <- sum(cc > 0.7, na.rm = TRUE)

  to_remove <- c()
  for (N in 1:n_pairs) {
    vars_remaining <- setdiff(colnames(covariates_clean), to_remove)
    form_iter <- paste("~", paste0(vars_remaining, collapse = " + "))

    if (!any(cc[vars_remaining, vars_remaining] > 0.7, na.rm = TRUE)) {
      break # We've removed all highly-correlated variables that can be removed, stop looping
    }

    varpart <- fitExtractVarPartModel(
      expr_norm[var_genes, ], form_iter, covariates_clean,
      BPPARAM = MulticoreParam(n_cores)
    ) |>
      sortCols()

    if (plot_var_explained) {
      print(plotVarPart(varpart, main = str_glue("{plot_title}: iteration {N}")))
    }

    cat("Mean variance:", "\n")
    print(colMeans(varpart))
    cat("\n", "\n")

    vp <- vp_summary(varpart)
    cor_pairs <- cc[vars_remaining, vars_remaining] |>
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
      ungroup() |>
      select(-lowest_pct)

    cat(str_glue("Iteration {N} pairs:"), "\n")
    print(as.data.frame(cor_pairs))
    cat("\n", "\n")

    pair <- c(cor_pairs$var1[1], cor_pairs$var2[1])

    # Only remove one of the two variables if neither has been removed yet
    if (!any(pair %in% to_remove)) {
      to_remove <- c(to_remove, cor_pairs$var_remove[1])
    }
  }

  cat(paste("Removing", length(to_remove), "variables:",
            paste(to_remove, collapse = ", ")), "\n")
  cat("\n", "\n")

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
