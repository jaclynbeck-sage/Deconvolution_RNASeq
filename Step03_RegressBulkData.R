# This script corrects bulk data counts for batch effects while leaving
# effects due to biological confounds in. Three different methods are used to
# correct the counts:
#   edgeR's glmQLFit fits a negative binomial model with fixed effects only
#   lme4's lmer fits a linear mixed model to log2-cpm counts, using batch as a
#     random variable. In the case of Mayo where batch was not included in the
#     model, a simple linear model is run instead.
#   Dream fits a linear mixed model to voom-normalized counts, using batch as a
#     random variable. It fits a fixed linear model for Mayo.
#
# In each case, the model is fit using the set of covariates that were
# determined to be significantly correlated with expression, and which were
# found to improve the model by Bayesian Inclusion Criteria. The set of
# covariates includes both biological confounds and technical confounds. After
# the model is fit, the estimated effects from technical confounds are
# subtracted from the log2(counts) (edgeR), or effects from biological
# confounds are added back to the residuals (lme4, Dream), leaving data which
# has only the technical confounds regressed out.

# All data is converted back to linear scale and rounded to produce integer
# counts.

library(SummarizedExperiment)
library(edgeR)
library(variancePartition)
library(dplyr)
library(lme4)
library(sageRNAUtils)
library(stringr)
library(parallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step03_Regression_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

# Assumes a ratio of 8 GB RAM per core
n_cores <- max(parallel::detectCores() / 2, 1)

for (dataset in datasets) {

  # Load data and covariates ---------------------------------------------------

  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)
  covariates <- Load_Covariates(dataset)


  # Per-tissue regression ------------------------------------------------------

  for (tissue in levels(bulk$tissue)) {
    message(str_glue("Regressing {dataset} / {tissue}"))
    bulk_tissue <- bulk[, bulk$tissue == tissue]

    covar_tissue <- Clean_BulkCovariates(colData(bulk_tissue), covariates,
                                         scale_numerical = TRUE)

    expr <- DGEList(assay(bulk_tissue, "counts"), samples = colData(bulk_tissue))
    expr$samples$norm.factors <- expr$samples$tmm_factors
    expr$samples$group <- expr$samples$diagnosis

    genes_keep <- edgeR::filterByExpr(expr, group = expr$samples$group)
    bulk_tissue <- bulk_tissue[genes_keep, ]

    # I'm not sure if there is a random component to any of these algorithms,
    # but just in case, we set a seed specific to the data set / tissue
    seed <- sageRNAUtils::string_to_seed(paste(dataset, tissue, "regression"))
    set.seed(seed)

    message("Winsorizing counts...")

    # Winsorize -- Replace values above 99th percentile for each gene with the
    # 99th percentile value

    # Normalize to CPM
    winsor <- sageRNAUtils::simple_cpm(assay(bulk_tissue, "counts"))

    pcts <- apply(winsor, 1, quantile, probs = 0.99)
    pcts <- matrix(rep(pcts, ncol(winsor)), nrow(winsor))

    winsor[winsor > pcts] <- pcts[winsor > pcts]

    # Convert back to counts
    winsor <- sageRNAUtils::cpm_to_counts(winsor,
                                          colSums(assay(bulk_tissue, "counts")))
    assay(bulk_tissue, "counts") <- winsor

    # Re-calculate norm factors
    bulk_tissue$tmm_factors <- edgeR::normLibSizes(winsor)

    rm(pcts, winsor, expr)


    # Load or calculate formulas for fixed/mixed models ------------------------

    #formulas <- Load_ModelFormulas(str_glue("{dataset}_{tissue}"))
    formulas <- NULL

    # Only run the stepwise regression if the formulas don't already exist. This
    # part takes a long time.
    if (is.null(formulas)) {
      message("Determining best models...")
      formulas <- Find_BestModel(dataset, tissue, bulk_tissue,
                                 covar_tissue, plot_var_explained = TRUE)
    }

    formula_fixed <- as.formula(formulas$formula_fixed)
    formula_mixed <- as.formula(formulas$formula_mixed)

    # edgeR glmQLfit -----------------------------------------------------------

    message("Performing edgeR glmQLfit...")
    expr <- DGEList(assay(bulk_tissue, "counts"), samples = colData(bulk_tissue))
    expr$samples$norm.factors <- expr$samples$tmm_factors

    mod <- model.matrix(formula_fixed, data = covar_tissue)

    expr <- estimateDisp(expr, mod)
    fit_edger <- glmQLFit(expr, mod)

    # edgeR uses unshrunk coefficients to calculate fitted values, so we use them too
    coefs_edger <- fit_edger$unshrunk.coefficients
    coefs_bio <- grepl("Intercept|diagnosis|sex|race|ageDeath|apoeGenotype",
                       colnames(coefs_edger))
    coefs_tech <- !coefs_bio

    # Since we don't have residuals from edgeR, we are subtracting technical
    # confounds from the data rather than adding the intercept back to the
    # residuals.
    adjust <- tcrossprod(coefs_edger[, coefs_tech], mod[, coefs_tech])

    # Note: There are multiple ways this value could be corrected. For a linear
    # model, it would just be exp(log(expr$counts) - adjust). For a GLM, it could
    # also be corrected on a linear scale by adjusting the pearson residuals:
    #   pearson = (expr$counts - mu_orig)/sqrt(variance_orig)
    #   y = mu_new + pearson * sqrt(variance_new)
    # where 'mu_orig' is the original fitted values (linear scale) and 'mu_new' is
    # exp(tcrossprod(coefs[,coefs_bio], mod[,coefs_bio]) + offset). Variance
    # is (mu + dispersion * mu^2) for edgeR. Since I'm not sure whether correcting
    # on the linear or log scale is "better" and the two methods produce very
    # similar results, I use log scale, because it's more intuitive and makes it
    # similar to how I correct with Dream and lme4.
    # [expr$counts / exp(adjust)] is functionally equivalent to
    # exp(log(expr$counts) - adjust) but avoids using pseudocount
    corrected_edger <- expr$counts / exp(adjust)

    # There shouldn't be any negative numbers but just in case
    corrected_edger[corrected_edger < 0] <- 0

    rm(expr, fit_edger, adjust)


    # Basic LM or LME fit ------------------------------------------------------

    message("Performing LME fit...")
    expr_norm <- sageRNAUtils::simple_log2norm(
      as.matrix(assay(bulk_tissue, "counts")),
      size_factors = bulk_tissue$tmm_factors,
      pseudocount = 0.5
    )

    # For Mayo TCX: No mixed effects in model, run a linear fixed-effects model
    if (formulas$formula_mixed == formulas$formula_fixed) {
      formula_lm <- paste("t(expr_norm)", formulas$formula_mixed)
      fits_lme <- lm(formula_lm, data = covar_tissue)

      coefs_lme <- t(coef(fits_lme))
      mod_lme <- model.matrix(formula_mixed, data = covar_tissue)

      resid_lme <- t(residuals(fits_lme))
    } else {
      fits_lme <- fitVarPartModel(expr_norm, formula_mixed, covar_tissue,
                                  BPPARAM = MulticoreParam(n_cores))

      coefs_lme <- t(sapply(fits_lme, fixef))
      coefs_lme <- coefs_lme[rownames(expr_norm), ] # Ensure rows are in the original order

      # Fixed-effects model matrix should be the same for all fits so just grab
      # the first one. "X" corresponds to fixed effects in this object.
      mod_lme <- getME(fits_lme[[1]], name = "X")

      resid_lme <- t(sapply(fits_lme, residuals))
      resid_lme <- resid_lme[rownames(expr_norm), ] # Ensure rows are in the original order
    }

    # Add biological covariates back to residuals
    coefs_bio <- grepl("Intercept|diagnosis|sex|race|ageDeath|apoeGenotype",
                       colnames(coefs_lme))
    adjust <- tcrossprod(coefs_lme[, coefs_bio], mod_lme[, coefs_bio])

    # Convert back to counts
    corrected_lme <- sageRNAUtils::log2_cpm_to_counts(
      data = resid_lme + adjust,
      library_size = colSums(assay(bulk_tissue, "counts")),
      size_factors = bulk_tissue$tmm_factors,
      pseudocount = 0.5
    )
    corrected_lme[corrected_lme < 0] <- 0

    rm(adjust, resid_lme, expr_norm, fits_lme)


    # ComBat batch adjustment --------------------------------------------------

    message("Performing ComBat batch adjustment...")

    # Biological covariates that might be in the formula. Diagnosis should
    # always be there.
    c_vars <- c("diagnosis",
                intersect(all.vars(formula_mixed),
                          c("sex", "ageDeath", "race", "apoeGenotype")))

    c_form <- paste("~", paste(c_vars, collapse = " + "))
    c_mod <- model.matrix(as.formula(c_form), data = covar_tissue)

    corrected_combat <- sva::ComBat_seq(as.matrix(assay(bulk_tissue, "counts")),
                                        batch = covar_tissue$batch,
                                        group = NULL,
                                        covar_mod = c_mod,
                                        full_mod = TRUE)


    # Create final bulk data object --------------------------------------------

    message("Saving final results...")

    assay(bulk_tissue, "corrected_edger") <- round(corrected_edger)
    assay(bulk_tissue, "corrected_lme") <- round(corrected_lme)
    assay(bulk_tissue, "corrected_combat") <- round(corrected_combat)

    bulk_tissue$tmm_factors_edger <- normLibSizes(round(corrected_edger))
    bulk_tissue$tmm_factors_lme <- normLibSizes(round(corrected_lme))
    bulk_tissue$tmm_factors_dream <- normLibSizes(round(corrected_combat))

    Save_BulkData(str_glue("{dataset}_{tissue}"), bulk_tissue)

    # Plot the effects of the corrections
    expr_norm <- simple_log2norm(assay(bulk_tissue, "counts"))
    var_genes <- rowVars(as.matrix(expr_norm)) |> sort(decreasing = TRUE) |> names()
    var_genes <- var_genes[1:5000]

    plts <- lapply(
      list("counts", "corrected_edger", "corrected_lme", "corrected_combat"),
      function(data) {
        expr_norm <- simple_log2norm(assay(bulk_tissue, data))
        pc <- prcomp(t(expr_norm[var_genes, ]))$x[, 1:2] |>
          merge(covar_tissue, by = "row.names")

        ggplot(pc, aes(x = PC1, y = PC2, color = batch)) +
          geom_point() +
          theme_bw() +
          ggtitle(data)
      }
    )

    print(patchwork::wrap_plots(plts))

    rm(bulk_tissue, corrected_edger, corrected_lme, corrected_combat)
    gc()
  }
}
