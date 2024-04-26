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
#
# The intercept coefficient changes depending on the scale of the other
# covariates, so to avoid having to figure out which covariates should be scaled
# and which shouldn't, we remove the intercept as a techincal covariate, which
# leaves the data approximately centered, and add back the mean log-counts (or
# log-cpm) for each gene. The mean values are either already in cpm with no
# adjustment needed per sample (lme4, Dream), or adjusted for library size
# (edgeR) prior to addition. Adding back the mean is necessary for comparing
# relative gene levels against the single cell reference or signature, as those
# values are un-centered. This produces better results from the algorithms than
# leaving the original intercept in and using all scaled or all unscaled
# covariates.

# All data is converted back to linear scale and rounded to produce integer
# counts.

library(SummarizedExperiment)
library(edgeR)
library(variancePartition)
library(scuttle)
library(dplyr)
library(lme4)
library(sageseqr)
library(stringr)
library(parallel)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step08_Regression_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

batch_vars <- list("Mayo" = "flowcell",
                   "MSBB" = "sequencingBatch",
                   "ROSMAP" = "final_batch")

for (dataset in datasets) {

  # Load data and covariates ---------------------------------------------------

  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)
  covariates <- Load_Covariates(dataset)
  covariates <- Clean_BulkCovariates(dataset, colData(bulk), covariates,
                                     scale_numerical = TRUE)

  # I'm not sure if there is a random component to any of these algorithms, but
  # just in case, we set a seed specific to the data set
  numerical_name <- as.numeric(charToRaw(dataset))
  set.seed(sum(numerical_name))


  # Load or calculate formulas for fixed/mixed models --------------------------

  formulas <- Load_ModelFormulas(dataset)

  # Only run the stepwise regression if the formulas don't already exist. This
  # part takes a long time.
  if (is.null(formulas)) {
    formulas <- Find_BestModel(dataset, bulk, covariates, batch_vars)
  }

  formula_fixed <- as.formula(formulas$formula_fixed)
  formula_mixed <- as.formula(formulas$formula_mixed)


  # edgeR glmQLfit -------------------------------------------------------------

  expr <- DGEList(assay(bulk, "counts"))
  expr$samples$norm.factors <- bulk$tmm_factors # these were pre-calculated

  mod <- model.matrix(formula_fixed, data = covariates)

  expr <- estimateDisp(expr, mod)
  fit_edger <- glmQLFit(expr, mod)

  # edgeR uses unshrunk coefficients to calculate fitted values, so use them too
  coefs_edger <- fit_edger$unshrunk.coefficients
  coefs_bio <- grepl("diagnosis|tissue|sex", colnames(coefs_edger))
  coefs_tech <- !coefs_bio

  # Since we don't have residuals from edgeR, we are subtracting technical
  # confounds from the data rather than adding biological confounds to the
  # residuals. The intercept is included in the technical confounds.
  adjust <- tcrossprod(coefs_edger[, coefs_tech], mod[, coefs_tech])

  # Edger uses an offset in its GLM (log(lib.size * norm.factors)), so the
  # formula is actually log(expr) ~ Intercept + covariates + offset. To replace
  # the intercept with mean log-counts, we have to subtract the offset from the
  # log-counts to get values on the same scale. i.e. if we had no covariates, we
  # can assume that mean(log(expr)) ~ Intercept + offset, so to subtract the
  # intercept in this equation we need to add (mean(log(expr)) - offset) in its
  # place. A pseudocount is added to avoid taking the log of 0. EdgeR values
  # are on the natural log scale.
  new_intercept <- rowMeans(log(expr$counts + 0.5) - fit_edger$offset)

  # Turn the vector into a matrix the same size as expr, where each column is
  # new_intercept, so it can be added to expr.
  new_intercept <- matrix(rep(new_intercept, ncol(expr)),
                          nrow = length(new_intercept))

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
  corrected_log <- log(expr$counts + 0.5) - adjust + new_intercept
  corrected_edger <- exp(corrected_log) - 0.5 # subtract the pseudo-count we just added
  corrected_edger[corrected_edger < 0] <- 0


  # Basic LM or LME fit --------------------------------------------------------

  expr_norm <- scuttle::calculateCPM(assay(bulk, "counts"))
  expr_norm <- scuttle::normalizeCounts(expr_norm,
                                        size.factors = bulk$tmm_factors,
                                        transform = "log",
                                        pseudo.count = 0.5,
                                        center.size.factors = FALSE)

  # For Mayo: No mixed effects in model, run a linear fixed-effects model
  if (formulas$formula_mixed == formulas$formula_fixed) {
    formula_lm <- paste("t(expr_norm)", formulas$formula_mixed)
    fits_lme <- lm(formula_lm, data = covariates)

    coefs_lme <- t(coef(fits_lme))
    mod_lme <- model.matrix(formula_mixed, data = covariates)

    resid_lme <- residuals(fits_lme)
  } else {
    formula_lme <- paste("row", formulas$formula_mixed)

    # LME4 can only do one gene at a time.
    # Hack to get row names back after mclapply, since I'm not sure the rows are
    # always returned in the same order
    genes <- rownames(expr_norm)
    names(genes) <- genes

    n_cores <- parallel::detectCores()-1

    fits_lme <- parallel::mclapply(genes, function(gene) {
      row <- expr_norm[gene, ]
      fit <- lme4::lmer(as.formula(formula_lme), data = covariates)
      return(fit)
    }, mc.cores = n_cores)

    coefs_lme <- t(sapply(fits_lme, fixef))
    coefs_lme <- coefs_lme[rownames(expr_norm), ] # Ensure rows are in the original order

    # Fixed-effects model matrix should be the same for all fits so just grab
    # the first one. "X" corresponds to fixed effects in this object.
    mod_lme <- getME(fits_lme[[1]], name = "X")

    resid_lme <- t(sapply(fits_lme, residuals))
    resid_lme <- resid_lme[rownames(expr_norm), ] # Ensure rows are in the original order
  }

  coefs_bio <- grepl("diagnosis|tissue|sex", colnames(coefs_lme))

  # Add biological covariates back to residuals -- excludes intercept
  adjust <- tcrossprod(coefs_lme[, coefs_bio], mod_lme[, coefs_bio])

  # No library size adjustment needed here, expr_norm is log2-cpm and there is
  # no offset in the formula
  new_intercept <- rowMeans(expr_norm)
  new_intercept <- matrix(rep(new_intercept, ncol(expr_norm)),
                          nrow = length(new_intercept))

  corrected_lme <- 2^(resid_lme + adjust + new_intercept) - 0.5 # subtract pseudocount

  # corrected_lme is in CPM. We need to reverse the CPM operation by multipling
  # by (library size / 1e6) to get counts again, in order to avoid rounding a
  # lot of CPM values < 0.5 to 0.
  size_adjust <- (colSums(assay(bulk, "counts")) * bulk$tmm_factors) / 1e6

  # Each column gets multiplied by its size_adjust entry
  corrected_lme <- sweep(corrected_lme, 2, size_adjust, "*")
  corrected_lme[corrected_lme < 0] <- 0


  # Dream fit ------------------------------------------------------------------

  # This is what sageseqr uses to fit a mixed effect model. The function
  # in sageseqr to calculate the fit + residuals requires CQN counts, which
  # we don't use, and requires jumping through extra hoops to get other args and
  # output into the right format, so I just call voomWithDreamWeights() and
  # dream() myself here.
  expr <- DGEList(assay(bulk, "counts"))
  expr$samples$norm.factors <- bulk$tmm_factors

  n_cores <- parallel::detectCores()-1
  voom_counts <- voomWithDreamWeights(counts = expr,
                                      formula = formula_mixed,
                                      data = covariates,
                                      BPPARAM = BiocParallel::SnowParam(n_cores))

  fit_dream <- dream(exprObj = voom_counts,
                     formula = formula_mixed,
                     data = covariates,
                     computeResiduals = TRUE,
                     BPPARAM = BiocParallel::SnowParam(n_cores))

  resid_dream <- fit_dream$residuals
  coefs_dream <- coef(fit_dream)
  mod_dream <- fit_dream$design

  # Add biological confounds back to the residuals -- excludes intercept.
  coefs_bio <- grepl("diagnosis|tissue|sex", colnames(coefs_dream))
  adjust <- tcrossprod(coefs_dream[, coefs_bio], mod_dream[, coefs_bio])

  # There is no offset in the formula and the voom counts are on the log2-cpm
  # scale, so we don't need to adjust for library size before taking the mean.
  new_intercept <- rowMeans(voom_counts$E)
  new_intercept <- matrix(rep(new_intercept, ncol(voom_counts$E)),
                          nrow = length(new_intercept))

  # Residuals are on the log2-cpm scale, so after converting to linear scale
  # we need to reverse the CPM operation to get counts.
  corrected_dream <- 2^(resid_dream + adjust + new_intercept)

  # Voom uses edgeR's cpm() function. This reverses how edgeR calculates cpm,
  # which involves scaling the prior count by mean library size.
  adj_lib <- expr$samples$lib.size * expr$samples$norm.factors
  adj_lib <- adj_lib / mean(adj_lib)
  prior <- 0.5 * adj_lib

  #voomWithDreamWeights uses lib_size + 1 so we add the 1 here too.
  size_adjust <- (expr$samples$lib.size * expr$samples$norm.factors + 1) / 1e6

  # Each column gets multiplied by its size_adjust entry and subtracts its
  # prior. EdgeR adds the prior before calculating CPM, so we have to subtract
  # the prior after reversing the CPM operation.
  corrected_dream <- sweep(corrected_dream, 2, size_adjust, "*")
  corrected_dream <- sweep(corrected_dream, 2, prior, "-")
  corrected_dream[corrected_dream < 0] <- 0


  # Save results/fits for further examination ----------------------------------

  saveRDS(list(fit_edger = fit_edger,
               corrected_edger = corrected_edger,
               fits_lme = fits_lme,
               corrected_lme = corrected_lme,
               fit_dream = fit_dream,
               corrected_dream = corrected_dream),
          file.path(dir_tmp, str_glue("{dataset}_models.rds")))


  # Create final bulk data object ----------------------------------------------

  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)

  assay(bulk, "corrected_edger") <- round(corrected_edger)
  assay(bulk, "corrected_lme") <- round(corrected_lme)
  assay(bulk, "corrected_dream") <- round(corrected_dream)

  bulk$tmm_factors_edger <- calcNormFactors(round(corrected_edger))
  bulk$tmm_factors_lme <- calcNormFactors(round(corrected_lme))
  bulk$tmm_factors_dream <- calcNormFactors(round(corrected_dream))

  Save_BulkData(dataset, bulk)
}
