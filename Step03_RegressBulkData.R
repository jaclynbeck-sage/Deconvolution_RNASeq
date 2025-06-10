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

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step03_Regression_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

for (dataset in datasets) {

  # Load data and covariates ---------------------------------------------------

  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)
  covariates <- Load_Covariates(dataset) |>
    # Remove some very skewed PCT metrics which result in over-fitting
    select(-starts_with("Alignment"), -RnaSeqMetrics__PCT_CORRECT_STRAND_READS,
           -RnaSeqMetrics__PCT_RIBOSOMAL_BASES, -RnaSeqMetrics__PCT_INTERGENIC_BASES)


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


    # Load or calculate formulas for fixed/mixed models ------------------------

    formulas <- Load_ModelFormulas(str_glue("{dataset}_{tissue}"))

    # Only run the stepwise regression if the formulas don't already exist. This
    # part takes a long time.
    if (is.null(formulas)) {
      message("Determining best models...")
      formulas <- Find_BestModel(dataset, tissue, bulk_tissue,
                                 covar_tissue)
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
    coefs_bio <- grepl("Intercept|diagnosis|sex|race|ethnicity|ageDeath",
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
    corrected_log <- log(expr$counts + 1e-8) - adjust # Use a very small pseudocount
    corrected_edger <- exp(corrected_log) - 1e-8 # subtract the pseudo-count we just added
    corrected_edger[corrected_edger < 0] <- 0


    # Basic LM or LME fit ------------------------------------------------------

    message("Performing LME fit...")
    expr_norm <- sageRNAUtils::simple_log2norm(as.matrix(assay(bulk_tissue, "counts")),
                                               size_factors = bulk_tissue$tmm_factors,
                                               pseudocount = 0.5)

    # For Mayo: No mixed effects in model, run a linear fixed-effects model
    if (formulas$formula_mixed == formulas$formula_fixed) {
      formula_lm <- paste("t(expr_norm)", formulas$formula_mixed)
      fits_lme <- lm(formula_lm, data = covar_tissue)

      coefs_lme <- t(coef(fits_lme))
      mod_lme <- model.matrix(formula_mixed, data = covar_tissue)

      resid_lme <- t(residuals(fits_lme))
    } else {
      formula_lme <- paste("row", formulas$formula_mixed)

      # LME4 can only do one gene at a time.
      # Hack to get row names back after mclapply, since I'm not sure the rows
      # are always returned in the same order
      genes <- rownames(expr_norm)
      names(genes) <- genes

      n_cores <- max(parallel::detectCores()/2, 1)

      fits_lme <- parallel::mclapply(genes, function(gene) {
        row <- expr_norm[gene, ]
        fit <- lme4::lmer(as.formula(formula_lme), data = covar_tissue)
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

    # Add biological covariates back to residuals
    coefs_bio <- grepl("Intercept|diagnosis|sex|race|ethnicity|ageDeath",
                       colnames(coefs_lme))
    adjust <- tcrossprod(coefs_lme[, coefs_bio], mod_lme[, coefs_bio])

    # Convert back to counts
    corrected_lme <- sageRNAUtils::log2_cpm_to_counts(
      data = resid_lme + adjust,
      library_size = colSums(assay(bulk_tissue, "counts")) * bulk_tissue$tmm_factors,
      pseudocount = 0.5
    )
    corrected_lme[corrected_lme < 0] <- 0


    # Dream fit ----------------------------------------------------------------

    message("Performing dream fit...")
    gc()
    # This is what sageseqr uses to fit a mixed effect model. The function
    # in sageseqr to calculate the fit + residuals requires CQN counts, which
    # we don't use, and requires jumping through extra hoops to get other args and
    # output into the right format, so I just call voomWithDreamWeights() and
    # dream() myself here.
    expr <- DGEList(assay(bulk_tissue, "counts"))
    expr$samples$norm.factors <- bulk_tissue$tmm_factors

    n_cores <- max(parallel::detectCores()/4, 1)
    voom_counts <- voomWithDreamWeights(counts = expr,
                                        formula = formula_mixed,
                                        data = covar_tissue,
                                        BPPARAM = BiocParallel::SnowParam(n_cores))

    fit_dream <- dream(exprObj = voom_counts,
                       formula = formula_mixed,
                       data = covar_tissue,
                       computeResiduals = TRUE,
                       BPPARAM = BiocParallel::SnowParam(n_cores))

    resid_dream <- fit_dream$residuals
    coefs_dream <- coef(fit_dream)
    mod_dream <- fit_dream$design

    # Add biological covariates back to residuals
    coefs_bio <- grepl("Intercept|diagnosis|sex|race|ethnicity|ageDeath",
                       colnames(coefs_lme))
    adjust <- tcrossprod(coefs_dream[, coefs_bio], mod_dream[, coefs_bio])

    if (nrow(resid_dream) != nrow(bulk_tissue)) {
      failed_genes <- setdiff(rownames(bulk_tissue), rownames(resid_dream))
      print(paste0("WARNING: Dream failed on genes [",
                   paste(failed_genes, collapse = ", "),
                   "]. They will be removed from the final data set."))
    }

    # voomWithDreamWeights uses edgeR to normalize counts so we reverse the
    # operation
    corrected_dream <- sageRNAUtils::edger_log2_cpm_to_counts(
      data = resid_dream + adjust,
      library_size = expr$samples$lib.size * expr$samples$norm.factors,
      prior_count = 0.5 # this is what voomWithDreamWeights uses
    )
    corrected_dream[corrected_dream < 0] <- 0


    # Save results/fits for further examination --------------------------------

    message("Saving final results...")
    #saveRDS(list(fit_edger = fit_edger,
    #             corrected_edger = corrected_edger,
    #             fits_lme = fits_lme,
    #             corrected_lme = corrected_lme,
    #             fit_dream = fit_dream,
    #             corrected_dream = corrected_dream),
    #        file.path(dir_tmp, str_glue("{dataset}_{tissue}_models.rds")))


    # Create final bulk data object --------------------------------------------

    # Dream can fail on specific genes and not return them in the matrix, so we
    # need to cut any failed genes out of the other matrices too
    genes_keep <- rownames(corrected_dream)
    corrected_edger <- corrected_edger[genes_keep, ]
    corrected_lme <- corrected_lme[genes_keep, ]

    bulk_tissue <- bulk_tissue[genes_keep, ] # Adjusts the counts array if necessary

    assay(bulk_tissue, "corrected_edger") <- round(corrected_edger)
    assay(bulk_tissue, "corrected_lme") <- round(corrected_lme)
    assay(bulk_tissue, "corrected_dream") <- round(corrected_dream)

    bulk_tissue$tmm_factors_edger <- normLibSizes(round(corrected_edger))
    bulk_tissue$tmm_factors_lme <- normLibSizes(round(corrected_lme))
    bulk_tissue$tmm_factors_dream <- normLibSizes(round(corrected_dream))

    Save_BulkData(str_glue("{dataset}_{tissue}"), bulk_tissue)

    rm(bulk_tissue, corrected_edger, corrected_lme, corrected_dream,
       fit_edger, fits_lme, fit_dream)
    gc()
  }
}
