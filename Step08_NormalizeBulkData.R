# This script corrects bulk data counts for batch effects while leaving
# effects due to biological confounds in. Three different methods are used to
# correct the counts:
#   edgeR's glmQLFit fits a negative binomial model with fixed effects only
#   DESeq2's DESeq also fits a negative binomial model with fixed effects only,
#     but uses different size factor and dispersion estimation than edgeR
#   sageseqr / dream fits a linear mixed model to voom-normalized counts, using
#     batch as a random variable
#
# In each case, the model is fit using the set of covariates that were
# determined to be significantly correlated with expression, and which were
# found to improve the model by Bayesian Inclusion Criteria (performed in
# Step02x). The set of covariates includes both biological confounds and
# technical confounds. After the model is fit, the estimated effects from
# biological confounds are added back to the residuals, leaving data which has
# only the technical confounds regressed out. All data is converted back to
# linear scale and rounded to produce integer counts.

library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(variancePartition)
library(scuttle)
library(dplyr)
library(lme4)

source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("Mayo", "MSBB", "ROSMAP")

designs <- list("Mayo" = paste0("~ 1 + diagnosis + tissue + percent_mito + ",
                                "RnaSeqMetrics_PCT_INTRONIC_BASES + ",
                                "RnaSeqMetrics_PCT_CODING_BASES + sex + ",
                                "RnaSeqMetrics_PCT_INTERGENIC_BASES + pmi"),

                "MSBB" = paste0("~ 1 + diagnosis + tissue + ",
                                "RnaSeqMetrics_PCT_INTRONIC_BASES + ",
                                "percent_mito + ",
                                "RnaSeqMetrics_PCT_CODING_BASES + ",
                                "RnaSeqMetrics_PCT_INTERGENIC_BASES + ",
                                "AlignmentSummaryMetrics_PCT_PF_READS_ALIGNED"),

                "ROSMAP" = paste0("~ 1 + diagnosis + tissue + percent_mito + ",
                                  "RnaSeqMetrics_PCT_INTRONIC_BASES + ",
                                  "sex + RIN + RIN2 + ",
                                  "RnaSeqMetrics_PCT_CODING_BASES + ",
                                  "RnaSeqMetrics_PCT_INTERGENIC_BASES"))

batch_vars <- list("Mayo" = "flowcell",
                   "MSBB" = c("individualID", "sequencingBatch"),
                   "ROSMAP" = "final_batch")

for (dataset in datasets) {

  ##### Create formulas for fixed/mixed models #####

  design_fixed <- designs[[dataset]]

  rand_effects <- paste0("(1|", batch_vars[[dataset]], ")")
  design_mixed <- paste(c(design_fixed, rand_effects), collapse = " + ")

  design_fixed <- as.formula(design_fixed)
  design_mixed <- as.formula(design_mixed)


  ##### Clean up covariates matrix #####

  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)
  covariates <- Load_Covariates(dataset)
  covariates <- subset(covariates, specimenID %in% colnames(bulk))

  # Ensure categorical covariates are factors
  for (batch in batch_vars[[dataset]]) {
    covariates[,batch] <- factor(covariates[,batch])
  }

  covariates$sex <- factor(covariates$sex)

  # Scale numerical covariates
  for (colname in colnames(covariates)) {
    if (is.numeric(covariates[,colname])) {
      covariates[,colname] <- scale(covariates[,colname])
    }
  }

  # Remove duplicate columns that already exist in colData
  covariates <- covariates %>% select(-diagnosis, -tissue)

  # Merge covariates into the metadata
  col_order <- rownames(colData(bulk))
  metadata <- merge(colData(bulk), covariates, by.x = "sample",
                    by.y = "specimenID", sort = FALSE)
  rownames(metadata) <- metadata$sample

  metadata <- data.frame(metadata[col_order,])


  ##### EdgeR glmQLfit #####

  expr <- DGEList(assay(bulk, "counts"))
  expr$samples$norm.factors <- bulk$tmm_factors # these were pre-calculated

  mod <- model.matrix(design_fixed, data = metadata)

  expr <- estimateDisp(expr, mod)
  fit_edger <- glmQLFit(expr, mod)

  # edgeR uses unshrunk coefficients to calculate fitted values, so I do too
  coefs_edger <- fit_edger$unshrunk.coefficients
  coefs_bio <- grepl("Intercept|diagnosis|tissue|sex",
                     colnames(fit_edger$unshrunk.coefficients))
  coefs_tech <- !coefs_bio

  # Since we don't have residuals from edgeR, we are subtracting technical
  # confounds from the data rather than adding biological confounds to the
  # residuals.
  adjust <- tcrossprod(coefs_edger[,coefs_tech], mod[,coefs_tech])

  # Correct for depth -- not used but keeping code in for reference
  #offset <- fit_edger$offset[1,]
  #offset <- offset - mean(offset)
  #adjust <- sweep(adjust, 2, offset, "+")

  # Note: There are multiple ways this value could be corrected. For a linear
  # model, it would just be exp(log(expr$counts - adjust)). For a GLM, it could
  # also be corrected on a linear scale by adjusting the pearson residuals:
  #   pearson = (expr$counts - mu_orig)/sqrt(variance_orig)
  #   y = mu_new + pearson * sqrt(variance_new)
  # where 'mu_orig' is the original fitted values (linear scale) and 'mu_new' is
  # exp(tcrossprod(coefs[,coefs_bio], mod[,coefs_bio]) + offset). Variance
  # could be (mu) for poisson, (mu + mu^2/theta) for negative binomial, or
  # something else, but edger doesn't provide theta or pearson residuals, so
  # I'm not sure which it is. Since I'm also not sure whether correcting on the
  # linear or log scale is "better", I use log scale, to make it similar to a
  # linear model. Adjusting pearson residuals using variance = mu produces a
  # distribution of counts similar to adjusting on the log scale, so using log
  # scale seems adequate.
  corrected_edger <- exp(log(expr$counts) - adjust) # edgeR is on the natural log scale
  tmm_edger <- calcNormFactors(corrected_edger)


  ##### DESeq2 fit #####

  dds <- DESeqDataSetFromMatrix(assay(bulk, "counts"),
                                metadata,
                                design = design_fixed)

  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  coefs_deseq <- coef(dds)
  mod_deseq <- attr(dds, "modelMatrix")

  coefs_bio <- grepl("Intercept|diagnosis|tissue|sex", colnames(coefs_deseq))
  coefs_tech <- !coefs_bio

  # Since we don't have residuals from DESeq2, we are subtracting technical
  # confounds from the data rather than adding biological confounds to the
  # residuals.
  adjust <- tcrossprod(coefs_deseq[,coefs_tech], mod_deseq[,coefs_tech])

  corrected_deseq2 <- 2^(log2(counts(dds)) - adjust) # deseq2 results are log2 scale
  tmm_deseq2 <- calcNormFactors(corrected_deseq2)


  ##### sageseqr/dream #####

  # This is what sageseqr uses to fit a mixed effect model. The function
  # in sageseqr to calculate the fit + residuals requires CQN counts, which
  # we don't use, and requires jumping through extra hoops to get other args and
  # output into the right format, so I just call voomWithDreamWeights() and
  # dream() myself here.
  n_cores <- parallel::detectCores()-1
  voom_counts <- voomWithDreamWeights(counts = expr,
                                      formula = design_mixed,
                                      data = metadata,
                                      BPPARAM = BiocParallel::SnowParam(n_cores))

  fit_dream <- dream(exprObj = voom_counts,
                     formula = design_mixed,
                     data = metadata,
                     computeResiduals = TRUE,
                     BPPARAM = BiocParallel::SnowParam(n_cores))

  resid_dream <- fit_dream$residuals
  coefs_dream <- coef(fit_dream)
  mod_dream <- fit_dream$design

  # Add biological confounds back to the residuals
  coefs_bio <- grepl("Intercept|diagnosis|tissue|sex", colnames(coefs_dream))
  adjust <- tcrossprod(coefs_dream[,coefs_bio], mod_dream[,coefs_bio])

  corrected_dream <- 2^(resid_dream + adjust) # log2 scale residuals
  tmm_dream <- calcNormFactors(corrected_dream)


  ##### Save results/fits for further examination #####

  saveRDS(list(fit_edger = fit_edger,
               corrected_edger = corrected_edger,
               corrected_deseq2 = corrected_deseq2,
               fit_dream = fit_dream,
               corrected_dream = corrected_dream),
          file.path(dir_tmp, str_glue("{dataset}_models.rds")))


  ##### Create final bulk data object #####

  bulk <- Load_PreprocessedData(dataset, remove_excluded = TRUE)
  assay(bulk, "corrected_edger") <- round(corrected_edger)
  assay(bulk, "corrected_deseq2") <- round(corrected_deseq2)
  assay(bulk, "corrected_dream") <- round(corrected_dream)

  bulk$tmm_factors_edger <- tmm_edger
  bulk$tmm_factors_deseq2 <- tmm_deseq2
  bulk$tmm_factors_dream <- tmm_dream

  Save_BulkData(dataset, bulk)
}
