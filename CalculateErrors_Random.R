library(Matrix)
library(SummarizedExperiment)
library(Metrics)
library(dplyr)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Error_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef") #, "seaAD")
granularity <- "broad_class"
bulk_dataset <- "Mayo"

random_types <- c("uniform", "mean_singlecell")

# The data that will be used in the LM (needs unadjusted log2(cpm))
bulk_se <- Load_BulkData(bulk_dataset, output_type = "log_cpm",
                         regression_method = "none")

meas_expr_log <- as.matrix(assay(bulk_se, "counts"))

# Get highly variable genes
var_genes <- rowVars(meas_expr_log, useNames = TRUE)
var_genes <- sort(var_genes, decreasing = TRUE)

highly_variable <- names(var_genes)[1:5000]
meas_expr_log <- meas_expr_log[highly_variable,]
bulk_metadata <- colData(bulk_se)

covariates <- Load_Covariates(bulk_dataset)
covariates <- Clean_BulkCovariates(bulk_dataset, bulk_metadata, covariates)

rm(bulk_se)
gc()

for (dataset in datasets) {
  data <- Load_AlgorithmInputData(dataset, bulk_dataset, granularity,
                                  reference_input_type = "signature",
                                  output_type = "cpm",
                                  regression_method = "none")

  signature <- data$reference
  bulk_cpm <- assay(data$test, "counts")

  errs_gof_mean <- list()
  errs_gof_subject <- list()
  params <- list()

  genes_use <- intersect(highly_variable, rownames(signature))
  bulk_filt <- bulk_cpm[genes_use,]
  sig_filt <- signature[genes_use,]

  for (rand_type in random_types) {
    param_id <- str_glue("{rand_type}_{dataset}_{bulk_dataset}_{granularity}")

    # Randomly generate proportions based on uniform distribution
    if (rand_type == "uniform") {
      set.seed(1234)
      est_pct <- matrix(runif(ncol(bulk_filt) * ncol(sig_filt)),
                        nrow = ncol(bulk_filt))
      est_pct <- sweep(est_pct, 1, rowSums(est_pct), "/") # rows sum to 1
    }
    # Use the mean of the single cell dataset plus gaussian noise
    else if (rand_type == "mean_singlecell") {
      pb <- Load_Pseudobulk(dataset, "sc_samples", granularity)
      pcts <- colMeans(pb@metadata[["pctRNA"]])

      est_pct <- matrix(rep(pcts/sum(pcts), ncol(bulk_filt)), nrow = ncol(bulk_filt),
                        byrow = TRUE)

      set.seed(5678)
      for (ct in 1:ncol(est_pct)) {
        # Generate noise with an SD relative to the mean percent of each cell type
        est_pct[,ct] <- est_pct[,ct] + rnorm(ncol(bulk_filt), mean = 0,
                                             sd = pcts[ct]/4)
      }

      # Ensure no negative numbers and that all rows sum to 1
      est_pct[est_pct < 0] <- 0
      est_pct <- sweep(est_pct, 1, rowSums(est_pct), "/")
    }

    rownames(est_pct) <- colnames(bulk_filt)
    colnames(est_pct) <- colnames(sig_filt)

    est_expr <- t(est_pct %*% t(sig_filt))

    gof_by_sample <- CalcGOF_BySample(bulk_filt, est_expr, param_id)
    gof_means <- CalcGOF_Means(gof_by_sample, bulk_metadata, param_id)

    gof_by_sample_lm <- CalcGOF_BySample_LM(bulk_dataset, covariates,
                                            meas_expr_log, est_pct, param_id)
    gof_means_lm <- CalcGOF_Means(gof_by_sample_lm, bulk_metadata,
                                  param_id)



    errs_gof_mean[[param_name]] <- gof_list[["means"]]
    errs_gof_subject[[param_name]] <- gof_list[["by_subject"]]

    params[[param_name]] <- data.frame(reference_data_name = dataset,
                                       test_data_name = bulk_dataset,
                                       granularity = granularity,
                                       random_type = rand_type)
  }

  err_list <- list("gof_mean" = do.call(rbind, errs_gof_mean),
                   "gof_subject" = errs_gof_subject,
                   "errs_ihc_props" = errs_ihc_props,
                   "errs_ihc_pct" = errs_ihc_pct,
                   "params" = params)

  Save_ErrorList(err_list, "random", dataset, bulk_dataset, granularity)
  print("Done")
}

