library(Matrix)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Error_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

cores <- 12
cl <- makeCluster(cores, type = "FORK", outfile = "errors_output.txt")
registerDoParallel(cl)

est_fields = list("Dtangle" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "HSPE" = "estimates",
                  "DeconRNASeq" = "out.all",
                  "DWLS" = "estimates",
                  "Baseline" = "estimates")

algorithms <- names(est_fields)

for (bulk_dataset in bulk_datasets) {
  dir_out <- switch(bulk_dataset,
                    "Mayo" = dir_mayo_output,
                    "MSBB" = dir_msbb_output,
                    "ROSMAP" = dir_rosmap_output)

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

  # Loop over each algorithm's results. Looping isn't strictly necessary but
  # is helpful for controlling what gets processed
  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {bulk_dataset}: {algorithm}"))
    est_field <- est_fields[[algorithm]]

    dir_alg <- file.path(dir_out, algorithm)
    res_files <- list.files(dir_alg, pattern = granularity, full.names = TRUE)
    if (length(res_files) == 0) {
      message(str_glue("No data for {algorithm} found. Skipping..."))
      next
    }

    # Process files
    for (file in res_files) {
      deconv_list <- readRDS(file)

      # If the file contains a null list, skip it
      if (length(deconv_list) == 0) {
        next
      }

      params <- do.call(rbind, lapply(deconv_list, "[[", "params"))

      # For backwards compatibility -- algorithms that input a signature didn't
      # originally have a 'reference_input_type' field so this puts it back if
      # it's missing.
      if (!("reference_input_type" %in% colnames(params))) {
        params$reference_input_type <- "signature"
      }

      # All entries in this file should have the same values for these parameters,
      # so we end up with a 1-row data frame
      params_data <- params %>%
                        select(reference_data_name, test_data_name, granularity,
                               reference_input_type, normalization, regression_method) %>%
                        distinct()

      # If the error file already exists, don't re-process
      tmp <- Load_ErrorList(algorithm, params_data)
      if (!is.null(tmp)) {
        msg <- paste("Error file for", algorithm,
                     paste(params_data, collapse = " "),
                     "found. Skipping...")
        message(msg)
        next
      }

      msg <- paste("Calculating errors for", algorithm,
                   paste(params_data, collapse = " "))
      message(msg)

      # Input data needs to be normalized by depth (cpm, tmm, or tpm). If the
      # original input was 'counts', normalize to CPM. If the original input was
      # on the log scale, normalize to linear scale.
      params_mod <- params_data
      if (params_mod$normalization == "counts") {
        params_mod$normalization <- "cpm"
      }
      params_mod$normalization <- str_replace(params_mod$normalization, "log_", "")

      # We need the signature matrix for error calculation, not the raw data
      params_mod$reference_input_type <- "signature"

      # The data that was used to generate the estimates
      data <- Load_AlgorithmInputData_FromParams(params_mod)
      signature <- data$reference
      bulk_cpm <- assay(data$test, "counts")

      rm(data)
      gc()

      # Needed for print output at the end
      total_length <- length(deconv_list)

      ##### Calculate error for each param set #####
      # Calculated in parallel
      deconv_list <- foreach (P = 1:length(deconv_list)) %dopar% {
        source(file.path("functions", "General_HelperFunctions.R"))
        source(file.path("functions", "Error_HelperFunctions.R"))

        param_id = names(deconv_list)[P]

        # If the error calculation for this param_id exists, don't re-calculate
        tmp <- Load_ErrorIntermediate(algorithm, deconv_list[[param_id]]$params)
        if (!is.null(tmp)) {
          message(paste("Using previously-calculated errors for", algorithm,
                        "/", paste(deconv_list[[param_id]]$params, collapse = " ")))
          return(tmp)
        }

        est_pct <- deconv_list[[param_id]][[est_field]]
        est_pct <- est_pct[colnames(bulk_cpm), colnames(signature)]

        if (any(is.na(est_pct))) {
          message(str_glue("Param set {param_id} has NA values. Skipping..."))
          return(NULL)
        }

        # For broad cell types, at least 25% of the samples need a non-zero
        # estimate for each cell type
        zeros <- colSums(est_pct == 0)
        zero_thresh <- 0.75 * nrow(est_pct) # 75% can be zeros

        if (granularity == "broad_class" && any(zeros > zero_thresh) && algorithm != "Baseline") {
          cts <- paste(colnames(est_pct)[which(zeros > zero_thresh)], collapse = ", ")
          msg <- str_glue(paste("Param set '{param_id}' has too many 0 estimates",
                                "for cell type(s) [{cts}]. Skipping..."))
          message(msg)
          return(NULL)
        }

        if (granularity == "broad_class" &
            any(est_pct[,"Inhibitory"] > est_pct[,"Excitatory"])) {
          num_ests <- sum(est_pct[,"Inhibitory"] > est_pct[,"Excitatory"])
          msg <- str_glue(paste("WARNING: Param set '{param_id}' has {num_ests}",
                                "estimates that contain more inhibitory than",
                                "excitatory neurons."))
          #message(msg)
        }

        params <- deconv_list[[param_id]]$params

        #genes_use <- rownames(signature)
        genes_use <- intersect(highly_variable, rownames(signature))

        sig_filt <- signature[genes_use,]
        bulk_cpm_filt <- bulk_cpm[genes_use,]

        est_expr <- t(est_pct %*% t(sig_filt))

        gof_by_sample <- CalcGOF_BySample(bulk_cpm_filt, est_expr, param_id)
        gof_means <- CalcGOF_Means(gof_by_sample, bulk_metadata, param_id)

        gof_by_sample_lm <- CalcGOF_BySample_LM(bulk_dataset, covariates,
                                                meas_expr_log, est_pct, param_id)
        gof_means_lm <- CalcGOF_Means(gof_by_sample_lm, bulk_metadata, param_id)

        decon_new <- deconv_list[[param_id]]
        decon_new$gof_by_sample <- gof_by_sample
        decon_new$gof_means_all <- gof_means$all_tissue
        decon_new$gof_means_by_tissue <- gof_means$by_tissue

        decon_new$gof_by_sample_lm <- gof_by_sample_lm
        decon_new$gof_means_all_lm <- gof_means_lm$all_tissue
        decon_new$gof_means_by_tissue_lm <- gof_means_lm$by_tissue

        decon_new$params$total_markers_used <- length(decon_new$markers)
        decon_new$estimates_df <- melt(est_pct)
        colnames(decon_new$estimates_df) <- c("sample", "celltype", "percent")

        Save_ErrorIntermediate(decon_new, algorithm)
        print(param_id)
        return(decon_new)
      }

      # Remove NULL (skipped) entries
      deconv_list <- deconv_list[lengths(deconv_list) > 0]

      gof_means_all <- lapply(deconv_list, "[[", "gof_means_all")
      gof_means_tissue <- lapply(deconv_list, "[[", "gof_means_by_tissue")
      gof_means_all_lm <- lapply(deconv_list, "[[", "gof_means_all_lm")
      gof_means_tissue_lm <- lapply(deconv_list, "[[", "gof_means_by_tissue_lm")
      params <- lapply(deconv_list, "[[", "params")
      all_ests <- do.call(rbind, lapply(deconv_list, "[[", "estimates_df"))

      est_stats <- CalcEstimateStats(all_ests, bulk_metadata)

      err_list <- list("means" = list("all_tissue" = do.call(rbind, gof_means_all),
                                      "by_tissue" = do.call(rbind, gof_means_tissue),
                                      "all_tissue_lm" = do.call(rbind, gof_means_all_lm),
                                      "by_tissue_lm" = do.call(rbind, gof_means_tissue_lm)),
                       "params" = do.call(rbind, params),
                       "by_sample" = lapply(deconv_list, "[[", "gof_by_sample"),
                       "by_sample_lm" = lapply(deconv_list, "[[", "gof_by_sample_lm"),
                       "estimate_stats" = est_stats)

      Save_ErrorList(bulk_dataset, err_list, algorithm, params_data)

      print(paste("Errors calculated for", nrow(err_list$params), "of",
                  total_length, "parameter sets."))
    }
  }
}
