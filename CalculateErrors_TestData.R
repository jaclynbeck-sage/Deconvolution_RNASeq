library(Matrix)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

source(file.path("functions", "General_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

cores <- 12
cl <- makeCluster(cores, type = "FORK", outfile = str_glue("errors_output.txt"))
registerDoParallel(cl)

est_fields = list("Dtangle" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "HSPE" = "estimates",
                  "DeconRNASeq" = "out.all",
                  "DWLS" = "estimates")

algorithms <- names(est_fields)

err_list <- list()

for (bulk_dataset in bulk_datasets) {
  dir_out <- switch(bulk_dataset,
                    "Mayo" = dir_mayo_output,
                    "MSBB" = dir_msbb_output,
                    "ROSMAP" = dir_rosmap_output)

  # The data that will be used in the LM (needs unadjusted log2(cpm))
  #bulk_se <- Load_BulkData(bulk_dataset, output_type = "log_cpm",
  #                         regression_method = "none")
  #covariates <- Load_Covariates(bulk_dataset)
  #covariates <- Clean_BulkCovariates(bulk_dataset, colData(bulk_se), covariates)

  # Loop over each algorithm's results. Looping isn't strictly necessary but
  # is helpful for controlling what gets processed
  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {bulk_dataset}: {algorithm}"))
    est_field <- est_fields[[algorithm]]

    dir_alg <- file.path(dir_out, algorithm)
    res_files <- list.files(dir_alg, pattern = algorithm, full.names = TRUE)
    if (length(res_files) == 0) {
      message(str_glue("No data for {algorithm} found. Skipping..."))
      next
    }

    # Process files
    for (F in res_files) {
      deconv_list <- readRDS(F)

      # If the file contains a null list, skip it
      if (length(deconv_list) == 0) {
        next
      }

      params <- do.call(rbind, lapply(deconv_list, "[[", "params"))

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

      # The data that was used to generate the estimate
      data <- Load_AlgorithmInputData_FromParams(params_mod)
      signature <- data$reference
      bulk_cpm <- assay(data$test, "counts")

      # The data that will be used in the GLM (needs unadjusted signature in CPM)
      #sig_for_glm <- Load_SignatureMatrix(params_data$reference_data_name,
      #                                    params_data$granularity,
      #                                    output_type = "cpm")

      #highly_variable <- var_genes[var_genes %in% rownames(signature)][1:5000]

      #marker_rows <- which(params$marker_type != "None")

      #common_markers <- c()
      #for (row in marker_rows) {
      #  common_markers <- c(common_markers, deconv_list[[row]]$markers)
      #}

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

        if (granularity == "broad_class" & any(zeros > zero_thresh)) {
          cts <- paste(colnames(est_pct)[which(zeros > zero_thresh)], collapse = ", ")
          msg <- str_glue(paste("Param set '{param_id}' has too many 0 estimates",
                                "for cell type(s) [{cts}]. Skipping..."))
          message(msg)
          return(NULL)
        }

        params <- deconv_list[[param_id]]$params

        markers <- rownames(signature)
        #markers <- highly_variable
        #markers <- common_markers

        sig_filt <- signature[markers,]
        bulk_cpm_filt <- bulk_cpm[markers,]
        #bulk_se_filt <- bulk_se[markers,]
        #sig_glm_filt <- sig_for_glm[markers,]

        est_expr <- t(est_pct %*% t(sig_filt))

        gof_by_sample <- CalcGOF_BySample(bulk_cpm_filt, est_expr, param_id)
        gof_means <- CalcGOF_Means(gof_by_sample, colData(data$test), param_id)

        #log_lib_size <- log(colSums(assay(bulk_se, "counts")))
        #gof_glm_by_sample <- CalcGOF_BySample_GLM(bulk_dataset, covariates,
        #                                          bulk_se_filt, est_pct,
        #                                          sig_glm_filt, log_lib_size,
        #                                          param_id)
        #gof_glm_means <- CalcGOF_Means(gof_glm_by_sample, colData(data$test),
        #                               param_id)

        decon_new <- deconv_list[[param_id]]
        decon_new$gof_by_sample <- gof_by_sample
        decon_new$gof_means_all <- gof_means$all_tissue
        decon_new$gof_means_by_tissue <- gof_means$by_tissue
        decon_new$params$total_markers_used <- length(decon_new$markers)
        #deconv_list[[param_id]]$gof_glm_by_sample <- gof_glm_by_sample
        #deconv_list[[param_id]]$gof_glm_means_all <- gof_glm_means$all_tissue
        #deconv_list[[param_id]]$gof_glm_means_by_tissue <- gof_glm_means$by_tissue

        Save_ErrorIntermediate(decon_new, algorithm)
        print(param_id)
        return(decon_new)
      }

      # Remove NULL (skipped) entries
      deconv_list <- deconv_list[lengths(deconv_list) > 0]

      gof_means_all <- lapply(deconv_list, "[[", "gof_means_all")
      gof_means_tissue <- lapply(deconv_list, "[[", "gof_means_by_tissue")
      params <- lapply(deconv_list, "[[", "params")

      err_list <- list("means" = list("all_tissue" = do.call(rbind, gof_means_all),
                                      "by_tissue" = do.call(rbind, gof_means_tissue)),
                       "params" = do.call(rbind, params),
                       "by_sample" = lapply(deconv_list, "[[", "gof_by_sample"))

      Save_ErrorList(bulk_dataset, err_list, algorithm, params_data)
      print(paste("Errors calculated for", nrow(err_list$params), "of",
                  length(deconv_list), "parameter sets."))
      print("Done")
    }
  }
}
