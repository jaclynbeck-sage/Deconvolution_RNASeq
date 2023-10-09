library(Matrix)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

source(file.path("functions", "General_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

cores <- 6
cl <- makeCluster(cores, type = "FORK", outfile = str_glue("errors_output.txt"))
registerDoParallel(cl)

est_fields = list("Dtangle" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "HSPE" = "estimates",
                  "DeconRNASeq" = "out.all",
                  "DWLS" = "estimates")

algorithms <- names(est_fields)
algorithms <- c("DeconRNASeq")

err_list <- list()

for (bulk_dataset in bulk_datasets) {
  dir_out <- switch(bulk_dataset,
                    "Mayo" = dir_mayo_output,
                    "MSBB" = dir_msbb_output,
                    "ROSMAP" = dir_rosmap_output)

  # The data that will be used in the GLM (needs unadjusted counts for bulk
  # data)
  bulk_se <- Load_BulkData(bulk_dataset, output_type = "counts",
                           regression_method = "none")
  covariates <- Load_Covariates(bulk_dataset)
  covariates <- Clean_BulkCovariates(bulk_dataset, colData(bulk_se), covariates)

  # The set of genes to use for error calculation -- highly variable genes as
  # determined by log_cpm values
  bulk_tmp <- Load_BulkData(bulk_dataset, output_type = "log_cpm")
  var_genes <- apply(assay(bulk_tmp, "counts"), 1, var)
  var_genes <- sort(var_genes, decreasing = TRUE)
  rm(bulk_tmp)

  # Loop over each algorithm's results. Looping isn't strictly necessary but
  # is helpful for controlling what gets processed
  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {bulk_dataset}: {algorithm}"))
    est_field <- est_fields[[algorithm]]

    res_files <- list.files(dir_out, pattern = algorithm, full.names = TRUE)

    # Process files in parallel
    foreach (F = 1:length(res_files)) %dopar% {
      source(file.path("functions", "General_HelperFunctions.R"))
      source(file.path("functions", "Error_HelperFunctions.R"))

      deconv_list <- readRDS(res_files[F])

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

      # Input data needs to be normalized by depth (cpm, tmm, or tpm). If the
      # original input was 'counts', normalize to CPM. If the original input was
      # on the log scale, normalize to linear scale.
      if (params_data$normalization == "counts") {
        params_data$normalization <- "cpm"
      }
      params_data$normalization <- str_replace(params_data$normalization, "log_", "")

      # We need the signature matrix for error calculation, not the raw data
      params_data$reference_input_type <- "signature"

      # The data that was used to generate the estimate
      data <- Load_AlgorithmInputData_FromParams(params_data)
      signature <- data$reference
      bulk_cpm <- assay(data$test, "counts")

      # The data that will be used in the GLM (needs unadjusted signature in CPM)
      sig_for_glm <- Load_SignatureMatrix(params_data$reference_data_name,
                                          params_data$granularity,
                                          output_type = "cpm")

      highly_variable <- var_genes[var_genes %in% rownames(signature)][1:5000]

      #marker_rows <- which(params$marker_type != "None")

      #common_markers <- c()
      #for (row in marker_rows) {
      #  common_markers <- c(common_markers, deconv_list[[row]]$markers)
      #}

      ##### Calculate error for each param set #####
      for (param_id in names(deconv_list)) {
        est_pct <- deconv_list[[param_id]][[est_field]]
        est_pct <- est_pct[colnames(bulk_cpm), colnames(signature)]

        if (any(is.na(est_pct))) {
          message(str_glue("Param set {param_id} has NA values. Skipping..."))
          next
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
          next
        }

        params <- deconv_list[[param_id]]$params

        #markers <- rownames(signature)
        markers <- highly_variable
        #markers <- common_markers

        sig_filt <- signature[markers,]
        bulk_cpm_filt <- bulk_cpm[markers,]
        bulk_se_filt <- bulk_se[markers,]
        sig_glm_filt <- sig_for_glm[markers,]

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

        deconv_list[[param_id]]$gof_by_sample <- gof_by_sample
        deconv_list[[param_id]]$gof_means_all <- gof_means$all_tissue
        deconv_list[[param_id]]$gof_means_by_tissue <- gof_means$by_tissue
        #deconv_list[[param_id]]$gof_glm_by_sample <- gof_glm_by_sample
        #deconv_list[[param_id]]$gof_glm_means_all <- gof_glm_means$all_tissue
        #deconv_list[[param_id]]$gof_glm_means_by_tissue <- gof_glm_means$by_tissue
      }

      gof_means_all <- lapply(deconv_list, "[[", "gof_means_all")
      gof_means_tissue <- lapply(deconv_list, "[[", "gof_means_by_tissue")
      params <- lapply(deconv_list, "[[", "params")

      err_list <- list("means" = list("all_tissue" = do.call(rbind, gof_means_all),
                                      "by_tissue" = do.call(rbind, gof_means_tissue)),
                       "params" = do.call(rbind, params),
                       "by_sample" = lapply(deconv_list, "[[", "gof_by_sample"))

      Save_ErrorList(bulk_dataset, err_list, algorithm, params_data)
      print(paste("Errors calculated for", length(err_list), "of",
                  length(deconv_list), "parameter sets."))
      print("Done")
    }
  }
}

