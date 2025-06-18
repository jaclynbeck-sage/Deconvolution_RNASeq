# Loops through estimate files for each bulk data set and calculates the error
# (correlation, rMSE, and mAPE) between [signature * estimated_percents] and the
# actual bulk data, for each estimate in the file. Estimates with too many zero
# values for major cell types are discarded.
library(Matrix)
library(dplyr)
library(stringr)
library(purrr)
library(parallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- all_bulk_datasets()
singlecell_datasets <- all_singlecell_datasets()

cores <- 12
cluster_type <- "FORK"
cluster_outfile <- "errors_output.txt"

# Which algorithms to calculate errors for
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

# Pre-load all signature matrices
signatures <- Get_Signatures(singlecell_datasets, granularity)

# Get all genes shared by all data sets so everything is judged on the same criteria
common_genes <- Get_CommonGenes(bulk_datasets, signatures)

# Filter the signatures to only the genes being used
filtered_signatures_list <- Get_FilteredSignatures(signatures, common_genes)
rm(signatures)

gc()

for (bulk_dataset in bulk_datasets) {
  # Loop over each algorithm's results
  for (algorithm in algorithms) {
    res_files <- list.files(file.path(dir_estimates, bulk_dataset, algorithm),
                            pattern = granularity, full.names = TRUE)

    if (length(res_files) == 0) {
      message(str_glue("No data for {algorithm} found. Skipping..."))
      next
    }

    all_possible <- 0
    all_valid <- 0

    # Process each file
    for (file in res_files) {
      deconv_list <- readRDS(file)

      # If the file contains a null list, skip it
      if (length(deconv_list) == 0) {
        next
      }

      # All entries in this file should have the same values for these parameters,
      # so we end up with a 1-row data frame
      file_params <- List_to_DF(deconv_list, "params") %>%
        FileParams_FromParams()

      # If the error file already exists, don't re-process
      tmp <- Load_ErrorList(file_params, top_params = FALSE)
      if (!is.null(tmp)) {
        msg <- paste("Error file for", paste(file_params, collapse = " "),
                     "found. Skipping...")
        message(msg)
        next
      }

      message(paste("Calculating errors for", paste(file_params, collapse = " ")))

      # Input data needs to be normalized by depth (cpm, tmm, or tpm). If the
      # original input was 'counts', normalize to CPM. If the original input was
      # on the log scale, normalize to linear scale.
      params_mod <- file_params %>%
        mutate(normalization = str_replace(normalization, "counts", "cpm"),
               normalization = str_replace(normalization, "log_", ""))

      # The data that was used to generate the estimates
      data <- Load_BulkData(params_mod$test_data_name,
                            output_type = params_mod$normalization,
                            regression_method = params_mod$regression_method)

      bulk_metadata <- as.data.frame(colData(data))

      # The "counts" assay is actually CPM, TMM, or TPM-normalized data
      bulk_cpm <- assay(data, "counts")

      rm(data)
      gc()

      filtered_signatures <- filtered_signatures_list[[params_mod$normalization]]

      # Needed for print output at the end
      total_length <- length(deconv_list)
      all_possible <- all_possible + total_length

      ## Calculate error for each parameter set --------------------------------

      cl <- makeCluster(cores, type = cluster_type, outfile = cluster_outfile)

      error_list <- parLapply(cl, deconv_list, function(deconv_result) {
        source(file.path("functions", "General_HelperFunctions.R"))
        source(file.path("functions", "Step11_Error_HelperFunctions.R"))

        param_id <- deconv_result$param_id
        params <- deconv_result$params
        est_pct <- deconv_result$estimates

        if (any(is.na(est_pct))) {
          message(str_glue("Param set {param_id} has NA values. Skipping..."))
          return(NULL)
        }

        # If the error calculation for this param_id exists, don't re-calculate
        tmp <- Load_ErrorIntermediate(params)
        if (!is.null(tmp)) {
          message(paste("Using previously-calculated errors for",
                        paste(params, collapse = " ")))
          return(tmp)
        }

        # At least 75% of the samples need a non-zero estimate for the most
        # abundant cell types. We assume that if more than 25% of estimates for
        # a supposedly-abundant cell type are *exactly* 0, that estimate is not
        # reliable.
        if (Too_Many_Zeros(est_pct, granularity, zero_thresh = 0.25)) {
          return(NULL)
        }

        bulk_cpm_filt <- bulk_cpm[common_genes, rownames(est_pct)]

        # Calculate error using each signature, since we don't know which, if
        # any, are more accurate
        gof_by_sample <- lapply(names(filtered_signatures), function(N) {
          est_expr <- filtered_signatures[[N]] %*% t(est_pct)

          gof <- CalcGOF_BySample(bulk_cpm_filt, est_expr, param_id)
          gof$signature <- N
          return(gof)
        })

        gof_means <- lapply(gof_by_sample, function(gof) {
          return(CalcGOF_Means(gof, bulk_metadata, param_id))
        })

        errors <- list("gof_by_sample" = List_to_DF(gof_by_sample),
                       "gof_means" = List_to_DF(gof_means),
                       "params" = params,
                       "param_id" = param_id)

        errors$params$total_markers_used <- length(deconv_result$markers)
        rownames(errors$params) <- param_id

        Save_ErrorIntermediate(errors)
        return(errors)
      })

      stopCluster(cl)

      # Remove NULL (skipped) entries
      error_list <- error_list[lengths(error_list) > 0]
      all_valid <- all_valid + length(error_list)

      # If there are no valid parameter sets, we just skip this file
      if (length(error_list) == 0) {
        print(paste("No valid parameter sets for", paste(file_params, collapse = " ")))
        next
      }

      # Combine all the separate results into single data frames (or vectors)
      errs_all <- list("means" = List_to_DF(error_list, "gof_means"),
                       "by_sample" = lapply(error_list, "[[", "gof_by_sample"),
                       "params" = List_to_DF(error_list, "params"),
                       "param_ids" = sapply(error_list, "[[", "param_id"))

      names(errs_all$by_sample) <- rownames(errs_all$params)

      Save_ErrorList(bulk_dataset, errs_all, algorithm, file_params,
                     top_params = FALSE)

      print(paste("Errors calculated for", length(error_list), "of",
                  total_length, "parameter sets."))
    }
    print(str_glue("{all_valid} of {all_possible}"))
  }
}
