# Loops through estimate files for each bulk data set and calculates the error
# (correlation, rMSE, and mAPE) between [signature * estimated_percents] and the
# actual bulk data, for each estimate in the file. Estimates with too many zero
# values for major cell types are discarded.
library(Matrix)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(doParallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

use_top_estimates <- FALSE

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")
singlecell_datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

cores <- 12
cl <- makeCluster(cores, type = "FORK", outfile = "errors_output.txt")
registerDoParallel(cl)

# Which algorithms to calculate errors for
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

# Pre-load all signature matrices
all_signatures_cpm <- sapply(singlecell_datasets, function(X) {
  Load_SignatureMatrix(X, granularity, "cpm")
})
all_signatures_tmm <- sapply(singlecell_datasets, function(X) {
  Load_SignatureMatrix(X, granularity, "tmm")
})

# Get all genes shared by all data sets so everything is judged on the same criteria
common_genes <- lapply(bulk_datasets, function(X) {
  bulk_se <- Load_BulkData(X, output_type = "log_cpm", regression_method = "none")
  rownames(bulk_se)
})
common_genes <- append(common_genes, lapply(all_signatures_cpm, rownames))
common_genes <- purrr::reduce(common_genes, intersect)

# Filter the signatures to only the genes being used
filtered_signatures_cpm <- lapply(all_signatures_cpm, function(X) {
  X[common_genes, ]
})
filtered_signatures_tmm <- lapply(all_signatures_tmm, function(X) {
  X[common_genes, ]
})

gc()

for (bulk_dataset in bulk_datasets) {
  # Loop over each algorithm's results
  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {bulk_dataset}: {algorithm}"))

    if (use_top_estimates) {
      dir_alg <- file.path(dir_top_estimates, bulk_dataset, algorithm)
    } else {
      dir_alg <- file.path(dir_estimates, bulk_dataset, algorithm)
    }

    res_files <- list.files(dir_alg, pattern = granularity, full.names = TRUE)

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

      params <- do.call(rbind, lapply(deconv_list, "[[", "params"))

      # All entries in this file should have the same values for these parameters,
      # so we end up with a 1-row data frame
      params_data <- params %>%
        select(reference_data_name, test_data_name, granularity,
               reference_input_type, normalization, regression_method) %>%
        distinct()

      if (use_top_estimates) {
        params_data$mode <- "best_estimates"
      }

      # If the error file already exists, don't re-process
      tmp <- Load_ErrorList(algorithm, params_data, top_params = use_top_estimates)
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

      # The data that was used to generate the estimates
      data <- Load_BulkData(params_mod$test_data_name,
                            output_type = params_mod$normalization,
                            regression_method = params_mod$regression_method)

      bulk_metadata <- colData(data)

      # The "counts" assay is actually CPM, TMM, or TPM-normalized data
      bulk_cpm <- assay(data, "counts")

      rm(data)
      gc()

      if (params_mod$normalization == "tmm") {
        filtered_signatures <- filtered_signatures_tmm
      } else {
        filtered_signatures <- filtered_signatures_cpm
      }

      # Needed for print output at the end
      total_length <- length(deconv_list)
      all_possible <- all_possible + total_length

      ## Calculate error for each parameter set --------------------------------

      deconv_list <- foreach(P = 1:length(deconv_list)) %dopar% {
        source(file.path("functions", "General_HelperFunctions.R"))
        source(file.path("functions", "Step11_Error_HelperFunctions.R"))

        param_id <- names(deconv_list)[[P]]

        if (use_top_estimates) {
          deconv_list[[param_id]]$params$mode <- "best_estimates"
        }

        # If the error calculation for this param_id exists, don't re-calculate
        tmp <- Load_ErrorIntermediate(algorithm, deconv_list[[param_id]]$params)
        if (!is.null(tmp)) {
          message(paste("Using previously-calculated errors for", algorithm, "/",
                        paste(deconv_list[[param_id]]$params, collapse = " ")))
          return(tmp)
        }

        est_pct <- deconv_list[[param_id]]$estimates

        if (any(is.na(est_pct))) {
          message(str_glue("Param set {param_id} has NA values. Skipping..."))
          return(NULL)
        }

        # At least 75% of the samples need a non-zero estimate for the most
        # abundant cell types (>10% of the population in the Cain dataset). We
        # assume that if more than 25% of estimates for a supposedly-abundant
        # cell type are *exactly* 0, that estimate is not reliable.
        if (granularity == "broad_class") {
          major_celltypes <- c("Astrocyte", "Excitatory", "Inhibitory", "Oligodendrocyte")

        } else {
          # sub_class uses cell types >5% of the population, plus the most abundant
          # excitatory and inhibitory subclasses
          major_celltypes <- c("Astrocyte", "Exc.1", "Inh.1", "Oligodendrocyte")
        }

        zeros <- colSums(est_pct == 0)[major_celltypes]
        zero_thresh <- 0.25 * nrow(est_pct) # 25% can be zeros

        if (algorithm != "Baseline" && any(zeros > zero_thresh)) {
          cts <- paste(names(which(zeros > zero_thresh)), collapse = ", ")
          msg <- str_glue(paste("Param set '{param_id}' has too many 0 estimates",
                                "for cell type(s) [{cts}]. Skipping..."))
          # message(msg)
          return(NULL)
        }

        params <- deconv_list[[param_id]]$params

        bulk_cpm_filt <- bulk_cpm[common_genes, rownames(est_pct)]

        # Calculate error using each signature, since we don't know which, if
        # any, are more accurate
        gof_by_sample <- lapply(names(filtered_signatures), function(N) {
          est_expr <- t(est_pct %*% t(filtered_signatures[[N]]))

          gof <- CalcGOF_BySample(bulk_cpm_filt, est_expr, param_id)
          gof$signature <- N
          return(gof)
        })
        names(gof_by_sample) <- names(filtered_signatures)

        gof_means <- lapply(gof_by_sample, function(gof) {
          return(CalcGOF_Means(gof, bulk_metadata, param_id))
        })
        gof_means <- do.call(rbind, gof_means)

        decon_new <- list("gof_by_sample" = do.call(rbind, gof_by_sample),
                          "gof_means" = gof_means,
                          "params" = deconv_list[[param_id]]$params,
                          "param_id" = param_id)

        decon_new$params$total_markers_used <- length(deconv_list[[param_id]]$markers)
        rownames(decon_new$params) <- param_id

        Save_ErrorIntermediate(decon_new, algorithm)
        # print(param_id)
        return(decon_new)
      }

      # Remove NULL (skipped) entries
      deconv_list <- deconv_list[lengths(deconv_list) > 0]
      all_valid <- all_valid + length(deconv_list)

      # If there are no valid parameter sets, we just skip this file
      if (length(deconv_list) == 0) {
        print(paste("No valid parameter sets for", algorithm,
                    paste(params_data, collapse = " ")))
        next
      }

      # Combine all the separate results into single data frames (or vectors)
      gof_means_all <- do.call(rbind, lapply(deconv_list, "[[", "gof_means"))
      params <- do.call(rbind, lapply(deconv_list, "[[", "params"))

      err_list <- list("means" = gof_means_all,
                       "by_sample" = lapply(deconv_list, "[[", "gof_by_sample"),
                       "params" = params,
                       "param_ids" = sapply(deconv_list, "[[", "param_id"))

      names(err_list$by_sample) <- rownames(err_list$params)

      Save_ErrorList(bulk_dataset, err_list, algorithm, params_data,
                     top_params = use_top_estimates)

      print(paste("Errors calculated for", length(deconv_list), "of",
                  total_length, "parameter sets."))
    }
    print(str_glue("{all_valid} of {all_possible}"))
  }
}

stopCluster(cl)
