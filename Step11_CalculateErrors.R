library(Matrix)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

cores <- 12
cl <- makeCluster(cores, type = "FORK", outfile = "errors_output.txt")
registerDoParallel(cl)

est_fields = list("CibersortX" = "estimates",
                  "DeconRNASeq" = "out.all",
                  "Dtangle" = "estimates",
                  "DWLS" = "estimates",
                  "HSPE" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "Scaden" = "estimates",
                  "Baseline" = "estimates")

algorithms <- names(est_fields)

all_signatures_cpm <- sapply(c("cain", "lau", "leng", "mathys", "seaRef"), function(X) {
  Load_SignatureMatrix(X, granularity, "cpm")
})
all_signatures_tmm <- sapply(c("cain", "lau", "leng", "mathys", "seaRef"), function(X) {
  Load_SignatureMatrix(X, granularity, "tmm")
})

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

  bulk_metadata <- colData(bulk_se)

  covariates <- Load_Covariates(bulk_dataset)
  covariates <- Clean_BulkCovariates(bulk_dataset, bulk_metadata, covariates)

  rm(bulk_se)
  gc()

  # Get the top 5000 most variable genes that exist in all signatures and in
  # the bulk dataset
  genes_use <- names(var_genes)
  for (name in names(all_signatures_cpm)) {
    genes_use <- intersect(genes_use, rownames(all_signatures_cpm[[name]]))
  }
  genes_use <- genes_use[1:5000]

  filtered_signatures_cpm <- lapply(all_signatures_cpm, function(X) {
    X[genes_use,]
  })
  filtered_signatures_tmm <- lapply(all_signatures_tmm, function(X) {
    X[genes_use,]
  })

  meas_expr_log <- meas_expr_log[genes_use,]

  # Loop over each algorithm's results
  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {bulk_dataset}: {algorithm}"))
    est_field <- est_fields[[algorithm]]

    dir_alg <- file.path(dir_out, algorithm)
    res_files <- list.files(dir_alg, pattern = granularity, full.names = TRUE)
    if (length(res_files) == 0) {
      message(str_glue("No data for {algorithm} found. Skipping..."))
      next
    }

    # Process files in parallel
    for (file in res_files) {
      source(file.path("functions", "General_HelperFunctions.R"))
      source(file.path("functions", "Step11_Error_HelperFunctions.R"))

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
      #if (!(params_mod$reference_input_type %in% c("signature", "cibersortx"))) {
      if (params_mod$reference_input_type != "signature") {
        params_mod$reference_input_type <- "signature"
      }

      # The data that was used to generate the estimates
      #data <- Load_AlgorithmInputData_FromParams(params_mod)
      data <- Load_BulkData(params_mod$test_data_name,
                            output_type = params_mod$normalization,
                            regression_method = params_mod$regression_method)
      #signature <- data$reference

      if (params_mod$normalization == "tmm") {
        filtered_signatures <- filtered_signatures_tmm
      }
      else {
        filtered_signatures <- filtered_signatures_cpm
      }

      # CibersortX doesn't subtract the 1 when calculating its signature
      #if (params_mod$reference_input_type == "cibersortx") {
      #  signature[signature == 1] <- 0
      #}

      bulk_cpm <- assay(data, "counts")

      rm(data)
      gc()

      # Needed for print output at the end
      total_length <- length(deconv_list)

      ##### Calculate error for each param set #####

      foreach (P = 1:length(deconv_list)) %dopar% {
      #deconv_list <- lapply(names(deconv_list), function(param_id) {
        param_id <- names(deconv_list)[[P]]

        # If the error calculation for this param_id exists, don't re-calculate
        tmp <- Load_ErrorIntermediate(algorithm, deconv_list[[param_id]]$params)
        if (!is.null(tmp)) {
          message(paste("Using previously-calculated errors for", algorithm,
                        "/", paste(deconv_list[[param_id]]$params, collapse = " ")))
          return(tmp)
        }

        est_pct <- deconv_list[[param_id]][[est_field]]
        est_pct <- est_pct[colnames(bulk_cpm), colnames(filtered_signatures[[1]])]

        if (any(is.na(est_pct))) {
          message(str_glue("Param set {param_id} has NA values. Skipping..."))
          return(NULL)
        }

        # At least 75% of the samples need a non-zero estimate for the most
        # abundant cell types (>10% of the population in the Cain dataset). We
        # assume that if more than 25% of estimates for a supposedly-abundant
        # cell type are *exactly* 0, that estimate is not reliable.
        if (algorithm != "Baseline") {
          if (granularity == "broad_class") {
            major_celltypes <- c("Astrocyte", "Excitatory", "Inhibitory", "Oligodendrocyte")
          }
          # sub_class uses cell types >5% of the population, plus the most abundant
          # inhibitory subclass (due to the difficulty predicting inhibitory percentages
          # in the broad_class case)
          else {
            major_celltypes <- c("Astrocyte", "Exc.1", "Exc.4", "Inh.1", "Oligodendrocyte")
          }

          zeros <- colSums(est_pct == 0)[major_celltypes]
          zero_thresh <- 0.25 * nrow(est_pct) # 25% can be zeros

          if (any(zeros > zero_thresh)) {
            cts <- paste(names(which(zeros > zero_thresh)), collapse = ", ")
            msg <- str_glue(paste("Param set '{param_id}' has too many 0 estimates",
                                  "for cell type(s) [{cts}]. Skipping..."))
            message(msg)
            return(NULL)
          }
        }

        # Check for if the algorithm estimated more inhibitory than excitatory
        # neurons, which generally means the estimate is bad
        pct_bad_inhibitory_ratio <- 0
        if (granularity == "broad_class") {
          num_ests <- sum(est_pct[,"Inhibitory"] > est_pct[,"Excitatory"])
        }
        else if (granularity == "sub_class") {
          excitatory_cols <- grepl("Exc", colnames(est_pct))
          inhibitory_cols <- grepl("Inh", colnames(est_pct))
          summed_estimates <- data.frame(Excitatory = rowSums(est_pct[,excitatory_cols]),
                                         Inhibitory = rowSums(est_pct[,inhibitory_cols]))
          num_ests <- sum(summed_estimates$Inhibitory > summed_estimates$Excitatory)
        }

        pct_bad_inhibitory_ratio <- num_ests / nrow(est_pct)

        params <- deconv_list[[param_id]]$params

        #sig_filt <- signature[intersect(genes_use, rownames(signature)),]
        bulk_cpm_filt <- bulk_cpm[genes_use,]

        #est_expr <- t(est_pct %*% t(sig_filt))
        #gof_by_sample <- CalcGOF_BySample(bulk_cpm_filt, est_expr, param_id)

        # Calculate error using all signatures
        gof_by_sample <- lapply(names(filtered_signatures), function(N) {
          est_expr <- t(est_pct %*% t(filtered_signatures[[N]]))
          gof <- CalcGOF_BySample(bulk_cpm_filt, est_expr, param_id)
          gof$signature <- N
          gof$solve_type <- "signature"
          return(gof)
        })
        names(gof_by_sample) <- names(filtered_signatures)

        #gof_means <- CalcGOF_Means(gof_by_sample, bulk_metadata, param_id)
        gof_means <- lapply(gof_by_sample, function(gof) {
          return(CalcGOF_Means(gof, bulk_metadata, param_id))
        })
        gof_means <- do.call(rbind, gof_means)

        # Calculate error using signature-less lm
        gof_by_sample_lm <- CalcGOF_BySample_LM(bulk_dataset, covariates,
                                                meas_expr_log, est_pct, param_id)
        gof_by_sample_lm$signature <- "none"
        gof_by_sample_lm$solve_type <- "lm"

        gof_means_lm <- CalcGOF_Means(gof_by_sample_lm, bulk_metadata, param_id)

        decon_new <- list("gof_by_sample" = do.call(rbind, gof_by_sample),
                          "gof_means" = gof_means,
                          "gof_by_sample_lm" = gof_by_sample_lm,
                          "gof_means_lm" = gof_means_lm,
                          "params" = deconv_list[[param_id]]$params,
                          "pct_bad_inhibitory_ratio" = pct_bad_inhibitory_ratio,
                          "estimates" = as.data.frame(est_pct))

        decon_new$params$total_markers_used <- length(deconv_list[[param_id]]$markers)
        rownames(decon_new$params) <- param_id

        decon_new$estimates$sample <- rownames(decon_new$estimates)
        decon_new$estimates$param_id <- param_id

        Save_ErrorIntermediate(decon_new, algorithm)
        #print(param_id)
        return(decon_new)
      }

      # Remove NULL (skipped) entries
      deconv_list <- deconv_list[lengths(deconv_list) > 0]

      if (length(deconv_list) == 0) {
        print(paste("No valid parameter sets for", algorithm,
                    paste(params_data, collapse = " ")))
        next
      }

      gof_means_all <- do.call(rbind, lapply(deconv_list, "[[", "gof_means_all"))
      gof_means_all_lm <- do.call(rbind, lapply(deconv_list, "[[", "gof_means_all_lm"))
      params <- do.call(rbind, lapply(deconv_list, "[[", "params"))
      pct_inh <- sapply(deconv_list, "[[", "pct_bad_inhibitory_ratio")
      names(pct_inh) <- rownames(params)
      all_ests <- do.call(rbind, lapply(deconv_list, "[[", "estimates"))
      all_ests <- melt(all_ests, variable.name = "celltype",
                       value.name = "percent", id.vars = c("param_id", "sample"))

      est_stats <- CalcEstimateStats(all_ests, bulk_metadata)

      err_list <- list("means" = list("all_signature" = gof_means_all,
                                      "all_lm" = gof_means_all_lm),
                       "params" = params,
                       "by_sample" = lapply(deconv_list, "[[", "gof_by_sample"),
                       "by_sample_lm" = lapply(deconv_list, "[[", "gof_by_sample_lm"),
                       "estimate_stats" = est_stats,
                       "pct_bad_inhibitory_ratio" = pct_inh,
                       "n_valid_results" = length(deconv_list),
                       "n_possible_results" = total_length,
                       "estimates" = all_ests)

      names(err_list$by_sample) <- rownames(err_list$params)
      names(err_list$by_sample_lm) <- rownames(err_list$params)

      Save_ErrorList(bulk_dataset, err_list, algorithm, params_data)

      print(paste("Errors calculated for", nrow(err_list$params), "of",
                  total_length, "parameter sets."))
    }
  }
}
