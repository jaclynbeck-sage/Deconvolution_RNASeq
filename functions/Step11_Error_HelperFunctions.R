library(dplyr)
library(Matrix)
library(Metrics)
library(lme4)

source(file.path("functions", "FileIO_HelperFunctions.R"))

# Goodness-of-fit error calculations -------------------------------------------

# Calculate goodness of fit (GOF) for each sample across 3 metrics (correlation,
# root mean squared error, and mean absolute percent error).
#
# Arguments:
#   meas_expr_cpm - the observed bulk data matrix, on the linear scale (CPM-,
#                   TMM-, or TPM-normalized values depending on the parameters)
#   est_expr - the estimated expression matrix calculated from multiplying a
#              signature matrix by the estimated cell type percentages, on the
#              linear scale
#   param_id - a string uniquely identifying this set of parameters
#
# Returns:
#   a data frame where rows are samples and columns are the error metrics plus
#   some metadata
CalcGOF_BySample <- function(meas_expr_cpm, est_expr, param_id) {
  meas_expr_cpm <- as.matrix(meas_expr_cpm)

  # Using log2 transformation for correlation so distributions are more gaussian
  cor_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    cor(log2(meas_expr_cpm[, sample] + 1), log2(est_expr[, sample] + 1),
        use = "na.or.complete")
  })

  rmse_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    rmse(meas_expr_cpm[, sample], est_expr[, sample])
  })

  mape_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    ok <- meas_expr_cpm[, sample] != 0 | est_expr[, sample] != 0
    smape(meas_expr_cpm[ok, sample], est_expr[ok, sample]) / 2
  })

  gof_by_sample <- data.frame(cor = cor_sample,
                              rMSE = rmse_sample,
                              mAPE = mape_sample,
                              param_id = param_id,
                              sample = names(cor_sample),
                              row.names = names(cor_sample))

  return(gof_by_sample)
}


# Calculate goodness of fit (GOF) for each sample across 3 metrics using a
# linear model and the formula for each bulk data set.
#
# Arguments:
#   bulk_dataset_name - the name of the bulk data set
#   covariates - a data frame of covariates for this data set
#   meas_expr_log - the observed bulk data matrix, on the log2-scale (log2 of
#                   CPM-normalized values)
#   est_pct - a matrix of estimated percentages where rows are samples and
#             columns are cell types
#   param_id - a string uniquely identifying this set of parameters
#
# Returns:
#   same data frame format as CalcGOF_BySample
CalcGOF_BySample_LM <- function(bulk_dataset_name, covariates, meas_expr_log,
                                est_pct, param_id) {
  covariates_ct <- cbind(covariates[colnames(meas_expr_log), ],
                         est_pct[colnames(meas_expr_log), ])

  formulas <- Load_ModelFormulas(bulk_dataset_name)

  # Replace all biological variables with the cell type proportions instead
  form <- str_replace(formulas$formula_mixed, "~ diagnosis \\+ tissue",
                      paste(colnames(est_pct), collapse = " + "))
  form <- str_replace(form, "sex \\+ ", "")

  # ROSMAP benefits from a linear mixed effect model
  # For speed we're using a fixed effect model for ROSMAP, but saving this code
  # just in case.
  if (bulk_dataset_name == "UNUSED") { # "ROSMAP") {
    form <- paste("expr ~ 0 +", form)

    lm_est <- sapply(rownames(meas_expr_log), function(gene) {
      expr <- meas_expr_log[gene, ]

      res <- lmer(as.formula(form), data = covariates_ct)
      return(fitted(res))
    })
  }
  # Mayo and MSBB can use fixed effects only -- using batch as a random effect
  # causes singular fit for a lot of genes, so we don't do that
  else {
    form <- str_replace(form, " \\+ \\(.*", "") # remove random effect term
    form <- paste("t(meas_expr_log) ~ 0 +", form)
    lm_est <- lm(as.formula(form), data = covariates_ct)
    lm_est <- fitted(lm_est)
  }

  # Reverse the log2(cpm+1) transform
  lm_est <- t(2^lm_est-1)

  return(CalcGOF_BySample(2^meas_expr_log-1, lm_est, param_id))
}


# Takes the mean goodness-of-fit across all samples, for each error metric
#
# Arguments:
#   gof_by_sample - a data frame where rows are samples and columns are the
#                   error metrics plus some metadata, as returned by CalcGOF_BySample
#   bulk_metadata - a data frame of metadata where rows are samples, and one
#                   column is "tissue"
#   param_id - a string uniquely identifying this set of parameters
#
# Returns:
#   a data frame containing one row for the means across all samples, plus one
#   row for each tissue containing the means across all samples in that tissue,
#   for each signature used in the error calculation
CalcGOF_Means <- function(gof_by_sample, bulk_metadata, param_id) {
  # take the mean of each numerical column across all samples
  gof_means_all <- gof_by_sample %>%
    summarise(across(where(is.numeric), mean),
              param_id = unique(param_id),
              .groups = "drop")

  gof_means_all$tissue <- "All"

  tissue <- bulk_metadata[rownames(gof_by_sample), "tissue"]

  gof_by_sample <- cbind(gof_by_sample, tissue)

  # mean of numerical columns for samples within a given tissue
  gof_means_tissue <- gof_by_sample %>%
    group_by(tissue) %>%
    summarise(across(where(is.numeric), mean),
              param_id = unique(param_id),
              .groups = "drop")

  gof_means <- rbind(gof_means_all, gof_means_tissue[, colnames(gof_means_all)])
  gof_means$signature <- unique(gof_by_sample$signature)
  gof_means$solve_type <- unique(gof_by_sample$solve_type)

  return(gof_means)
}


# Calculates statistics about goodness-of-fit data:
#   - the mean and SD of estimated percentages for each cell type for a given
#     sample across all estimates in the file
#   - the same mean and SD except only across the top 10-scoring estimates in
#     the file
#   - the mean of these means across all samples
#
# Arguments:
#   all_ests - a melted data frame with columns "param_id", "sample", "celltype",
#              and "percent", representing the estimates for each cell type for
#              each sample
#   bulk_metadata - a data frame of metadata, where rows are bulk samples and
#                   one column is "tissue"
#   gof_means_all - a data frame with N rows per parameter set (where N = [1 +
#                   the number of tissues] * [the number of signatures]), and
#                   columns are the mean error metrics as returned by CalcGOF_Means
#
# Returns:
#   a list of statistics, containing items "by_sample", "all_tissue",
#   "by_tissue" (all data frames), and "top_10_stats" (a list, one entry for
#   each error metric)
CalcEstimateStats <- function(all_ests, bulk_metadata, gof_means_all) {
  # Mean of each cell type percentage over all estimates in the file, per
  # celltype per sample
  estimate_stats_sample <- all_ests %>%
    group_by(celltype, sample) %>%
    summarize(mean_pct = mean(percent),
              sd_pct = sd(percent),
              rel_sd_pct = sd_pct / mean_pct,
              .groups = "drop_last")

  estimate_stats_sample$tissue <- bulk_metadata[estimate_stats_sample$sample, "tissue"]

  # Top 10 scoring parameter sets per tissue, per signature used to
  # calculate the error. Each error metric may have different top 10 lists. If
  # there are less than 10 parameter sets, this just selects all of them.
  top_cor <- gof_means_all %>% group_by(tissue, signature) %>% top_n(10, wt = cor)
  top_rMSE <- gof_means_all %>% group_by(tissue, signature) %>% top_n(10, wt = -rMSE)
  top_mAPE <- gof_means_all %>% group_by(tissue, signature) %>% top_n(10, wt = -mAPE)

  top_10 <- list(cor = top_cor,
                 rMSE = top_rMSE,
                 mAPE = top_mAPE)

  ests_tmp <- all_ests
  ests_tmp$tissue <- bulk_metadata[ests_tmp$sample, "tissue"]

  # Stats for each of the error metrics on the top 10 scoring parameter sets
  # for each tissue
  top_10_stats <- lapply(top_10, function(top_errs) {
    # Calculate stats for each tissue (which may be a specific tissue or "All")
    ests_top_10 <- lapply(unique(top_errs$tissue), function(tiss) {
      top_tmp <- subset(top_errs, tissue == tiss)

      top_stats <- lapply(unique(top_tmp$signature), function(sig) {
        top_sub <- subset(top_tmp, signature == sig)
        ests_filt <- subset(ests_tmp, param_id %in% top_sub$param_id)

        if (tiss != "All") {
          ests_filt <- subset(ests_filt, tissue == tiss)
        }

        # Mean cell type percent for each cell type for each sample, across the
        # top 10 scoring parameter sets
        ests_stats <- ests_filt %>%
          group_by(celltype, sample) %>%
          summarize(mean_pct = mean(percent),
                    sd_pct = sd(percent),
                    rel_sd_pct = sd_pct / mean_pct,
                    .groups = "drop_last")

        ests_stats$tissue <- tiss
        ests_stats$signature <- sig
        return(ests_stats)
      })

      top_stats <- do.call(rbind, top_stats)
    })

    ests_top_10 <- do.call(rbind, ests_top_10)
    return(ests_top_10)
  })

  # average and sd of the *mean* values for each sample, per cell type
  estimate_stats_all <- estimate_stats_sample %>%
    group_by(celltype) %>%
    summarize(mean_pct = mean(mean_pct),
              mean_sd_pct = mean(sd_pct),
              mean_rel_sd_pct = mean(rel_sd_pct),
              .groups = "drop_last")

  # Same but broken down by tissue
  estimate_stats_tissue <- estimate_stats_sample %>%
    group_by(tissue, celltype) %>%
    summarize(mean_pct = mean(mean_pct),
              mean_sd_pct = mean(sd_pct),
              mean_rel_sd_pct = mean(rel_sd_pct),
              .groups = "drop_last")

  return(list("by_sample" = estimate_stats_sample,
              "all_tissue" = estimate_stats_all,
              "by_tissue" = estimate_stats_tissue,
              "top_10_stats" = top_10_stats))
}


# IHC error functions ----------------------------------------------------------

# NOTE: This function has not been used or tested in a long time so it probably
# doesn't work.
CalcError_MeanIHC_AllTissues <- function(err_list) {
  # IHC props and percent errors are per individual, get the mean for each
  # error type across all individuals
  errs_ihc_props <- lapply(err_list[["errs_ihc_props"]], function(err) {
    err %>% summarize(across(colnames(err), ~ mean(.x, na.rm = TRUE)))
  })

  errs_ihc_props <- do.call(rbind, errs_ihc_props) %>%
    dplyr::rename(cor_props = cor,
                  rMSE_props = rMSE,
                  mAPE_props = mAPE)

  errs_ihc_pct <- lapply(err_list[["errs_ihc_pct"]], function(err) {
    err %>% summarize(across(colnames(err), ~ mean(.x, na.rm = TRUE)))
  })

  errs_ihc_pct <- do.call(rbind, errs_ihc_pct) %>%
    dplyr::rename(cor_pct = cor,
                  rMSE_pct = rMSE,
                  mAPE_pct = mAPE)

  return(cbind(errs_ihc_props, errs_ihc_pct))
}


# NOTE: This function has not been used or tested in a long time so it probably
# doesn't work.
CalcError_MeanIHC_ByTissue <- function(err_list, params, bulk_meta) {
  errs_ihc_props <- err_list[["errs_ihc_props"]][names(params)]
  errs_ihc_props <- lapply(names(errs_ihc_props), function(N) {
    err <- errs_ihc_props[[N]]
    err <- cbind(err, tissue = bulk_meta[rownames(err), "tissue"])

    df <- err %>%
      group_by(tissue) %>%
      summarize(across(!contains("tissue"), ~ mean(.x, na.rm = TRUE)))
    df$name <- N

    return(df)
  })

  errs_ihc_props <- do.call(rbind, errs_ihc_props) %>%
    dplyr::rename(cor_props = cor,
                  rMSE_props = rMSE,
                  mAPE_props = mAPE)

  errs_ihc_pct <- err_list[["errs_ihc_pct"]][names(params)]

  errs_ihc_pct <- lapply(names(errs_ihc_pct), function(N) {
    err <- errs_ihc_pct[[N]]
    err <- cbind(err, tissue = bulk_meta[rownames(err), "tissue"])

    df <- err %>%
      group_by(tissue) %>%
      summarize(across(!contains("tissue"), ~ mean(.x, na.rm = TRUE)))
    df$name <- N

    return(df)
  })

  errs_ihc_pct <- do.call(rbind, errs_ihc_pct) %>%
    dplyr::rename(cor_pct = cor,
                  rMSE_pct = rMSE,
                  mAPE_pct = mAPE)

  errs_ihc <- merge(errs_ihc_props, errs_ihc_pct, by = c("name", "tissue"))
  return(errs_ihc)
}


# Utility functions for loading error files ------------------------------------

# UNUSED -- delete if Get_AllErrorsAsDf is deleted
CalcError_MeanByTissue <- function(bulk_dataset, err_list, params, bulk_meta) {
  # TODO add metadata to error df in calculate errors function
  errs_by_sample <- err_list[["gof_sample"]][names(params)]
  errs_by_sample <- lapply(errs_by_sample, function(df) {
    cbind(df, tissue = bulk_meta[rownames(df), "tissue"])
  })

  errs_by_tissue <- lapply(names(errs_by_sample), function(E) {
    df <- errs_by_sample[[E]] %>%
      group_by(tissue) %>%
      summarize(across(c("cor", "rMSE", "mAPE"), ~ mean(.x, na.rm = TRUE)))
    df$name <- E
    return(df)
  })
  errs_by_tissue <- do.call(rbind, errs_by_tissue)

  if (bulk_dataset == "ROSMAP" & "errs_ihc_props" %in% names(err_list)) {
    errs_ihc <- CalcError_MeanIHC_ByTissue(err_list, params, bulk_meta)
    errs_by_tissue <- merge(errs_by_tissue, errs_ihc, by = c("name", "tissue"))
  }

  return(errs_by_tissue)
}


# UNUSED -- delete if Get_AllErrorsAsDf is deleted
Filter_Params <- function(err_list, best_params) {
  params <- lapply(err_list[["params"]], function(X) {
    as.list(X %>% select(-test_data_name))
  })

  params_use <- params %in% best_params

  return(params[params_use])
}


# Loads all the "best_errors" files into one big list. The errors can be
# concatenated into one data frame but the params lists have to be coded into
# strings because they differ between algorithms so they have differing numbers
# of columns.
#
# Arguments:
#   bulk_datasets - a list of bulk dataset names to load ("Mayo", "MSBB", and/or
#                   "ROSMAP")
#   granularity - either "broad_class" or "sub_class"
#
# Returns:
#   a list containing three entries: "errors" is a data frame with all errors
#   from all files concatenated together, "quality stats" is a data frame
#   with all statistics from all files concatenated together, and "estimates"
#   is a data frame with all estimates from all files concatenated together
Get_AllBestErrorsAsDf <- function(bulk_datasets, granularity) {
  # List of "best errors" files matching the bulk data set names and granularity
  file_list <- list.files(dir_best_errors,
                          pattern = paste0("(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*", granularity),
                          full.names = TRUE,
                          recursive = TRUE)

  best_err_list <- lapply(file_list, function(file) {
    data <- readRDS(file)
    print(basename(file))

    # Encode the parameters as strings. We don't need the reference_data_name
    # or test_data_name because they're the same for every entry in this file,
    # and total_markers_used isn't really a parameter
    param_strings <- format(data$params, trim = TRUE) %>%
      select(-reference_data_name, -test_data_name, -total_markers_used) %>%
      apply(1, paste, collapse = "_")

    param_strings <- data.frame(param_id = names(param_strings),
                                param_string = as.character(param_strings))

    # Add the parameter strings to the main errors data frame
    errs_df <- merge(data$means$all_signature, param_strings, by = "param_id")

    # Merge a subset of the actual parameters into the data frame. These
    # parameters exist in every error file for every algorithm
    errs_df <- merge(errs_df,
                     select(data$params, reference_data_name, test_data_name,
                            reference_input_type, normalization,
                            regression_method),
                     by.x = "param_id", by.y = "row.names")

    errs_df$algorithm <- str_replace(errs_df$param_id, "_.*", "")

    qual_stats <- Get_QualityStatsAsDf(data)

    return(list("errors" = errs_df,
                #"estimates" = data$estimates,
                "quality_stats" = qual_stats))
  })

  return(list("errors" = do.call(rbind, lapply(best_err_list, "[[", "errors")),
              #"estimates" = do.call(rbind, lapply(best_err_list, "[[", "estimates")),
              "quality_stats" = do.call(rbind, lapply(best_err_list, "[[", "quality_stats"))))
}


Get_AllBestEstimatesAsDf <- function(bulk_datasets, granularity, metadata, best_params = NULL) {
  # List of "best errors" files matching the bulk data set names and granularity
  file_list <- list.files(dir_top_estimates,
                          pattern = paste0("(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*", granularity),
                          full.names = TRUE,
                          recursive = TRUE)

  best_ests_list <- lapply(file_list, function(file) {
    data <- readRDS(file)
    print(basename(file))

    data <- data[names(data) %in% best_params$param_id]

    for (N in names(data)) {
      tissues <- subset(best_params, param_id == N)$tissue
      samples <- subset(metadata, tissue %in% tissues)$sample

      new_ests <- as.data.frame(data[[N]]$estimates[as.character(samples),])
      new_ests$sample <- rownames(new_ests)

      new_ests <- new_ests %>%
        melt(id.vars = c("sample"),
             variable.name = "celltype",
             value.name = "percent") %>%
        merge(combined_metadata, by = "sample", all.y = FALSE)

      new_ests$param_id <- N
      new_ests$algorithm <- str_replace(N, "_.*", "")

      data[[N]]$estimates <- new_ests
    }

    ests_df <- do.call(rbind, lapply(data, "[[", "estimates"))
    params <- do.call(rbind, lapply(data, "[[", "params"))

    return(list("estimates" = ests_df,
                "params" = params))
  })

  return(do.call(rbind, lapply(best_ests_list, "[[", "estimates")))
}


# UNUSED -- possibly delete
Get_AllBestParamsAsDf <- function(reference_datasets, granularity) {
  best_params_all <- lapply(reference_datasets, function(ref_dataset) {
    data <- readRDS(file.path(dir_output, str_glue("best_params_{ref_dataset}_{granularity}.rds")))
    data$reference_data_name <- ref_dataset
    return(data)
  })
  best_params_all <- do.call(rbind, best_params_all)
  return(best_params_all)
}


# Helper function to extract statistics from an error file and concatenate them
# into a data frame
#
# Arguments:
#   data - a named list of the best errors from each file as generated by
#          Step12_Get_TopParamSets
#
# Returns:
#   a data frame containing columns for parameters and for percent bad inhibitory
#   ratio data and percent of valid estimates data from each file
Get_QualityStatsAsDf <- function(data) {
  # Get the basic parameters that define each file
  params_mod <- data$params %>%
    select(reference_data_name, test_data_name, granularity,
           normalization, regression_method) %>%
    distinct()

  rownames(params_mod) <- str_replace(rownames(params_mod), "_[0-9].*", "")

  params_mod$mean_pct_bad_inhibitory_ratio <- mean(data$pct_bad_inhibitory_ratio)
  params_mod$pct_valid_results <- data$n_valid_results / data$n_possible_results

  return(params_mod)
}


# UNUSED -- possibly delete
Get_AllErrorsAsDf <- function(bulk_datasets, reference_datasets, algorithms, granularity, best_params) {
  # Read in all errors for each dataset & data type,
  # put in one big dataframe
  errs_all <- lapply(bulk_datasets, function(bulk_dataset) {
    bulk_se <- Load_BulkData(bulk_dataset)
    bulk_meta <- colData(bulk_se)

    errs_r <- lapply(reference_datasets, function(ref_dataset) {
      errs_a <- lapply(algorithms, function(alg) {
        err_list <- Load_ErrorList(alg, ref_dataset, bulk_dataset, granularity)

        if (length(err_list) == 0) {
          return(NULL)
        }

        params <- Filter_Params(err_list, best_params)
        param_ids <- sapply(params, paste, collapse = " ")

        errs <- err_list[["gof_mean"]]

        if (bulk_dataset == "ROSMAP" & "errs_ihc_props" %in% names(err_list)) {
          errs_ihc <- CalcError_MeanIHC_AllTissues(err_list)
          errs <- cbind(errs, errs_ihc)
        }

        errs <- errs[names(params), ]
        errs$tissue <- "All"
        errs$name <- rownames(errs)

        errs_by_tissue <- CalcError_MeanByTissue(bulk_dataset, err_list, params, bulk_meta)

        errs <- rbind(errs, errs_by_tissue)

        if (!("cor_props" %in% colnames(errs))) {
          errs$cor_props <- NA
          errs$rMSE_props <- NA
          errs$mAPE_props <- NA
          errs$cor_pct <- NA
          errs$rMSE_pct <- NA
          errs$mAPE_pct <- NA
        }

        errs$method <- alg
        errs$reference_data_name <- ref_dataset
        errs$test_data_name <- bulk_dataset
        errs$param_id <- param_ids[errs$name]

        return(errs)
      })

      if (all(sapply(errs_a, is.null)) || length(errs_a) == 0) {
        return(NULL)
      }

      return(do.call(rbind, errs_a))
    })

    errs_r <- do.call(rbind, errs_r)
    return(errs_r)
  })

  errs_all <- do.call(rbind, errs_all)

  return(errs_all)
}


# Load all estimates from all files matching the arguments to this function,
# into a single data frame
#
# Arguments:
#   ref_dataset - the name of the reference (single cell) dataset
#   bulk_dataset - the name of the bulk data set
#   algorithms - a list of algorithms to load data for
#   granularity - either "broad_class" or "sub_class"
#   best_errors - the "errors" data frame from the list output by
#                 Get_AllBestErrorsAsDf
#   est_fields - a named list describing what field corresponds to cell type
#                estimates for each algorithm
#
# Returns:
#   a melted data frame with columns for cell type, estimated percent, sample,
#   param id, and algorithm
Get_AllEstimatesAsDf <- function(ref_dataset, bulk_dataset, algorithms,
                                 granularity, best_errors, est_fields) {
  ests_alg <- list()

  for (algorithm in algorithms) {
    if (algorithm == "Baseline") {
      next
    }
    est_field <- est_fields[[algorithm]]

    # Find the base parameters for each file we want to load
    alg_name <- algorithm
    file_params <- best_errors %>%
      subset(reference_data_name == ref_dataset &
               test_data_name == bulk_dataset &
               algorithm == alg_name) %>%
      select(reference_data_name, test_data_name, algorithm,
             reference_input_type, normalization, regression_method) %>%
      distinct()

    if (nrow(file_params) == 0) {
      print(str_glue(paste0("No valid data found for {algorithm} / ",
                            "{bulk_dataset} / {ref_dataset} / {granularity}. ",
                            "Skipping...")))
      next
    }

    # Each row of file_params describes a single estimate file -- load that
    # file, filter to the best parameter sets, and add some extra metadata to
    # the data frame for each parameter set
    ests_list <- lapply(1:nrow(file_params), function(P) {
      ests <- Load_AlgorithmOutputList(
        algorithm,
        ref_dataset,
        bulk_dataset,
        granularity,
        reference_input_type = file_params$reference_input_type[P],
        normalization = file_params$normalization[P],
        regression_method = file_params$regression_method[P]
      )

      # Filter to the best parameter sets and get the estimates data frames
      ests <- ests[names(ests) %in% best_errors$param_id]
      ests <- lapply(ests, "[[", est_field)

      # We need to add "sample" and "param_id" columns to the estimates data
      # frames so they can be identified when concatenated
      ests <- lapply(names(ests), FUN = function(X) {
        ests[[X]] <- as.data.frame(ests[[X]])
        ests[[X]]$sample <- rownames(ests[[X]])
        ests[[X]]$param_id <- X
        return(ests[[X]])
      })

      ests_df <- do.call(rbind, ests)

      return(ests_df)
    })

    # Concatenate all the estimates together
    ests_df <- do.call(rbind, ests_list)

    if (length(ests_df) == 0) {
      next
    }

    ests_df <- melt(ests_df, variable.name = "celltype", value.name = "pct_est")
    ests_df$algorithm <- algorithm

    ests_alg[[algorithm]] <- ests_df
  }

  # Concatenate all data from all algorithms together
  ests_alg <- do.call(rbind, ests_alg)

  return(ests_alg)
}


# Takes a nested list of errors and "flattens" it so that it is no longer nested
# (i.e. all entries are data frames instead of lists of data frames)
#
# Arguments:
#   err_list - a nested list, where each entry is a list of data frames or vectors
#
# Returns:
#   a list where each entry is either a data frame or vector
Flatten_ErrorList <- function(err_list) {
  fields <- names(err_list[[1]])
  err_list_new <- lapply(fields, function(field) {
    if (is(err_list[[1]][[field]], "data.frame")) {
      return(do.call(rbind, lapply(err_list, "[[", field))) # concat data frames
    } else {
      return(unlist(lapply(err_list, "[[", field))) # unlist vectors
    }
  })
  names(err_list_new) <- fields
  return(err_list_new)
}
