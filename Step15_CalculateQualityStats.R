# ...
# This script also calculates some statistics about the estimates in each file:
#   - how many samples in each estimate have a "bad" inhibitory:excitatory ratio
#   - the mean and SD of estimated percentages for each cell type for a given
#     sample across all estimates in the file
#   - the same mean and SD except only across the top 10-scoring estimates in
#     the file
#   - the mean of these means across all samples

library(dplyr)
library(stringr)
library(parallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step14_Analysis_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

cores <- 12

# Which algorithms to calculate stats for
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

for (bulk_dataset in bulk_datasets) {
  # Get bulk metadata
  bulk <- Load_BulkData(bulk_dataset)
  bulk_metadata <- colData(bulk) %>%
    as.data.frame() %>%
    select(sample, tissue)
  rm(bulk)
  gc()

  alg_qstats_all <- lapply(algorithms, function(algorithm) {
    est_files <- list.files(file.path(dir_estimates, bulk_dataset, algorithm),
                            pattern = granularity,
                            recursive = TRUE, full.names = TRUE)

    if (length(est_files) == 0) {
      message(str_glue("No data found for {bulk_dataset}/{algorithm}/{granularity}. Skipping..."))
      next
    }

    print(paste(bulk_dataset, algorithm))

    # Collect some per-file stats ----------------------------------------------

    # Results from this loop will be combined at the end to calculate top-level
    # stats.
    file_qstats <- mclapply(est_files, function(est_f) {
      file_id <- str_replace(basename(est_f), "estimates_", "") %>%
        str_replace(".rds", "")

      est_list <- readRDS(est_f)
      file_params <- FileParams_FromParams(est_list[[1]]$params)

      err_f <- list.files(file.path(dir_errors, bulk_dataset, algorithm),
                          pattern = file_id, full.names = TRUE)

      best_est_f <- list.files(file.path(dir_top_estimates, bulk_dataset, algorithm),
                               pattern = file_id, full.names = TRUE)

      best_err_f <- list.files(file.path(dir_best_errors, bulk_dataset, algorithm),
                               pattern = file_id, full.names = TRUE)

      # If there weren't any valid results for this file, the error file won't
      # exist. Create an abbreviated quality stats file for the sole purpose of
      # keeping track of number of failures.
      if (length(err_f) != 1 || length(best_est_f) != 1 || length(best_err_f) != 1) {
        print(str_glue(
          paste("No valid error, best estimates, and/or best errors files found",
                "for {basename(est_f)}. Returning abbreviated stats data.")))

        n_valid_results <- data.frame(n_valid = 0,
                                      n_possible = length(est_list),
                                      file_id = file_id)
        n_valid_results <- cbind(n_valid_results, file_params)

        return(list("n_valid_results" = n_valid_results))
      }

      err_list <- readRDS(err_f)
      best_err_list <- readRDS(best_err_f)
      best_est_list <- readRDS(best_est_f)


      ## Valid vs possible results ---------------------------------------------

      # Must be calculated before subsetting to valid results
      n_valid_results <- data.frame(n_valid = length(err_list$param_ids),
                                    n_possible = length(est_list),
                                    file_id = file_id)
      n_valid_results <- cbind(n_valid_results, file_params)


      ## Fix Baseline params columns -------------------------------------------

      # Baseline params have a "trial" column in the middle of the file-specific
      # parameters, but are missing a "reference_input_type" column or have it
      # in the wrong place, so we remove trial and add a reference_input_type
      # column in the right place

      if (algorithm == "Baseline") {
        if ("trial" %in% colnames(best_err_list$params)) {
          best_err_list$params <- select(best_err_list$params, -trial)
        }

        if (!("reference_input_type" %in% colnames(best_err_list$params))) {
          best_err_list$params <- mutate(best_err_list$params,
                                         reference_input_type = "signature",
                                         .after = granularity)
        } else {
          best_err_list$params <- relocate(best_err_list$params,
                                           reference_input_type,
                                           .after = granularity)
        }
      }

      ## Subset to valid results -----------------------------------------------

      # Backward compatibility: Only keep errors that were calculated against
      # the same signature as the reference data. We don't subset Baseline
      # data since it doesn't use the single cell references as input
      if (algorithm != "Baseline") {
        ref_data <- unique(best_err_list$params$reference_data_name)
        best_err_list$means <- subset(best_err_list$means,
                                      signature == ref_data)
      } else if (unique(best_err_list$params$reference_data_name) == "zeros") {
        # For the "zeros" Baseline data, errors for all signatures are the same,
        # so we just subset to "cain".
        best_err_list$means <- best_err_list$means %>%
          mutate(signature = NA) %>%
          distinct()
      }

      # Backward compatibility: Getting rid of "best_mean_rank" errors -- we
      # only want top cor, RMSE, or MAPE
      errors_weighted <- best_err_list$means %>%
        subset(tissue != "All") %>%
        Get_TopRanked("tissue", n_top = 1, with_mean_rank = FALSE)

      best_param_ids <- intersect(best_err_list$param_ids,
                                  errors_weighted$param_id)
      best_params <- best_err_list$params[best_param_ids, ] %>%
        mutate(param_id = best_param_ids)

      # Subset to valid estimates only. For the "best estimates", we got rid of
      # any best errors that were for mean_rank only, so we need to subset the
      # best estimates list too
      est_list <- est_list[err_list$param_ids]
      best_est_list <- best_est_list[best_param_ids]

      # Get all best estimates as one data frame
      est_pcts <- lapply(best_est_list, function(est_item) {
        estimates <- as.data.frame(est_item$estimates)
        estimates$sample <- rownames(estimates)

        estimates %>% merge(bulk_metadata, by = "sample", all = FALSE) %>%
          mutate(param_id = est_item$param_id)
      })
      est_pcts <- List_to_DF(est_pcts)

      # Create the same kind of duplication in the estimates df for estimate
      # stats calculations
      weights <- select(errors_weighted, param_id, tissue, type)

      est_pcts_weighted <- merge(est_pcts, weights,
                                 by = c("tissue", "param_id"),
                                 all = FALSE) %>%
        select(-type)

      # Use non-duplicated data for exc:inh ratio and number of zeros
      est_pcts <- distinct(est_pcts_weighted)


      ## Ratio of excitatory to inhibitory cells -------------------------------
      # For best estimates

      # This cbind works for both broad and sub classes
      neuron_ests <- cbind(
        select(est_pcts, sample, tissue, param_id),
        data.frame(
          Excitatory = rowSums(select(est_pcts, starts_with("Exc"))),
          Inhibitory = rowSums(select(est_pcts, starts_with("Inh")))
        )) %>%
        mutate(is_bad_estimate = Inhibitory > Excitatory,
               exc_inh_ratio = Excitatory / Inhibitory)

      # Change these to NA so they get removed when doing mean/median below
      neuron_ests$exc_inh_ratio[is.infinite(neuron_ests$exc_inh_ratio)] <- NA

      exc_inh_ratio <- neuron_ests %>%
        group_by(tissue, param_id) %>%
        dplyr::summarize(n_bad = sum(is_bad_estimate),
                         count = n(),
                         pct_bad_inhibitory_ratio = n_bad / count,
                         mean_exh_inh_ratio = mean(exc_inh_ratio, na.rm = TRUE),
                         median_exh_inh_ratio = median(exc_inh_ratio, na.rm = TRUE),
                         param_id = unique(param_id),
                         .groups = "drop") %>%
        select(-n_bad, -count)


      ## Number of 0 guesses for each cell type --------------------------------
      # For best estimates

      num_zeros <- est_pcts %>%
        group_by(tissue, param_id) %>%
        summarize(across(where(is.numeric), ~ sum(.x == 0)), .groups = "drop")


      ## Error and estimate stats ----------------------------------------------
      # For best errors/estimates

      err_stats <- errors_weighted %>%
        Calculate_ErrorStats(group_cols = "tissue") %>%
        mutate(file_id = file_id) %>%
        cbind(file_params)

      est_stats <- est_pcts_weighted %>%
        Calculate_EstimateStats(group_cols = c("sample", "tissue")) %>%
        mutate(file_id = file_id) %>%
        cbind(file_params)


      # Remove duplication for top level analysis
      errs_sub <- select(errors_weighted, -type) %>% distinct()

      return(list("n_valid_results" = n_valid_results,
                  "exc_inh_ratio" = exc_inh_ratio,
                  "n_zero_guesses" = num_zeros,
                  "error_stats" = err_stats,
                  "estimate_stats" = est_stats,
                  "best_estimates" = est_pcts,
                  "best_errors" = errs_sub,
                  "best_params" = best_params))
    }, mc.cores = cores)


    # Top errors for the algorithm --------------------------------------------

    best_errors <- List_to_DF(file_qstats, "best_errors")
    best_ests <- List_to_DF(file_qstats, "best_estimates")
    best_params <- List_to_DF(file_qstats, "best_params")

    # Merge some params in for top-level analysis
    best_errors <- merge(best_errors,
                         best_params[, c("param_id", Get_ParameterColumnNames())],
                         by = "param_id")

    group_cols <- c("tissue", Get_ParameterColumnNames()) %>%
      setdiff(c("reference_data_name", "reference_input_type"))

    # Only keep error information for the best parameters. Calling this function
    # may produce duplicated rows where a param_id is the best for multiple
    # error metrics, which is what we want to get a weighted mean: if a
    # parameter set is the best for multiple metrics, it should count more than
    # a parameter set that is only the best for one metric.
    best_errs_toplevel <- Get_TopRanked(best_errors, group_cols, n_top = 1,
                                        with_mean_rank = FALSE)

    # Create the same kind of duplication in the estimates df. Also pull in
    # parameters listed in group_cols for this calculation
    weights <- best_errs_toplevel %>%
      select_at(c("param_id", group_cols, "type"))

    best_ests_toplevel <- merge(best_ests, weights,
                                by = c("tissue", "param_id"), all = FALSE)

    err_stats_toplevel <- Calculate_ErrorStats(best_errs_toplevel, group_cols)
    est_stats_toplevel <- select(best_ests_toplevel, -type) %>%
      Calculate_EstimateStats(c("sample", group_cols))

    best_params_toplevel <- subset(best_params,
                                   param_id %in% best_errs_toplevel$param_id)

    if (algorithm != "Baseline") {
    # Check for algorithm-specific parameters and if they exist calculate how
    # many times each value for those parameters shows up in the top 3 errors
    # for each error metric
      alg_specific <- select(best_params,
                             -all_of(Get_ParameterColumnNames()), # remove file params
                             -contains("filter_level"), -contains("marker"), # remove general marker params
                             -mode)
      if (ncol(alg_specific) > 1) {
        top3 <- Get_TopRanked(best_errors, "tissue", n_top = 3, with_mean_rank = FALSE)
        top3 <- merge(top3, alg_specific, by = "param_id")

        params_test <- setdiff(colnames(alg_specific), "param_id")
        param_frequency <- lapply(params_test, function(param_col) {
          Count_ParamFrequency(top3, groups = c("tissue", param_col),
                               pivot_column = param_col,
                               algorithm_specific = algorithm)
        })
        param_frequency <- List_to_DF(param_frequency)

      } else {
        param_frequency <- NULL
      }
    } else {
      param_frequency <- NULL
    }

    # Write statistics specific to this algorithm, which are not yet compared
    # to baseline or other algorithms
    items_save <- c("n_valid_results", "exc_inh_ratio", "n_zero_guesses",
                    "error_stats", "estimate_stats", "best_errors",
                    "best_estimates", "best_params")
    alg_qstats <- lapply(items_save, function(item) {
      List_to_DF(file_qstats, item)
    })
    names(alg_qstats) <- items_save

    alg_qstats$best_errors_toplevel <- select(best_errs_toplevel, -type,
                                              -all_of(Get_ParameterColumnNames()))
    alg_qstats$best_estimates_toplevel <- select(best_ests_toplevel, -type,
                                                 -any_of(Get_ParameterColumnNames()))
    alg_qstats$best_params_toplevel <- best_params_toplevel

    alg_qstats$param_frequency <- param_frequency

    saveRDS(alg_qstats,
            file.path(dir_analysis,
                      str_glue("quality_stats_{algorithm}_{bulk_dataset}_{granularity}.rds")))

    # Pass these variables forward for the purpose of comparing all algorithms
    return(alg_qstats)
  })
  names(alg_qstats_all) <- algorithms

  n_valid <- List_to_DF(alg_qstats_all, "n_valid_results")
  n_valid_by_norm <- n_valid %>%
    group_by(algorithm, test_data_name, granularity, normalization,
             regression_method) %>%
    summarize(n_valid = sum(n_valid),
              n_possible = sum(n_possible),
              pct_valid = n_valid / n_possible,
              .groups = "drop")

  n_valid_by_algorithm <- n_valid %>%
    group_by(algorithm) %>%
    summarize(n_valid = sum(n_valid),
              n_possible = sum(n_possible),
              pct_valid = n_valid / n_possible,
              .groups = "drop")


  # Calculate frequency of parameters ------------------------------------------
  # In the top 3 estimates per tissue. We ignore Baseline results for this
  # calculation since they aren't real estimates.

  best_errors <- List_to_DF(alg_qstats_all, "best_errors")

  # Non-algorithm-specific parameters we are interested in
  cols_keep <- c("algorithm", "reference_data_name", "normalization",
                 "regression_method", "marker_subtype", "marker_type",
                 "marker_input_type", "marker_order", "param_id")

  # Concat all algorithm params into one data frame -- ignoring Baseline
  best_params <- lapply(alg_qstats_all[setdiff(algorithms, "Baseline")],
                        "[[", "best_params") %>%
    lapply(select_at, cols_keep) %>%  # subset to only cols_keep columns
    List_to_DF()

  top3 <- best_errors %>%
    merge(best_params, by = "param_id") %>%
    Get_TopRanked(group_cols = "tissue", n_top = 3, with_mean_rank = FALSE) %>%
    # These three marker specifications are not independent of each other so we
    # combine them into one variable
    mutate(marker_algorithm = paste(marker_type, marker_subtype, marker_input_type),
           marker_algorithm = str_replace_all(marker_algorithm, " None", ""))

  # Subbing in marker_algorithm for the other 3 marker fields, removing param_id
  param_cols <- c("algorithm", "reference_data_name", "normalization",
                  "regression_method", "marker_algorithm", "marker_order")

  # Creates a data frame of parameter values vs counts in long format (columns
  # for tissue, parameter name, parameter value, and count) so the different
  # paramter dfs can be concatenated.
  param_frequency_list <- lapply(param_cols, function(p_col) {
    Count_ParamFrequency(top3, groups = c("tissue", p_col),
                         pivot_column = p_col, algorithm_specific = "None")
  })

  # Combine this list with any algorithm-specific data frames
  param_frequency <- List_to_DF(param_frequency_list) %>%
    rbind(List_to_DF(alg_qstats_all, "param_frequency"))


  # Percent of errors better than baseline -------------------------------------
  # Using top-level data

  best_errors_toplevel <- List_to_DF(alg_qstats_all, "best_errors_toplevel") %>%
    distinct()

  best_params_toplevel <- lapply(alg_qstats_all[setdiff(algorithms, "Baseline")],
                                 "[[", "best_params_toplevel") %>%
    lapply(select_at, cols_keep) %>%  # subset to only cols_keep columns
    List_to_DF()

  baseline_bests <- best_errors_toplevel %>%
    merge(alg_qstats_all$Baseline$best_params_toplevel) %>%
    subset(reference_data_name != "zeros") %>%
    group_by(tissue, normalization, regression_method) %>%
    summarize(cor = max(cor),
              rMSE = min(rMSE),
              mAPE = min(mAPE),
              .groups = "drop") %>%
    pivot_longer(cols = c(cor, rMSE, mAPE), names_to = "error_metric",
                 values_to = "baseline_value")


  # TODO this isn't the same calculation being done in Step 16 ?
  better_than_baseline <- best_errors_toplevel %>%
    merge(best_params_toplevel) %>%
    pivot_longer(cols = c(cor, rMSE, mAPE), names_to = "error_metric") %>%
    # "counts" should be in the same category as "cpm", all "log_X" normalizations
    # should be in the same category as their linear counterparts
    mutate(normalization = str_replace(normalization, "counts", "cpm"),
           normalization = str_replace(normalization, "log_", "")) %>%
    merge(baseline_bests) %>%
    group_by(error_metric) %>%
    mutate(better = if(unique(error_metric) == "cor")
      (value > baseline_value) else (value < baseline_value)) %>%
    group_by(tissue, algorithm) %>%
    dplyr::summarize(count = n(),
                     pct_better_than_baseline = sum(better) / n(),
                     .groups = "drop")


  # Best data transform per tissue

  # TODO below
  # Get ranked errors that will be used in analysis
  # best data transform
  # N markers vs error

  saveRDS(list("n_valid" = n_valid,
               "n_valid_by_norm" = n_valid_by_norm,
               "n_valid_by_algorithm" = n_valid_by_algorithm,
               "best_errors_all" = best_errors,
               "best_errors_toplevel" = best_errors_toplevel,
               "best_params" = best_params,
               "best_params_toplevel" = best_params_toplevel,
               "top3_errors" = top3,
               "parameter_frequency" = param_frequency,
               "best_baselines" = baseline_bests,
               "better_than_baseline" = better_than_baseline),
          file.path(dir_analysis,
                    str_glue("quality_stats_all_{bulk_dataset}_{granularity}.rds")))
}
