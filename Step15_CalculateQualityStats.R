# ...
# This script also calculates some statistics about the estimates in each file:
#   - how many samples in each estimate have a "bad" inhibitory:excitatory ratio
#   - the mean and SD of estimated percentages for each cell type for a given
#     sample across all estimates in the file
#   - the same mean and SD except only across the top 10-scoring estimates in
#     the file
#   - the mean of these means across all samples
#
# TODO:
# * Make plotting functions use the pre-calculated values from this script.

library(dplyr)
library(stringr)
library(parallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step15_Analysis_HelperFunctions.R"))

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

n_cores <- 12

# Which algorithms to calculate stats for
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

# Find the top errors along each error metric ----------------------------------

# Get the full set of all errors from Step 14. Subset to errors where the
# signature used to calculate the error matches the reference data used to
# generate the estimate, except for Baseline data which doesn't use reference
# data. Keep only errors for individual tissues.
best_errors_step14 <- Get_AllBestErrorsAsDf(bulk_datasets, granularity, n_cores) %>%
  subset(signature == reference_data_name | algorithm == "Baseline") %>%
  subset(tissue != "All")

# Group by all parameters common to all algorithms (ignoring input type), and
# then create a top-level grouping that also ignores which single cell reference
# was used.
group_cols <- c("tissue", Get_ParameterColumnNames())  %>%
  setdiff("reference_input_type") # Ignore input type
group_cols_toplevel <- setdiff(group_cols, "reference_data_name")

top_errors <- list(
  "all" = Get_TopErrors(best_errors_step14, group_cols, n_cores, with_mean_rank = FALSE),
  "toplevel" = Get_TopErrors(best_errors_step14, group_cols_toplevel, n_cores, with_mean_rank = FALSE)
)

# This doesn't necessarily need to be true but I want to know if it happens
stopifnot(all(top_errors$toplevel$errors$param_id %in% top_errors$all$errors$param_id))

saveRDS(top_errors, file.path(dir_analysis, str_glue("best_errors_{granularity}.rds")))


for (bulk_dataset in bulk_datasets) {
  # Get bulk metadata
  bulk <- Load_BulkData(bulk_dataset)
  bulk_metadata <- colData(bulk) %>%
    as.data.frame() %>%
    select(sample, tissue, diagnosis)
  rm(bulk)
  gc()

  alg_qstats_all <- lapply(algorithms, function(algorithm) {
    est_files_step09 <- list.files(file.path(dir_estimates, bulk_dataset, algorithm),
                                   pattern = granularity,
                                   recursive = TRUE, full.names = TRUE)

    if (length(est_files_step09) == 0) {
      message(str_glue("No data found for {bulk_dataset}/{algorithm}/{granularity}. Skipping..."))
      next
    }

    print(paste(bulk_dataset, algorithm))

    # Collect some per-file stats ----------------------------------------------

    # Results from this loop will be combined at the end to calculate top-level
    # stats.
    file_qstats <- mclapply(est_files_step09, function(est_f) {

      ## Load all error and estimate files related to est_f --------------------
      file_id <- str_replace(basename(est_f), "estimates_", "") %>%
        str_replace(".rds", "")

      est_list_step09 <- readRDS(est_f)
      file_params <- FileParams_FromParams(est_list_step09[[1]]$params)

      err_f_step11 <- list.files(file.path(dir_errors, bulk_dataset, algorithm),
                                 pattern = file_id, full.names = TRUE)

      best_est_f_step13 <- list.files(file.path(dir_top_estimates, bulk_dataset,
                                                algorithm),
                               pattern = file_id, full.names = TRUE)

      best_err_f_step14 <- list.files(file.path(dir_best_errors, bulk_dataset,
                                                algorithm),
                                      pattern = file_id, full.names = TRUE)

      # If there weren't any valid results for this file, the error file won't
      # exist. Create an abbreviated quality stats file for the sole purpose of
      # keeping track of number of failures.
      if (length(err_f_step11) != 1 || length(best_est_f_step13) != 1) {
        print(str_glue(
          paste("No valid error or best estimates files found for",
                "{basename(est_f)}. Returning abbreviated stats data.")))

        n_valid_results <- data.frame(n_valid = 0,
                                      n_possible = length(est_list_step09),
                                      file_id = file_id)
        n_valid_results <- cbind(n_valid_results, file_params)

        return(list("n_valid_results" = n_valid_results))
      }

      err_list_step11 <- readRDS(err_f_step11)
      best_est_list_step13 <- readRDS(best_est_f_step13)


      ## Valid vs possible results ---------------------------------------------
      # Using all errors from step 11

      # Must be calculated before subsetting to valid results
      n_valid_results <- data.frame(n_valid = length(err_list_step11$param_ids),
                                    n_possible = length(est_list_step09),
                                    file_id = file_id)
      n_valid_results <- cbind(n_valid_results, file_params)

      ## Subset to valid results -----------------------------------------------

      best_param_ids <- intersect(names(best_est_list_step13),
                                  top_errors$all$errors$param_id)

      if (length(best_param_ids) == 0) {
        print(str_glue(
          paste("No estimates from {basename(est_f)} were included in the best",
                "error set. Returning abbreviated stats data.")))
        return(list("n_valid_results" = n_valid_results))
      }

      # Subset to the best parameter IDs.
      best_est_list_step13 <- best_est_list_step13[best_param_ids]

      best_params <- List_to_DF(best_est_list_step13, "params") %>%
        mutate(param_id = rownames(.))

      # Get all best estimates as one data frame
      est_pcts_step13 <- lapply(best_est_list_step13, function(est_item) {
        estimates <- as.data.frame(est_item$estimates) %>%
          mutate(sample = rownames(.),
                 param_id = est_item$param_id) %>%
          merge(bulk_metadata)
      })
      est_pcts_step13 <- List_to_DF(est_pcts_step13)
      rownames(est_pcts_step13) <- NULL

      # Create the same kind of duplication in the estimates df for estimate
      # stats calculations, using only the top errors
      weights <- top_errors$all$ranks %>%
        subset(param_id %in% best_param_ids) %>%
        select(param_id, tissue, type)

      est_pcts_weighted_step_13 <- merge(est_pcts_step13, weights) %>%
        #subset(type != "best_mean") %>% # TODO should we exclude?
        select(-type)

      # Use non-duplicated data for exc:inh ratio and number of zeros
      est_pcts_step13 <- distinct(est_pcts_weighted_step_13)


      ## Ratio of excitatory to inhibitory cells -------------------------------
      # For best estimates

      # This cbind works for both broad and sub classes
      neuron_ests <- cbind(
        select(est_pcts_step13, sample, tissue, param_id),
        data.frame(
          Excitatory = rowSums(select(est_pcts_step13, starts_with("Exc"))),
          Inhibitory = rowSums(select(est_pcts_step13, starts_with("Inh")))
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

      num_zeros <- est_pcts_step13 %>%
        group_by(tissue, param_id) %>%
        summarize(across(where(is.numeric), ~ sum(.x == 0)), .groups = "drop")

      return(list("n_valid_results" = n_valid_results,
                  "exc_inh_ratio" = exc_inh_ratio,
                  "n_zero_guesses" = num_zeros,
                  "best_estimates" = est_pcts_step13,
                  "best_params" = best_params))
    }, mc.cores = n_cores)


    # Top estimates for the algorithm ------------------------------------------

    best_ests <- List_to_DF(file_qstats, "best_estimates")
    best_params <- List_to_DF(file_qstats, "best_params")

    ## Estimate stats ----------------------------------------------
    # For best estimates

    # Create the same kind of duplication in the estimates df as we used to
    # calculate best errors, so that each estimate is weighted by the number of
    # times it shows up as a "best" among the 4 error metrics. Also pull in
    # parameters listed in group_cols for this calculation.
    weights_all <- top_errors$all$ranks %>%
      merge(top_errors$all$errors) %>%
      select(param_id, all_of(group_cols), type)

    weights_toplevel <- top_errors$toplevel$ranks %>%
      merge(top_errors$toplevel$errors) %>%
      select(param_id, all_of(group_cols), type)

    est_stats_all <- best_ests %>%
      merge(weights_all) %>%
      select(param_id, sample, all_of(group_cols), where(is.numeric)) %>%
      Calculate_EstimateStats(c("sample", group_cols))

    best_ests_toplevel <- merge(best_ests, weights_toplevel)

    est_stats_toplevel <- best_ests_toplevel %>%
      select(param_id, sample, all_of(group_cols_toplevel), where(is.numeric)) %>%
      Calculate_EstimateStats(c("sample", group_cols_toplevel))

    best_params_toplevel <- subset(best_params,
                                   param_id %in% top_errors$toplevel$errors$param_id)

    if (algorithm != "Baseline") {
      # Check for algorithm-specific parameters and if they exist calculate how
      # many times each value for those parameters shows up in the top 3 errors
      # for each error metric
      alg_specific <- select(best_params,
                             -all_of(Get_ParameterColumnNames()), # remove file params
                             -contains("filter_level"), -contains("marker"), # remove general marker params
                             -mode)

      if (ncol(alg_specific) > 1) {
        top3 <- top_errors$all$errors %>%
          subset(test_data_name == bulk_dataset &
                   algorithm == unique(best_params$algorithm)) %>%
          Get_TopRanked("tissue", n_top = 3, with_mean_rank = TRUE)
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
                    "best_estimates", "best_params")
    alg_qstats <- lapply(items_save, function(item) {
      List_to_DF(file_qstats, item)
    })
    names(alg_qstats) <- items_save

    alg_qstats$best_estimates_toplevel <- best_ests_toplevel %>%
      select(-type, -any_of(Get_ParameterColumnNames())) %>%
      distinct()

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
    group_by(algorithm, granularity) %>%
    summarize(n_valid = sum(n_valid),
              n_possible = sum(n_possible),
              pct_valid = n_valid / n_possible,
              .groups = "drop")


  # Calculate frequency of parameters ------------------------------------------
  # In the top 3 estimates per tissue. We ignore Baseline results for this
  # calculation since they aren't real estimates.

  # Non-algorithm-specific parameters we are interested in
  cols_keep <- c("algorithm", "reference_data_name", "normalization",
                 "regression_method", "marker_subtype", "marker_type",
                 "marker_input_type", "marker_order", "param_id")

  # Concat all algorithm params into one data frame -- ignoring Baseline
  best_params <- lapply(alg_qstats_all[setdiff(algorithms, "Baseline")],
                        "[[", "best_params") %>%
    lapply(select_at, cols_keep) %>%  # subset to only cols_keep columns
    List_to_DF()

  top3 <- top_errors$all$errors %>%
    merge(best_params) %>%
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

  best_params_toplevel <- lapply(alg_qstats_all[setdiff(algorithms, "Baseline")],
                                 "[[", "best_params_toplevel") %>%
    lapply(select_at, cols_keep) %>%  # subset to only cols_keep columns
    List_to_DF()

  baseline_bests <- top_errors$toplevel$errors %>%
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
  better_than_baseline <- top_errors$toplevel$errors %>%
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


  # Significance calculations --------------------------------------------------

  best_errors_tmp <- subset(top_errors$all$errors, test_data_name == bulk_dataset) %>%
    merge(top_errors$all$ranks)
  best_errors_tmp$avg_id <- unlist(apply(best_errors_tmp, 1, function(row) {
    paste(row[group_cols], collapse = "_")
  }))

  best_errors_toplevel_tmp <- subset(top_errors$toplevel$errors, test_data_name == bulk_dataset) %>%
    merge(top_errors$toplevel$ranks)
  best_errors_toplevel_tmp$avg_id <- unlist(apply(best_errors_toplevel_tmp, 1, function(row) {
    paste(row[group_cols_toplevel], collapse = "_")
  }))

  best_estimates <- List_to_DF(alg_qstats_all, "best_estimates")
  best_estimates_toplevel <- List_to_DF(alg_qstats_all, "best_estimates_toplevel")

  # Average the estimates corresponding to best correlation, best rMSE, and best
  # mAPE together for each data input type
  avg_list <- Create_AveragesList(best_errors_tmp,
                                  best_estimates,
                                  group_cols,
                                  n_cores,
                                  with_mean_rank = FALSE)
  avg_list_toplevel <- Create_AveragesList(best_errors_toplevel_tmp,
                                           best_estimates_toplevel,
                                           group_cols_toplevel,
                                           n_cores,
                                           with_mean_rank = FALSE)

  # Calculate significance of cell type differences on a tissue-by-tissue basis
  mean_props_all <- mclapply(avg_list,
                             Get_MeanProps_Significance,
                             group_cols = group_cols,
                             mc.cores = n_cores)
  mean_props_toplevel <- mclapply(avg_list_toplevel,
                                  Get_MeanProps_Significance,
                                  group_cols = group_cols_toplevel,
                                  mc.cores = n_cores)

  # Best data transform per tissue

  # TODO below
  # best data transform
  # N markers vs error

  saveRDS(list("n_valid" = n_valid,
               "n_valid_by_norm" = n_valid_by_norm,
               "n_valid_by_algorithm" = n_valid_by_algorithm,
               "best_params" = best_params,
               "best_params_toplevel" = best_params_toplevel,
               "top3_errors" = top3,
               "parameter_frequency" = param_frequency,
               "best_baselines" = baseline_bests,
               "better_than_baseline" = better_than_baseline),
          file.path(dir_analysis,
                    str_glue("quality_stats_all_{bulk_dataset}_{granularity}.rds")))
}
