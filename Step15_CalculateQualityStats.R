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

granularities <- c("broad_class", "sub_class")
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
best_errors_step14 <- Get_AllBestErrorsAsDf(bulk_datasets, granularities, n_cores) %>%
  subset(signature == reference_data_name | algorithm == "Baseline") %>%
  subset(tissue != "All")

top3_by_tissue <- best_errors_step14 %>%
  subset(algorithm != "Baseline") %>%
  Get_TopRanked(group_cols = c("tissue", "granularity"),
                n_top = 3,
                with_mean_rank = TRUE)

top3_by_algorithm <- best_errors_step14 %>%
  subset(algorithm != "Baseline") %>%
  Get_TopRanked(group_cols = c("tissue", "granularity", "algorithm"),
                n_top = 3,
                with_mean_rank = TRUE)

best_dt <- Get_BestDataTransform(Standardize_DataTransform(top3_by_tissue),
                                 algorithms)

# Group by all parameters common to all algorithms (ignoring input type), and
# then create a top-level grouping that also ignores which single cell reference
# was used.
group_cols <- c("tissue", Get_ParameterColumnNames())  %>%
  setdiff("reference_input_type") # Ignore input type
group_cols_toplevel <- setdiff(group_cols, "reference_data_name")

top_errors <- list(
  "by_reference" = Get_TopErrors(best_errors_step14, group_cols, n_cores, with_mean_rank = FALSE),
  "by_algorithm" = Get_TopErrors(best_errors_step14, group_cols_toplevel, n_cores, with_mean_rank = FALSE),
  "top3_by_tissue" = top3_by_tissue,
  "top3_by_algorithm" = top3_by_algorithm,
  "best_data_transform" = best_dt
)

top_errors$param_ids <- unique(c(top_errors$by_reference$param_ids,
                                 top_errors$by_algorithm$param_ids,
                                 top3_by_tissue$param_id,
                                 top3_by_algorithm$param_id))

# This doesn't necessarily need to be true but I want to know if it happens
stopifnot(all(top_errors$by_algorithm$param_ids %in% top_errors$by_reference$param_ids))

saveRDS(top_errors, file.path(dir_analysis, "top_errors.rds"))

iter_vars <- expand.grid(bulk_dataset = bulk_datasets,
                         algorithm = algorithms,
                         granularity = granularities,
                         stringsAsFactors = FALSE)

combined_metadata <- lapply(bulk_datasets, Get_BulkMetadata,
                            columns = c("sample", "tissue", "diagnosis")) %>%
  List_to_DF()

qstats_all <- lapply(1:nrow(iter_vars), function(iter_row) {
  bulk_dataset <- iter_vars$bulk_dataset[iter_row]
  algorithm <- iter_vars$algorithm[iter_row]
  granularity <- iter_vars$granularity[iter_row]

  est_files_step09 <- list.files(file.path(dir_estimates, bulk_dataset, algorithm),
                                 pattern = granularity,
                                 recursive = TRUE, full.names = TRUE)

  if (length(est_files_step09) == 0) {
    message(str_glue("No data found for {bulk_dataset} / {algorithm} / ",
                     "{granularity}. Skipping..."))
    next
  }

  print(paste(bulk_dataset, algorithm, granularity))

  # Collect some per-file stats ----------------------------------------------

  # Results from this loop will be combined at the end to calculate top-level
  # stats.
  file_qstats <- mclapply(est_files_step09, function(est_f) {

    ## Load all error and estimate files related to est_f --------------------
    file_id <- str_replace(basename(est_f), "estimates_", "") %>%
      str_replace(".rds", "")

    est_list_step09 <- readRDS(est_f)
    file_params <- FileParams_FromParams(est_list_step09[[1]]$params)

    err_f_step11 <- Find_ErrorFiles(bulk_dataset, algorithm, file_id)
    best_est_f_step13 <- Find_BestEstimateFiles(bulk_dataset, algorithm, file_id)

    # If there weren't any valid results for this file, the error file won't
    # exist. Create an abbreviated quality stats file for the sole purpose of
    # keeping track of number of failures.
    if (length(err_f_step11) != 1 || length(best_est_f_step13) != 1) {
      print(str_glue("No valid error or best estimates files found for ",
                     "{basename(est_f)}. Returning abbreviated stats data."))

      n_valid_results <- data.frame(n_valid = 0,
                                    n_possible = length(est_list_step09))
      n_valid_results <- cbind(n_valid_results, file_params)

      return(list("n_valid_results" = n_valid_results))
    }

    err_list_step11 <- readRDS(err_f_step11)
    best_est_list_step13 <- readRDS(best_est_f_step13)


    ## Valid vs possible results ---------------------------------------------
    # Using all errors from step 11

    # Must be calculated before subsetting to valid results
    n_valid_results <- data.frame(n_valid = length(err_list_step11$param_ids),
                                  n_possible = length(est_list_step09))
    n_valid_results <- cbind(n_valid_results, file_params)


    ## Subset to results from top errors -------------------------------------

    best_param_ids <- intersect(names(best_est_list_step13),
                                top_errors$param_ids)

    if (length(best_param_ids) == 0) {
      print(str_glue("No estimates from {basename(est_f)} were included in ",
                     "the best error set. Returning abbreviated stats data."))
      return(list("n_valid_results" = n_valid_results))
    }

    best_params <- List_to_DF(best_est_list_step13, "params") %>%
      mutate(param_id = rownames(.)) %>%
      subset(param_id %in% best_param_ids) %>%
      tibble::remove_rownames()

    est_pcts_step13 <- Subset_BestEstimates(best_param_ids,
                                            best_est_list_step13,
                                            combined_metadata)


    ## Number of 0 guesses for each cell type --------------------------------
    # For best estimates

    num_zeros <- est_pcts_step13 %>%
      group_by(tissue, param_id) %>%
      summarize(across(where(is.numeric), ~ sum(.x == 0)),
                .groups = "drop") %>%
      as.data.frame()

    return(list("n_valid_results" = n_valid_results,
                "n_zero_guesses" = num_zeros,
                "best_estimates" = est_pcts_step13,
                "best_params" = best_params))
  }, mc.cores = n_cores)


  # Top estimates for the algorithm ------------------------------------------

  best_ests_all <- List_to_DF(file_qstats, "best_estimates")
  best_params_all <- List_to_DF(file_qstats, "best_params")

  if (is.null(best_ests_all) || is.null(best_params_all)) {
    print(str_glue("No valid data found for {bulk_dataset} / {algorithm} / ",
                   "{granularity}. Returning abbreviated stats data."))
    # Add some identifying information
    info <- data.frame("bulk_dataset" = bulk_dataset,
                       "algorithm" = algorithm,
                       "granularity" = granularity)
    return(list("n_valid_results" = List_to_DF(file_qstats, "n_valid_results"),
                "info" = info))
  }

  ## Ratio of excitatory to inhibitory cells -----------------------------------
  # For best estimates

  # Limit to combinations of parameters/tissues that appear in the top-level
  # best errors list
  valid <- top_errors$by_algorithm$errors %>%
    subset(param_id %in% best_ests_all$param_id) %>%
    select(param_id, tissue) %>%
    distinct()

  ests_filt <- merge(best_ests_all, valid)
  exc_inh_ratio <- Calculate_ExcInhRatio(ests_filt, best_params_all)


  ## Parameter stats -----------------------------------------------------------
  # For best estimates

  # Count parameters that show up the most in the top 3 errors for this algorithm
  if (algorithm != "Baseline") {
    best_params_top3 <- subset(best_params_all,
                               param_id %in% top_errors$top3_by_algorithm$param_id)

    param_frequency <- Count_AlgSpecificParameters(best_params_top3,
                                                   top_errors$top3_by_algorithm)
  } else {
    param_frequency <- NULL
  }

  # Flatten the file_qstats list from a list of lists of data frames to a list
  # of data frames
  items_save <- c("n_valid_results", "n_zero_guesses", "best_estimates",
                  "best_params")
  alg_qstats <- lapply(items_save, function(item) {
    List_to_DF(file_qstats, item)
  })
  names(alg_qstats) <- items_save

  # Add the new variables calculated here
  alg_qstats$param_frequency <- param_frequency
  alg_qstats$exc_inh_ratio <- exc_inh_ratio

  # Add some identifying information
  alg_qstats$info <- data.frame("bulk_dataset" = bulk_dataset,
                                "algorithm" = algorithm,
                                "granularity" = granularity)

  saveRDS(alg_qstats,
          file.path(dir_analysis,
                    str_glue("quality_stats_{algorithm}_{bulk_dataset}_{granularity}.rds")))

  # Pass these variables forward for the purpose of comparing all algorithms
  return(alg_qstats)
})

# Top-level stats for all data -------------------------------------------------

qstats_info <- List_to_DF(qstats_all, "info")

n_valid <- List_to_DF(qstats_all, "n_valid_results") %>%
  subset(algorithm != "Baseline")

n_valid_by_norm <- n_valid %>%
  group_by(algorithm, granularity, normalization,
           regression_method) %>%
  summarize(n_valid = sum(n_valid),
            n_possible = sum(n_possible),
            pct_valid = n_valid / n_possible,
            .groups = "drop") %>%
  as.data.frame()

n_valid_by_algorithm <- n_valid %>%
  group_by(algorithm, granularity) %>%
  summarize(n_valid = sum(n_valid),
            n_possible = sum(n_possible),
            pct_valid = n_valid / n_possible,
            .groups = "drop") %>%
  as.data.frame()

n_zeros <- lapply(qstats_all, function(qstat) {
  if ("n_zero_guesses" %in% names(qstat)) {
    qstat$n_zero_guesses %>%
      subset(param_id %in% top_errors$by_algorithm$param_ids) %>%
      merge(select(top_errors$by_algorithm$errors, tissue, param_id, algorithm)) %>%
      distinct() %>%
      pivot_longer(cols = where(is.numeric),
                   names_to = "celltype",
                   values_to = "count")
  } else {
    NULL
  }
}) %>%
  List_to_DF() %>%
  as.data.frame()

exc_inh_ratio <- List_to_DF(qstats_all, "exc_inh_ratio") %>%
  subset(algorithm != "Baseline") %>%
  # Best data transforms only
  Standardize_DataTransform() %>%
  merge(top_errors$best_data_transform) %>%
  # Weight by number of apperances
  merge(top_errors$by_algorithm$ranks) %>%
  group_by(tissue, algorithm, granularity) %>%
  summarize(median_exc_inh_ratio = median(mean_exc_inh_ratio),
            mean_exc_inh_ratio = mean(mean_exc_inh_ratio),
            .groups = "drop")


# Calculate frequency of parameters --------------------------------------------
# In the top 3 estimates per tissue. We ignore Baseline results for this
# calculation since they aren't real estimates.

# Non-algorithm-specific parameters we are interested in
cols_keep <- c("algorithm", "reference_data_name", "granularity",
               "normalization", "regression_method", "marker_subtype",
               "marker_type", "marker_input_type", "marker_order", "param_id")

# Concat all algorithm params into one data frame. Requires two steps because
# some algorithms don't have data for a specific dataset and the NULLs need
# to be filtered out before calling select()
best_params <- lapply(qstats_all, "[[", "best_params")

# Remove any algorithms with NULL entries, which can happen if no output was
# valid. Also remove baseline data.
ok <- (lengths(best_params) > 0) & (qstats_info$algorithm != "Baseline")

best_params <- best_params[ok] %>%
  lapply(select_at, cols_keep) %>%  # subset to only cols_keep columns
  List_to_DF() %>%
  subset(algorithm != "Baseline")

# Limit to only parameter sets that appear in the top 3 errors
params_top3 <- top_errors$top3_by_tissue %>%
  merge(best_params) %>%
  # These three marker specifications are not independent of each other so we
  # combine them into one variable
  mutate(marker_algorithm = paste(marker_type, marker_subtype, marker_input_type),
         marker_algorithm = str_replace_all(marker_algorithm, " None", ""))

# Columns to count
param_cols <- c("algorithm", "reference_data_name", "normalization",
                "regression_method", "marker_algorithm", "marker_order")

# Creates a data frame of parameter values vs counts in long format (columns
# for tissue, parameter name, parameter value, and count) so the different
# paramter dfs can be concatenated.
param_frequency_list <- lapply(param_cols, function(p_col) {
  Count_ParamFrequency(params_top3, groups = c("tissue", "granularity", p_col),
                       pivot_column = p_col, algorithm_specific = "None")
})

# Combine this list with any algorithm-specific data frames
param_frequency <- List_to_DF(param_frequency_list) %>%
  rbind(List_to_DF(qstats_all, "param_frequency")) %>%
  as.data.frame()


# Percent of errors better than baseline ---------------------------------------
# Using top-level data

baseline_bests <- top_errors$by_algorithm$errors %>%
  subset(algorithm == "Baseline" & reference_data_name != "zeros") %>%
  Standardize_DataTransform() %>%
  group_by(tissue, granularity, test_data_name, normalization, regression_method) %>%
  summarize(cor = max(cor),
            rMSE = min(rMSE),
            mAPE = min(mAPE),
            .groups = "drop") %>%
  pivot_longer(cols = c(cor, rMSE, mAPE), names_to = "error_metric",
               values_to = "baseline_value")

# Sub class is missing some algorithms, this fills them in
alg_info <- expand.grid(tissue = unique(best_errors_step14$tissue),
                        algorithm = unique(best_errors_step14$algorithm),
                        granularity = unique(best_errors_step14$granularity),
                        error_metric = unique(baseline_bests$error_metric)) %>%
  subset(algorithm != "Baseline")

better_than_baseline <- top_errors$by_algorithm$errors %>%
  subset(algorithm != "Baseline") %>%
  Standardize_DataTransform() %>%
  pivot_longer(cols = c(cor, rMSE, mAPE), names_to = "error_metric",
               values_to = "value") %>%
  merge(baseline_bests) %>%
  distinct() %>%
  merge(top_errors$best_data_transform) %>%
  merge(alg_info, all = TRUE) %>%
  mutate(better = case_when(is.na(value) ~ FALSE,
                            error_metric == "cor" ~ value > baseline_value,
                            TRUE ~ value < baseline_value))

better_than_baseline_by_tissue <- better_than_baseline %>%
  group_by(tissue, granularity, test_data_name, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / count,
                   .groups = "drop") %>%
  as.data.frame()

better_than_baseline_by_algorithm <- better_than_baseline %>%
  group_by(granularity, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / count,
                   .groups = "drop") %>%
  as.data.frame()


# Significance calculations --------------------------------------------------
# For top-level errors by algorithm

best_errors_by_alg <- merge(top_errors$by_algorithm$errors,
                            top_errors$by_algorithm$ranks) %>%
  subset(reference_data_name != "zeros") %>%
  Standardize_DataTransform()

best_errors_by_alg$avg_id <- unlist(apply(best_errors_by_alg, 1, function(row) {
  paste(row[group_cols_toplevel], collapse = "_")
}))

valid <- best_errors_by_alg %>%
  select(param_id, tissue) %>%
  distinct()

# Separate by granularity because of differing numbers of cell types
sig_stats <- lapply(granularities, function(granularity) {
  best_errors_tmp <- best_errors_by_alg[best_errors_by_alg$granularity == granularity, ]

  best_estimates_all <- lapply(qstats_all, function(qstat) {
    if (("best_estimates" %in% names(qstat)) &&
        (qstat$info$granularity == granularity)) {
      ests <- qstat$best_estimates %>%
        subset(param_id %in% best_errors_tmp$param_id) %>%
        merge(valid)
      return(ests)
    } else {
      return(NULL)
    }
  }) %>%
    List_to_DF()

  # Average the estimates corresponding to best correlation, best rMSE, and best
  # mAPE together for each data input type
  avg_list <- Create_AveragesList(best_errors_tmp,
                                  best_estimates_all,
                                  group_cols_toplevel,
                                  n_cores,
                                  with_mean_rank = FALSE)

  avg_ests <- List_to_DF(avg_list, "avg_estimates")

  significance <- Get_MeanProps_Significance(avg_list,
                                              group_cols_toplevel,
                                              n_cores = n_cores)

  return(list("significance" = significance, "avg_ests" = avg_ests))
})


# TODO below
# N markers vs error

results <- list("n_valid_results" = n_valid,
                "n_valid_by_norm" = n_valid_by_norm,
                "n_valid_by_algorithm" = n_valid_by_algorithm,
                "exc_inh_ratio" = exc_inh_ratio,
                "n_zero_guesses" = n_zeros,
                "parameter_frequency" = param_frequency,
                "better_than_baseline_by_tissue" = better_than_baseline_by_tissue,
                "better_than_baseline_by_algorithm" = better_than_baseline_by_algorithm,
                "significance" = List_to_DF(sig_stats, "significance"),
                "avg_estimates_by_algorithm" = List_to_DF(sig_stats, "avg_ests"))

saveRDS(results, file.path(dir_analysis, str_glue("quality_stats_all.rds")))
