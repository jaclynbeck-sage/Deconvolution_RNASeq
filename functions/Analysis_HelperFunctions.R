library(parallel)
source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

# Loads all the "best_errors" files matching specific parameters into one big
# list.
#
# Arguments:
#   bulk_datasets - a vector of bulk dataset names to load ("Mayo", "MSBB",
#                   and/or "ROSMAP")
#   granularities - a vector of granularities to load ("broad_class" and/or
#                   "sub_class")
#
# Returns:
#   a data frame with all errors from all files concatenated together
Get_AllBestErrorsAsDf <- function(bulk_datasets, granularities, n_cores = 2) {
  # List of "best errors" files matching the bulk data set names and granularity
  file_list <- list.files(dir_best_errors,
                          pattern = paste0("(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*",
                                           "(",
                                           paste(granularities, collapse = "|"),
                                           ")"),
                          full.names = TRUE,
                          recursive = TRUE)

  best_err_list <- mclapply(file_list, function(file) {
    data <- readRDS(file)

    # Merge a subset of the parameters into the mean errors data frame. These
    # parameters exist in every error file for every algorithm
    errs_df <- cbind(data$means, FileParams_FromParams(data$params))
    return(errs_df)
  }, mc.cores = n_cores)

  return(do.call(rbind, best_err_list))
}


Get_AllBestEstimatesAsDf <- function(bulk_datasets, granularity, metadata,
                                     best_params = NULL, n_cores = 2) {
  # List of "best estimates" files matching the bulk data set names and granularity
  file_list <- list.files(dir_top_estimates,
                          pattern = paste0("(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*", granularity),
                          full.names = TRUE,
                          recursive = TRUE)

  best_ests_list <- mclapply(file_list, function(file) {
    data <- readRDS(file)
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
        merge(metadata, by = "sample", all.y = FALSE)

      new_ests$param_id <- N
      new_ests$algorithm <- data[[N]]$params$algorithm

      data[[N]]$estimates <- new_ests
    }

    ests_df <- do.call(rbind, lapply(data, "[[", "estimates"))
    params <- do.call(rbind, lapply(data, "[[", "params"))

    return(list("estimates" = ests_df,
                "params" = params))
  }, mc.cores = n_cores)

  return(do.call(rbind, lapply(best_ests_list, "[[", "estimates")))
}


# Extracts general quality stats from the top parameters objects. Currently,
# we're only interested in n_valid_results and n_possible_results, though this
# function could be expanded.
Get_AllQualityStatsAsDf <- function(bulk_datasets, granularity, n_cores = 2) {
  file_list <- list.files(dir_analysis,
                          pattern = paste0("quality_stats_.*(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*", granularity),
                          full.names = TRUE,
                          recursive = TRUE)

  file_list_per_alg <- file_list[!grepl("all", file_list)]
  file_list_all <- setdiff(file_list, file_list_per_alg)

  qstats_list <- mclapply(file_list_per_alg, function(file) {
    data <- readRDS(file)

    file_id <- str_replace(basename(file), "quality_stats_", "") %>%
      str_replace(".rds", "")

    # TODO

    return(n_valid)
  }, mc.cores = n_cores)

  return(do.call(rbind, qstats_list))
}


Get_BulkMetadata <- function(bulk_dataset, columns) {
  bulk <- Load_BulkData(bulk_dataset)
  bulk_metadata <- colData(bulk) %>%
    as.data.frame() %>%
    select(all_of(columns))
  return(bulk_metadata)
}


Find_ErrorFiles <- function(bulk_dataset, algorithm, file_id) {
  list.files(file.path(dir_errors, bulk_dataset, algorithm),
             pattern = file_id, full.names = TRUE)
}


Find_BestEstimateFiles <- function(bulk_dataset, algorithm, file_id) {
  list.files(file.path(dir_top_estimates, bulk_dataset, algorithm),
             pattern = file_id, full.names = TRUE)
}


Find_BestSignature <- function(errs_df) {
  ranks <- errs_df %>%
    subset(algorithm != "Baseline") %>%
    Rank_Errors(group_cols = c("tissue", "data_transform")) %>%
    Get_TopRank(n_top = 1) %>%
    subset(cor_rank == 1 | rMSE_rank == 1 | mAPE_rank == 1) # drop top mean_rank

  get_best_signature <- function(signature) {
    tmp <- table(signature)
    names(tmp)[which.max(tmp)]
  }

  best_signatures <- ranks_sub %>%
    group_by(tissue, data_transform) %>%
    dplyr::summarize(best_signature = get_best_signature(signature),
                     .groups = "drop")

  return(best_signatures)
}


Standardize_DataTransform <- function(data) {
  data %>%
    mutate(normalization = str_replace(normalization, "counts", "cpm"),
           normalization = str_replace(normalization, "log_", ""),
           data_transform = paste(normalization, regression_method, sep = " + "))
}


Get_BestDataTransform <- function(ranked_df, algorithms) {
  params <- ranked_df %>%
    select(tissue, normalization, regression_method, data_transform) %>%
    distinct()

  best_dt <- ranked_df %>%
    group_by(tissue, data_transform) %>%
    dplyr::summarize(count = n(),
                     mean_rank = mean(mean_rank),
                     .groups = "drop_last") %>%
    # First, pick the data transform(s) that show up the most in the top3
    slice_max(order_by = count, with_ties = TRUE) %>%
    # If there's a tie, use the transform with the lowest mean rank
    slice_min(order_by = mean_rank, with_ties = FALSE) %>%
    merge(params)

  best_dt <- tidyr::expand_grid(best_dt,
                                algorithm = algorithms)

  # MuSiC always has to use CPM (counts). CibersortX uses CPM when the
  # normalization is TMM since TMM isn't a valid normalization in CibersortX.
  best_dt$normalization[best_dt$algorithm == "Music"] <- "cpm"
  best_dt$normalization[best_dt$normalization == "tmm" &
                          best_dt$algorithm == "CibersortX"] <- "cpm"

  best_dt <- best_dt %>%
    # Fix the data_transform field for MuSiC and CibersortX
    mutate(data_transform = paste(normalization, "+", regression_method)) %>%
    select(-count, -mean_rank) %>%
    as.data.frame()

  return(best_dt)
}


Rank_Errors <- function(errs_df, group_cols) {
  ranks <- errs_df %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::mutate(cor_rank = rank(-cor),
                  rMSE_rank = rank(rMSE),
                  mAPE_rank = rank(mAPE)) %>%
    dplyr::ungroup()
  ranks$mean_rank <- rowMeans(ranks[, c("cor_rank", "rMSE_rank", "mAPE_rank")])
  return(ranks)
}


Get_TopRanked <- function(df, group_cols, n_top = 1, with_mean_rank = TRUE) {
  df_ranked <- Rank_Errors(df, group_cols) %>%
    group_by_at(group_cols)

  top_ranked <- do.call(rbind, list(
    mutate(dplyr::slice_min(df_ranked, order_by = cor_rank, n = n_top, with_ties = FALSE),
           type = "best_cor"),
    mutate(dplyr::slice_min(df_ranked, order_by = rMSE_rank, n = n_top, with_ties = FALSE),
           type = "best_rMSE"),
    mutate(dplyr::slice_min(df_ranked, order_by = mAPE_rank, n = n_top, with_ties = FALSE),
           type = "best_mAPE")
  ))

  if (with_mean_rank) {
    top_ranked <- rbind(
      top_ranked,
      mutate(dplyr::slice_min(df_ranked, order_by = mean_rank, n = n_top, with_ties = FALSE),
             type = "best_mean"))
  }

  top_ranked %>%
    ungroup() %>%
    as.data.frame()
}


Get_TopErrors <- function(errors_df, group_cols, n_cores, with_mean_rank = TRUE) {
  # Special case: The signature doesn't matter when all percentages are zeros, so
  # for the "zeros" baseline data we just subset to have one unique param_id per
  # tissue/data transform. The zeros data also needs to be held out separately
  # from the other baseline data and not combined with it at the top level.
  best_zeros <- subset(errors_df, reference_data_name == "zeros") %>%
    mutate(signature = NA) %>%
    distinct() %>%
    Get_TopRanked(group_cols, n_top = 1, with_mean_rank = with_mean_rank)

  # Remove the "zeros" data from the best errors df for ranking
  errors_df <- subset(errors_df, reference_data_name != "zeros")

  top_errors <- errors_df %>%
    Get_TopRanked(group_cols, n_top = 1, with_mean_rank = with_mean_rank) %>%
    rbind(best_zeros)

  error_stats <- top_errors %>%
    Calculate_ErrorStats(group_cols)

  # Break into smaller data frames for storage. These dfs currently have one "best"
  # param ID per error metric, which results in a lot of duplicate rows with the
  # same data where the param ID is the "best" for multiple error metrics. We
  # save the info for mapping param ID -> best metric as one data frame, and the actual
  # errors for each param ID/tissue in another
  error_ranks <- top_errors %>%
    select(param_id, tissue, signature, ends_with("rank"), type)

  top_errors <- top_errors %>%
    select(-type, -ends_with("rank")) %>%
    distinct()

  return(list(
    "ranks" = error_ranks,
    "errors" = top_errors,
    "stats" = error_stats,
    "param_ids" = top_errors$param_id
  ))
}


Find_BestParameters <- function(errs_df, group_cols, with_mean_rank = TRUE) {
  ranks <- Get_TopRanked(errs_df, group_cols, n_top = 1,
                         with_mean_rank = with_mean_rank)

  # Some parameter IDs will be duplicated across multiple error metrics, this
  # creates a data frame with the list of unique parameter IDs associated with
  # each tissue.
  best_params <- ranks %>%
    dplyr::select(tissue, param_id) %>%
    dplyr::distinct()

  return(best_params)
}


Calculate_ErrorStats <- function(errs_df, group_cols) {
  err_stats <- errs_df %>%
    melt(variable.name = "metric",
         id.vars = c("param_id", group_cols),
         measure.vars = c("cor", "rMSE", "mAPE")) %>%
    group_by_at(c(group_cols, "metric")) %>%
    summarize(mean_err = mean(value),
              sd_err = sd(value),
              rel_sd_err = sd_err / mean_err,
              .groups = "drop")

  return(err_stats)
}


Calculate_EstimateStats <- function(best_ests, top_errors_list, group_cols) {
  # Create the same kind of duplication in the estimates df as we used to
  # calculate best errors, so that each estimate is weighted by the number of
  # times it shows up as a "best" among the 4 error metrics. Also pull in
  # parameters listed in group_cols for this calculation.
  weights <- top_errors_list$ranks %>%
    merge(top_errors_list$errors) %>%
    select(param_id, all_of(group_cols), type)

  est_stats <- best_ests %>%
    merge(weights) %>%
    select(param_id, sample, all_of(group_cols), where(is.numeric)) %>%
    melt(variable.name = "celltype",
         value.name = "percent",
         id.vars = c("param_id", "sample", group_cols)) %>%
    group_by_at(c("sample", group_cols, "celltype")) %>%
    summarize(mean_pct = mean(percent),
              sd_pct = sd(percent),
              rel_sd_pct = sd_pct / mean_pct,
              .groups = "drop")

  return(est_stats)
}


Subset_BestEstimates <- function(best_param_ids, best_est_list, bulk_metadata) {
  # Subset to the best parameter IDs.
  best_est_list <- best_est_list[best_param_ids]

  # Get all best estimates as one data frame
  est_pcts <- lapply(best_est_list, function(est_item) {
    estimates <- as.data.frame(est_item$estimates) %>%
      mutate(sample = rownames(.),
             param_id = est_item$param_id) %>%
      merge(bulk_metadata)
  })
  est_pcts <- List_to_DF(est_pcts) %>%
    tibble::remove_rownames()

  return(est_pcts)
}


Calculate_ExcInhRatio <- function(est_pcts, params) {
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
                     mean_exc_inh_ratio = mean(exc_inh_ratio, na.rm = TRUE),
                     median_exc_inh_ratio = median(exc_inh_ratio, na.rm = TRUE),
                     param_id = unique(param_id),
                     .groups = "drop") %>%
    select(-n_bad, -count) %>%
    as.data.frame()

  params <- params %>%
    select(param_id, all_of(Get_ParameterColumnNames())) %>%
    subset(param_id %in% exc_inh_ratio$param_id)

  return(merge(exc_inh_ratio, params))
}


Count_AlgSpecificParameters <- function(best_params, top3_errors_df) {
  alg <- unique(best_params$algorithm)
  gran <- unique(best_params$granularity)

  # Check for algorithm-specific parameters and if they exist calculate how
  # many times each value for those parameters shows up in the top 3 errors
  # for each error metric
  alg_specific <- select(best_params,
                         # remove file params
                         -all_of(Get_ParameterColumnNames()),
                         # remove general marker params
                         -contains("filter_level"), -contains("marker"),
                         # remove unneeded "mode" indicator
                         -mode)

  if (ncol(alg_specific) > 1) {
    top3 <- merge(top3_errors_df, alg_specific)

    params_test <- setdiff(colnames(alg_specific), "param_id")
    param_frequency <- lapply(params_test, function(param_col) {
      Count_ParamFrequency(top3, groups = c("tissue", param_col),
                           pivot_column = param_col,
                           algorithm_specific = alg)
    })
    param_frequency <- List_to_DF(param_frequency) %>%
      mutate(granularity = gran)

  } else {
    param_frequency <- NULL
  }

  return(param_frequency)
}


Count_ParamFrequency <- function(df, groups, pivot_column, algorithm_specific = FALSE) {
  df %>%
    group_by_at(groups) %>%
    dplyr::summarize(count = n(), .groups = "drop") %>%
    pivot_longer(cols = all_of(pivot_column), names_to = "parameter_name",
                 values_to = "parameter_value") %>%
    mutate(algorithm_specific = algorithm_specific)
}


Get_AverageStats <- function(avg_id, errs_df, ests_df, group_cols, with_mean_rank = TRUE) {
  # There can be up to 4 parameter sets per combination of data input
  # parameters, but sometimes a single parameter set was the best for multiple
  # error metrics and is only represented once in the df. When this happens, it
  # should be weighted higher when taking the average. Weights only matter for
  # when there is more than one parameter set.
  errs_df <- errs_df[errs_df$avg_id == avg_id, ]

  if (!with_mean_rank) {
    errs_df <- subset(errs_df, type != "best_mean")
  }

  avg_err <- errs_df %>%
    dplyr::group_by_at(c(group_cols, "avg_id")) %>%
    dplyr::summarize(across(c(cor, rMSE, mAPE),
                            list("mean" = ~ mean(.x),
                                 "sd" = ~ sd(.x))),
                     .groups = "drop") %>%
    as.data.frame()

  weights <- as.data.frame(table(errs_df$param_id))
  colnames(weights) <- c("param_id", "weight")

  ests_df <- subset(ests_df, param_id %in% errs_df$param_id &
                      tissue == unique(errs_df$tissue)) %>%
    tidyr::pivot_longer(cols = where(is.numeric),
                        names_to = "celltype",
                        values_to = "percent") %>%
    dplyr::mutate(avg_id = avg_id,
                  algorithm = unique(errs_df$algorithm),
                  granularity = unique(errs_df$granularity))

  cols_keep_ests <- setdiff(colnames(ests_df), c("percent", "param_id"))

  ests_df <- merge(ests_df, weights, by = "param_id")

  # Speed-up optimization: wtd.mean takes longer on data with 1 row and throws a
  # lot of warnings in that case. It's also unnecessary to take the mean/sd when
  # there is only 1 row. With the number of times this function gets called, it
  # can take a while to run on the full data set, and separating it out like
  # this speeds it up considerably even though it's messier.
  if (nrow(weights) == 1) {
    mean_fun <- function(value, weight) { value }
    sd_fun <- function(value, weight) { 0 }
  } else if (nrow(weights) == nrow(errs_df)) {
    mean_fun <- function(values, weights) { mean(values) }
    sd_fun <- function(values, weights) { sd(values) }
  } else { # any other number of rows
    mean_fun <- Hmisc::wtd.mean
    sd_fun <- function(value, weights) { sqrt(Hmisc::wtd.var(value, weights)) }
  }

  avg_est <- ests_df %>%
    dplyr::group_by_at(cols_keep_ests) %>%
    dplyr::summarize(percent_mean = mean_fun(percent, weight),
                     percent_sd = sd_fun(percent, weight),
                     .groups = "drop") %>%
    as.data.frame()

  return(list("avg_error" = avg_err, "avg_estimates" = avg_est))
}


Create_AveragesList <- function(errs_df, ests_df, group_cols, n_cores = 2, with_mean_rank = TRUE) {
  cl <- parallel::makeCluster(n_cores, outfile = "")
  parallel::clusterEvalQ(cl, library(magrittr, include.only = c("%>%")))

  # Although we could do this in one lapply statement, subsetting ests_df
  # takes a non-trivial amount of time because it's so large, so we do this to
  # break it into smaller pieces before having to call subset every time
  avg_list <- lapply(unique(errs_df$tissue), function(tiss) {
    errs_tmp <- subset(errs_df, tissue == tiss)
    ests_tmp <- subset(ests_df, tissue == tiss)

    avg_tmp <- parallel::parLapply(cl,
                                   X = unique(errs_tmp$avg_id),
                                   fun = Get_AverageStats,
                                   errs_df = errs_tmp,
                                   ests_df = ests_tmp,
                                   group_cols = group_cols,
                                   with_mean_rank = with_mean_rank)

    return(list("avg_errors" = List_to_DF(avg_tmp, "avg_error"),
                "avg_estimates" = List_to_DF(avg_tmp, "avg_estimates")))
  })

  names(avg_list) <- unique(errs_df$tissue)

  parallel::stopCluster(cl)

  return(avg_list)
}


Calculate_Significance <- function(avg_id, ests_df) {
  ests_param <- ests_df[ests_df$avg_id == avg_id, ]

  anov <- aov(percent_mean ~ diagnosis*celltype, data = ests_param)
  summ <- summary(anov)[[1]]
  tuk <- TukeyHSD(anov, "diagnosis:celltype")

  comparisons <- paste0("CT:", unique(ests_param$celltype),
                        "-AD:", unique(ests_param$celltype))

  tuk <- as.data.frame(tuk[[1]][comparisons,])
  tuk$p_adj <- tuk$`p adj`
  tuk$p_adj[is.na(tuk$p_adj)] <- 1 # Baseline zeros case causes this

  tuk$celltype <- stringr::str_split(rownames(tuk),
                                     pattern = ":",
                                     simplify = TRUE)[,3]

  tuk$anova_pval <- summ["diagnosis:celltype", "Pr(>F)"]
  tuk$avg_id <- avg_id

  return(tuk)
}


Get_MeanProps_Significance <- function(avg_list, group_cols, n_cores = 2) {
  cl <- parallel::makeCluster(n_cores, outfile = "")

  sig_list <- lapply(avg_list, function(avg_item) {
    errs_tmp <- avg_item$avg_errors
    tissue <- unique(errs_tmp$tissue)
    # Tissues are unique to a single bulk dataset so this works
    bulk_dataset <- unique(errs_tmp$test_data_name)

    print(paste(bulk_dataset, tissue))

    ests_ad <- avg_item$avg_estimates

    # Calculate significance for estimates from each parameter ID separately
    significant <- parallel::parLapply(cl,
                                       X = unique(ests_ad$avg_id),
                                       fun = Calculate_Significance,
                                       ests_df = ests_ad)

    significant <- List_to_DF(significant) %>%
      dplyr::mutate(tissue = tissue,
                    significant_05 = (p_adj <= 0.05),
                    significant_01 = (p_adj <= 0.01)) %>%
      dplyr::select(celltype, tissue, p_adj, significant_05, significant_01,
                    anova_pval, avg_id)

    # cap minimum p to avoid log(0) in downstream analysis
    significant$p_adj_thresh <- significant$p_adj
    significant$p_adj_thresh[significant$p_adj < 1e-8] <- 1e-8

    # Average cell type percentages across all AD or all CT samples in a given
    # parameter set
    mean_props <- ests_ad %>% subset(diagnosis %in% c("CT", "AD")) %>%
      group_by(avg_id, diagnosis, tissue, celltype, algorithm) %>%
      dplyr::summarize(mean_pct = mean(percent_mean),
                       sd_pct = sd(percent_mean),
                       rel_sd_pct = sd_pct / mean_pct,
                       count = n(),
                       .groups = "drop") %>%
      # Make one column for AD and one for CT for each of mean_pct, sd_pct,
      # rel_sd_pct, and count, calculate fold-change between AD and CT means
      pivot_wider(names_from = "diagnosis",
                  values_from = c("mean_pct", "sd_pct", "rel_sd_pct", "count")) %>%
      dplyr::mutate(fc = mean_pct_AD / mean_pct_CT,
                    log2_fc = log2(fc), # equivalent to log2(AD)-log2(CT)
                    # If one or both of the means is 0, rel_sd_pct values and the
                    # FC would be NA due to divide by zero, so they need to be
                    # changed to 0
                    fc = ifelse(is.na(fc), 0, fc),
                    log2_fc = ifelse(is.na(log2_fc), 0, log2_fc),
                    rel_sd_pct_AD = ifelse(is.na(rel_sd_pct_AD), 0, rel_sd_pct_AD),
                    rel_sd_pct_CT = ifelse(is.na(rel_sd_pct_CT), 0, rel_sd_pct_CT))

    mean_props <- merge(mean_props, significant) %>%
      dplyr::mutate(log_p = log(p_adj_thresh)) %>%
      # pull in missing parameter fields from errs_tmp
      merge(errs_tmp[, c("avg_id", group_cols)]) %>%
      as.data.frame()

    return(mean_props)
  })

  parallel::stopCluster(cl)

  return(List_to_DF(sig_list))
}


# TODO outdated
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
