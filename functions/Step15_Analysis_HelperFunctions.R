library(parallel)
source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

# Loads all the "best_errors" files matching specific parameters into one big
# list.
#
# Arguments:
#   bulk_datasets - a vector of bulk dataset names to load ("Mayo", "MSBB",
#                   and/or "ROSMAP")
#   granularity - either "broad_class" or "sub_class"
#
# Returns:
#   a data frame with all errors from all files concatenated together
Get_AllBestErrorsAsDf <- function(bulk_datasets, granularity, n_cores = 2) {
  # List of "best errors" files matching the bulk data set names and granularity
  file_list <- list.files(dir_best_errors,
                          pattern = paste0("(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*", granularity),
                          full.names = TRUE,
                          recursive = TRUE)

  best_err_list <- mclapply(file_list, function(file) {
    data <- readRDS(file)
    print(basename(file))

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
    print(basename(file))

    file_id <- str_replace(basename(file), "quality_stats_", "") %>%
      str_replace(".rds", "")

    # TODO

    return(n_valid)
  }, mc.cores = n_cores)

  return(do.call(rbind, qstats_list))
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


Get_BestDataTransform <- function(ranked_df) {
  best_dt <- Count_ParamFrequency(ranked_df, c("tissue", "data_transform")) %>%
    group_by(tissue) %>%
    slice_max(order_by = count, with_ties = FALSE) %>%
    mutate(normalization = str_replace(data_transform, " \\+.*", ""),
           regression_method = str_replace(data_transform, ".*\\+ ", "")) %>%
    as.data.frame()

  best_dt <- tidyr::expand_grid(best_dt,
                                algorithm = c(unique(ranked_df$algorithm), "Baseline"))
  best_dt$normalization[best_dt$algorithm == "MuSiC"] <- "CPM"
  best_dt$normalization[best_dt$normalization == "TMM" &
                          best_dt$algorithm == "CibersortX"] <- "CPM"

  best_dt <- best_dt %>%
    mutate(data_transform = paste(normalization, "+", regression_method)) %>%
    select(tissue, algorithm, data_transform) %>%
    distinct()

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

  return(ungroup(top_ranked))
}


Get_TopErrors <- function(errors_df, group_cols, n_cores, with_mean_rank = TRUE) {
  # Special case: The signature doesn't matter when all percentages are zeros, so
  # for the "zeros" baseline data we just subset to have one unique param_id per
  # tissue/data transform. The zeros data also needs to be held out separately
  # from the other baseline data and not combined with it at the top level.
  best_zeros <- subset(errors_df, reference_data_name == "zeros") %>%
    Get_TopRanked(group_cols, n_top = 1, with_mean_rank = with_mean_rank) %>%
    mutate(signature = NA)

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
    "stats" = error_stats
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


Calculate_EstimateStats <- function(est_pcts, group_cols) {
  est_stats <- est_pcts %>%
    melt(variable.name = "celltype",
         value.name = "percent",
         id.vars = c("param_id", group_cols)) %>%
    group_by_at(c(group_cols, "celltype")) %>%
    summarize(mean_pct = mean(percent),
              sd_pct = sd(percent),
              rel_sd_pct = sd_pct / mean_pct,
              .groups = "drop")

  return(est_stats)
}


Count_ParamFrequency <- function(df, groups, pivot_column, algorithm_specific = FALSE) {
  df %>%
    group_by_at(groups) %>%
    dplyr::summarize(count = n(), .groups = "drop") %>%
    pivot_longer(cols = all_of(pivot_column), names_to = "parameter_name",
                 values_to = "parameter_value") %>%
    mutate(algorithm_specific = algorithm_specific)
}


Get_AverageStats <- function(errs_df, ests_df, group_cols, with_mean_rank = TRUE) {
  # There can be up to 4 parameter sets per combination of data input
  # parameters, but sometimes a single parameter set was the best for multiple
  # error metrics and is only represented once in the df. When this happens, it
  # should be weighted higher when taking the average. Weights only matter for
  # when there is more than one parameter set.

  avg_id <- unique(errs_df$avg_id)

  if (!with_mean_rank) {
    errs_df <- subset(errs_df, type != "best_mean")
  }

  avg_err <- errs_df %>%
    group_by_at(c(group_cols, "avg_id")) %>%
    dplyr::summarize(across(c(cor, rMSE, mAPE),
                            list("mean" = ~ mean(.x),
                                 "sd" = ~ sd(.x))),
                     .groups = "drop") %>%
    as.data.frame()

  weights <- as.data.frame(table(errs_df$param_id))
  colnames(weights) <- c("param_id", "weight")

  ests_df <- subset(ests_df, param_id %in% errs_df$param_id &
                      tissue == unique(errs_df$tissue)) %>%
    melt(id.vars = c("param_id", "tissue", "sample", "diagnosis"),
         variable.name = "celltype",
         value.name = "percent") %>%
    dplyr::mutate(avg_id = avg_id,
                  algorithm = unique(errs_df$algorithm))

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
    group_by_at(cols_keep_ests) %>%
    dplyr::summarize(percent_mean = mean_fun(percent, weight),
                     percent_sd = sd_fun(percent, weight),
                     .groups = "drop") %>%
    as.data.frame()

  return(list("avg_error" = avg_err, "avg_estimates" = avg_est))
}


Create_AveragesList <- function(errs_df, ests_df, group_cols, n_cores = 2, with_mean_rank = TRUE) {
  # Although we could do this in one lapply statement, subsetting ests_df
  # takes a non-trivial amount of time because it's so large, so we do this to
  # break it into smaller pieces before having to call subset every time
  avg_list <- lapply(unique(errs_df$tissue), function(tiss) {
    errs_tmp <- subset(errs_df, tissue == tiss)
    ests_tmp <- subset(ests_df, tissue == tiss)

    avg_tmp <- mclapply(unique(errs_tmp$avg_id), function(a_id) {
      print(a_id)
      errs_filt <- subset(errs_tmp, avg_id == a_id)

      return(Get_AverageStats(errs_filt, ests_tmp, group_cols, with_mean_rank))
    }, mc.cores = n_cores)

    return(list("avg_errors" = do.call(rbind, lapply(avg_tmp, "[[", "avg_error")),
                "avg_estimates" = do.call(rbind, lapply(avg_tmp, "[[", "avg_estimates"))))
  })

  names(avg_list) <- unique(errs_df$tissue)

  return(avg_list)
}


Calculate_Significance <- function(ests_df, tissue) {
  significant <- lapply(unique(ests_df$avg_id), function(a_id) {
    ests_param <- subset(ests_df, avg_id == a_id)

    anov <- aov(percent_mean ~ diagnosis*celltype, data = ests_param)
    summ <- summary(anov)[[1]]
    tuk <- TukeyHSD(anov, "diagnosis:celltype")

    comparisons <- paste0("CT:", levels(ests_param$celltype),
                          "-AD:", levels(ests_param$celltype))

    tuk <- as.data.frame(tuk[[1]][comparisons,])
    tuk$p_adj <- tuk$`p adj`
    tuk$p_adj[is.na(tuk$p_adj)] <- 1 # Baseline zeros case causes this

    tuk$celltype <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,3]
    tuk$tissue <- tissue
    tuk$significant_05 <- tuk$p_adj <= 0.05
    tuk$significant_01 <- tuk$p_adj <= 0.01
    tuk$anova_pval <- summ["diagnosis:celltype", "Pr(>F)"]

    tuk$avg_id <- a_id
    return(tuk)
  })

  significant <- do.call(rbind, significant) %>%
    dplyr::select(celltype, tissue, p_adj, significant_05, significant_01,
                  anova_pval, avg_id)

  # cap minimum p to avoid log(0) further down
  significant$p_adj_thresh <- significant$p_adj
  significant$p_adj_thresh[significant$p_adj < 1e-8] <- 1e-8

  return(significant)
}


Get_MeanProps_Significance <- function(avg_item, group_cols) {
  errs_tmp <- avg_item$avg_errors
  tissue <- unique(errs_tmp$tissue)
  bulk_dataset <- unique(errs_tmp$test_data_name) # Tissues are unique to a single bulk dataset so this works

  print(paste(bulk_dataset, tissue))

  ests_ad <- avg_item$avg_estimates

  # Calculate significance for estimates from each parameter ID separately
  significant <- Calculate_Significance(ests_ad, tissue)

  # Average cell type percentages across all AD or all CT samples in a given
  # parameter set
  mean_props <- ests_ad %>% subset(diagnosis %in% c("CT", "AD")) %>%
    group_by(avg_id, diagnosis, tissue, celltype, algorithm) %>%
    dplyr::summarize(mean_pct = mean(percent_mean),
                     sd_pct = sd(percent_mean),
                     rel_sd_pct = sd_pct / mean_pct,
                     count = n(),
                     .groups = "drop")

  # Make one column for AD and one for CT for each of mean_pct, sd_pct,
  # rel_sd_pct, and count, calculate fold-change between AD and CT means
  mean_props <- pivot_wider(mean_props, names_from = "diagnosis",
                            values_from = c("mean_pct", "sd_pct",
                                            "rel_sd_pct", "count")) %>%
    dplyr::mutate(fc = mean_pct_AD / mean_pct_CT,
                  log2_fc = log2(fc)) # equivalent to log2(AD)-log2(CT)

  # If one or both of the means is 0, rel_sd_pct values and the FC would be NA
  # due to divide by zero, so they need to be changed to 0
  mean_props$fc[is.na(mean_props$fc)] <- 0
  mean_props$log2_fc[is.na(mean_props$log2_fc)] <- 0
  mean_props$rel_sd_pct_AD[is.na(mean_props$rel_sd_pct_AD)] <- 0
  mean_props$rel_sd_pct_CT[is.na(mean_props$rel_sd_pct_CT)] <- 0

  mean_props <- merge(mean_props, significant,
                      by = c("avg_id", "celltype", "tissue")) %>%
    dplyr::mutate(log_p = log(p_adj_thresh)) %>%
    # pull in missing parameter fields from errs_tmp
    merge(errs_tmp[, c("avg_id", group_cols)]) %>%
    as.data.frame()

  return(mean_props)
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
