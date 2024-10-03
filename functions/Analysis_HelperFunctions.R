library(parallel)
source(file.path("functions", "FileIO_HelperFunctions.R"))

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
    errs_df <- merge(data$means,
                     dplyr::select(data$params, reference_data_name,
                                   test_data_name, reference_input_type,
                                   normalization, regression_method),
                     by.x = "param_id", by.y = "row.names") %>%
      dplyr::mutate(algorithm = str_replace(param_id, "_.*", "")) %>%
      merge(data$exc_inh_ratio,
            by = c("tissue", "param_id"), all = TRUE)

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
      new_ests$algorithm <- str_replace(N, "_.*", "")

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
  file_list <- list.files(dir_top_parameters,
                          pattern = paste0("(",
                                           paste(bulk_datasets, collapse = "|"),
                                           ").*", granularity),
                          full.names = TRUE,
                          recursive = TRUE)

  qstats_list <- mclapply(file_list, function(file) {
    data <- readRDS(file)
    print(basename(file))

    file_id <- str_replace(basename(file), "top_parameters_", "") %>%
      str_replace(".rds", "")

    n_valid <- data.frame(n_valid_results = data$n_valid_results,
                          n_possible_results = data$n_possible_results,
                          algorithm = str_replace(file_id, "_.*", ""))

    # Add file parameters to the data frame
    if ("params" %in% names(data)) {
      params <- data$params %>%
        select(reference_data_name, test_data_name, granularity,
               reference_input_type, normalization, regression_method) %>%
        distinct()
    } else {
      # Extract parameters from the file_id string. These str_replaces are so
      # that normalizations like "log_cpm" and the granularity don't get split
      # up by str_split.
      file_id <- str_replace(file_id, "log_", "log.") %>%
        str_replace("_class", ".class")

      params <- as.data.frame(str_split(file_id, "_", simplify = TRUE))
      colnames(params) <- c("algorithm", "reference_data_name", "test_data_name",
                            "granularity", "reference_input_type", "normalization",
                            "regression_method")
      params$granularity <- str_replace(params$granularity, "\\.", "_")
      params$normalization <- str_replace(params$normalization, "\\.", "_")

      params <- params[, -1]
    }

    n_valid <- cbind(n_valid, params)

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


Rank_Errors <- function(errs_df, group_cols) {
  ranks <- errs_df %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::mutate(cor_rank = rank(-cor),
                  rMSE_rank = rank(rMSE),
                  mAPE_rank = rank(mAPE))
  ranks$mean_rank <- rowMeans(ranks[, c("cor_rank", "rMSE_rank", "mAPE_rank")])
  return(ranks)
}


Get_TopRanked <- function(df, n_top = 1) {
  top_ranked <- do.call(rbind, list(
    mutate(dplyr::slice_min(df, order_by = cor_rank, n = n_top, with_ties = FALSE),
           type = "best_cor"),
    mutate(dplyr::slice_min(df, order_by = rMSE_rank, n = n_top, with_ties = FALSE),
           type = "best_rMSE"),
    mutate(dplyr::slice_min(df, order_by = mAPE_rank, n = n_top, with_ties = FALSE),
           type = "best_mAPE"),
    mutate(dplyr::slice_min(df, order_by = mean_rank, n = n_top, with_ties = FALSE),
           type = "best_mean")
  ))

  return(top_ranked)
}


Find_BestParameters <- function(errs_df, group_cols, with_mean_rank = TRUE) {
  ranks <- errs_df %>%
    Rank_Errors(group_cols) %>%
    group_by_at(group_cols) %>%
    Get_TopRanked(n_top = 1) %>%
    ungroup()

  if (!with_mean_rank) {
    ranks <- subset(ranks, type != "best_mean")
  }

  # Some parameter IDs will be duplicated across multiple error metrics, this
  # creates a data frame with the list of unique parameter IDs associated with
  # each tissue.
  best_params <- ranks %>%
    dplyr::select(tissue, param_id) %>%
    dplyr::distinct()

  return(best_params)
}


Get_AverageStats <- function(errs_df, ests_df) {
  # There can be up to 4 parameter sets per combination of data input
  # parameters, but sometimes a single parameter set was the best for multiple
  # error metrics and is only represented once in the df. When this happens, it
  # should be weighted higher when taking the average. Weights only matter for
  # when there are 2 parameter sets, if there is 1 or 3 the average doesn't need
  # to be weighted.

  avg_id <- unique(errs_df$avg_id)

  # Ensures that param_ids that are the best for multiple error metrics are
  # represented that many times in the data frame. We also don't include the
  # param set with the top mean_rank here
  ranked <- Get_TopRanked(errs_df, n_top = 1) %>%
    subset(cor_rank == 1 | mAPE_rank == 1 | rMSE_rank == 1)

  avg_err <- ranked %>%
    dplyr::summarize(across(c(cor, rMSE, mAPE),
                            list("mean" = ~ mean(.x),
                                 "sd" = ~ sd(.x))),
                     .groups = "drop") %>%
    mutate(avg_id = avg_id)
  #cbind(errs_df[1, cols_keep_errs])

  weights <- as.data.frame(table(ranked$param_id))
  colnames(weights) <- c("param_id", "weight")

  ests_df <- subset(ests_df, param_id %in% ranked$param_id &
                      tissue == unique(ranked$tissue)) %>%
    dplyr::mutate(avg_id = avg_id)

  #cols_keep_errs <- setdiff(colnames(errs_df),
  #                          c("param_id", "cor", "rMSE", "mAPE"))

  data_keep_ests <- ests_df %>%
    dplyr::select(-percent, -param_id) %>%
    dplyr::distinct()

  ests_df <- merge(ests_df, weights, by = "param_id")

  # Speed-up optimization: wtd.mean takes longer on data with 1 row and throws a
  # lot of warnings in that case. It's also unnecessary to take the mean/sd when
  # there is only 1 row. With the number of times this function gets called, it
  # can take a while to run on the full data set, and separating it out like
  # this speeds it up considerably even though it's messier.
  if (nrow(ranked) == 1) {
    mean_fun <- function(value, weight) { value }
    sd_fun <- function(value, weight) { 0 }
  } else if (nrow(ranked) == 3) {
    mean_fun <- function(values, weights) { mean(values) }
    sd_fun <- function(values, weights) { sd(values) }
  } else { # any other number of rows
    mean_fun <- wtd.mean
    sd_fun <- function(value, weights) { sqrt(wtd.var(value, weights)) }
  }

  avg_est <- ests_df %>%
    group_by(sample, celltype) %>%
    dplyr::summarize(percent_mean = mean_fun(percent, weight),
                     percent_sd = sd_fun(percent, weight),
                     .groups = "drop") %>%
    merge(data_keep_ests, by = c("sample", "celltype"))

  return(list("avg_error" = avg_err, "avg_estimates" = avg_est))
}


Create_AveragesList <- function(errs_df, ests_df, n_cores = 2) {
  # Although we could do this in one lapply statement, subsetting best_errors
  # takes a non-trivial amount of time because it's so large, so we do this to
  # break it into smaller pieces before having to call subset every time
  avg_list <- lapply(unique(errs_df$tissue), function(tiss) {
    errs_tmp <- subset(errs_df, tissue == tiss)
    ests_tmp <- subset(ests_df, tissue == tiss)

    avg_tmp <- mclapply(unique(errs_tmp$avg_id), function(a_id) {
      print(a_id)
      errs_filt <- subset(errs_tmp, avg_id == a_id)
      return(Get_AverageStats(errs_filt, ests_tmp))
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


Get_MeanProps_Significance <- function(avg_list, n_cores = 2) {
  mean_props_all <- mclapply(avg_list, function(avg_item) {
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
      merge(select(errs_tmp, test_data_name, normalization, regression_method, avg_id),
            by = "avg_id")

    return(mean_props)
  }, mc.cores = n_cores)

  return(mean_props_all)
}
