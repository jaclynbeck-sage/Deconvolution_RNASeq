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

    # Merge a subset of the parameters into the mean errors data frame. These
    # parameters exist in every error file for every algorithm
    errs_df <- merge(data$means$all_signature,
                     select(data$params, reference_data_name, test_data_name,
                            reference_input_type, normalization,
                            regression_method),
                     by.x = "param_id", by.y = "row.names") %>%
      mutate(algorithm = str_replace(param_id, "_.*", ""),
             pct_valid_results = data$n_valid_results / data$n_possible_results)

    inh_ratio <- data.frame(param_id = names(data$pct_bad_inhibitory_ratio),
                            pct_bad_inhibitory_ratio = data$pct_bad_inhibitory_ratio)

    errs_df <- merge(errs_df, inh_ratio, by = "param_id")

    return(errs_df)
  })

  return(do.call(rbind, best_err_list))
}


Get_AllBestEstimatesAsDf <- function(bulk_datasets, granularity, metadata, best_params = NULL) {
  # List of "best estimates" files matching the bulk data set names and granularity
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
        merge(metadata, by = "sample", all.y = FALSE)

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


Find_BestSignature <- function(errs_df) {
  # Estimates for each param_id were scored against all 5 signatures. This
  # selects the signature that gave the best score along each error metric,
  # for each tissue separately
  best_signatures <- errs_df %>%
    subset(algorithm != "Baseline") %>%
    group_by(tissue, param_id, normalization, regression_method) %>%
    dplyr::summarize(best_cor = signature[which.max(cor)],
                     best_rMSE = signature[which.min(rMSE)],
                     best_mAPE = signature[which.min(mAPE)],
                     .groups = "drop")

  get_best_signature <- function(best_cor, best_rMSE, best_mAPE) {
    tmp <- table(c(best_cor, best_rMSE, best_mAPE))
    names(tmp)[which.max(tmp)]
  }

  best_signatures <- best_signatures %>%
    group_by(tissue, normalization, regression_method) %>%
    dplyr::summarize(best_signature = get_best_signature(best_cor, best_rMSE, best_mAPE),
                     .groups = "drop")

  return(best_signatures)
}


Find_BestParameters <- function(errs_df, group_cols) {
  # Because we scored against all 5 signature matrices, there will be parameter
  # ids in the errors df where the errors were the best score from a different
  # signature, not the single signature we are using going forward. This grabs
  # all parameter IDs associated with a specific set of parameters and picks the
  # one with the best score along each error metric, to get a single parameter
  # ID per error metric per set of parameters.
  # NOTE: The baseline "zeros" entries have NA correlation so they need to
  # return NA for best_cor, since otherwise calling which.max would return 0
  # entries and cause an error.
  best_params <- errs_df %>%
    group_by_at(group_cols) %>%
    dplyr::summarize(best_cor = ifelse(all(is.na(cor)), NA, param_id[which.max(cor)]),
                     best_rMSE = param_id[which.min(rMSE)],
                     best_mAPE = param_id[which.min(mAPE)],
                     .groups = "drop")

  # Some parameter IDs will be duplicated across multiple error metrics, this
  # creates a data frame with the list of unique parameter IDs associated with
  # each tissue.
  best_params <- best_params %>%
    melt(measure.vars = c("best_cor", "best_rMSE", "best_mAPE"),
         value.name = "param_id") %>%
    select(tissue, param_id) %>%
    subset(!is.na(param_id)) %>%
    distinct()

  return(best_params)
}


Get_AverageStats <- function(errs_df, ests_df) {
  # There can be up to 3 parameter sets per combination of data input
  # parameters, but sometimes a single parameter set was the best for multiple
  # error metrics and is only represented once in the df. When this happens, it
  # should be weighted higher when taking the average. Weights only matter for
  # when there are 2 parameter sets, if there is 1 or 3 the average doesn't need
  # to be weighted.

  avg_id <- unique(errs_df$avg_id)

  ests_df <- subset(ests_df, param_id %in% errs_df$param_id &
                      tissue == unique(errs_df$tissue)) %>%
    mutate(avg_id = avg_id)

  cols_keep_errs <- setdiff(colnames(errs_df),
                            c("param_id", "cor", "rMSE", "mAPE"))

  data_keep_ests <- ests_df %>%
    select(-percent, -param_id) %>%
    distinct()

  # Speed-up optimization: wtd.mean takes longer on data with 1 row and throws a
  # lot of warnings in that case. It's also unnecessary to take the mean/sd when
  # there is only 1 row, and unnecessary to do a weighted mean/sd when there are
  # 3 rows. With the number of times this function gets called, it can take a
  # while to run on the full data set, and separating it out like this speeds it
  # up considerably even though it's messier.
  if (nrow(errs_df) == 1) {
    mean_fun <- function(value, weight) { value }
    sd_fun <- function(value, weight) { 0 }

    errs_df$weight <- 1
    ests_df$weight <- 1

  } else if (nrow(errs_df) == 2) {
      bests <- c(which.max(errs_df$cor), which.min(errs_df$rMSE), which.min(errs_df$mAPE))
      weights <- as.data.frame(table(errs_df$param_id[bests]))
      colnames(weights) <- c("param_id", "weight")

      errs_df <- merge(errs_df, weights, by = "param_id")
      ests_df <- merge(ests_df, weights, by = "param_id")

      mean_fun <- wtd.mean
      sd_fun <- wtd.var

  } else { # nrow = 3
      mean_fun <- function(values, weights) { mean(values) }
      sd_fun <- function(values, weights) { sd(values) }

      errs_df$weight <- 1
      ests_df$weight <- 1
  }

  avg_err <- errs_df %>%
    dplyr::summarize(across(c(cor, rMSE, mAPE),
                            list("mean" = ~ mean_fun(.x, weight),
                                 "sd" = ~ sqrt(sd_fun(.x, weight))))) %>%
    cbind(errs_df[1, cols_keep_errs])

  avg_est <- ests_df %>%
    group_by(sample, celltype) %>%
    dplyr::summarize(percent_mean = mean_fun(percent, weight),
                     percent_sd = sqrt(sd_fun(percent, weight)),
                     .groups = "drop") %>%
    merge(data_keep_ests, by = c("sample", "celltype"))

  return(list("avg_error" = avg_err, "avg_estimates" = avg_est))
}


Create_AveragesList <- function(errs_df, ests_df) {
  # Although we could do this in one lapply statement, subsetting best_errors
  # takes a non-trivial amount of time because it's so large, so we do this to
  # break it into smaller pieces before having to call subset every time
  avg_list <- lapply(unique(errs_df$tissue), function(tiss) {
    errs_tmp <- subset(errs_df, tissue == tiss)
    ests_tmp <- subset(ests_df, tissue == tiss)

    avg_tmp <- lapply(unique(errs_tmp$avg_id), function(a_id) {
      print(a_id)
      errs_filt <- subset(errs_tmp, avg_id == a_id)
      return(Get_AverageStats(errs_filt, ests_tmp))
    })

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

    comparisons <- paste0("CT:", levels(ests_df$celltype),
                          "-AD:", levels(ests_df$celltype))

    tuk <- as.data.frame(tuk[[1]][comparisons,])
    tuk$p_adj <- tuk$`p adj`
    tuk$p_adj[is.na(tuk$p_adj)] <- 1 # Baseline zeros case causes this

    tuk$celltype <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,3]
    tuk$tissue <- tissue
    tuk$significant <- tuk$p_adj <= 0.05
    tuk$anova_significant <- summ["diagnosis:celltype", "Pr(>F)"] <= 0.05
    tuk$anova_significant[is.na(tuk$anova_significant)] <- FALSE # for Baseline zeros

    tuk$avg_id <- a_id
    return(tuk)
  })

  significant <- do.call(rbind, significant) %>%
    select(celltype, tissue, p_adj, significant, anova_significant, avg_id)

  # cap minimum p to avoid log(0) further down
  significant$p_adj_thresh <- significant$p_adj
  significant$p_adj_thresh[significant$p_adj < 1e-8] <- 1e-8

  return(significant)
}


Get_MeanProps_Significance <- function(avg_list) {
  mean_props_all <- lapply(avg_list, function(avg_item) {
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
      summarise(mean_pct = mean(percent_mean),
                sd_pct = sd(percent_mean),
                rel_sd_pct = sd_pct / mean_pct,
                count = n(),
                .groups = "drop")

    # Make one column for AD and one for CT for each of mean_pct, sd_pct,
    # rel_sd_pct, and count, calculate fold-change between AD and CT means
    mean_props <- pivot_wider(mean_props, names_from = "diagnosis",
                              values_from = c("mean_pct", "sd_pct",
                                              "rel_sd_pct", "count")) %>%
      mutate(fc = mean_pct_AD / mean_pct_CT,
             log2_fc = log2(fc)) # equivalent to log2(AD)-log2(CT)

    # If one or both of the means is 0, rel_sd_pct values and the FC would be NA
    # due to divide by zero, so they need to be changed to 0
    mean_props$fc[is.na(mean_props$fc)] <- 0
    mean_props$log2_fc[is.na(mean_props$log2_fc)] <- 0
    mean_props$rel_sd_pct_AD[is.na(mean_props$rel_sd_pct_AD)] <- 0
    mean_props$rel_sd_pct_CT[is.na(mean_props$rel_sd_pct_CT)] <- 0

    mean_props <- merge(mean_props, significant,
                        by = c("avg_id", "celltype", "tissue")) %>%
      mutate(log_p = log(p_adj_thresh))

    return(mean_props)
  })

  return(mean_props_all)
}
