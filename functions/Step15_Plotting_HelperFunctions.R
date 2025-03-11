
MakePropsDataframe <- function(props_truth, props_est, est_field) {
  props_melt <- as.data.frame(props_truth)
  props_melt$subject <- rownames(props_melt)
  props_melt <- melt(props_melt)
  colnames(props_melt) <- c("subject", "celltype", "prop_truth")

  ests_melt <- lapply(props_est, "[[", est_field)
  ests_melt <- lapply(names(ests_melt), FUN = function(X) {
    ests_melt[[X]] <- as.data.frame(ests_melt[[X]])
    ests_melt[[X]]$name <- X
    ests_melt[[X]]$subject <- rownames(ests_melt[[X]])
    ests_melt[[X]]
  })
  ests_melt <- do.call(rbind, ests_melt)
  ests_melt <- melt(ests_melt)
  colnames(ests_melt) <- c("name", "subject", "celltype", "prop_est")

  ests_df <- merge(props_melt, ests_melt, by = c("subject", "celltype"))
  return(ests_df)
}


# Each facet is an error type + algorithm, with datasets on the x-axis.
# It is assumed that 'error_names' is either a single name of a single error,
# or a vector of errors that should all have the same axis constraints,
# since the axes are fixed in the figure.
Plot_ErrsByAlgorithm <- function(errors_df, error_names,
                                 x_axis = "reference_data_name",
                                 facet_vars = "algorithm",
                                 fill_colors = NULL, ...) {
  errs <- subset(errors_df, error_type %in% error_names)

  plt <- Plot_FacetBoxPlot(errs, x_axis = x_axis, facet_vars = facet_vars,
                           fill_colors = fill_colors, ...)

  return(plt)
}


# Each facet is an error type + dataset, with algorithms on the x-axis
Plot_ErrsByDataset <- function(errors_df, params, error_names,
                               test_data_names = c(),
                               x_axis = "algorithm",
                               facet_var = "reference_data_name",
                               fill_colors = NULL, ...) {
  errs <- subset(errors_df, error_type %in% error_names &
                   normalization %in% params[["normalization"]] &
                   regression_method %in% params[["regression_method"]])

  if (!is.null(test_data_names) & length(test_data_names) > 0) {
    errs <- subset(errs, test_data_name %in% test_data_names)
  }

  plt <- Plot_FacetBoxPlot(errs, x_axis = x_axis,
                           facet_var = facet_var,
                           fill_colors = fill_colors, ...)
  return(plt)
}


Plot_BasicBoxPlot <- function(errors_df, x_axis,
                              fill_colors = NULL, ...) {
  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) {
    if (is.character(X)) sym(X) else X
  })

  plt <- ggplot(errors_df, aes(x = !!sym(x_axis), !!!aes_opts,
                               ymin = min_val, ymax = max_val,
                               lower = lower_quartile, middle = median_val,
                               upper = upper_quartile))

  plt <- plt + geom_boxplot(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(fill_colors)) {
    plt <- plt + scale_fill_manual(values = fill_colors)
  }
  return(plt)
}


# A generic faceted box plot function that plots data as a facet_grid, with
# arguments:
#   errors_df - dataframe that must at least have an 'error_type' column
#   x_axis - string, name of variable on the x axis
#   facet_vars - vector or string, names of variables to divide the facet by,
#                after it is divided by error_type. If facet_vars contains two
#                elements, the plot will use facet_grid, and the first element
#                should be facet rows, second element should be facet columns.
#                For any other number of elements, the plot will use facet_wrap.
#   ... - User can specify additional arguments to be passed to aes, like:
#           fill = "test_data_name", color = "tissue"
Plot_FacetBoxPlot <- function(errors_df, x_axis, facet_vars,
                              fill_colors = NULL, ...) {
  plt <- Plot_BasicBoxPlot(errors_df, x_axis = x_axis,
                           fill_colors = fill_colors, ...)

  # The lines below will auto-inject "..." into the aes() statement
  #aes_opts <- lapply(list(...), function(X) { if (is.character(X)) sym(X) else X })

  #plt <- ggplot(errors_df, aes(x = !!sym(x_axis), !!!aes_opts,
  #                             ymin = min_val, ymax = max_val,
  #                             lower = lower_quartile, middle = median_val,
  #                             upper = upper_quartile))

  #plt <- plt + geom_boxplot(stat = "identity") +
  #  theme_bw() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (length(facet_vars) == 2) {
    plt <- plt + facet_grid(as.formula(paste(facet_vars, collapse = "~")),
                            scales = "fixed")
  } else {
    if (length(facet_vars) == 1) {
      facet_cols = NULL
    } else {
      # Number of columns matches the number of unique items in the last facet column
      facet_cols <- length(unique(errors_df[,facet_vars[length(facet_vars)]]))
    }
    plt <- plt + facet_wrap(facet_vars, scales = "fixed", ncol = facet_cols)
  }

  #if (!is.null(fill_colors)) {
  #  plt <- plt + scale_fill_manual(values = fill_colors)
  #}
  return(plt)
}


Plot_BasicViolinPlot <- function(errors_df, x_axis,
                                 fill_colors = NULL, ...) {
  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) {
    if (is.character(X)) sym(X) else X
  })

  plt <- ggplot(errors_df, aes(x = !!sym(x_axis), !!!aes_opts,
                               ymin = min_val, ymax = max_val,
                               lower = lower_quartile, middle = median_val,
                               upper = upper_quartile))

  plt <- plt + geom_boxplot(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(fill_colors)) {
    plt <- plt + scale_fill_manual(values = fill_colors)
  }
  return(plt)
}


Plot_FacetViolinPlot <- function(errors_df, x_axis, facet_vars,
                                 fill_colors = NULL, ...) {

  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) { if (is.character(X)) sym(X) else X })

  plt <- ggplot(errors_df, aes(x = !!sym(x_axis), y = value, !!!aes_opts))

  plt <- plt + geom_violin(drop = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (length(facet_vars) == 2) {
    plt <- plt + facet_grid(as.formula(paste(facet_vars, collapse = "~")),
                            scales = "fixed")
  } else {
    if (length(facet_vars) == 1) {
      facet_cols = NULL
    } else {
      # Number of columns matches the number of unique items in the last facet column
      facet_cols <- length(unique(errors_df[,facet_vars[length(facet_vars)]]))
    }
    plt <- plt + facet_wrap(facet_vars, scales = "fixed", ncol = facet_cols)
  }

  if (!is.null(fill_colors)) {
    plt <- plt + scale_fill_manual(values = fill_colors)
  }
  return(plt)
}



Create_BoxStats <- function(errs_df, grouping_cols) {
  errs_df %>%
    group_by_at(grouping_cols) %>%
    dplyr::summarize(max_val = max(value, na.rm = TRUE),
                     min_val = min(value, na.rm = TRUE),
                     median_val = median(value, na.rm = TRUE),
                     upper_quartile = quantile(value, probs = 0.75, na.rm = TRUE),
                     lower_quartile = quantile(value, probs = 0.25, na.rm = TRUE),
                     .groups = "drop") %>%
    as.data.frame()
}


Extract_MeanErrors <- function(datasets, datatypes, best_params, dir_output) {
  errs_means <- lapply(datasets, function(dataset) {
    errs_m <- lapply(datatypes, function(datatype) {
      errs <- readRDS(file = file.path(dir_output, paste0("errors_", dataset, "_",
                                                          datatype, "_broad.rds")))
      # Pull out the errors in the "means" list
      # Pull out the errors in the "means" list
      errs <- lapply(errs, '[[', "means")

      # Add the algorithm name ("method") and the params name to the data frame
      # so we can tell all the algorithms apart
      errs <- lapply(names(errs), FUN = function(X) {
        errs[[X]]$method <- X
        errs[[X]]$name <- rownames(errs[[X]])
        errs[[X]] <- subset(errs[[X]], name %in% best_params$name)
        return(errs[[X]])
      })

      errs <- melt(do.call(rbind, errs))
      errs$datatype <- datatype
      return(errs)
    })

    # Resulting data frame has fields for method, name, variable (error type),
    # value (error value), and datatype
    return(do.call(rbind, errs_m))
  })

  # All params names start with the dataset name. The mutate statement extracts
  # the dataset name from this and adds it as a field to the data frame
  errs_means <- do.call(rbind, errs_means) %>% mutate(dataset = str_replace(name, "_.*", ""))
  errs_means <- errs_means %>% rename(error_type = "variable")
  return(errs_means)
}


Extract_ErrorsByCelltype <- function(errs_list, algorithm, best_params) {
  errs_by_celltype <- errs_list[[algorithm]][["by_celltype"]]
  errs_by_celltype <- lapply(names(errs_by_celltype), FUN = function(X) {
    errs_by_celltype[[X]]$name <- X
    errs_by_celltype[[X]]$celltype <- rownames(errs_by_celltype[[X]])
    errs_by_celltype[[X]]
  })
  errs_df <- melt(do.call(rbind, errs_by_celltype)) %>% rename(error_type = "variable")
  errs_df <- subset(errs_df, name %in% best_params$name)
  return(errs_df)
}

Extract_ErrorsBySubject <- function(errs_list, algorithm, best_params) {
  errs_by_subject <- errs_list[[algorithm]][["by_subject"]]
  errs_by_subject <- lapply(names(errs_by_subject), FUN = function(X) {
    errs_by_subject[[X]]$name <- X
    errs_by_subject[[X]]$subject <- rownames(errs_by_subject[[X]])
    errs_by_subject[[X]]
  })
  errs_df <- melt(do.call(rbind, errs_by_subject)) %>% rename(error_type = "variable")
  errs_df <- subset(errs_df, name %in% best_params$name)
  return(errs_df)
}


Plot_TruthVsEstimates_Dots <- function(ests_df, titles, n_col = 4, axis_limits = c(0, 1)) {
  ggplot(ests_df, aes(x = prop_truth, y = prop_est, color = name)) +
    geom_jitter(size = 0.5) + geom_abline(slope = 1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom", aspect.ratio = 1) +
    lims(x = axis_limits, y = axis_limits) +
    scale_color_discrete(labels = titles) +
    #scale_color_viridis(labels = titles, option = "turbo", discrete = TRUE) + # TODO try better colors like this
    facet_wrap(~celltype, scales = "free", ncol = n_col)
}

Plot_TruthVsEstimates_Lines <- function(ests_df, titles, n_col = 4, axis_limits = c(0,1)) {
  ggplot(ests_df, aes(x = prop_truth, y = avg_est, color = name)) +
    geom_line() + geom_abline(slope = 1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none", aspect.ratio = 1) +
    lims(x = axis_limits, y = axis_limits) +
    scale_color_discrete(labels = titles) +
    facet_wrap(~celltype, scales = "free", ncol = n_col)
}

Paper_Renames <- function(df) {
  dataset_renames <- c("cain" = "Cain 2020", "lau" = "Lau 2020",
                       "leng" = "Leng 2021", "mathys" = "Mathys 2019",
                       "seaRef" = "Gabitto 2023",
                       "random_biased" = "Random (biased)",
                       "random_educated" = "Random (educated)",
                       "random_uniform" = "Random (uniform)",
                       "zeros" = "All zeros")
  regression_renames <- c("none" = "No regression", "edger" = "edgeR",
                          "lme" = "LME", "dream" = "dream")

  if (all(c("normalization", "regression_method") %in% colnames(df))) {
    df <- df %>%
      mutate(regression_method = regression_renames[regression_method],
             normalization = str_replace(normalization, "counts", "cpm"),
             normalization = str_replace(normalization, "log_", ""),
             normalization = toupper(normalization),
             data_transform = paste(normalization, "+", regression_method),)
  }

  if ("algorithm" %in% colnames(df)) {
    df <- df %>%
      mutate(algorithm = str_replace(algorithm, "Music", "MuSiC"))
  }

  if ("tissue" %in% colnames(df)) {
    df <- df %>%
      mutate(tissue = factor(tissue,
                             levels = c("CBE", "TCX", "FP", "IFG", "PHG",
                                        "STG", "ACC", "DLPFC", "PCC"))
    )

    if ("test_data_name" %in% colnames(df)) {
      df <- df %>%
        mutate(tissue_full = paste(test_data_name, as.character(tissue)))
    }
  }

  if ("cor" %in% colnames(df)) {
    df <- df %>% dplyr::rename(Correlation = cor, RMSE = rMSE, MAPE = mAPE)
  }

  if ("reference_data_name" %in% colnames(df)) {
    df <- df %>%
      mutate(reference_data_name = dataset_renames[reference_data_name])
  }

  if ("type" %in% colnames(df)) {
    df <- df %>%
      mutate(error_type = str_replace(type, "best_", ""),
             error_type = case_when(error_type == "cor" ~ "Correlation",
                                    TRUE ~ toupper(error_type)))
  }

  if ("error_metric" %in% colnames(df)) {
    df <- df %>%
      mutate(error_metric = case_when(error_metric == "cor" ~ "Correlation",
                                      TRUE ~ toupper(error_metric)))
  }

  return(df)
}


Load_BestErrors <- function(granularity) {
  best_errors_list <- readRDS(file.path(dir_analysis,
                                        str_glue("best_errors_{granularity}.rds")))
  best_errors <- Paper_Renames(best_errors_list$all$errors)
  best_errors_top <- Paper_Renames(best_errors_list$toplevel$errors)

  #quality_stats <- Paper_Renames(best_errors_list$all$stats)
  #quality_stats_top <- Paper_Renames(best_errors_list$toplevel$stats)

  top3_by_tissue <- best_errors_list$top3_by_tissue %>%
    Paper_Renames() %>%
    subset(type != "best_mean") # drop mean rank if it exists TODO do we want to drop this?

  # There are up to 3 param_ids per set of input parameters. Get the best of each
  # error metric for each set
  best_errors <- best_errors  %>%
    group_by(tissue, tissue_full, reference_data_name, test_data_name, algorithm,
             normalization, regression_method, data_transform) %>%
    dplyr::summarize(Correlation = max(Correlation),
                     RMSE = min(RMSE),
                     MAPE = min(MAPE),
                     .groups = "drop") %>%
    melt(variable.name = "error_metric")

  best_errors_top <- best_errors_top  %>%
    group_by(tissue, tissue_full, test_data_name, algorithm,
             normalization, regression_method, data_transform) %>%
    dplyr::summarize(Correlation = max(Correlation),
                     RMSE = min(RMSE),
                     MAPE = min(MAPE),
                     .groups = "drop") %>%
    melt(variable.name = "error_metric")

  baselines <- subset(best_errors, algorithm == "Baseline")
  best_errors <- subset(best_errors, algorithm != "Baseline")

  baselines_top <- subset(best_errors_top, algorithm == "Baseline")
  best_errors_top <- subset(best_errors_top, algorithm != "Baseline")

  best_dt <- Get_BestDataTransform(top3_by_tissue,
                                   algorithms = c(unique(best_errors$algorithm), "Baseline"))

  items <- list(best_errors, best_errors_top, baselines, baselines_top,
    top3_by_tissue, best_dt) #, quality_stats, quality_stats_top)
  names(items) <- paste(c("best_errors", "best_errors_top", "baselines",
                          "baselines_top","top3_by_tissue", "best_dt"),
                          #"quality_stats", "quality_stats_top"),
                        str_replace(granularity, "_class", ""),
                        sep = "_")

  return(items)
}


Load_QualityStats <- function(bulk_datasets, granularity) {
  qstats_list <- lapply(bulk_datasets, function(bulk) {
    readRDS(file.path(dir_analysis,
                      str_glue("quality_stats_all_{bulk}_{granularity}.rds")))
  })

  # qstats_list will be a list of lists, so we combine inner lists of the same
  # name from each qstats_list item
  qstats_list_final <- lapply(names(qstats_list[[1]]), function(item_name) {
    List_to_DF(qstats_list, item_name) %>%
      Paper_Renames()
  })

  names(qstats_list_final) <- names(qstats_list[[1]])

  return(qstats_list_final)
}


Count_Grouped <- function(ranked_df, group_cols) {
  ranked_df %>%
    group_by_at(group_cols) %>%
    dplyr::count() %>%
    dplyr::rename(count = n)
}


Load_Significance <- function(significance_df, best_dt, p_sig = 0.01, log2_cap = 1) {
  # Fill in missing data with NA
  params <- expand.grid(celltype = unique(significance_df$celltype),
                        tissue = unique(significance_df$tissue),
                        algorithm = unique(significance_df$algorithm),
                        normalization = unique(significance_df$normalization),
                        regression_method = unique(significance_df$regression_method))

  significance_df <- merge(significance_df, params,
                           by = colnames(params),
                           all = TRUE)
  significance_df$p_adj_thresh[is.na(significance_df$p_adj_thresh)] <- 1

  # Fill in any NA data transforms created by missing data
  significance_df$data_transform <- paste(significance_df$normalization, "+",
                                          significance_df$regression_method)

  if (granularity == "sub_class") {
    ct_order <- c("Astrocyte", "Endothelial", paste0("Exc.", 1:10),
                  paste0("Inh.", 1:7), "Microglia", "Oligodendrocyte", "OPC",
                  "Pericyte", "VLMC")
    significance_df$celltype <- factor(significance_df$celltype, levels = ct_order)
  }

  sig_final <- merge(significance_df, best_dt)

  # Cap log2_fc values to +/- 1 so color scaling is better. Need to account for
  # Inf and -Inf values. Set non-significant log2 values to NA.
  sig_final$log2_fc[sig_final$fc == 0] <- NA
  sig_final$log2_fc[is.infinite(sig_final$log2_fc)] <- log2_cap
  sig_final$log2_fc[sig_final$log2_fc > 1] <- log2_cap
  sig_final$log2_fc[sig_final$log2_fc < -1] <- -log2_cap

  sig_final$log2_fc[sig_final$p_adj_thresh >= p_sig] <- NA

  return(sig_final)
}
