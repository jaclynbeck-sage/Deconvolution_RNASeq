
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
                   solve_type %in% params[["solve_type"]] &
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

  df <- df %>% mutate(
    tissue_full = paste(test_data_name, tissue),
    tissue = factor(tissue, levels = c("CBE", "TCX", "FP", "IFG", "PHG",
                                       "STG", "ACC", "DLPFC", "PCC")),
    regression_method = regression_renames[regression_method],
    normalization = str_replace(normalization, "counts", "cpm"),
    normalization = str_replace(normalization, "log_", ""),
    normalization = toupper(normalization),
    algorithm = str_replace(algorithm, "Music", "MuSiC"),
    data_transform = paste(normalization, "+", regression_method),
  )

  if ("cor" %in% colnames(df)) {
    df <- df %>% dplyr::rename(Correlation = cor, RMSE = rMSE, MAPE = mAPE)
  }

  if ("reference_data_name" %in% colnames(df)) {
    df <- df %>%
      mutate(reference_data_name = dataset_renames[reference_data_name])
  }

  if ("type" %in% colnames(df)) {
    df <- df %>%
      mutate(error_type = str_replace(type, "best_", "")) %>%
      group_by(error_type) %>%
      mutate(error_type = if (unique(error_type) == "cor") "Correlation"
             else toupper(error_type)) %>%
      ungroup()
  }

  return(df)
}

Count_Grouped <- function(df, groups) {
  df %>%
    group_by_at(groups) %>%
    dplyr::summarize(count = n(), .groups = "drop")
}


Get_BestDataTransform <- function(ranked_df) {
  best_dt <- Count_Grouped(ranked_df, c("tissue", "data_transform")) %>%
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

Load_BestErrors <- function(granularity) {
  best_errors_list <- readRDS(file.path(dir_analysis,
                                        str_glue("best_errors_{granularity}.rds")))
  best_errors <- Paper_Renames(best_errors_list$best_errors_all)
  best_errors_top <- Paper_Renames(best_errors_list$best_errors_toplevel)

  quality_stats <- subset(best_errors, algorithm != "Baseline") %>%
    select(-Correlation, -RMSE, -MAPE)

  quality_stats_top <- subset(best_errors_top, algorithm != "Baseline") %>%
    select(-Correlation, -RMSE, -MAPE)

  # There are up to 3 param_ids per set of input parameters. Get the best of each
  # error metric for each set
  best_errors <- best_errors  %>%
    group_by(tissue, tissue_full, reference_data_name, test_data_name, algorithm,
             normalization, regression_method) %>%
    dplyr::summarize(Correlation = max(Correlation),
                     RMSE = min(RMSE),
                     MAPE = min(MAPE),
                     .groups = "drop") %>%
    melt(variable.name = "error_type")

  best_errors_top <- best_errors_top  %>%
    group_by(tissue, tissue_full, reference_data_name, test_data_name, algorithm,
             normalization, regression_method) %>%
    dplyr::summarize(Correlation = max(Correlation),
                     RMSE = min(RMSE),
                     MAPE = min(MAPE),
                     .groups = "drop") %>%
    melt(variable.name = "error_type")

  baselines <- subset(best_errors, algorithm == "Baseline")
  best_errors <- subset(best_errors, algorithm != "Baseline")

  baselines_top <- subset(best_errors_top, algorithm == "Baseline")
  best_errors_top <- subset(best_errors_top, algorithm != "Baseline")

  ranked <- readRDS(file.path(dir_analysis, str_glue("ranked_errors_{granularity}.rds")))

  ranked_df <- Paper_Renames(ranked$ranked_errors_best_signatures) %>%
    subset(type != "best_mean") # drop mean_rank

  best_dt <- Get_BestDataTransform(ranked_df)

  items <- list(best_errors, best_errors_top, baselines, baselines_top,
    ranked_df, best_dt, quality_stats, quality_stats_top)
  names(items) <- paste(c("best_errors", "best_errors_top", "baselines",
                          "baselines_top","ranked_df", "best_dt",
                          "quality_stats", "quality_stats_top"),
                        str_replace(granularity, "_class", ""),
                        sep = "_")

  return(items)
}
