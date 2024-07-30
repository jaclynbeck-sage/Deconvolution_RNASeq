
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
Plot_ErrsByAlgorithm <- function(errors_df, params, error_names,
                                 test_data_names = c(),
                                 x_axis = "reference_data_name",
                                 facet_var = "algorithm",
                                 fill_colors = NULL, ...) {
  errs <- subset(errors_df, error_type %in% error_names &
                   solve_type %in% params[["solve_type"]] &
                   normalization %in% params[["normalization"]] &
                   regression_method %in% params[["regression_method"]])

  if (!is.null(test_data_names) & length(test_data_names) > 0) {
    errs <- subset(errs, test_data_name %in% test_data_names)
  }

  plt <- Plot_FacetBoxPlot(errs, x_axis = x_axis, facet_var = facet_var,
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

# A generic faceted box plot function with arguments:
#   errors_df - dataframe that must at least have an 'error_type' column
#   x_axis - string, name of variable on the x axis
#   facet_var - string, name of variable to divide the facet by, after it is
#               divided by error_type
#   ... - User can specify additional arguments to be passed to aes, like:
#           fill = "test_data_name", color = "tissue"
Plot_FacetBoxPlot <- function(errors_df, x_axis, facet_var,
                              fill_colors = NULL, ...) {
  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) { sym(X) })

  plt <- ggplot(errors_df, aes(x = !!sym(x_axis), y = value, !!!aes_opts))

  # Number of columns matches the number of unique items in the last facet column
  facet_cols <- length(unique(errors_df[,facet_var[length(facet_var)]]))

  plt <- plt + geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(facet_var, scales = "fixed",
               ncol = facet_cols)

  if (!is.null(fill_colors)) {
    plt <- plt + scale_fill_manual(values = fill_colors)
  }
  return(plt)
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
