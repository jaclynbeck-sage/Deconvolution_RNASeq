library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(viridis)


source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE", "Music", "Scaden")

do_plot <- TRUE

#algorithm <- "Dtangle"
#bulk_dataset <- "Mayo"
#dataset <- "leng"

plot_6_errs <- function(data, x_aes, color_aes = NULL, shape_aes = NULL) {
  plt1 <- ggplot(data,
                 aes(x = .data[[x_aes]], y = .data[["cor"]],
                     color = .data[[color_aes]], shape = .data[[shape_aes]])) +
    geom_point() + scale_x_log10() + theme_bw()

  plt2 <- ggplot(data,
                 aes(x = .data[[x_aes]], y = .data[["rMSE"]],
                     color = .data[[color_aes]], shape = .data[[shape_aes]])) +
    geom_point() + scale_x_log10() + theme_bw()

  plt3 <- ggplot(data,
                 aes(x = .data[[x_aes]], y = .data[["mAPE"]],
                     color = .data[[color_aes]], shape = .data[[shape_aes]])) +
    geom_point() + scale_x_log10() + theme_bw()

  plt4 <- ggplot(data,
                 aes(x = .data[[x_aes]], y = .data[["cor_lm"]],
                     color = .data[[color_aes]], shape = .data[[shape_aes]])) +
    geom_point() + scale_x_log10() + theme_bw()

  plt5 <- ggplot(data,
                 aes(x = .data[[x_aes]], y = .data[["rMSE_lm"]],
                     color = .data[[color_aes]], shape = .data[[shape_aes]])) +
    geom_point() + scale_x_log10() + theme_bw()

  plt6 <- ggplot(data,
                 aes(x = .data[[x_aes]], y = .data[["mAPE_lm"]],
                     color = .data[[color_aes]], shape = .data[[shape_aes]])) +
    geom_point() + scale_x_log10() + theme_bw()

  print((plt1 + plt2 + plt3) / (plt4 + plt5 + plt6))
}

plot_marker_combos <- function(data, target_patterns, error_metric, dataset_title,
                               color_aes = NULL, shape_aes = NULL, linetype_aes = NULL) {
  color_field <- sort(unique(data[[color_aes]]))
  all_colors <- viridis(length(color_field), option = "H")
  names(all_colors) <- color_field

  plt_list <- lapply(target_patterns, function(pattern) {
    colors <- rep("#DDDDDD88", length(color_field))
    names(colors) <- color_field

    target <- grepl(pattern, color_field)
    colors[target] <- all_colors[target]

    plt <- ggplot(data,
                  aes(x = .data[["total_markers_used"]], y = .data[[error_metric]],
                      color = .data[[color_aes]], shape = .data[[shape_aes]],
                      linetype = .data[[linetype_aes]])) +
      geom_point() + scale_x_log10() + theme_bw() + geom_line() +
      scale_color_manual(values = colors)

    return(plt)
  })

  dataset_title <- paste(dataset_title, error_metric, sep = " / ")
  print(wrap_plots(plt_list, nrow = 4) + plot_annotation(title = dataset_title))
}


plot_algorithm_args <- function(data, error_metric, dataset_title,
                                color_aes, shape_aes, linetype_aes, group_aes) {
  plt <- ggplot(data,
                aes(x = .data[["total_markers_used"]], y = .data[[error_metric]],
                    color = .data[[color_aes]], shape = .data[[shape_aes]],
                    linetype = .data[[linetype_aes]])) +
    geom_point() + scale_x_log10() + theme_bw() + geom_line() +
    scale_color_viridis(option = "H", discrete = TRUE, begin = 0.1, end = 0.9) +
    facet_wrap(group_aes, nrow = 4, scales = "fixed")

  dataset_title <- paste(dataset_title, error_metric, sep = " / ")
  print(plt + plot_annotation(title = dataset_title))
}


error_cols <- c("cor", "rMSE", "mAPE")

for (bulk_dataset in bulk_datasets) {
  for (algorithm in algorithms) {
    summary_alg <- lapply(datasets, function(dataset) {
      err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity, dataset)

      if (length(err_files) == 0) {
        message(str_glue("No data found for {bulk_dataset}/{dataset}/{algorithm}. Skipping..."))
        next
      }

      summary_dataset <- lapply(err_files, function(EF) {
        err_list <- readRDS(EF)

        if (length(err_list) == 0) {
          next
        }

        errs_all <- rbind(err_list$means$all_signature,
                          err_list$means$all_lm)
        errs_all <- merge(errs_all, err_list$params,
                          by.x = "param_id", by.y = "row.names")

        errs_all$marker_combo <- paste(errs_all$marker_type,
                                       errs_all$marker_subtype,
                                       errs_all$marker_input_type)

        errs_all$marker_order[errs_all$marker_order == "None"] <- "distance"

        errs_all$n_marker_type <- sapply(errs_all$n_markers, function(N) {
          if (N <= 1) {
            return("percent")
          }
          return("fixed")
        })

        target_patterns <- sort(unique(paste(errs_all$marker_type,
                                             errs_all$marker_subtype)))
        # Adding a pattern with just the wildcard will force the function to
        # color-code every marker combo, not just a few
        target_patterns <- c(paste0("^", target_patterns), "*")

        # backwards compatibility
        if (!("reference_input_type" %in% colnames(errs_all))) {
          errs_all$reference_input_type <- "signature"
        }

        file_params <- errs_all %>%
          select(reference_data_name, test_data_name, granularity,
                 reference_input_type, normalization, regression_method) %>%
          distinct()

        display_title <- paste(algorithm, paste(file_params, collapse = " "))
        file_title <- paste(algorithm, paste(file_params, collapse = "_"), sep = "_")

        pdf_dir <- file.path(dir_figures, bulk_dataset, algorithm)
        if (!dir.exists(pdf_dir)) {
          dir.create(pdf_dir, recursive = TRUE)
        }

        if (do_plot) {
          pdf(file.path(pdf_dir, paste0("param_set_graphs_", file_title, ".pdf")),
              width=24, height = 24)

          subsets_plot <- expand.grid(tissue = unique(errs_all$tissue),
                                      solve_type = unique(errs_all$solve_type),
                                      stringsAsFactors = FALSE)

          #subsets_plot <- subset(subsets_plot,
          #                       (solve_type == "signature" & signature != "none") |
          #                         (solve_type == "lm" & signature == "none"))

          errs_tmp = errs_all %>%
                        group_by(tissue, solve_type, marker_combo, marker_order,
                                 n_marker_type, total_markers_used) %>%
                        summarize(cor = max(cor),
                                  rMSE = min(rMSE),
                                  mAPE = min(mAPE),
                                  .groups = "drop")

          for (R in 1:nrow(subsets_plot)) {
            errs_plot <- subset(errs_tmp, tissue == subsets_plot$tissue[R] &
                                  solve_type == subsets_plot$solve_type[R])

            display_title2 <- str_glue(paste0("{display_title} / ",
                                              "tissue {subsets_plot$tissue[R]}, ",
                                              "using {subsets_plot$solve_type[R]}"))

            # total_markers_used vs marker_combo
            for (error_metric in error_cols) {
              plot_marker_combos(errs_plot, target_patterns, error_metric,
                                 dataset_title = display_title2,
                                 color_aes = "marker_combo",
                                 shape_aes = "marker_order",
                                 linetype_aes = "n_marker_type")
            }

            # These algorithms have algorithm-specific parameters to graph
            if (algorithm %in% c("DeconRNASeq", "DWLS", "Music")) {
              if (algorithm == "DeconRNASeq") {
                algorithm_arg <- "use_scale"
              }
              if (algorithm == "DWLS") {
                algorithm_arg <- "solver_type"
              }
              if (algorithm == "Music") {
                rename_centered <- c("TRUE" = "centered", "FALSE" = "not centered")
                rename_normalize <- c("TRUE" = "normalized", "FALSE" = "not normalized")

                errs_all$music_arg <- paste(rename_centered[as.character(errs_all$centered)],
                                            rename_normalize[as.character(errs_all$normalize)],
                                            sep = " / ")
                algorithm_arg <- "music_arg"
              }

              cols_group <- c(setdiff(colnames(errs_plot), error_cols),
                              algorithm_arg)

              errs_plot2 <- errs_all %>%
                              subset(tissue == subsets_plot$tissue[R] &
                                       solve_type == subsets_plot$solve_type[R]) %>%
                              group_by(across(all_of(cols_group))) %>%
                                summarize(cor = max(cor),
                                          rMSE = min(rMSE),
                                          mAPE = min(mAPE),
                                          .groups = "drop")

              for (error_metric in error_cols) {
                plot_algorithm_args(errs_plot2, error_metric, display_title2,
                                    color_aes = algorithm_arg,
                                    shape_aes = "marker_order",
                                    linetype_aes = "n_marker_type",
                                    group_aes = "marker_combo")
              }
            } # end if (algorithm %in% c("DeconRNASeq", "DWLS", "Music"))
          } # end for (R in 1:nrow(subsets_plot))
        } # end if (do_plot)

        for (col in colnames(errs_all)) {
          if (!is.numeric(errs_all[[col]])) {
            errs_all[[col]] <- factor(errs_all[[col]])
          }
        }

        # Stats on the top 20 parameters for each error metric
        best_params <- lapply(error_cols, function(error_metric) {
          if (grepl("cor", error_metric)) {
            errs_all %>% group_by(tissue, solve_type) %>%
              slice_max(order_by = .data[[error_metric]], n = 20)
          }
          else {
            errs_all %>% group_by(tissue, solve_type) %>%
              slice_min(order_by = .data[[error_metric]], n = 20)
          }
        })
        names(best_params) <- error_cols

        fields <- c("marker_combo", "marker_order", "n_marker_type", "total_markers_used")
        if (algorithm %in% c("DeconRNASeq", "DWLS", "Music")) {
          fields <- c(fields, algorithm_arg)
        }

        summary_stats <- lapply(error_cols, function(error_metric) {
          bp <- best_params[[error_metric]]
          best_totals <- lapply(fields, function(field) {
            df <- if (is.numeric(bp[[field]])) {
              bp %>% group_by(tissue, solve_type) %>%
                      summarize("best" = .data[[field]][1],
                                "min" = min(.data[[field]]),
                                "max" = max(.data[[field]]),
                                "mean" = mean(.data[[field]]),
                                "sd" = sd(.data[[field]]),
                                .groups = "drop")
            } else {
              tab <- table(bp$tissue, bp$solve_type, bp[[field]],
                           dnn = c("tissue", "solve_type", field))
              tab <- dcast(as.data.frame(tab, stringsAsFactors = FALSE),
                           as.formula(paste("tissue + solve_type ~", field)),
                           value.var = "Freq", stringsAsFactors = FALSE,
                           drop = FALSE)
              tab
            }
            colnames(df)[3:ncol(df)] <- paste0(field, ".", colnames(df)[3:ncol(df)])
            return(df)
          })
          names(best_totals) <- fields
          return(Reduce(function(d1, d2) {
            merge(d1, d2, by = c("tissue", "solve_type"))
          }, best_totals))
        })
        names(summary_stats) <- error_cols

        summary_stats <- do.call(rbind, summary_stats)
        summary_stats$error_metric <- str_replace(rownames(summary_stats), "\\..*", "")

        if (do_plot) {
          for_plot <- select(summary_stats, !starts_with("total_markers"))

          for_plot <- melt(for_plot, variable.name = "parameter", value.name = "count")
          for_plot$parameter <- str_replace(for_plot$parameter, "marker_combo.", "")
          for_plot[for_plot == 0] <- NA

          # Force the y-axis to have labels in a specific order. Order is reversed
          # so the first label is at the top instead of the bottom.
          col_order <- unique(for_plot$parameter)
          col_order <- col_order[length(col_order):1]

          plot_title <- "Count of appearance in top 20 results per error metric"

          plt <- ggplot(for_plot, aes(x = tissue, y = parameter, color = count, size = count)) +
            geom_count() + facet_wrap(solve_type ~ error_metric, ncol = 6) +
            scale_color_viridis(option = "plasma", direction = -1) + theme_bw() +
            scale_y_discrete(limits = col_order) +
            ggtitle(display_title, subtitle = plot_title) +
            theme(text = element_text(size = 20),
                  plot.margin = margin(t = 1, r = 5, b = 14, l = 5,
                                       unit = "in"))

          print(plt)

          dev.off()
        }

        # UNUSED
        if (FALSE) {
          plot_6_errs(errs_all, x_aes = "total_markers_used",
                      color_aes = "marker_type", shape_aes = "marker_order")

          plot_6_errs(errs_all, x_aes = "total_markers_used",
                      color_aes = "marker_order", shape_aes = "marker_type")

          plot_6_errs(errs_all, x_aes = "total_markers_used",
                      color_aes = "marker_subtype", shape_aes = "marker_type")

          plot_6_errs(errs_all, x_aes = "total_markers_used",
                      color_aes = "marker_input_type", shape_aes = "marker_type")
        }

        # UNUSED
        if (FALSE) {
          more_inhibitory <- more_inhibitory[errs_all$param_id]
          errs_all$more_inhibitory <- err_list$pct_bad_inhibitory_ratio[errs_all$param_id]

          plot_6_errs(errs_all, x_aes = "total_markers_used",
                      color_aes = "more_inhibitory")
        }


        #errs_all$algorithm <- algorithm
        #return(errs_all)
        summary_stats <- as.data.frame(summary_stats)
        summary_stats$algorithm <- algorithm
        file_params <- do.call(rbind, rep(list(file_params), nrow(summary_stats)))
        summary_stats <- cbind(summary_stats, file_params)

        print(display_title)

        return(summary_stats)
      })

      # Sometimes some factors don't show up at all in a file, so that summary
      # doesn't have a column for those. This ensures that all summaries have
      # the same columns in the same order.
      cols <- unique(unlist(sapply(summary_dataset, colnames)))
      file_fields <- c("reference_data_name", "test_data_name", "granularity",
                       "reference_input_type", "normalization", "regression_method")
      cols <- c(sort(setdiff(cols, file_fields)), file_fields)
      summary_dataset <- lapply(summary_dataset, function(S) {
        missing <- setdiff(cols, colnames(S))
        if (length(missing) > 0) {
          S[,missing] <- 0
        }
        return(S[,cols])
      })

      return(do.call(rbind, summary_dataset))
    })

    saveRDS(summary_alg, file.path(dir_tmp, "summary_alg.rds"))

    cols <- unique(unlist(sapply(summary_alg, colnames)))
    file_fields <- c("reference_data_name", "test_data_name", "granularity",
                     "reference_input_type", "normalization", "regression_method")
    cols <- c(sort(setdiff(cols, file_fields)), file_fields)
    summary_alg <- lapply(summary_alg, function(S) {
      missing <- setdiff(cols, colnames(S))
      if (length(missing) > 0) {
        S[,missing] <- 0
      }
      return(S[,cols])
    })
    summary_alg <- do.call(rbind, summary_alg)
    write.csv(summary_alg, file.path(dir_figures,
                                     paste0("param_set_stats_", bulk_dataset,
                                            "_", algorithm, "_", granularity, ".csv")))
  }
}
