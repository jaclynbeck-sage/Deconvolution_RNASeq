library(dplyr)
library(stringr)
library(reshape2)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef") #, "seaAD")

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE",
                "Music", "Scaden", "Baseline")

for (bulk_dataset in bulk_datasets) {
  for (algorithm in algorithms) {
    errs_alg <- lapply(datasets, function(dataset) {
      print(c(bulk_dataset, algorithm, dataset))
      err_files <- Get_ErrorFiles(algorithm, dataset, bulk_dataset, granularity)

      if (length(err_files) == 0) {
        message(str_glue("No data found for {bulk_dataset}/{dataset}/{algorithm}. Skipping..."))
        return(NULL)
      }

      errs_dataset <- lapply(err_files, function(EF) {
        err_list <- readRDS(EF)

        if (length(err_list) == 0) {
          next
        }

        err_list$means$all_tissue$tissue_assignments = "All"
        err_list$means$all_tissue_lm$tissue_assignments = "All"
        err_list$means$all_tissue$type = "signature"
        err_list$means$all_tissue_lm$type = "lm"
        err_list$means$by_tissue$type = "signature"
        err_list$means$by_tissue_lm$type = "lm"

        rownames(err_list$params) <- err_list$means$all_tissue$param_id
        names(err_list$by_sample) <- rownames(err_list$params)
        names(err_list$by_sample_lm) <- rownames(err_list$params)
        names(err_list$pct_bad_inhibitory_ratio) <- rownames(err_list$params)

        cols <- colnames(err_list$means$all_tissue)

        errs_all <- do.call(rbind, lapply(err_list$means, function(X) { X[,cols]}))

        cols_test <- setdiff(colnames(errs_all), c("param_id", "tissue_assignments", "type"))

        get_best_vals <- function(errs_df, col_name) {
          if (length(grep("cor", col_name)) > 0) {
            top_ind <- which.max(errs_df[[col_name]])
          }
          else {
            top_ind <- which.min(errs_df[[col_name]])
          }
          return(top_ind)
        }

        best_vals <- lapply(cols_test, function(col_name) {
          df <- errs_all %>% group_by(tissue_assignments, type) %>%
                  summarize(best_inds = get_best_vals(.data, col_name),
                            param_id = .data$param_id[best_inds],
                            .groups = "drop")
          df$best_group <- col_name
          return(df)
        })
        best_vals <- do.call(rbind, best_vals)

        errs_sub_all <- best_vals %>% group_by(param_id, tissue_assignments, type) %>%
                          summarize(metrics = str_c(best_group, collapse = ", "),
                                    .groups = "drop")

        errs_sub_all <- merge(errs_sub_all, errs_all,
                              by = c("param_id", "tissue_assignments", "type"),
                              all.x = TRUE, all.y = FALSE)

        best_param_ids <- unique(errs_sub_all$param_id)

        err_stats <- melt(errs_sub_all, variable.name = "metric") %>%
                        group_by(tissue_assignments, type, metric) %>%
                        summarize(mean_err = mean(value),
                                  sd_err = sd(value),
                                  rel_sd_err = sd_err / mean_err,
                                  .groups = "drop")

        err_list_new <- list(means = errs_sub_all,
                             params = err_list$params[best_param_ids,],
                             by_sample = do.call(rbind, err_list$by_sample[best_param_ids]),
                             by_sample_lm = do.call(rbind, err_list$by_sample_lm[best_param_ids]),
                             pct_bad_inhibitory_ratio = err_list$pct_bad_inhibitory_ratio[best_param_ids],
                             n_valid_results = err_list$n_valid_results,
                             n_possible_results = err_list$n_possible_results,
                             error_stats = err_stats)

        return(err_list_new)
      })

      return(Flatten_ErrorList(errs_dataset))
    })

    errs_alg <- errs_alg[lengths(errs_alg) > 0]

    best_errors <- Flatten_ErrorList(errs_alg)

    if (length(best_errors) == 0) {
      message(str_glue("No data found for {bulk_dataset}/{algorithm}. Skipping..."))
      next
    }

    # Save some information about how many times each parameter set came up as
    # the best, ignoring dataset and normalization/regression
    param_strings <- format(best_errors$params, trim = TRUE) %>%
                        select(-reference_data_name, -test_data_name,
                               -normalization, -regression_method,
                               -total_markers_used) %>%
                        apply(1, paste, collapse = "_")

    param_strings <- data.frame(param_id = names(param_strings),
                                param_string = as.character(param_strings))
    param_strings <- merge(param_strings, best_errors$means, by = "param_id")

    param_stats <- param_strings %>% group_by(param_string) %>%
                      summarize(metrics = str_c(unique(metrics), collapse = ", "),
                                param_ids = str_c(unique(param_id), collapse = ", "),
                                count = n())
    best_errors$param_stats <- param_stats
    best_errors$params$algorithm <- algorithm

    saveRDS(best_errors,
            file.path(dir_best_errors,
                      str_glue("best_errors_{bulk_dataset}_{algorithm}_{granularity}.rds")))
  }
}
