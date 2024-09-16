library(dplyr)
library(stringr)
library(reshape2)
library(foreach)
library(doParallel)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

options(scipen = 999)

cores <- 12
cl <- makeCluster(cores, type = "FORK", outfile = "top_params_output.txt")
registerDoParallel(cl)

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE",
                "Music", "Scaden", "Baseline")

for (bulk_dataset in bulk_datasets) {
  # TODO tmp
  if (FALSE) {
  # Dtangle only -- get top 10 scoring parameter sets for each file for use with HSPE
  algorithm <- "Dtangle"
  err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity,
                              reference_dataset = NULL)
  if (length(err_files) == 0) {
    message(str_glue("No data found for {bulk_dataset}/{algorithm}/{granularity}. Skipping..."))
    next
  }

  for (EF in err_files) {
    err_list <- readRDS(EF)

    if (length(err_list) == 0) {
      next
    }

    # Only use the errors across all tissues, not tissue-specific errors
    errs_all <- err_list$means$all_signature
    errs_all <- subset(errs_all, tissue == "All")

    cols_test <- colnames(select_if(errs_all, is.numeric))

    get_top10_vals <- function(errs_df, col_name) {
      if (length(grep("cor", col_name)) > 0) {
        slice_fun <- slice_max
      } else {
        slice_fun <- slice_min
      }

      # We calculated errors against 5 signatures for each parameter set, so
      # first we pick the best error among the 5 for each parameter set and
      # then pick the top 10 errors from that data
      top_vals <- errs_df %>%
        group_by(param_id) %>%
        slice_fun(.data[[col_name]], n = 1) %>% # top error for each param set
        group_by(tissue) %>%
        slice_fun(.data[[col_name]], n = 10) # top 10

      return(top_vals)
    }

    best_vals <- lapply(cols_test, function(col_name) {
      df <- errs_all %>%
        get_top10_vals(col_name)
      return(df$param_id)
    })
    best_vals <- unique(unlist(best_vals))
    best_params <- err_list$params[best_vals,]

    file_params <- best_params %>%
      select(reference_data_name, test_data_name, granularity,
             reference_input_type, normalization, regression_method) %>%
      distinct()

    saveRDS(best_params,
            file = file.path(dir_hspe_params,
                             paste0("hspe_params_",
                                    paste(file_params, collapse = "_"),
                                    ".rds")))
    print(paste(file_params, collapse = " "))
  }
  }

  # Top errors for each error metric, for each algorithm, for further analysis
  for (algorithm in algorithms) {
    err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity,
                                reference_dataset = NULL) # get all error files

    print(c(bulk_dataset, algorithm, granularity))

    if (length(err_files) == 0) {
      message(str_glue("No data found for {bulk_dataset}/{algorithm}/{granularity}. Skipping..."))
      next
    }

    foreach(EF = err_files) %dopar%  {
      err_list <- readRDS(EF)

      if (length(err_list) == 0) {
        next
      }

      #errs_all <- rbind(err_list$means$all_signature, err_list$means$all_lm)
      errs_all <- err_list$means$all_signature
      #errs_all <- subset(errs_all, signature == unique(err_list$params$reference_data_name))

      cols_test <- colnames(select_if(errs_all, is.numeric))

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
        df <- errs_all %>% group_by(tissue, solve_type, signature) %>%
          summarize(best_inds = get_best_vals(.data, col_name),
                    param_id = .data$param_id[best_inds],
                    .groups = "drop") %>%
          select(-best_inds)

        df$best_metric <- col_name
        return(df)
      })
      best_vals <- do.call(rbind, best_vals)

      errs_sub_all <- merge(best_vals, errs_all,
                            by = c("param_id", "tissue", "solve_type", "signature"),
                            all.x = TRUE, all.y = FALSE)

      best_param_ids <- unique(errs_sub_all$param_id)

      file_id <- str_replace(best_param_ids[1], "_[0-9].*", "")

      err_stats <- melt(errs_sub_all, variable.name = "metric",
                        id.vars = c("param_id", "tissue", "solve_type",
                                    "signature", "best_metric")) %>%
        group_by(tissue, solve_type, signature, metric) %>%
        summarize(mean_err = mean(value),
                  sd_err = sd(value),
                  rel_sd_err = sd_err / mean_err,
                  .groups = "drop")
      err_stats$file_id <- file_id

      est_stats <- subset(err_list$estimates, param_id %in% best_param_ids) %>%
        merge(errs_sub_all, by = "param_id") %>%
        select(-best_metric) %>%
        group_by(tissue, celltype, signature, sample) %>%
        summarize(mean_pct = mean(percent),
                  sd_pct = sd(percent),
                  rel_sd_pct = sd_pct / mean_pct,
                  .groups = "drop")


      best_errors <- list(means = errs_sub_all,
                          params = err_list$params[best_param_ids,],
                          by_sample = do.call(rbind, err_list$by_sample[best_param_ids]),
                          by_sample_lm = do.call(rbind, err_list$by_sample_lm[best_param_ids]),
                          pct_bad_inhibitory_ratio = subset(err_list$pct_bad_inhibitory_ratio, param_id %in% best_param_ids),
                          n_valid_results = err_list$n_valid_results,
                          n_possible_results = err_list$n_possible_results,
                          error_stats = err_stats,
                          estimate_stats = est_stats)

      names(best_errors$n_valid_results) <- file_id
      names(best_errors$n_possible_results) <- file_id

      print(basename(EF))

      # Save some information about how many times each parameter set came up as
      # the best. We ignore normalization/regression/etc since they're all the
      # same for this file
      param_strings <- format(best_errors$params, trim = TRUE) %>%
        select(-reference_data_name, -test_data_name,
               -normalization, -regression_method,
               -total_markers_used) %>%
        apply(1, paste, collapse = "_")

      param_strings <- data.frame(param_id = names(param_strings),
                                  param_string = as.character(param_strings))
      param_strings <- merge(param_strings, best_errors$means, by = "param_id")

      param_stats <- param_strings %>% group_by(param_string) %>%
        summarize(metrics = str_c(unique(best_metric), collapse = ", "),
                  param_ids = str_c(unique(param_id), collapse = ", "),
                  count = n())

      best_errors$param_stats <- param_stats
      best_errors$params$algorithm <- algorithm

      dir_top_params_alg <- file.path(dir_top_parameters, bulk_dataset, algorithm)
      dir.create(dir_top_params_alg, recursive = TRUE, showWarnings = FALSE)

      saveRDS(best_errors,
              file.path(dir_top_params_alg,
                        str_glue("top_parameters_{file_id}.rds")))
    }
  }
}

stopCluster(cl)
print("Finished")
