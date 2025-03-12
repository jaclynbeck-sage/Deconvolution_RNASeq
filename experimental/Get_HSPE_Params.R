library(dplyr)
library(stringr)
library(reshape2)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))

options(scipen = 999)

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

for (bulk_dataset in bulk_datasets) {
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

    # Only use the errors across all tissues, not tissue-specific errors
    errs_all <- err_list$means
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

stopCluster(cl)
print("Finished")
