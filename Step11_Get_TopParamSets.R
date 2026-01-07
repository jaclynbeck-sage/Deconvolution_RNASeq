library(dplyr)
library(stringr)
library(parallel)

source(file.path("functions", "Analysis_HelperFunctions.R"))

cores <- max(parallel::detectCores() - 1, 1)
cluster_type <- "FORK"
cluster_outfile <- "top_params_output.txt"

granularities <- c("broad_class", "sub_class")
bulk_datasets <- all_bulk_datasets()
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

for (granularity in granularities) {
  for (bulk_dataset in bulk_datasets) {
    # Top errors for each error metric, for each algorithm, for further analysis
    for (algorithm in algorithms) {
      err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity,
                                  reference_dataset = NULL) # get all error files

      print(paste(bulk_dataset, algorithm, granularity))

      if (length(err_files) == 0) {
        message(str_glue("No data found for {bulk_dataset}/{algorithm}/{granularity}. Skipping..."))
        next
      }

      cl <- makeCluster(cores, type = cluster_type, outfile = cluster_outfile)

      parLapply(cl, err_files, function(EF) {
        err_list <- readRDS(EF)
        errs_all <- err_list$means

        best_params <- Find_BestParameters(errs_all,
                                           group_cols = c("tissue", "signature"),
                                           with_mean_rank = TRUE)

        best_param_ids <- unique(best_params$param_id)

        file_id <- str_replace(basename(EF), "errors_", "") %>%
          str_replace(".rds", "")

        best_errors <- list(param_ids = best_param_ids,
                            params = err_list$params[best_param_ids,],
                            file_id = file_id)

        print(basename(EF))

        # Save the list of top parameters
        dir_top_params_alg <- file.path(dir_top_parameters, bulk_dataset, algorithm)
        dir.create(dir_top_params_alg, recursive = TRUE, showWarnings = FALSE)

        saveRDS(best_errors,
                file.path(dir_top_params_alg,
                          str_glue("top_parameters_{file_id}.rds")))

        # Save the top errors and top estimates into new files to avoid having
        # to open many large files to get a few items

        # Errors
        best_error_list <- err_list
        best_error_list$means <- subset(best_error_list$means, param_id %in% best_param_ids)
        best_error_list$by_sample <- best_error_list$by_sample[best_param_ids]
        best_error_list$params <- best_error_list$params[best_param_ids, ]
        best_error_list$param_ids <- best_param_ids

        file_params <- FileParams_FromParams(best_error_list$params) |>
          mutate(mode = "best_estimates")

        Save_ErrorList(bulk_dataset, best_error_list, algorithm, file_params,
                       top_params = TRUE)

        # Estimates
        ests_list <- Load_AlgorithmOutputList(algorithm,
                                              file_params$reference_data_name,
                                              file_params$test_data_name,
                                              file_params$granularity,
                                              file_params$reference_input_type,
                                              file_params$normalization,
                                              file_params$regression_method)

        ests_list <- ests_list[best_param_ids]

        Save_AlgorithmOutputList(ests_list, algorithm,
                                 test_dataset = file_params$test_data_name,
                                 name_base = paste(file_params, collapse = "_"),
                                 top_params = TRUE)
      })

      stopCluster(cl)
    }
  }
}

print("Finished")
