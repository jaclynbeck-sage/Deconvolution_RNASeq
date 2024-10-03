library(dplyr)
library(stringr)
library(reshape2)
library(foreach)
library(doParallel)

source(file.path("functions", "Analysis_HelperFunctions.R"))

options(scipen = 999)

cores <- 12
cl <- makeCluster(cores, type = "FORK", outfile = "top_params_output.txt")
registerDoParallel(cl)

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

for (bulk_dataset in bulk_datasets) {
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
