# This script copies the top-performing Baseline estimates into the top
# estimates directory, as if they were run through step 13. We don't run step 13
# on Baseline estimates directly because we don't use a Baseline config file
# or an inner loop script on Baseline data.
library(plyr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

top_param_files <- list.files(dir_top_parameters, recursive = TRUE,
                              full.names = TRUE)
top_param_files <- top_param_files[grepl("Baseline", top_param_files)]

for (param_file in top_param_files) {
  top_params <- readRDS(param_file)

  params <- top_params$params

  file_info <- FileParams_FromParams(params)

  prev_res <- Load_AlgorithmOutputList("Baseline",
                                       file_info$reference_data_name,
                                       file_info$test_data_name,
                                       file_info$granularity,
                                       file_info$reference_input_type,
                                       file_info$normalization,
                                       file_info$regression_method)

  results_list <- prev_res[top_params$param_ids]

  message(paste("Found Baseline results containing",
                nrow(results_list[[1]]$estimates),
                "samples."))

  # Save the completed list
  name_base <- paste(file_info, collapse = "_")
  print(str_glue("Saving final list for {name_base}..."))
  Save_AlgorithmOutputList(results_list, "Baseline",
                           test_dataset = file_info$test_data_name,
                           name_base = name_base,
                           top_params = TRUE)
}
