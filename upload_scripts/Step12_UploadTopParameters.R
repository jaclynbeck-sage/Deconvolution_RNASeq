library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
top_params_folder <- Folder("12_top_parameters", parent = "syn68238853")
top_params_folder <- synStore(top_params_folder, forceVersion = FALSE)

# Provenance TODO

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step12_Get_TopParamSets.R"

# Loop over each data set and algorithm
for (test_dataset in c("Mayo", "MSBB", "ROSMAP")) {
  bulk_params_folder <- Folder(test_dataset, parent = top_params_folder)
  bulk_params_folder <- synStore(bulk_params_folder)

  algorithms <- list.files(file.path(dir_top_parameters, test_dataset))

  for (algorithm in algorithms) {
    algorithm_params_folder <- Folder(algorithm, parent = bulk_params_folder)
    algorithm_params_folder <- synStore(algorithm_params_folder)

    dir_alg <- file.path(dir_top_parameters, test_dataset, algorithm)

    files <- list.files(dir_alg, full.names = TRUE)

    # Parameter output lists
    for (filename in files) {
      # The name of the param file should match the name of the algorithm output file
      #param_name_base <- str_replace(filename, "errors_", "")
      #param_name_base <- str_replace(param_name_base, ".rds", "")

      #provenance <- subset(params_list_df, name_base == error_name_base)

      UploadFile(filename,
                 parent_folder = algorithm_params_folder,
                 provenance = list("used" = c(), "executed" = github))
    }
  }
}
