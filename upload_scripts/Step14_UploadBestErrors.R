library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
errors_folder <- Folder("14_best_error_calculations", parent = "syn68238853")
errors_folder <- synStore(errors_folder, forceVersion = FALSE)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step14_CalculateBestErrors.R"

# Loop over each data set and algorithm
for (test_dataset in c("Mayo", "MSBB", "ROSMAP")) {
  bulk_errors_folder <- Folder(test_dataset, parent = errors_folder)
  bulk_errors_folder <- synStore(bulk_errors_folder)

  algorithms <- list.files(file.path(dir_errors, test_dataset))

  for (algorithm in algorithms) {
    algorithm_errors_folder <- Folder(algorithm, parent = bulk_errors_folder)
    algorithm_errors_folder <- synStore(algorithm_errors_folder)

    dir_alg <- file.path(dir_best_errors, test_dataset, algorithm)

    files <- list.files(dir_alg, full.names = TRUE)

    for (filename in files) {
      UploadFile(filename,
                 parent_folder = algorithm_errors_folder,
                 provenance = list("used" = c(), "executed" = github))
    }
  }
}
