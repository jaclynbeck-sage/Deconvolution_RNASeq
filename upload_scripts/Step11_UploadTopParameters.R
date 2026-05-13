library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
top_params_folder <- Folder(basename(dir_top_parameters),
                            parent = config::get("upload_synid")) |>
  synStore()

top_errors_folder <- Folder(basename(dir_top_errors),
                            parent = config::get("upload_synid")) |>
  synStore()

top_estimates_folder <- Folder(basename(dir_top_estimates),
                               parent = config::get("upload_synid")) |>
  synStore()

# Provenance TODO

github <- paste0(config::get("github_repo_url"), "Step11_Get_TopParamSets.R")

# Loop over each data set and algorithm
for (test_dataset in all_bulk_datasets()) {
  RecursiveUpload(folder_name = file.path(dir_top_parameters, test_dataset),
                  parent = top_params_folder,
                  provenance = list("executed" = github))

  RecursiveUpload(folder_name = file.path(dir_top_errors, test_dataset),
                  parent = top_errors_folder,
                  provenance = list("executed" = github))

  RecursiveUpload(folder_name = file.path(dir_top_estimates, test_dataset),
                  parent = top_estimates_folder,
                  provenance = list("executed" = github))
}
