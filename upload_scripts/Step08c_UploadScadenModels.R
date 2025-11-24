library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

models_folder <- Folder("08_algorithm_models", parent = config::get("upload_synid"))
models_folder <- synStore(models_folder, forceVersion = FALSE)

scaden_folder <- Folder("scaden_models", parent = models_folder)
scaden_folder <- synStore(scaden_folder, forceVersion = FALSE)

# Provenance
main_folders <- GetMainFolderIds()
sc_df <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)
bulk_df <- GetChildrenAsDf(main_folders[[basename(dir_bulk)]]$id)

provenance_df <- rbind(sc_df, bulk_df)

# Loop over each model and upload -- exclude the "tmp" directory
models <- list.files(dir_scaden_models, pattern = "[^tmp]", full.names = TRUE)

github <- paste0(config::get("github_repo_url"), "Step08_RunDeconvolutionAlgorithms.R")

for (model_folder_name in models) {
  datasets <- ExtractDatasetName(model_folder_name)
  provenance_sub <- subset(provenance_df, dataset %in% datasets)
  RecursiveUpload(model_folder_name,
                  parent_folder = scaden_folder,
                  provenance = list("used" = provenance_sub$id,
                                    "executed" = github))
}
