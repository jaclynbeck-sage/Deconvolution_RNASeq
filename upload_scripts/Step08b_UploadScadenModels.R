library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

models_folder <- Folder(basename(dir_algorithm_models), parent = config::get("upload_synid"))
models_folder <- synStore(models_folder, forceVersion = FALSE)

scaden_folder <- Folder(basename(dir_scaden_models), parent = models_folder)
scaden_folder <- synStore(scaden_folder, forceVersion = FALSE)

# Provenance
main_folders <- GetMainFolderIds()
bulk_df <- GetChildrenAsDf(main_folders[[basename(dir_bulk)]]$id)

# The scaden pseudobulk folder is a sub-folder of the main pseudobulk folder
pseudo_df <- GetChildrenAsDf(main_folders[[basename(dir_pseudobulk)]]$id,
                             types = list("folder"))
scaden_df <- GetChildrenAsDf(pseudo_df$id[pseudo_df$name == "scaden_simulated"])


scaden_df$normalization <- str_replace(scaden_df$name, ".*_", "") |>
  str_replace("\\.rds", "")

# Loop over each model and upload -- exclude the "tmp" directory
models <- list.files(dir_scaden_models, pattern = "[^tmp]", full.names = TRUE)

github <- paste0(config::get("github_repo_url"), "Step08_RunDeconvolutionAlgorithms.R")

for (model_folder_name in models) {
  datasets <- ExtractDatasetName(model_folder_name)
  norm <- ifelse(grepl("tmm", model_folder_name), "tmm", "cpm")
  gran <- ifelse(grepl("broad_class", model_folder_name), "broad_class", "sub_class")

  prov_bulk <- subset(bulk_df, dataset %in% datasets)
  prov_scaden <- subset(scaden_df, dataset %in% datasets & granularity == gran &
                          normalization == norm)

  RecursiveUpload(model_folder_name,
                  parent_folder = scaden_folder,
                  provenance = list("used" = c(prov_bulk$id, prov_scaden$id),
                                    "executed" = github))
}
