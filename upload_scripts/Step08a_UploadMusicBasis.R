library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

models_folder <- Folder(basename(dir_algorithm_models), parent = config::get("upload_synid"))
models_folder <- synStore(models_folder, forceVersion = FALSE)

music_folder <- Folder(basename(dir_music_basis), parent = models_folder)
music_folder <- synStore(music_folder, forceVersion = FALSE)

# Provenance
main_folders <- GetMainFolderIds()
sc_df <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)

github <- paste0(config::get("github_repo_url"), "Step08_RunDeconvolutionAlgorithms.R")
datasets <- unique(sc_df$dataset)

for (dataset in datasets) {
  files <- list.files(dir_music_basis, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename,
               parent_folder = music_folder,
               provenance = list("used" = sc_df$id[sc_df$dataset == dataset],
                                 "executed" = github))
  }
}
