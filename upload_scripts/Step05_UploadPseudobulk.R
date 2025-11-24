source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
pseudo_folder <- Folder(basename(dir_pseudobulk), parent = config::get("upload_synid"))
pseudo_folder <- synStore(pseudo_folder, forceVersion = FALSE)

main_folders <- GetMainFolderIds()
provenance <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)

github <- paste0(config::get("github_repo_url"), "Step05_GeneratePseudobulk.R")

for (dataset in provenance$dataset) {
  files <- list.files(dir_pseudobulk, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, pseudo_folder,
               list("used" = provenance$id[provenance$dataset == dataset],
                    "executed" = github))
  }
}
