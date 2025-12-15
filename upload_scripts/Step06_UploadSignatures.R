source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
sig_folder <- Folder(basename(dir_signatures), parent = config::get("upload_synid"))
sig_folder <- synStore(sig_folder, forceVersion = FALSE)

main_folders <- GetMainFolderIds()
provenance <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)

github <- paste0(config::get("github_repo_url"), "Step06_CalculateSignatures.R")

for (dataset in provenance$dataset) {
  files <- list.files(dir_signatures,
                      pattern = str_glue("{dataset}.*rds"),
                      full.names = TRUE)

  for (filename in files) {
    UploadFile(filename, sig_folder,
               list("used" = provenance$id[provenance$dataset == dataset],
                    "executed" = github))
  }
}
