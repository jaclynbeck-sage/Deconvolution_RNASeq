source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
sig_folder <- Folder(basename(dir_signatures), parent = config::get("upload_synid"))
sig_folder <- synStore(sig_folder, forceVersion = FALSE)

main_folders <- synGetChildren(config::get("upload_synid"))$asList()
main_folders <- do.call(rbind, main_folders) |> as.data.frame()
singlecell_folder <- main_folders$id[main_folders$name == basename(dir_singlecell)] |> unlist()

provenance <- GetChildrenAsDf(singlecell_folder)

github <- paste0(config::get("github_repo_url"), "Step06_CalculateSignatures.R")

for (dataset in provenance$dataset) {
  files <- list.files(dir_signatures, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, sig_folder,
               list("used" = provenance$id[provenance$dataset == dataset],
                    "executed" = github))
  }
}
