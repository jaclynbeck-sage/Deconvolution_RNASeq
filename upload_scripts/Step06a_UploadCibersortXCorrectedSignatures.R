library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

sig_folder <- Folder(basename(dir_signatures), parent = config::get("upload_synid"))
sig_folder <- synStore(sig_folder, forceVersion = FALSE)

cx_folder <- Folder(basename(dir_cibersort_corrected_signatures), parent = sig_folder)
cx_folder <- synStore(cx_folder, forceVersion = FALSE)

# Provenance
main_folders <- GetMainFolderIds()
sc_df <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)
bulk_df <- GetChildrenAsDf(main_folders[[basename(dir_bulk)]]$id)

github <- paste0(config::get("github_repo_url"), "Step06_CalculateSignatures.R")
ref_datasets <- unique(sc_df$dataset)
bulk_datasets <- unique(bulk_df$dataset)

for (ref_dataset in ref_datasets) {
  for (bulk_dataset in bulk_datasets) {
    files <- list.files(dir_cibersort_corrected_signatures,
                        pattern = str_glue("{ref_dataset}_{bulk_dataset}"),
                        full.names = TRUE)

    if (length(files) == 0) {
      next
    }

    used <- c(sc_df$id[sc_df$dataset == ref_dataset],
              bulk_df$id[bulk_df$dataset == bulk_dataset])

    for (filename in files) {
      UploadFile(filename,
                 parent_folder = cx_folder,
                 provenance = list("used" = used,
                                   "executed" = github))
    }
  }
}
