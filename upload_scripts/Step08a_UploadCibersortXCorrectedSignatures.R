library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

cx_folder <- Folder("cibersortx_batch_corrected", parent = config::get("upload_synid"))
cx_folder <- synStore(cx_folder, forceVersion = FALSE)

# Provenance
main_folders <- GetMainFolderIds()
sc_df <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)
bulk_df <- GetChildrenAsDf(main_folders[[basename(dir_bulk)]]$id)

github <- paste0(config::get("github_repo_url"), "Step08_RunDeconvolutionAlgorithms.R")
ref_datasets <- unique(sc_df$dataset)
bulk_datasets <- unique(bulk_df$dataset)

for (ref_dataset in ref_datasets) {
  for (bulk_dataset in bulk_datasets) {
    files <- list.files(file.path(dir_signatures, "cibersortx_batch_corrected"),
                        pattern = str_glue("{ref_dataset}_{bulk_dataset}"),
                        full.names = TRUE)
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
