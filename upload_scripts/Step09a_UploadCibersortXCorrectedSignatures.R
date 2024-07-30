library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

cx_folder <- Folder("cibersortx_batch_corrected", parent = "syn59480278")
cx_folder <- synStore(cx_folder, forceVersion = FALSE)

# Provenance
sc_df <- GetChildrenAsDf("syn58807549")
bulk_df <- GetChildrenAsDf("syn58841972")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step09_RunDeconvolutionAlgorithms.R"
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
