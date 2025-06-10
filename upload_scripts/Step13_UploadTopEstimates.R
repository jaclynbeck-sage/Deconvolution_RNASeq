library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# TODO provenance

# Deconvolution WG Synapse space
output_folder <- Folder("13_top_estimates", parent = "syn68238853")
output_folder <- synStore(output_folder, forceVersion = FALSE)

algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE",
                "Music", "Scaden", "Baseline")

test_datasets <- c("Mayo", "MSBB", "ROSMAP")

# Parameter output lists
for (bulk in test_datasets) {
  bulk_folder <- Folder(bulk, parent = output_folder)
  bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

  dir_bulk <- file.path(dir_top_estimates, bulk)

  for (alg in algorithms) {
    alg_folder <- Folder(alg, parent = bulk_folder)
    alg_folder <- synStore(alg_folder, forceVersion = FALSE)

    files_alg <- list.files(file.path(dir_bulk, alg), full.names = TRUE)

    github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step13_RunDeconvolution_TopParameters.R"

    for (filename in files_alg) {
      UploadFile(filename,
                 parent_folder = alg_folder,
                 provenance = list("used" = NULL,
                                   "executed" = github))
    }
  }
}
