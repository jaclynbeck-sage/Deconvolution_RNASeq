library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

hspe_folder <- Folder("hspe_params", parent = "syn59489760")
hspe_folder <- synStore(hspe_folder, forceVersion = FALSE)

# Provenance
github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step12_Get_TopParamSets.R"

files <- list.files(dir_hspe_params, full.names = TRUE)
for (filename in files) {
  UploadFile(filename,
             parent_folder = hspe_folder,
             provenance = list("executed" = github))
}

