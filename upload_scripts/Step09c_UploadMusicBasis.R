library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

music_folder <- Folder("music_basis", parent = "syn59489760")
music_folder <- synStore(music_folder, forceVersion = FALSE)

# Provenance
sc_df <- GetChildrenAsDf("syn58807549")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step09_RunDeconvolutionAlgorithms.R"
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
