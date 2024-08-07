library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
scaden_folder <- Folder("scaden_models", parent = "syn59489760")
scaden_folder <- synStore(scaden_folder, forceVersion = FALSE)

# Provenance
sc_df <- GetChildrenAsDf("syn58807549")
bulk_df <- GetChildrenAsDf("syn58841972")
provenance_df <- rbind(sc_df, bulk_df)

# Loop over each model and upload
models <- list.files(dir_scaden_models, pattern = ".*_.*", full.names = TRUE)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step09_RunDeconvolutionAlgorithms.R"

for (model_folder_name in models) {
  datasets <- str_split(basename(model_folder_name), "_", simplify = TRUE)[1:2]
  provenance_sub <- subset(provenance_df, dataset %in% datasets)
  RecursiveUpload(model_folder_name,
                  parent_folder = scaden_folder,
                  provenance = list("used" = provenance_sub$id,
                                    "executed" = github))
}
