source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
pseudo_folder <- Folder("05_pseudobulk", parent = "syn68238853")
pseudo_folder <- synStore(pseudo_folder, forceVersion = FALSE)

provenance <- list("cain" = "syn68239290.1",
                   "lau" = "syn68239291.1",
                   "leng" = "syn68239293.1",
                   "mathys" = "syn68239294.1",
                   "seaRef" = "syn68239299.1")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step05_GeneratePseudobulk.R"

for (dataset in names(provenance)) {
  files <- list.files(dir_pseudobulk, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, pseudo_folder,
               list("used" = provenance[[dataset]],
                    "executed" = github))
  }
}
