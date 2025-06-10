source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
sig_folder <- Folder("07_signatures", parent = "syn68238853")
sig_folder <- synStore(sig_folder, forceVersion = FALSE)

provenance <- GetChildrenAsDf("syn58807549")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step07_CalculateSignatures.R"

for (dataset in provenance$dataset) {
  files <- list.files(dir_signatures, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, sig_folder,
               list("used" = provenance$id[provenance$dataset == dataset],
                    "executed" = github))
  }
}
