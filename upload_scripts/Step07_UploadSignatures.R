source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
sig_folder <- Folder("07_signatures", parent = "syn58802522")
sig_folder <- synStore(sig_folder, forceVersion = FALSE)

# TODO
provenance <- list("cain" = "",
                   "lau" = "",
                   "leng" = "",
                   "mathys" = "",
                   "seaRef" = "")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step07_CalculateSignatures.R"

for (dataset in names(provenance)) {
  files <- list.files(dir_signatures, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, sig_folder,
               list("used" = provenance[[dataset]],
                    "executed" = github))
  }
}
