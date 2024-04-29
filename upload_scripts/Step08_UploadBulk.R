source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
bulk_folder <- Folder("08_bulk", parent = "syn58802522")
bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

provenance <- list("Mayo" = c("syn58803962.1", "syn29855549.2", "syn23277389.6"),
                   "MSBB" = c("syn58803964.1", "syn29855570.2", "syn6101474.9"),
                   "ROSMAP" = c("syn58803965.1", "syn29855598.2", "syn3191087.11"))

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step08_RegressBulkData.R"

# Processed bulk datasets
for (dataset in datasets) {
  files <- list.files(dir_input, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, bulk_folder,
               list("used" = provenance[[dataset]],
                    "executed" = github))
  }
}
