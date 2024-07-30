source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
input_folder <- Folder("02_preprocessed", parent = "syn58802522")
input_folder <- synStore(input_folder, forceVersion = FALSE)

provenance <- list("cain" = c("syn23554294.5", "syn21323366.17", "syn3191087.11",
                              "syn23554292.3", "syn23554293.2"),
                   "lau" = c("syn52308080.1", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157827&format=file"),
                   "leng" = c("syn22722817.1", "syn22722860.1"),
                   "mathys" = c("syn3191087.11", "syn18686372.1",
                                "syn18686381.1", "syn18686382.1"),
                   "seaRef" = c("syn31149116.5", url_searef_h5),
                   "Mayo" = c("syn29855549.2", "syn23277389.6", "syn27024951.1"),
                   "MSBB" = c("syn29855570.2", "syn6101474.9", "syn27068754.1"),
                   "ROSMAP" = c("syn29855598.2", "syn3191087.11", "syn26967451.1"))

genes <- "syn58803308.1"
github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step02_Preprocess_RNASeqData.R"

# Processed single cell datasets
for (dataset in names(provenance)) {
  files <- list.files(dir_preprocessed, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, input_folder,
               list("used" = c(provenance[[dataset]], genes),
                    executed = github))
  }
}



