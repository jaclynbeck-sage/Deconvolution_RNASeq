source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
input_folder <- Folder("02_preprocessed", parent = "syn68238853")
input_folder <- synStore(input_folder, forceVersion = FALSE)

provenance <- list(
  "cain" = c("syn23554294.5", "syn21323366.17", "syn3191087.11",
             "syn23554292.3", "syn23554293.2"),
  "lau" = c("syn52308080.1", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157827&format=file"),
  "leng" = c("syn22722817.1", "syn22722860.1"),
  "mathys" = c("syn3191087.11", "syn18686383.1", "syn18687958.1",
               "syn18687959.1"),
  "seaRef" = c("syn31149116.5", url_searef_h5),
  "Mayo" = c("syn66639062.1", "syn20827192.13", "syn20827193.4",
             "syn21544637.1", "syn21544635.1"),
  "MSBB" = c("syn66639063.2", "syn21893059.14", "syn22447899.6",
             "syn21544666.1", "syn21544664.1"),
  "ROSMAP" = c("syn66639064.2", "syn21323366.19", "syn21088596.5",
               "syn22283384.4", "syn22301603.4", "syn22314232.4",
               "syn25817663.4", "syn22283382.4", "syn22301601.4",
               "syn22314230.4", "syn25817661.4")
)

genes <- "syn68238911.1"
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



