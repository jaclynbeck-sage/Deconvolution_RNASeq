source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
input_folder <- Folder("04_single_cell", parent = "syn68238853")
input_folder <- synStore(input_folder, forceVersion = FALSE)

provenance <- list("cain" = c("syn68239073.1", "syn68239068.1"),
                   "lau" = c("syn68239080.1", "syn68239074.1"),
                   "leng" = c("syn68239086.1", "syn68239081.1"),
                   "mathys" = c("syn68239089.1", "syn68239087.1"),
                   "seaRef" = c("syn68239090.1", url_searef_h5))

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step04_TransferCellTypeLabels.R"

# Processed single cell datasets
for (dataset in names(provenance)) {
  files <- list.files(dir_input, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, input_folder,
               list("used" = c(provenance[[dataset]]),
                    executed = github))
  }
}
