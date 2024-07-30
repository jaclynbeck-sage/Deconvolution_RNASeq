source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
input_folder <- Folder("04_single_cell", parent = "syn58802522")
input_folder <- synStore(input_folder, forceVersion = FALSE)

provenance <- list("cain" = "syn58803769.1",
                   "lau" = "syn58803777.1",
                   "leng" = "syn58803786.1",
                   "mathys" = "syn58803789.1",
                   "seaRef" = c("syn58807209.1", url_searef_h5))

maps <- c("syn58803447.1", "syn58803667.1")
github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step04_MapDatasets.R"

# Processed single cell datasets
for (dataset in names(provenance)) {
  files <- list.files(dir_input, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, input_folder,
               list("used" = c(provenance[[dataset]], maps),
                    executed = github))
  }
}
