source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
map_folder <- Folder("03_map_reference", parent = "syn58802522")
map_folder <- synStore(map_folder, forceVersion = FALSE)

github_create_map <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step03_CreateMappingReference.R"

provenance <- list("used" = "syn58803769.1", # Cain pre-processed data
                   "executed" = github_create_map)

files <- list.files(dir_map_reference, full.names = TRUE)
for (filename in files) {
  UploadFile(filename, map_folder, provenance)
}
