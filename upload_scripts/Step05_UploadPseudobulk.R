source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
pseudo_folder <- Folder("05_pseudobulk", parent = "syn58802522")
pseudo_folder <- synStore(pseudo_folder, forceVersion = FALSE)

provenance <- list("cain" = "syn58807553.1",
                   "lau" = "syn58807566.1",
                   "leng" = "syn58807572.1",
                   "mathys" = "syn58808591.1",
                   "seaRef" = "syn58808855.1")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step05_GeneratePseudobulk.R"

for (dataset in names(provenance)) {
  files <- list.files(dir_pseudobulk, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, pseudo_folder,
               list("used" = provenance[[dataset]],
                    "executed" = github))
  }
}
