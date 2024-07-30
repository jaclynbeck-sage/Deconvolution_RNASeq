library(stringr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space

meta_folder <- Folder("01_metadata", parent = "syn58802522")
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

bulk_folder <- Folder("08_bulk", parent = "syn58802522")
bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

provenance <- list("Mayo" = c("syn58803962.1", "syn29855549.2", "syn23277389.6"),
                   "MSBB" = c("syn58803964.1", "syn29855570.2", "syn6101474.9"),
                   "ROSMAP" = c("syn58803965.1", "syn29855598.2", "syn3191087.11"))

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step08_RegressBulkData.R"

# Processed bulk datasets
for (dataset in names(provenance)) {
  filename <- file.path(dir_input, str_glue("{dataset}_se.rds"))
  UploadFile(filename, bulk_folder,
             list("used" = provenance[[dataset]],
                  "executed" = github))

  filename <- file.path(dir_metadata, str_glue("model_formulas_{dataset}.rds"))
  UploadFile(filename, meta_folder,
             list("used" = provenance[[dataset]],
                  "executed" = github))
}
