library(stringr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space

meta_folder <- Folder("01_metadata", parent = "syn68238853")
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

bulk_folder <- Folder("03_bulk", parent = "syn68238853")
bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

provenance <- list(
  "Mayo" = c("syn68239093.1", "syn66639062.1", "syn21544637.1"),
  "MSBB" = c("syn68239156.1", "syn66639063.2", "syn21544666.1"),
  "ROSMAP" = c("syn68239095.1", "syn66639064.2", "syn22283384.4",
               "syn22301603.4", "syn22314232.4", "syn25817663.4")
)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step03_RegressBulkData.R"

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
