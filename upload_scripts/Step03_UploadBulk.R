library(stringr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space

meta_folder <- Folder("01_metadata", parent = config::get("upload_synid"))
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

bulk_folder <- Folder("03_bulk", parent = config::get("upload_synid"))
bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

cfg <- config::get("step02_preprocess")

provenance <- list(
  "Mayo" = unlist(cfg$mayo),
  "MSBB" = unlist(cfg$msbb),
  "ROSMAP" = unlist(cfg$rosmap)
)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step03_RegressBulkData.R"

# Processed bulk datasets
for (dataset in names(provenance)) {
  # Each dataset has multiple RDS files, one per tissue
  all_rds <- list.files(dir_bulk, pattern = dataset, full.names = TRUE)

  lapply(all_rds, UploadFile,
         parent_folder = bulk_folder,
         provenance = list("used" = provenance[[dataset]],
                           "executed" = github))

  # Each dataset has multiple formula files, one per tissue
  all_forms <- list.files(dir_metadata,
                          pattern = str_glue("model_formulas_{dataset}"),
                          full.names = TRUE)

  lapply(all_forms, UploadFile,
         parent_folder = meta_folder,
         provenance = list("used" = provenance[[dataset]],
                           "executed" = github))
}
