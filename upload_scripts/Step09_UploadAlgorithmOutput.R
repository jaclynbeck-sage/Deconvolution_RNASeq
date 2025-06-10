library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# TODO add markers to provenance?

# Deconvolution WG Synapse space
output_folder <- Folder("09_estimates", parent = "syn68238853")
output_folder <- synStore(output_folder, forceVersion = FALSE)

# Get provenance IDs
input_ids <- list(sce = "syn58807549",
                  pseudobulk = "syn58808874",
                  signatures = "syn59480278",
                  bulk = "syn58841972")

dfs <- lapply(input_ids, GetChildrenAsDf)

dfs$pseudobulk <- subset(dfs$pseudobulk, grepl("puresamples", name))
dfs$signatures <- subset(dfs$signatures, grepl("signature", name))

colnames(dfs$sce) <- c("sc_name", "sc_id", "dataset")
colnames(dfs$pseudobulk) <- c("pb_name", "pb_id", "dataset", "granularity")
colnames(dfs$signatures) <- c("sig_name", "sig_id", "dataset")

provenance_df <- merge(dfs$sce, dfs$pseudobulk, by = "dataset")
provenance_df <- merge(provenance_df, dfs$signatures, by = "dataset")

reference_datasets <- unique(provenance_df$dataset)
test_datasets <- unique(dfs$bulk$dataset)
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE",
                "Music", "Scaden", "Baseline")

# Parameter output lists
for (bulk in test_datasets) {
  bulk_folder <- Folder(bulk, parent = output_folder)
  bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

  dir_bulk <- file.path(dir_estimates, bulk)

  for (alg in algorithms) {
    alg_folder <- Folder(alg, parent = bulk_folder)
    alg_folder <- synStore(alg_folder, forceVersion = FALSE)

    files_alg <- list.files(file.path(dir_bulk, alg), full.names = TRUE)

    if (alg == "Baseline") {
      github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step10_GenerateBaselineData.R"
      for (filename in files_alg) {
        UploadFile(filename,
                   parent_folder = alg_folder,
                   provenance = list("used" = c(), "executed" = github))
      }
    } else {
      github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step09_RunDeconvolutionAlgorithms.R"

      for (ref in reference_datasets) {
        files <- grep(ref, files_alg, value = TRUE)

        for (filename in files) {
          broadfine <- if (grepl("broad", filename)) "broad_class" else "sub_class"
          provenance <- subset(provenance_df, dataset == ref &
                                 granularity == broadfine)
          provenance_bulk <- subset(dfs$bulk, dataset == bulk)

          if (alg == "CibersortX") {
            used_ids <- c(provenance$sig_id, provenance$sc_id)
          } else if (alg %in% c("DeconRNASeq", "DWLS")) {
            used_ids <- provenance$sig_id
          } else if (alg %in% c("Dtangle", "HSPE")) {
            used_ids <- provenance$pb_id
          } else { # Baseline, Music, and Scaden
            used_ids <- provenance$sc_id
          }

          UploadFile(filename,
                     parent_folder = alg_folder,
                     provenance = list("used" = c(used_ids, provenance_bulk$id),
                                       "executed" = github))
        }
      }
    }
  }
}
