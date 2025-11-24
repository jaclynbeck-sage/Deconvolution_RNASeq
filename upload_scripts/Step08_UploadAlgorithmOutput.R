library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# TODO add markers to provenance?

# Deconvolution WG Synapse space
output_folder <- Folder("08_estimates", parent = config::get("upload_synid"))
output_folder <- synStore(output_folder, forceVersion = FALSE)

# Get provenance IDs
main_folders <- synGetChildren(config::get("upload_synid"))$asList()
main_folders <- do.call(rbind, main_folders) |> as.data.frame()
singlecell_folder <- main_folders$id[main_folders$name == basename(dir_singlecell)] |> unlist()
pseudo_folder <- main_folders$id[main_folders$name == basename(dir_pseudobulk)] |> unlist()
sig_folder <- main_folders$id[main_folders$name == basename(dir_signatures)] |> unlist()
bulk_folder <- main_folders$id[main_folders$name == basename(dir_bulk)] |> unlist()

input_ids <- list(sce = singlecell_folder,
                  pseudobulk = pseudo_folder,
                  signatures = sig_folder,
                  bulk = bulk_folder)

dfs <- lapply(input_ids, GetChildrenAsDf)

dfs$pseudobulk <- subset(dfs$pseudobulk, grepl("puresamples", name))
dfs$signatures <- subset(dfs$signatures, grepl("signature", name))
dfs$bulk$dataset <- str_replace(dfs$bulk$name, "_se.rds", "")

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
      github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step09_GenerateBaselineData.R"
      for (filename in files_alg) {
        UploadFile(filename,
                   parent_folder = alg_folder,
                   provenance = list("used" = c(), "executed" = github))
      }
    } else {
      github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step08_RunDeconvolutionAlgorithms.R"

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
