library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# TODO add markers to provenance?

# Deconvolution WG Synapse space
output_folder <- Folder("08_estimates", parent = config::get("upload_synid"))
output_folder <- synStore(output_folder, forceVersion = FALSE)

# Get provenance IDs
main_folders <- GetMainFolderIds()

input_ids <- list(sce = main_folders[[basename(dir_singlecell)]]$id,
                  pseudobulk = main_folders[[basename(dir_pseudobulk)]]$id,
                  signatures = main_folders[[basename(dir_signatures)]]$id,
                  bulk = main_folders[[basename(dir_bulk)]]$id)

# Scaden pseudobulk files are a sub-folder of the main pseudobulk folder
input_ids$scaden <- GetChildrenAsDf(input_ids$pseudobulk, types = list("folder")) |>
  pull(id)

dfs <- lapply(input_ids, GetChildrenAsDf)

dfs$pseudobulk <- subset(dfs$pseudobulk, grepl("puresamples", name))
dfs$signatures <- subset(dfs$signatures, grepl("signature", name))
dfs$scaden$norm <- str_replace(dfs$scaden$name, ".*_", "") |>
  str_replace("\\.rds", "")

colnames(dfs$sce)[1:2] <- c("sc_name", "sc_id")
colnames(dfs$pseudobulk)[1:2] <- c("pb_name", "pb_id")
colnames(dfs$signatures)[1:2] <- c("sig_name", "sig_id")

provenance_df <- purrr::reduce(dfs[c("sce", "pseudobulk", "signatures")],
                               merge)

reference_datasets <- unique(provenance_df$dataset)
test_datasets <- unique(dfs$bulk$dataset)
algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music",
                "Scaden", "Baseline")

# Parameter output lists
for (bulk in test_datasets) {
  bulk_folder <- Folder(bulk, parent = output_folder)
  bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

  dir_bulk_dataset <- file.path(dir_estimates, bulk)

  for (alg in algorithms) {
    alg_folder <- Folder(alg, parent = bulk_folder)
    alg_folder <- synStore(alg_folder, forceVersion = FALSE)

    files_alg <- list.files(file.path(dir_bulk_dataset, alg), full.names = TRUE)

    if (alg == "Baseline") {
      github <- paste0(config::get("github_repo_url"), "Step09_GenerateBaselineData.R")
      for (filename in files_alg) {
        UploadFile(filename,
                   parent_folder = alg_folder,
                   provenance = list("used" = c(), "executed" = github))
      }
    } else {
      github <- paste0(config::get("github_repo_url"), "Step08_RunDeconvolutionAlgorithms.R")

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

          } else if (alg == "Dtangle") {
            used_ids <- provenance$pb_id

          } else if (alg == "Scaden") {
            prov_scaden <- subset(dfs$scaden, dataset == ref &
                                    granularity == broadfine)
            norm <- ifelse(grepl("tmm", filename), "tmm", "cpm")
            used_ids <- prov_scaden$id[prov_scaden$norm == norm]

          } else { # Baseline, Music
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
