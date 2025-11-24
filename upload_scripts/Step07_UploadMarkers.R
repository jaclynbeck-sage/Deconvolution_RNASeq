library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
markers_folder <- Folder(basename(dir_markers), parent = config::get("upload_synid"))
markers_folder <- synStore(markers_folder, forceVersion = FALSE)

meta_folder <- Folder(basename(dir_metadata), parent = config::get("upload_synid"))
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

# Get provenance IDs
main_folders <- GetMainFolderIds()

sce_df <- GetChildrenAsDf(main_folders[[basename(dir_singlecell)]]$id)
pseudo_df <- GetChildrenAsDf(main_folders[[basename(dir_pseudobulk)]]$id)

colnames(sce_df) <- c("sce_name", "sce_id", "dataset")

pseudo_df <- subset(pseudo_df, grepl("puresamples", name))
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "granularity")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

datasets <- unique(provenance_df$dataset)

github <- paste0(config::get("github_repo_url"), "Step07_FindMarkers.R")

for (dset in datasets) {
  # Marker files
  files <- list.files(dir_markers, pattern = dset, full.names = TRUE)

  for (filename in files) {
    broadfine <- if (grepl("broad", filename)) "broad_class" else "sub_class"
    provenance <- subset(provenance_df, dataset == dset & granularity == broadfine)

    if (grepl("deseq2", filename) || grepl("dtangle", filename)) {
      provenance <- provenance$pseudobulk_id
    } else {
      provenance <- provenance$sce_id
    }

    UploadFile(filename, markers_folder,
               list("used" = provenance,
                    "executed" = github))
  }

  # Excluded genes go in the metadata folder
  filename <- file.path(dir_metadata, str_glue("{dset}_excluded_genes.rds"))
  if (file.exists(filename)) {
    provenance <- sce_df$id[sce_df$dataset == dset]

    UploadFile(filename, meta_folder,
               list("used" = provenance,
                    "executed" = github))
  }
}
