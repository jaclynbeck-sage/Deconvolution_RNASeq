library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
markers_folder <- Folder("06_markers", parent = "syn68238853")
markers_folder <- synStore(markers_folder, forceVersion = FALSE)

meta_folder <- Folder("01_metadata", parent = "syn68238853")
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

# Get provenance IDs
sce_df <- GetChildrenAsDf("syn58807549")
pseudo_df <- GetChildrenAsDf("syn58808874")

colnames(sce_df) <- c("sce_name", "sce_id", "dataset")

pseudo_df <- subset(pseudo_df, grepl("puresamples", name))
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "granularity")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

datasets <- unique(provenance_df$dataset)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step06_FindMarkers.R"

# Markers
for (dset in datasets) {
  files <- list.files(dir_markers, pattern = dset, full.names = TRUE)

  for (filename in files) {
    broadfine <- if (grepl("broad", filename)) "broad_class" else "sub_class"
    provenance <- subset(provenance_df, dataset == dset & granularity == broadfine)

    if (grepl("pseudobulk", filename)) {
      provenance <- provenance$pseudobulk_id
    } else {
      provenance <- provenance$sce_id
    }

    UploadFile(filename, markers_folder,
               list("used" = provenance,
                    "executed" = github))
  }
}

# Excluded genes go in the metadata folder
for (dset in datasets) {
  filename <- file.path(dir_metadata, str_glue("{dset}_excluded_genes.rds"))
  if (file.exists(filename)) {
    provenance <- sce_df$id[sce_df$dataset == dset]

    UploadFile(filename, meta_folder,
               list("used" = provenance,
                    "executed" = github))
  }
}

