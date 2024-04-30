library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
markers_folder <- Folder("06_markers", parent = "syn58802522")
markers_folder <- synStore(markers_folder, forceVersion = FALSE)

# Get provenance IDs
sce_list <- as.list(synGetChildren("syn58807549"))
pseudo_list <- as.list(synGetChildren("syn58808874"))

sce_df <- data.frame(do.call(rbind, sce_list))
pseudo_df <- data.frame(do.call(rbind, pseudo_list))

sce_df <- select(sce_df, name, id)
sce_df$dataset <- str_replace(sce_df$name, "_.*", "")
colnames(sce_df) <- c("sce_name", "sce_id", "dataset")

pseudo_df <- select(pseudo_df, name, id)
pseudo_df <- subset(pseudo_df, grepl("puresamples", name))
pseudo_df$dataset <- str_split(pseudo_df$name, "_", simplify = TRUE)[,2]
pseudo_df$granularity <- str_split(pseudo_df$name, "_", simplify = TRUE)[,4]
pseudo_df$granularity <- paste0(pseudo_df$granularity, "_class")
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "granularity")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

datasets <- unique(provenance_df$dataset)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step06_FindMarkers.R"

# Markers
for (dset in datasets) {
  files <- list.files(dir_markers, pattern = dset, full.names = TRUE)

  for (filename in files) {
    broadfine <- str_split(filename, "_", simplify = TRUE)[4]
    broadfine <- paste0(broadfine, "_class")
    provenance <- subset(provenance_df, dataset == dset & granularity == broadfine)

    if (grepl("pseudobulk", filename)) {
      provenance <- unique(provenance$pseudobulk_id[[1]])
    } else {
      provenance <- unique(provenance$sce_id[[1]])
    }

    UploadFile(filename, markers_folder,
               list("used" = provenance,
                    "executed" = github))
  }
}

