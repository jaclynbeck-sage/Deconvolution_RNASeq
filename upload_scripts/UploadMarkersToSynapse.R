library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
markers.folder <- Folder("markers", parent = "syn49332774")
markers.folder <- synStore(markers.folder, forceVersion = FALSE)

# Get provenance IDs
sce_list <- as.list(synGetChildren("syn52569484")) # sce files
pseudo_list <- as.list(synGetChildren("syn52570274")) # pseudobulk files

sce_df <- data.frame(do.call(rbind, sce_list))
pseudo_df <- data.frame(do.call(rbind, pseudo_list))

sce_df <- select(sce_df, name, id)
sce_df$dataset <- str_replace(sce_df$name, "_.*", "")
colnames(sce_df) <- c("sce_name", "sce_id", "dataset")
sce_df <- subset(sce_df, grepl("_sce.rds", sce_name))

pseudo_df <- select(pseudo_df, name, id)
pseudo_df <- subset(pseudo_df, grepl("puresamples", name))
pseudo_df$dataset <- str_split(pseudo_df$name, "_", simplify = TRUE)[,2]
pseudo_df$granularity <- str_split(pseudo_df$name, "_", simplify = TRUE)[,4]
pseudo_df$granularity <- paste0(pseudo_df$granularity, "_class")
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "granularity")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

datasets <- unique(provenance_df$dataset)

# Markers
for (dset in datasets) {
  files <- list.files(dir_markers, pattern = dset, full.names = TRUE)

  for (filename in files) {
    broadfine <- str_split(filename, "_", simplify = TRUE)[4]
    broadfine <- paste0(broadfine, "_class")
    input_type <- str_split(filename, "_", simplify = TRUE)[7]
    provenance <- subset(provenance_df, dataset == dset & granularity == broadfine)

    if (is.na(input_type) | input_type == "singlecell") {
      provenance <- unique(provenance$sce_id[[1]])
    }
    else {
      provenance <- unique(provenance$pseudobulk_id[[1]])
    }

    file <- File(path = filename, parent = markers.folder)
    file <- synStore(file, used = provenance, forceVersion = FALSE)
    print(paste0(file$properties$id, ": ", file$properties$name))
  }
}

