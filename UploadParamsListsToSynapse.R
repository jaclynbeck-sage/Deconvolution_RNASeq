library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# My private workspace
params.folder <- Folder("params_lists", parent = "syn49332774")
params.folder <- synStore(params.folder, forceVersion = FALSE)

markers.folder <- Folder("markers", parent = "syn49332774")
markers.folder <- synStore(markers.folder, forceVersion = FALSE)

# Get provenance IDs
sce_list <- as.list(synGetChildren("syn49340746")) # sce files
pseudo_list <- as.list(synGetChildren("syn49334026")) # pseudobulk files

sce_df <- data.frame(do.call(rbind, sce_list))
pseudo_df <- data.frame(do.call(rbind, pseudo_list))

sce_df <- select(sce_df, name, id)
sce_df$dataset <- str_replace(sce_df$name, "_.*", "")
colnames(sce_df) <- c("sce_name", "sce_id", "dataset")

pseudo_df <- select(pseudo_df, name, id)
pseudo_df <- subset(pseudo_df, !grepl("puresamples", name) & !grepl("finecelltypes", name))
pseudo_df$dataset <- str_split(pseudo_df$name, "_", simplify = TRUE)[,2]
pseudo_df$test.type <- str_split(pseudo_df$name, "_", simplify = TRUE)[,3]
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "test.type")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

datasets <- unique(provenance_df$dataset)

# Parameter output lists
for (dset in datasets) {
  for (test in c("donors", "training")) {
    patt <- paste(dset, test, sep = "_")
    files <- list.files(dir_params_lists, pattern = patt, full.names = TRUE)
    for (filename in files) {
      provenance <- subset(provenance_df, dataset == dset & test.type == test)
      file <- File(path = filename, parent = params.folder)
      file <- synStore(file,
                       used = c(provenance$sce_id[[1]], provenance$pseudobulk_id[[1]]),
                       forceVersion = FALSE)
      print(paste0(file$properties$id, ": ", file$properties$name))
    }
  }
}

# Markers
for (dset in datasets) {
  files <- list.files(dir_markers, pattern = dset, full.names = TRUE)
  provenance <- subset(sce_df, dataset == dset)

  for (filename in files) {
    file <- File(path = filename, parent = markers.folder)
    file <- synStore(file, used = provenance$sce_id[[1]], forceVersion = FALSE)
  }
}
