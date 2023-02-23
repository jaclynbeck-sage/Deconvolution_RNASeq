library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
params.folder <- Folder("params_lists", parent = "syn51119398")
params.folder <- synStore(params.folder, forceVersion = FALSE)

# Get provenance IDs
sce_list <- as.list(synGetChildren("syn51119405")) # sce files
pseudo_list <- as.list(synGetChildren("syn51119989")) # pseudobulk files

sce_df <- data.frame(do.call(rbind, sce_list))
pseudo_df <- data.frame(do.call(rbind, pseudo_list))

sce_df <- select(sce_df, name, id)
sce_df$dataset <- str_replace(sce_df$name, "_.*", "")
colnames(sce_df) <- c("sce_name", "sce_id", "dataset")
sce_df <- subset(sce_df, grepl("_sce.rds", sce_name))

pseudo_df <- select(pseudo_df, name, id)
pseudo_df$dataset <- str_split(pseudo_df$name, "_", simplify = TRUE)[,2]
pseudo_df$granularity <- str_split(pseudo_df$name, "_", simplify = TRUE)[,4]
pseudo_df$granularity <- str_replace(pseudo_df$granularity, ".rds", "")
pseudo_df$datatype <- str_split(pseudo_df$name, "_", simplify = TRUE)[,3]
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "granularity", "datatype")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

datasets <- unique(provenance_df$dataset)

# Parameter output lists
for (dset in datasets) {
  files <- list.files(dir_params_lists, pattern = dset, full.names = TRUE)

  for (filename in files) {
    broadfine <- str_split(filename, "_", simplify = TRUE)[6]
    broadfine <- str_replace(broadfine, ".rds", "")
    dtype <- str_split(filename, "_", simplify = TRUE)[5]
    provenance <- subset(provenance_df, dataset == dset &
                           granularity == broadfine &
                           datatype %in% c(dtype, "puresamples"))

    file <- File(path = filename, parent = params.folder)
    file <- synStore(file,
                     used = c(provenance$sce_id[[1]], unlist(provenance$pseudobulk_id)),
                     forceVersion = FALSE)
    print(paste0(file$properties$id, ": ", file$properties$name))
  }
}
