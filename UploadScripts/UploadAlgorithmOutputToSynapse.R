library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# TODO add markers to provenance?

# Deconvolution WG Synapse space
output_folder <- Folder("algorithm_output", parent = "syn49332774")
output_folder <- synStore(output_folder, forceVersion = FALSE)

# Get provenance IDs
sce_list <- as.list(synGetChildren("syn52245349")) # sce files
pseudo_list <- as.list(synGetChildren("syn52245380")) # pseudobulk files

input_df <- data.frame(do.call(rbind, sce_list))
pseudo_df <- data.frame(do.call(rbind, pseudo_list))

input_df <- select(input_df, name, id)
input_df$dataset <- str_replace(input_df$name, "_.*", "")
colnames(input_df) <- c("se_name", "se_id", "dataset")
bulk_df <- subset(input_df, grepl("_se.rds", se_name))
sce_df <- subset(input_df, grepl("_sce.rds", se_name))

pseudo_df <- select(pseudo_df, name, id)
pseudo_df$dataset <- str_split(pseudo_df$name, "_", simplify = TRUE)[,2]
pseudo_df$granularity <- str_split(pseudo_df$name, "_", simplify = TRUE)[,4]
pseudo_df$granularity <- str_replace(pseudo_df$granularity, ".rds", "")
pseudo_df$datatype <- str_split(pseudo_df$name, "_", simplify = TRUE)[,3]
colnames(pseudo_df) <- c("pseudobulk_name", "pseudobulk_id", "dataset", "granularity", "datatype")
pseudo_df <- subset(pseudo_df, datatype == "puresamples")

provenance_df <- merge(sce_df, pseudo_df, by = "dataset")

reference_datasets <- unique(provenance_df$dataset)
test_datasets <- unique(bulk_df$dataset)

# Parameter output lists
for (bulk in test_datasets) {
  bulk_folder <- Folder(bulk, parent = output_folder)
  bulk_folder <- synStore(bulk_folder, forceVersion = FALSE)

  dir_bulk <- switch(bulk,
                     "Mayo" = dir_mayo_output,
                     "MSBB" = dir_msbb_output,
                     "ROSMAP" = dir_rosmap_output,
                     dir_params_lists)

  for (ref in reference_datasets) {
    files <- list.files(dir_bulk, pattern = ref, full.names = TRUE)

    for (filename in files) {
      broadfine <- str_split(filename, "_", simplify = TRUE)[5]
      provenance <- subset(provenance_df, dataset == ref &
                             granularity == broadfine)
      provenance_bulk <- subset(bulk_df, dataset == bulk)

      file <- File(path = filename, parent = bulk_folder)
      file <- synStore(file,
                       used = c(unlist(provenance$se_id),
                                unlist(provenance$pseudobulk_id),
                                unlist(provenance_bulk$se_id)),
                       forceVersion = FALSE)
      print(paste0(file$properties$id, ": ", file$properties$name))
    }
  }
}
