library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
scaden_folder <- Folder("scaden_models", parent = "syn58802522")
scaden_folder <- synStore(scaden_folder, forceVersion = FALSE)

# Get provenance IDs from input data sets TODO
input_data <- as.list(synGetChildren("syn52569484"))
provenance_df <- do.call(rbind, lapply(input_data, as.data.frame))
provenance_df <- provenance_df[grepl("_sc?e.rds", provenance_df$name),] %>%
                  select(name, id)
provenance_df$dataset <- str_replace(provenance_df$name, "_.*", "")

# Loop over each model and upload
models <- list.files(dir_scaden_models, pattern = ".*_.*", full.names = TRUE)

recursive_upload <- function(folder_name, parent_folder, provenance) {
  folder <- Folder(basename(folder_name), parent = parent_folder)
  folder <- synStore(folder)

  files <- list.files(folder_name, full.names = TRUE, include.dirs = FALSE)
  for (filename in files) {
    # Recursive upload directories
    if (dir.exists(filename)) {
      recursive_upload(filename, folder, provenance)
    }
    # Directly upload files that are not directories
    else {
      file <- File(path = filename, parent = folder)
      file <- synStore(file,
                       used = unlist(provenance$id),
                       forceVersion = FALSE)
      print(paste0(file$properties$id, ": ", filename))
    }
  }
}

for (model_folder_name in models) {
  datasets <- str_split(basename(model_folder_name), "_", simplify = TRUE)[1:2]
  provenance <- subset(provenance_df, dataset %in% datasets)
  recursive_upload(model_folder_name, scaden_folder, provenance)
}
