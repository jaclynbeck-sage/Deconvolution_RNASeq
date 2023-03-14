library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
errors.folder <- Folder("error_calculations", parent = "syn51119398")
errors.folder <- synStore(errors.folder, forceVersion = FALSE)

# Get provenance IDs
params_lists <- as.list(synGetChildren("syn51124454")) # params lists

params_list_df <- data.frame(do.call(rbind, params_lists))
params_list_df <- subset(params_list_df, !grepl("delete", params_list_df$name))

fields <- str_split(params_list_df$name, "_", simplify = TRUE)
params_list_df$dataset <- fields[,3]
params_list_df$algorithm <- fields[,1]
params_list_df$datatype <- fields[,4]
params_list_df$granularity <- str_replace(fields[,5], ".rds", "")

files <- list.files(dir_errors)

# Parameter output lists
for (filename in files) {
  fields <- str_split(filename, "_", simplify = TRUE)
  provenance <- subset(params_list_df, dataset == fields[3] &
                         datatype == fields[4] &
                         granularity == str_replace(fields[5], ".rds", "") &
                         algorithm == fields[2])

  file <- File(path = file.path(dir_errors, filename), parent = errors.folder)
  file <- synStore(file,
                   used = unlist(provenance$id),
                   forceVersion = FALSE)
  print(paste0(file$properties$id, ": ", file$properties$name))
}
