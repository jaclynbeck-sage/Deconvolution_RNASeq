library(synapser)
library(stringr)
library(dplyr)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
errors.folder <- Folder("error_calculations", parent = "syn49332774")
errors.folder <- synStore(errors.folder, forceVersion = FALSE)

# Get provenance IDs
bulk_folders <- as.list(synGetChildren("syn52245555"))
params_lists <- sapply(bulk_folders, function(X) {
  as.list(synGetChildren(X$id))
})

params_list_df <- data.frame(do.call(rbind, params_lists)) %>%
  select(name, id)

fields <- str_split(params_list_df$name, "_", simplify = TRUE)
params_list_df$reference_dataset <- fields[,3]
params_list_df$algorithm <- fields[,2]
params_list_df$test_dataset <- fields[,4]
params_list_df$granularity <- fields[,5]

files <- list.files(dir_errors)

# Parameter output lists
for (filename in files) {
  fields <- unlist(str_split(filename, "_"))
  provenance <- subset(params_list_df, reference_dataset == fields[3] &
                         test_dataset == fields[4] &
                         granularity == str_replace(fields[5], "\\.rds", "") &
                         algorithm == fields[2])

  file <- File(path = file.path(dir_errors, filename), parent = errors.folder)
  file <- synStore(file,
                   used = unlist(provenance$id),
                   forceVersion = FALSE)
  print(paste0(file$properties$id, ": ", file$properties$name))
}
