library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
errors_folder <- Folder("11_error_calculations", parent = "syn58802522")
errors_folder <- synStore(errors_folder, forceVersion = FALSE)

# Get provenance IDs from algorithm output folder TODO
bulk_folders <- as.list(synGetChildren("syn52245555"))

# Folder structure on Synapse is <bulk_dataset>/<algorithm>/<params_file>
params_list_df <- lapply(bulk_folders, function(B) {
  algorithms <- as.list(synGetChildren(B$id))

  # Get the files under each algorithm folder
  params_files <- lapply(algorithms, function(A) {
    params_lists <- as.list(synGetChildren(A$id))
    if (length(params_lists) == 0) {
      return(NULL)
    }
    df <- data.frame(do.call(rbind, params_lists)) %>% select(name, id)
    df$test_dataset <- B$name
    return(df)
  })

  params_files <- do.call(rbind, params_files)
})

params_list_df <- do.call(rbind, params_list_df)

params_list_df$name_base <- str_replace(unlist(params_list_df$name), "estimates_", "")
params_list_df$name_base <- str_replace(params_list_df$name_base, ".rds", "")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step11_CalculateErrors.R"

# Loop over each data set and algorithm
for (test_dataset in c("Mayo", "MSBB", "ROSMAP")) { #unique(params_list_df$test_dataset)) {
  bulk_errors_folder <- Folder(test_dataset, parent = errors_folder)
  bulk_errors_folder <- synStore(bulk_errors_folder)

  algorithms <- list.files(file.path(dir_errors, test_dataset))

  for (algorithm in algorithms) {
    algorithm_errors_folder <- Folder(algorithm, parent = bulk_errors_folder)
    algorithm_errors_folder <- synStore(algorithm_errors_folder)

    dir_alg <- file.path(dir_errors, test_dataset, algorithm)

    files <- list.files(dir_alg, full.names = TRUE)

    # Parameter output lists
    for (filename in files) {
      # The name of the error file should match the name of the algorithm output file
      #error_name_base <- str_replace(filename, "errors_", "")
      #error_name_base <- str_replace(error_name_base, ".rds", "")

      #provenance <- subset(params_list_df, name_base == error_name_base)

      UploadFile(filename,
                 parent_folder = algorithm_errors_folder,
                 provenance = list("used" = c(), "executed" = github))
    }
  }
}
