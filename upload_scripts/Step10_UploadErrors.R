library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
errors_folder <- Folder(basename(dir_errors), parent = config::get("upload_synid"))
errors_folder <- synStore(errors_folder, forceVersion = FALSE)

github <- paste0(config::get("github_repo_url"), "Step10_CalculateErrors.R")

# Get provenance IDs
main_folders <- GetMainFolderIds()

input_ids <- list(estimates = main_folders[[basename(dir_estimates)]]$id,
                  bulk = main_folders[[basename(dir_bulk)]]$id)

bulk_datasets <- GetChildrenAsDf(input_ids$bulk) |>
  dplyr::rename(data_id = id) |>
  select(-name)
bulk_folders <- GetChildrenAsDf(input_ids$estimates, types = list("folder")) |>
  merge(bulk_datasets)

# Folder structure on Synapse is <bulk_dataset>/<algorithm>/<params_file>
for (B in 1:nrow(bulk_folders)) {
  b_folder <- bulk_folders[B, ]
  algorithms <- as.list(synGetChildren(b_folder$id))

  provenance_df <- lapply(algorithms, function(A) {
    # Helper function to figure out which of the two dataset names in each
    # "dataset" column entry belongs to reference or test
    split_dataset_names <- function(d_names, d_type) {
      if (d_type == "reference") {
        possible_datasets <- c(all_singlecell_datasets(), "random_biased",
                               "random_educated", "random_uniform", "zeros")
      } else {
        possible_datasets <- all_bulk_datasets()
      }

      d_names[which(d_names %in% possible_datasets)]
    }

    df <- GetChildrenAsDf(A$id)
    if (!is.null(df)) {
      df <- df |>
        mutate(reference_dataset = sapply(dataset, split_dataset_names, d_type = "reference"),
               test_dataset = sapply(dataset, split_dataset_names, d_type = "test")) |>
        select(-dataset)
    }
    return(df)
  })

  provenance_df <- do.call(rbind, provenance_df)

  provenance_df$name_base <- str_replace(provenance_df$name, "estimates_", "") |>
    str_replace(".rds", "")

  bulk_errors_folder <- Folder(b_folder$name, parent = errors_folder)
  bulk_errors_folder <- synStore(bulk_errors_folder)

  algorithms <- list.files(file.path(dir_errors, b_folder$dataset))

  for (algorithm in algorithms) {
    algorithm_errors_folder <- Folder(algorithm, parent = bulk_errors_folder)
    algorithm_errors_folder <- synStore(algorithm_errors_folder)

    dir_alg <- file.path(dir_errors, b_folder$name, algorithm)

    files <- list.files(dir_alg, full.names = TRUE)

    # Parameter output lists
    for (filename in files) {
      # The name of the error file should match the name of the algorithm output file
      error_name_base <- str_replace(basename(filename), "errors_", "") |>
        str_replace(".rds", "")

      provenance <- subset(provenance_df, name_base == error_name_base)

      UploadFile(filename,
                 parent_folder = algorithm_errors_folder,
                 provenance = list("used" = c(provenance$id, b_folder$data_id),
                                   "executed" = github))
    }
  }
}
