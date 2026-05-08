library(dplyr)
library(stringr)
library(synapser)
source(file.path("functions", "General_HelperFunctions.R"))

synLogin(silent = TRUE)

UploadFile <- function(filename, parent_folder, provenance) {
  if (!(file.info(filename)$isdir)) { # Exclude directories
    file <- File(path = filename, parent = parent_folder)
    file <- synStore(file,
                     used = provenance$used,
                     executed = provenance$executed,
                     forceVersion = FALSE)
    print(paste0(file$properties$id, ": ", file$properties$name))
  }
}


RecursiveUpload <- function(folder_name, parent_folder, provenance) {
  folder <- Folder(basename(folder_name), parent = parent_folder)
  folder <- synStore(folder)

  files <- list.files(folder_name, full.names = TRUE, include.dirs = FALSE)
  for (filename in files) {
    # Recursive upload directories
    if (dir.exists(filename)) {
      RecursiveUpload(filename, folder, provenance)
    }
    # Directly upload files that are not directories
    else {
      UploadFile(filename,
                 parent_folder = folder,
                 provenance = provenance)
    }
  }
}

ExtractDatasetName <- function(filename) {
  all_datasets <- c(all_singlecell_datasets(), all_bulk_datasets(),
                    "random_biased", "random_educated", "random_uniform", "zeros")

  # Can return multiple dataset names (i.e. bulk + single cell)
  names <- all_datasets[sapply(all_datasets, grepl, x = filename)]

  if (length(names) == 0) {
    return(NA)
  } else if (length(names) > 1) {
    return(list(names))
  }

  return(names)
}

GetChildrenAsDf <- function(syn_id, types = list("file")) {
  syn_list <- as.list(synGetChildren(syn_id, includeTypes = types))
  if (length(syn_list) == 0) {
    return(NULL)
  }

  syn_df <- data.frame(do.call(rbind, syn_list)) |>
    select(name, id) |>
    mutate(across(everything(), unlist)) |>
    # Find which dataset and granularity (if applicable) this is
    rowwise() |>
    mutate(
      dataset = ExtractDatasetName(name),
      granularity = case_when(
        grepl("broad_class", name) ~ "broad_class",
        grepl("sub_class", name) ~ "sub_class",
        .default = NA
      )
    ) |>
    ungroup()

  if (all(is.na(syn_df$granularity))) {
    syn_df <- select(syn_df, -granularity)
  }

  return(syn_df)
}

GetMainFolderIds <- function() {
  main_folders <- synGetChildren(config::get("upload_synid"))$asList()
  names(main_folders) <- lapply(main_folders, "[[", "name")
  return(main_folders)
}
