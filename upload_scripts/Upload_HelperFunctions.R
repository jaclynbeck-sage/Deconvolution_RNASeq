library(dplyr)
library(stringr)
library(synapser)
synLogin()

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


GetChildrenAsDf <- function(syn_id) {
  syn_list <- as.list(synGetChildren(syn_id))
  syn_df <- data.frame(do.call(rbind, syn_list)) %>%
      select(name, id) %>%
      mutate(name = unlist(name),
             id = unlist(id))

  syn_df$dataset <- str_replace(str_replace(syn_df$name, "pseudobulk_", ""),
                                "_.*", "")

  if (any(grepl("[broad|sub]_class", syn_df$name))) {
    syn_df$granularity <- "broad_class"
    syn_df$granularity[grepl("sub", syn_df$name)] <- "sub_class"
  }

  return(syn_df)
}
