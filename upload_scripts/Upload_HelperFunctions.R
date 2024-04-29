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
