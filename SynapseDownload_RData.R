library(synapser)

source("Filenames.R")

# First run only -- run with appropriate values
# For generating an auth token, see: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
#synLogin(<username>, authToken = <auth token>, rememberMe=True)

synLogin(silent = TRUE)

# List of folders to download to disk. Names on Synapse mirror the folder names
# locally.
download_folders <- c(dir_metadata, dir_bulk, dir_singlecell, dir_pseudobulk,
                      dir_signatures, dir_markers, dir_estimates, dir_errors,
                      dir_top_parameters, dir_top_estimates, dir_top_errors)
# uncomment only if these files are needed, they take a lot of disk space
#dir_music_basis, dir_scaden_models)
names(download_folders) <- basename(download_folders)

# Get which of these folders are available on Synapse
synapse_folders <- synGetChildren(config::get("upload_synid"))$asList()
names(synapse_folders) <- sapply(synapse_folders, "[[", "name")

synapse_folders <- synapse_folders[intersect(names(synapse_folders), names(download_folders))]


# Helper function for recursive download
download <- function(synID, downloadLocation) {
  parent <- synGetChildren(synID)
  parent <- parent$asList()

  for (child in parent) {
    res <- synGet(child$id, downloadLocation = downloadLocation,
                  ifcollision = "overwrite.local")

    if (is(res, "synapseclient.entity.Folder")) {
      new_folder <- file.path(downloadLocation, res$name)
      dir.create(new_folder, recursive = TRUE, showWarnings = FALSE)
      download(res$id, new_folder)
    }
    else {
      print(file.path(downloadLocation, res$name))
    }
  }
}

# Download everything in each folder
for (folder in synapse_folders) {
  download(folder$id, download_folders[folder$name])
}
