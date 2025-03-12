library(synapser)

source("Filenames.R")

# First run only -- run with appropriate values
# For generating an auth token, see: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
#synLogin(<username>, authToken = <auth token>, rememberMe=True)

synLogin()

folders <- list("metadata" = list("synID" = "syn58803307", "downloadLocation" = dir_metadata),
                "singlecell" = list("synID" = "syn58807549", "downloadLocation" = dir_input),
                "pseudobulk" = list("synID" = "syn58808874", "downloadLocation" = dir_pseudobulk),
                "bulk" = list("synID" = "syn58841972", "downloadLocation" = dir_input),
                "markers" = list("synID" = "syn58842534", "downloadLocation" = dir_markers),
                "signatures" = list("synID" = "syn59480278", "downloadLocation" = dir_signatures),
                "estimates" = list("synID" = "syn59489760", "downloadLocation" = dir_estimates),
                "errors" = list("synID" = "syn60969684", "downloadLocation" = dir_errors),
                "top_parameters" = list("synID" = "syn61917690", "downloadLocation" = dir_top_parameters),
                "top_estimates" = list("synID" = "syn61922322", "downloadLocation" = dir_top_estimates),
                "top_errors" = list("synID" = "syn63660315", "downloadLocation" = dir_best_errors)
                # uncomment only if these files are needed
                #"hspe_params" = list("synID" = "syn61709305", "downloadLocation" = dir_hspe_params),
                #"music_basis" = list("synID" = "syn59500238", "downloadLocation" = dir_music_basis),
                # This takes a lot of disk space
                #"scaden_models" = list("synID" = "syn63151233", "downloadLocation" = dir_scaden_models)
)


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
for (folder in folders) {
  download(folder$synID, folder$downloadLocation)
}

rm(list=ls())
gc()
