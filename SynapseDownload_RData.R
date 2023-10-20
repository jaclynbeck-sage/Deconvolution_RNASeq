library(synapser)

source("Filenames.R")

# First run only -- run with appropriate values
# For generating an auth token, see: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
#synLogin(<username>, authToken = <auth token>, rememberMe=True)

synLogin()

folders <- list("singlecell" = list("synID" = "syn52569484", "downloadLocation" = dir_input),
                "pseudobulk" = list("synID" = "syn52570274", "downloadLocation" = dir_pseudobulk),
                "markers" = list("synID" = "syn52570296", "downloadLocation" = dir_markers),
                "errors" = list("synID" = "syn52245653", "downloadLocation" = dir_errors),
                "external_metadata" = list("synID" = "syn52539337", "downloadLocation" = dir_metadata),
                "Mayo" = list("synID" = "syn52587834", "downloadLocation" = dir_mayo_output),
                "MSBB" = list("synID" = "syn52587544", "downloadLocation" = dir_msbb_output),
                "ROSMAP" = list("synID" = "syn52587984", "downloadLocation" = dir_rosmap_output))


# Helper function for recursive download
download <- function(synID, downloadLocation) {
  parent <- synGetChildren(synID)
  parent <- parent$asList()

  for (child in parent) {
    res = synGet(child$id, downloadLocation = downloadLocation,
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
