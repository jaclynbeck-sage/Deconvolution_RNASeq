library(synapser)

source("Filenames.R")

# First run only -- run with appropriate values
# For generating an auth token, see: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
#synLogin(<username>, authToken = <auth token>, rememberMe=True)

synLogin()

folders <- list("singlecell" = list("synID" = "syn52245349", "downloadLocation" = dir_input),
                "pseudobulk" = list("synID" = "syn52245380", "downloadLocation" = dir_pseudobulk),
                "markers" = list("synID" = "syn52245416", "downloadLocation" = dir_markers),
                "errors" = list("synID" = "syn52245653", "downloadLocation" = dir_errors),
                "external_metadata" = list("synID" = "syn52539337", "downloadLocation" = dir_metadata),
                "Mayo" = list("synID" = "syn52245606", "downloadLocation" = dir_mayo_output),
                "MSBB" = list("synID" = "syn52245570", "downloadLocation" = dir_msbb_output),
                "ROSMAP" = list("synID" = "syn52245624", "downloadLocation" = dir_rosmap_output))

# Download everything in each folder
for (folder in folders) {
  parent <- synGetChildren(folder$synID)
  parent <- parent$asList()

  for (child in parent) {
    synGet(child$id, downloadLocation = folder$downloadLocation)
  }
}

rm(list=ls())
gc()
