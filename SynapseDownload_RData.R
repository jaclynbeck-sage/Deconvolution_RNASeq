library(synapser)

source("Filenames.R")

# First run only -- run with appropriate values
# For generating an auth token, see: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
#synLogin(<username>, authToken = <auth token>, rememberMe=True)

synLogin()

folders <- list("singlecell" = list("synID" = "syn51119405", "downloadLocation" = dir_input),
                "pseudobulk" = list("synID" = "syn51119989", "downloadLocation" = dir_pseudobulk),
                "markers" = list("synID" = "syn51120961", "downloadLocation" = dir_markers),
                "params" = list("synID" = "syn51124454", "downloadLocation" = dir_params_lists),
                "errors" = list("synID" = "syn51188500", "downloadLocation" = dir_errors))

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
