library(synapser)

# First run only -- run with appropriate values
# For generating an auth token, see: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
#synLogin(<username>, authToken = <auth token>, rememberMe=True)

synLogin()

synapseDir <- file.path("DeconvolutionData", "input")
dir.create(synapseDir, recursive = TRUE)

# Download everything in Starting_data_and_models_20220930
# TODO download original data sets from original Synapse IDs
parent <- synGetChildren("syn38559392")
parent <- parent$asList()

for (child in parent) {
  synGet(child$id, downloadLocation = synapseDir)
}

rm(list=ls())
gc()


