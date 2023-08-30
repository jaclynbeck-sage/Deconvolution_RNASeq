# Helper script to get covariate data if the full Step02 was not run on this
# machine. This is necessary because while I upload pre-processed data sets
# to Synapse, these do not include the clinical metadata that may be protected
# health information, in order to ensure that I don't accidentally give access
# to data that otherwise requires a DUC. Covariates are only necessary for
# integrating single cell data sets and regressing bulk data sets.

source(file.path("functions", "Step02_Preprocess_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", #"morabito",
              "seaRef", #"seaAD", # Single cell
              "Mayo", "MSBB", "ROSMAP") # Bulk

for (dataset in datasets) {
  files <- DownloadData(dataset, metadata_only = TRUE)
  covariates <- ReadCovariates(dataset, files)
  saveRDS(covariates, file.path(dir_covariates,
                                str_glue("{dataset}_covariates.rds")))
}
