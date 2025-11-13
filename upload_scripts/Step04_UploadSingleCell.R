source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
input_folder <- Folder(basename(dir_singlecell), parent = config::get("upload_synid"))
input_folder <- synStore(input_folder, forceVersion = FALSE)

cfg1 <- config::get("step02_preprocess")
cfg2 <- config::get("step04_label_transfer")

cfg1$lau$geo_accession <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=",
                                 cfg1$lau$geo_accession)

provenance <- list(
  cain = c(unlist(cfg1$cain), cfg2$cain_annotations),
  lau = c(unlist(cfg1$lau), cfg2$lau_annotations),
  mathys = c(unlist(cfg1$mathys), cfg2$mathys_annotations),
  seaRef = unlist(cfg1$searef) # no annotations for seaRef
)

github <- paste0(config::get("github_repo_url"), "Step04_TransferCellTypeLabels.R")

# Processed single cell datasets
for (dataset in names(provenance)) {
  files <- list.files(dir_singlecell, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, input_folder,
               list("used" = c(provenance[[dataset]]),
                    "executed" = github))
  }
}
