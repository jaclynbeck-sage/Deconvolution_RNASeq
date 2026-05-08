source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
input_folder <- Folder(basename(dir_preprocessed), parent = config::get("upload_synid"))
input_folder <- synStore(input_folder, forceVersion = FALSE)

provenance <- config::get("step02_preprocess")
provenance$lau$geo_accession <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=",
                                       provenance$lau$geo_accession)

# Unlist any lists of IDs inside each dataset item
provenance <- lapply(provenance, unlist)

genes <- config::get("gene_metadata_synid")
github <- paste0(config::get("github_repo_url"), "Step02_Preprocess_RNASeqData.R")

# Processed single cell datasets
for (dataset in names(provenance)) {
  files <- list.files(dir_preprocessed, pattern = dataset,
                      ignore.case = TRUE, full.names = TRUE)
  for (filename in files) {
    UploadFile(filename, input_folder,
               list("used" = c(as.character(provenance[[dataset]]), genes),
                    executed = github))
  }
}



