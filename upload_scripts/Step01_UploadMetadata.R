source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
meta_folder <- Folder("01_metadata", parent = config::get("upload_synid"))
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

cfg <- config::get("step01_gene_metadata")

provenance <- list(
  "genes" = c(
    cfg$gtf_bulk,
    cfg$fasta_bulk,
    cfg$gtf_cain,
    cfg$gtf_sea_ad,
    cfg$genes_mathys,
    paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", cfg$geo_lau)
  ),
  "ihc" = cfg$git_ihc)

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step01_Preprocess_ExternalMetadata.R"

UploadFile(file_gene_list, meta_folder,
           list("used" = provenance$genes,
                "executed" = github))

UploadFile(file_rosmap_ihc_proportions, meta_folder,
           list("used" = provenance$ihc,
                "executed" = github))
