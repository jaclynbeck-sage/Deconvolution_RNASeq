source("Filenames.R")
source(file.path("upload_scripts", "Upload_HelperFunctions.R"))

# Deconvolution WG Synapse space
meta_folder <- Folder("01_metadata", parent = "syn58802522")
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

provenance <- list("genes" = c("syn26967452.1", # ROSMAP, Mayo, MSBB gene conversions
                               "syn18687959.1", # Mathys gene conversions
                               "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/67/39/67390730-a684-47a5-b9f4-89c47cd4e3fc/genesgtf.gz", # seaRef gene conversions
                               "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz",
                               "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz",
                               "https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"),
                   "ihc" = "https://github.com/ellispatrick/CortexCellDeconv.git")

github <- "https://github.com/jaclynbeck-sage/Deconvolution_RNASeq/blob/main/Step01_Preprocess_ExternalMetadata.R"

UploadFile(file_gene_list, meta_folder,
           list("used" = provenance$genes,
                "executed" = github))

UploadFile(file_rosmap_ihc_proportions, meta_folder,
           list("used" = provenance$ihc,
                "executed" = github))
