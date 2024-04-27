library(synapser)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
meta_folder <- Folder("metadata", parent = "syn49332774")
meta_folder <- synStore(meta_folder, forceVersion = FALSE)

provenance <- list("genes" = c("syn26967452", # ROSMAP, Mayo, MSBB gene conversions
                               "syn18687959", # Mathys gene conversions
                               "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/67/39/67390730-a684-47a5-b9f4-89c47cd4e3fc/genesgtf.gz", # seaRef gene conversions
                               "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz",
                               "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz",
                               "https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"),
                   "ihc" = "https://github.com/ellispatrick/CortexCellDeconv.git")

file1 <- File(path = file_gene_list, parent = meta_folder)
file1 <- synStore(file1, used = provenance[["genes"]], forceVersion = FALSE)
print(paste0(file1$properties$id, ": ", file1$properties$name))

file2 <- File(path = file_rosmap_ihc_proportions, parent = meta_folder)
file2 <- synStore(file2, used = provenance[["ihc"]], forceVersion = FALSE)
print(paste0(file2$properties$id, ": ", file2$properties$name))
