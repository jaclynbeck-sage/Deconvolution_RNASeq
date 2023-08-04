library(synapser)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
pseudo.folder <- Folder("pseudobulk", parent = "syn49332774")
pseudo.folder <- synStore(pseudo.folder, forceVersion = FALSE)

provenance <- list("cain" = c("syn23554294", "syn23554292", "syn23554293"),
                   "lau" = c("syn38609700", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157827&format=file"),
                   "leng" = c("syn22722817", "syn22722860"),
                   "mathys" = c("syn3191087", "syn18686372",
                                "syn18686381", "syn18686382"),
                   "morabito" = c("syn22130806", "syn24978676", "syn24978737",
                                  "syn22130805"),
                   "seaRef" = c("syn31149116", url_searef_h5),
                   "seaAD" = c("syn31149116", url_seaad_h5))

datasets <- c("cain", "lau", "leng", "mathys", "morabito", "seaRef", "seaAD")

for (dataset in datasets) {
  files <- list.files(dir_pseudobulk, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    file <- File(path = filename, parent = pseudo.folder)
    file <- synStore(file, used = provenance[[dataset]], forceVersion = FALSE)
    print(paste0(file$properties$id, ": ", file$properties$name))
  }
}
