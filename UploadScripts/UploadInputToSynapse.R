library(synapser)
source("Filenames.R")

synLogin()

# Deconvolution WG Synapse space
input_folder <- Folder("input", parent = "syn49332774")
input_folder <- synStore(input_folder, forceVersion = FALSE)

provenance <- list("cain" = c("syn23554294", "syn23554292", "syn23554293"),
                   "lau" = c("syn38609700", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157827&format=file"),
                   "leng" = c("syn22722817", "syn22722860"),
                   "mathys" = c("syn3191087", "syn18686372",
                                "syn18686381", "syn18686382"),
                   "morabito" = c("syn22130806", "syn24978676", "syn24978737",
                                  "syn22130805"),
                   "seaRef" = c("syn31149116", url_searef_h5),
                   "seaAD" = c("syn31149116", url_seaad_h5),
                   "Mayo" = c("syn29855549", "syn27024951", "syn27024965",
                              "syn27024967", "syn20827192"),
                   "MSBB" = c("syn29855570", "syn27068754", "syn27068756",
                              "syn27068760", "syn21893059"),
                   "ROSMAP" = c("syn29855598", "syn26967451", "syn26967453",
                                "syn26967455", "syn21323366", "syn3191087"))

datasets <- c("cain", "lau", "leng", "mathys", "morabito",
              "seaRef", "seaAD", "Mayo", "MSBB", "ROSMAP")

for (dataset in datasets) {
  files <- list.files(dir_input, pattern = dataset, full.names = TRUE)
  for (filename in files) {
    if (!(file.info(filename)$isdir)) { # Exclude directories
      file <- File(path = filename, parent = input_folder)
      file <- synStore(file, used = provenance[[dataset]], forceVersion = FALSE)
      print(paste0(file$properties$id, ": ", file$properties$name))
    }
  }
}
