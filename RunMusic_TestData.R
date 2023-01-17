library(MuSiC)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito") #,
              #"seaRef") #, "seaAD")

for (sndata in datasets) {
  sce <- readRDS(file.path(dir_input, paste(sndata, "sce.rds", sep = "_")))

  # Convert gene names to Ensembl IDs
  genes <- rowData(sce)
  rownames(sce) <- genes[rownames(sce), "Ensembl.ID"]

  # Specific to ROSMAP for now
  bulk <- read.table(file_rosmap, header = TRUE, row.names = 1)

  keepgene <- intersect(rownames(sce), rownames(bulk))

  bulk <- as.matrix(bulk[keepgene, ]) # Needs to be a matrix, not data.frame
  sce <- sce[keepgene, ]

  best_params <- readRDS(file.path(dir_output, paste0("best_params_", sndata,
                                                      "_", cellclasstype, ".rds")))
  best_params <- subset(best_params, algorithm == "music")
  params <- best_params %>%
    extract(name, c("ct.cov", "centered", "normalize"),
            "ct.cov_([A-Z]+)_centered_([A-Z]+)_normalize_([A-Z]+)",
            convert = TRUE)

  music_list <- list()

  for (R in 1:nrow(params)) {
    samples = "donor"
    ct.cov <- params$ct.cov[R]
    centered <- params$centered[R]
    normalize <- params$normalize[R]

    name <- paste(sndata, cellclasstype,
                  "samples", samples,
                  "ct.cov", ct.cov,
                  "centered", centered,
                  "normalize", normalize,
                  "counts", sep = "_")

    music_list[[name]] = music_prop(bulk.mtx = bulk, sc.sce = sce,
                                    clusters = 'broadcelltype',
                                    samples = samples, verbose = TRUE,
                                    ct.cov = ct.cov, centered = centered,
                                    normalize = normalize)

    gc()
    print(name)
  } # End params loop

  print("Saving final list...")
  saveRDS(music_list, file = file.path(dir_output,
                                       paste0("music_list_", sndata, "_",
                                              cellclasstype, "_ROSMAP.rds")))

  rm(music_list, bulk, sce)
  gc()
}




