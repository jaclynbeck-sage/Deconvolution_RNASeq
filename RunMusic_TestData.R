library(MuSiC)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

# HACK to MuSiC to account for the case where some weights are < 0, which leads
# to error-causing NaNs later. seaRef is the only data set this has happened on,
# so I need to figure out what's actually going on here.
weight.cal.ct.orig <- MuSiC::weight.cal.ct
weight.cal.ct <- function (...) {
  weight = weight.cal.ct.orig(...)
  weight[weight < 0] = 0
  return(weight)
}

R.utils::reassignInPackage("weight.cal.ct", pkgName="MuSiC", weight.cal.ct)

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

  # The seaRef dataset will fit in memory all at once, so this converts it
  # to a sparse matrix. The seaAD data set will NOT fit so this won't work on it.
  if (is(counts(sce), "DelayedArray")) {
    counts(sce) <- as(counts(sce), "dgCMatrix")
  }

  best_params <- readRDS(file.path(dir_output, paste0("best_params_", sndata,
                                                      "_", cellclasstype, ".rds")))
  best_params <- subset(best_params, algorithm == "music_nnls" | algorithm == "music_wt")
  params <- best_params %>%
    extract(name, c("ct.cov", "centered", "normalize"),
            "ct.cov_([A-Z]+)_centered_([A-Z]+)_normalize_([A-Z]+)",
            convert = TRUE) %>%
    select(ct.cov, centered, normalize) %>% distinct()


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

    result <- music_prop(bulk.mtx = bulk, sc.sce = sce,
                         clusters = 'broadcelltype',
                         samples = samples, verbose = TRUE,
                         ct.cov = ct.cov, centered = centered,
                         normalize = normalize)

    # Remove "Weight.gene", "r.squared.full", and "Var.prop". "Weight.gene"
    # especially is a very large array and is unneeded, so this reduces
    # output size.
    music_list[[name]] <- result[1:2]

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




