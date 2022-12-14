library(MuSiC)
library(SingleCellExperiment)
library(SummarizedExperiment)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatypes <- list("donors", "training")

for (sndata in datasets) {
  for (datatype in datatypes) {
    if (cellclasstype == "fine") {
      load(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_finecelltypes_30percentlimit.rda")))
    }
    if (cellclasstype=="broad") {
      pseudobulk <- readRDS(file.path(dir_pseudobulk,
                                      paste0("pseudobulk_", sndata, "_",
                                             datatype, "_broadcelltypes.rds")))
    }

    sce <- readRDS(file.path(dir_input, paste(sndata, "sce.rds", sep = "_")))

    pseudobulk <- assays(pseudobulk)[["counts"]]

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(sce), rownames(pseudobulk))
    pseudobulk <- as.matrix(pseudobulk[keepgene, ])
    sce <- sce[keepgene, ]

    # Every combination of parameters being tested
    params <- expand.grid(ct.cov = c(TRUE, FALSE),
                          centered = c(TRUE, FALSE),
                          normalize = c(TRUE, FALSE))

    music_list <- list()

    # Test different combinations of parameters
    for (samples in c('donor')) { #, 'cellid')) { # cellid is way too slow for quick testing
      for (R in 1:nrow(params)) {
        ct.cov <- params$ct.cov[R]
        centered <- params$centered[R]
        normalize <- params$normalize[R]

        name <- paste(sndata, cellclasstype,
                      "samples", samples,
                      "ct.cov", ct.cov,
                      "centered", centered,
                      "normalize", normalize,
                      "counts", sep = "_")

        music_list[[name]] = music_prop(bulk.mtx = pseudobulk, sc.sce = sce,
                                        clusters = 'broadcelltype',
                                        samples = samples, verbose = TRUE,
                                        ct.cov = ct.cov, centered = centered,
                                        normalize = normalize)

        gc()
        print(name)
      } # End params loop

      # Periodically save the list, in case of crashes
      print("Saving list checkpoint...")
      saveRDS(music_list, file = file.path(dir_output,
                                           paste0("music_list_", sndata, "_",
                                                  datatype, "_",
                                                  cellclasstype, ".rds")))
    } # End samples loop

    print("Saving final list...")
    saveRDS(music_list, file = file.path(dir_output,
                                         paste0("music_list_", sndata, "_",
                                                datatype, "_",
                                                cellclasstype, ".rds")))

    rm(music_list, pseudobulk, sce)
    gc()
  } # End datatypes loop
}




