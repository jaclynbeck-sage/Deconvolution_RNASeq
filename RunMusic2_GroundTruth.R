library(MuSiC)
library(SingleCellExperiment)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- list("mathys")#,"cain","lau","morabito","lengSFG","lengEC")

for (sndata in datasets) {
  if (cellclasstype == "fine") {
    load(paste0("pseudobulk_", sndata, "_finecelltypes_30percentlimit.rda"))
  }
  if (cellclasstype=="broad") {
    #se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_",sndata,"_broadcelltypes.rds")))
    se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_bydonor_broadcelltypes.rds")))
  }

  load(file.path(dir_input, paste0(sndata,"_counts.rda")))

  # TODO: metadata processing should be done in a function since multiple files do this
  meta <- read.csv(file.path(dir_input, paste0(sndata,"_metadata.csv")), as.is=T)
  rownames(meta) <- meta$cellid

  meta$broadcelltype <- factor(meta$broadcelltype)
  meta$subcluster <- factor(meta$subcluster)

  keep <- intersect(meta$cellid, colnames(fullmat))
  meta <- meta[keep,]
  fullmat <- fullmat[,keep]

  pseudobulk <- assays(se)[["counts"]]

  # These SHOULD have the same rownames, but just in case.
  keepgene <- intersect(rownames(fullmat),rownames(pseudobulk))
  pseudobulk <- as.matrix(pseudobulk[keepgene,])

  controls <- unique(meta$donor[meta$diagnosis == "Control"])
  disease <- unique(metadata$donor[meta$diagnosis == "??"])

  sce <- SingleCellExperiment(assays = list(counts = fullmat[keepgene,]),
                              colData = meta)

  # Clear up some memory
  rm(fullmat, se)
  gc()

  music_list <- list()

  # Test different combinations of parameters
  for (samples in c('donor', 'cellid')) {
    for (music2type in c("default", "t_statistics", "TOAST")) {

      # The TOAST function doesn't have ct.cov, centered, or normalize parameters
      if (music2type == "TOAST") {
        name <- paste(sndata, cellclasstype,
                      "samples", samples,
                      "music2Type", music2Type,
                      "counts", sep = "_")
        music_list[[name]] = music2_toast(bulk.control.mtx = pseudobulk[, controls],
                                          bulk.case.mtx = pseudobulk[, disease],
                                          sc.sce = sce,
                                          clusters = 'broadcelltype',
                                          samples = samples)
        next # Skip iterating over parameters below
      }

      for (ct.cov in c(TRUE, FALSE)) {
        for (centered in c(TRUE, FALSE)) {
          for (normalize in c(TRUE, FALSE)) {
            name <- paste(sndata, cellclasstype,
                          "samples", samples,
                          "music2Type", music2Type,
                          "ct.cov", ct.cov,
                          "centered", centered,
                          "normalize", normalize,
                          "counts", sep = "_")

            if (music2type == "default") {
              music_list[[name]] = music2_prop(bulk.control.mtx = pseudobulk[, controls],
                                               bulk.case.mtx = pseudobulk[, disease],
                                               sc.sce = sce,
                                               clusters = 'broadcelltype',
                                               samples = samples, verbose = TRUE,
                                               ct.cov = ct.cov, centered = centered,
                                               normalize = normalize)
            }
            else if (music2type == "t_statistics") {
              music_list[[name]] = music2_prop_t_statistics(
                                               bulk.control.mtx = pseudobulk[, controls],
                                               bulk.case.mtx = pseudobulk[, disease],
                                               sc.sce = sce,
                                               clusters = 'broadcelltype',
                                               samples = samples, verbose = TRUE,
                                               ct.cov = ct.cov, centered = centered,
                                               normalize = normalize)
            }

            gc()
            print(name)
          } # End normalize loop
        } # End centered loop
      } # End ct.cov loop

      # Periodically save the list, in case of crashes
      print("Saving list checkpoint...")
      saveRDS(music_list, file = file.path(dir_output,
                                           paste0("music2_list_", sndata, "_donors_",
                                                  cellclasstype, ".rds")))
    } # End music2type loop
  } # End samples loop

  print("Saving final list...")
  saveRDS(music_list, file = file.path(dir_output,
                                       paste0("music2_list_", sndata, "_donors_",
                                              cellclasstype, ".rds")))
}




