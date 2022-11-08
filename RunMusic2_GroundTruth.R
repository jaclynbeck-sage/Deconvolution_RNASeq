library(MuSiC)
library(SingleCellExperiment)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- list("mathys")#,"cain","lau","morabito","lengSFG","lengEC")

# Workaround for a bug in music2_prop_t_statistics
exprs <- function(X) {X}

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

  controls <- paste0("donor", unique(meta$donor[meta$diagnosis == "Control"]))
  disease <- paste0("donor", unique(meta$donor[meta$diagnosis == "AD"]))

  sce <- SingleCellExperiment(assays = list(counts = fullmat[keepgene,]),
                              colData = meta)

  # Clear up some memory
  rm(fullmat, se)
  gc()

  # Every combination of parameters being tested
  params <- expand.grid(ct.cov = c(TRUE, FALSE),
                        centered = c(TRUE, FALSE),
                        normalize = c(TRUE, FALSE))

  music_list <- list()

  # Test different combinations of parameters
  for (samples in c('donor')) { #}, 'cellid')) { # cellid is way too slow for quick testing
    for (music2type in c("default", "t_statistics")) { #, "TOAST")) { # TOAST function is broken

      # The TOAST function doesn't have ct.cov, centered, or normalize parameters
      if (music2type == "TOAST") {
        name <- paste(sndata, cellclasstype,
                      "samples", samples,
                      "music2type", music2type,
                      "counts", sep = "_")

        music_list[[name]] = music2_prop_toast_fix(bulk.control.mtx = pseudobulk[, controls],
                                               bulk.case.mtx = pseudobulk[, disease],
                                               sc.sce = sce, select.ct = levels(meta$broadcelltype),
                                               clusters = 'broadcelltype',
                                               samples = samples)

        gc()
        print(name)

        next # Skip iterating over parameters below
      }

      for (R in 1:nrow(params)) {
        ct.cov <- params$ct.cov[R]
        centered <- params$centered[R]
        normalize <- params$normalize[R]

        name <- paste(sndata, cellclasstype,
                      "samples", samples,
                      "music2type", music2type,
                      "ct.cov", ct.cov,
                      "centered", centered,
                      "normalize", normalize,
                      "counts", sep = "_")

        if (music2type == "default") {
          music_list[[name]] = music2_prop(bulk.control.mtx = pseudobulk[, controls],
                                           bulk.case.mtx = pseudobulk[, disease],
                                           sc.sce = sce,
                                           clusters = 'broadcelltype',
                                           samples = samples,
                                           select.ct = levels(meta$broadcelltype),
                                           ct.cov = ct.cov, centered = centered,
                                           normalize = normalize)
        }
        else if (music2type == "t_statistics") {
          music_list[[name]] = music2_prop_t_statistics(
                                          bulk.control.mtx = pseudobulk[, controls],
                                          bulk.case.mtx = pseudobulk[, disease],
                                          sc.sce = sce,
                                          clusters = 'broadcelltype',
                                          samples = samples, select.ct = levels(meta$broadcelltype),
                                          ct.cov = ct.cov, centered = centered,
                                          normalize = normalize)
        }

        gc()
        print(name)
      } # End params loop

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




