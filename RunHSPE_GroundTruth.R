library(hspe)
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- list("mathys")#,"cain","lau","morabito","lengSFG","lengEC")

datatypes <- list("donors", "training")

###load bulk and snRNA-seq data###
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
    meta <- colData(sce)

    sce_cpm <- calculateCPM(counts(sce))

    sc_mats <- list()

    sc_mats[["logcpm"]] <- sce_cpm
    sc_mats[["logcpm"]]@x <- log2(sce_cpm@x + 1)

    sc_mats[["logcpm_scuttle"]] <- normalizeCounts(counts(sce), log = TRUE, pseudo.count = 1)

    sc_mats[["logcounts"]] <- counts(sce)
    sc_mats[["logcounts"]]@x <- log2(counts(sce)@x + 1)

    if (cellclasstype == "broad") {
      broadtypes <- levels(meta$broadcelltype)
      pure_samples <- list()
      for (ii in broadtypes) {
        pure_samples[[ii]] <- which(meta$broadcelltype == ii)
      }
    }

    if (cellclasstype == "fine") {
      finetypes <- levels(meta$subcluster)
      pure_samples <- list()
      for (ii in finetypes) {
        pure_samples[[ii]] <- which(meta$subcluster == ii)
      }
    }

    pseudobulk <- assays(pseudobulk)[["counts"]]

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(sce), rownames(pseudobulk))

    pseudobulk_cpm <- calculateCPM(pseudobulk)

    pb_mats <- list()

    pb_mats[["logcpm"]] <- log2(pseudobulk_cpm + 1)
    pb_mats[["logcpm_scuttle"]] <- normalizeCounts(pseudobulk, log = TRUE, pseudo.count = 1)
    pb_mats[["logcounts"]] <- log2(pseudobulk + 1)

    # Clear up as much memory as possible
    rm(sce, sce_cpm, pseudobulk, pseudobulk_cpm)
    gc()

    hspe_list <- list()

    ###test hspe parameters on pseudo-bulk combinations

    for (normtype in c("logcpm")) { #, "logcpm_scuttle", "logcounts")) {

      # New list for each normtype. Otherwise the list takes up a lot of memory
      #hspe_list <- list()

      # Pre-combine matrices so this isn't repeatedly done on every hspe call.
      # Single cell data must be first so indices in pure_samples are correct
      # TODO consider quantile normalization of the matrix
      Y <- as.matrix(t(cbind(sc_mats[[normtype]][keepgene,],
                             pb_mats[[normtype]][keepgene,])))
      gc()

      # We don't need to re-calculate markers every type we change n_markers
      # or loss_fn, so do it outside those loops
      # Note: "regression" is also an option for marker method but it needs more
      # than 64 GB of memory so I haven't run it.
      for (marker_meth in c("ratio", "diff", "p.value")) {
        markers <- find_markers(Y = Y, pure_samples = pure_samples,
                                marker_method = marker_meth)
        gc()

        markers_use <- list(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, # Percent of markers in each cell type
                            lengths(markers$L)) # All markers for each cell type

        for (loss_fn in c("var", "L2")) {
          for (n_markers in markers_use) {
            n_markers_name <- n_markers
            if (is(n_markers, "integer")) { # the list of lengths
              n_markers_name <- "all"
            }

            name <- paste(sndata, cellclasstype,
                          "method", marker_meth,
                          "lossfn", loss_fn,
                          "nmarkers", n_markers_name,
                          "normalization", normtype, sep = "_")

            result <- hspe(Y = Y,
                           pure_samples = pure_samples,
                           n_markers = n_markers,
                           markers = markers,
                           loss_fn = loss_fn,
                           seed = 12345)

            # Only keep results for pseudobulk samples
            result$estimates <- result$estimates[colnames(pb_mats[[normtype]]), ]

            hspe_list[[name]] <- result

            gc()
            print(name)
          } # End n_markers loop
        } # End loss_fn loop

        # Periodically save the list in case of crash
        print("Saving list checkpoint...")
        saveRDS(hspe_list, file = file.path(dir_output,
                                            paste0("hspe_list_", sndata,
                                                   "_", datatype,
                                                   "_", cellclasstype, ".rds")))
      } # End marker_meth loop

      # Next iteration will start with new data, remove the old data
      rm(Y)
      gc()
    } # End normtype loop

    print("Saving final list...")
    saveRDS(hspe_list, file = file.path(dir_output,
                                        paste0("hspe_list_", sndata,
                                               "_", datatype,
                                               "_", cellclasstype, ".rds")))
  } # End datatype loop
}
