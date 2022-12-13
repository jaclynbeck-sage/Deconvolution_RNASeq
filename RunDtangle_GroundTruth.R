library(dtangle)
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
    rm(sce_cpm, pseudobulk, pseudobulk_cpm)
    gc()

    dtangle_list <- list()

    ###test dtangle parameters on pseudo-bulk combinations

    for (normtype in c("logcpm")) { #, "logcpm_scuttle", "logcounts")) {

      # New list for each normtype. Otherwise the list takes up a lot of memory
      #dtangle_list <- list()

      # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
      # Single cell data must be first so indices in pure_samples are correct
      # TODO consider quantile normalization of the matrix
      Y <- as.matrix(t(cbind(sc_mats[[normtype]][keepgene,],
                             pb_mats[[normtype]][keepgene,])))
      gc()

      # We don't need to re-calculate markers every type we change n_markers
      # or summary_fn, so do it outside those loops
      # Note: "regression" is also an option for marker method but it needs more
      # than 64 GB of memory so I haven't run it.
      for (marker_meth in c("ratio", "diff", "p.value")) {

        for (gamma_name in list("auto", 1)) {
          gamma <- gamma_name
          if (gamma_name == "auto") {
            gamma = NULL
          }

          markers <- find_markers(Y = Y, pure_samples = pure_samples,
                                  data_type = "rna-seq",
                                  gamma = gamma,
                                  marker_method = marker_meth)
          gc()

          markers_use <- list(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, # Percent of markers in each cell type
                              lengths(markers$L)) # All markers for each cell type

          for (sum_fn_type in c("mean", "median")) {
            sum_fn <- mean
            if (sum_fn_type == "median") {
              sum_fn <- median
            }

            for (n_markers in markers_use) {
              n_markers_name <- n_markers
              if (is(n_markers, "integer")) { # the list of lengths
                n_markers_name <- "all"
              }

              name <- paste(sndata, cellclasstype,
                            "method", marker_meth,
                            "gamma", gamma_name,
                            "summaryfn", sum_fn_type,
                            "nmarkers", n_markers_name,
                            "normalization", normtype, sep = "_")

              result <- dtangle(Y = Y,
                                pure_samples = pure_samples,
                                data_type = "rna-seq",
                                gamma = gamma, # If gamma is not NULL, it will override data_type argument
                                n_markers = n_markers,
                                markers = markers,
                                summary_fn = sum_fn)

              # Only keep results for pseudobulk samples
              result$estimates <- result$estimates[colnames(pb_mats[[normtype]]), ]

              dtangle_list[[name]] <- result

              gc()
              print(name)
            } # End n_markers loop
          } # End sum_fn_type loop
        } # End gamma_name loop

        # Periodically save the list in case of crash
        print("Saving list checkpoint...")
        saveRDS(dtangle_list, file = file.path(dir_output,
                                               paste0("dtangle_list_", sndata,
                                                     "_", datatype,
                                                     "_", cellclasstype, ".rds")))
      } # End marker_meth loop

      # Next iteration will start with new data, remove the old data
      rm(Y)
      gc()
    } # End normtype loop

    # Save the completed list
    print("Saving final list...")
    saveRDS(dtangle_list, file = file.path(dir_output,
                                           paste0("dtangle_list_", sndata,
                                                  "_", datatype,
                                                  "_", cellclasstype, ".rds")))
  } # end datatypes loop
}


# Test code -- out of date
se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_bydonor_broadcelltypes.rds")))
propval <- as.matrix(colData(se))

ses = list()
for (N in names(dtangle_list)) {
  tmp = dtangle_list[[N]]$estimates
  ses[[N]] = ((tmp[rownames(propval), colnames(propval)] - propval)/(propval + 1e-10))^2 # Proportional to cell proportion -- TODO this biases against ones where propval is 0
  mse = mean(ses[[N]])
  rmse = sqrt(mse)
  print(paste(N, ": ", rmse / 1000000))
}

tmp <- sapply(ses, mean)
best <- min(tmp[!is.na(tmp)])
best1 <- names(tmp)[which(tmp == best)]


ses2 = list()
for (N in names(dtangle_list)) {
  tmp = dtangle_list[[N]]$estimates
  ses2[[N]] = ((tmp[rownames(propval), colnames(propval)] - propval))^2
  mse = mean(ses2[[N]])
  rmse = sqrt(mse)
  print(paste(N, ": ", rmse))
}

tmp <- sapply(ses2, mean)
best <- min(tmp[!is.na(tmp)])
best2 <- names(tmp)[which(tmp == best)]


ses3 = list()
for (N in names(dtangle_list)) {
  tmp = dtangle_list[[N]]$estimates
  ses3[[N]] = abs(tmp[rownames(propval), colnames(propval)] - propval)
  me = mean(ses3[[N]])
  print(paste(N, ": ", me))
}

tmp <- sapply(ses3, mean)
best <- min(tmp[!is.na(tmp)])
best3 <- names(tmp)[which(tmp == best)]


matplot(propval, dtangle_list[[best1]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
        xlab="Truth", ylab="Estimates", main = best1)
abline(a = 0, b = 1)

matplot(propval, dtangle_list[[best2]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
        xlab="Truth", ylab="Estimates", main = best2)
abline(a = 0, b = 1)

matplot(propval, dtangle_list[[best3]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
        xlab="Truth", ylab="Estimates", main = best3)
abline(a = 0, b = 1)
