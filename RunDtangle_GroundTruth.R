library(dtangleSparse) # dtangle with my mods to make it sparse matrix-friendly
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(stringr)

source("Filenames.R")

marker_file_format <- paste0("dtangle_markers_{sndata}_{cellclasstype}",
                             "_normalization_{normtype}_method_{marker_meth}",
                             "_gamma_{gamma_name}.rds")

list_name_format <- paste0("{sndata}_{cellclasstype}_method_{marker_meth}",
                           "_gamma_{gamma_name}_summaryfn_{sum_fn_type}",
                           "_nmarkers_{n_markers_name}_normalization_{normtype}")

list_file_format <- "dtangle_list_{sndata}_{datatype}_{cellclasstype}.rds"

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito") #,
              #"seaRef") #, "seaAD")

datatypes <- list("donors", "training")

###load bulk and snRNA-seq data###
for (sndata in datasets) {
  sce <- readRDS(file.path(dir_input, paste(sndata, "sce.rds", sep = "_")))
  meta <- colData(sce)

  sce_cpm <- calculateCPM(counts(sce))

  sc_mats <- list()

  sc_mats[["logcpm"]] <- sce_cpm
  sc_mats[["logcpm"]]@x <- log2(sce_cpm@x + 1)

  #sc_mats[["logcpm_scuttle"]] <- normalizeCounts(counts(sce), log = TRUE, pseudo.count = 1)

  #sc_mats[["logcounts"]] <- counts(sce)
  #sc_mats[["logcounts"]]@x <- log2(counts(sce)@x + 1)

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

  # Clear up as much memory as possible
  rm(meta, sce, sce_cpm)
  gc()

  params <- expand.grid(normtype = c("logcpm"), # room for other types of normalizations
                        marker_meth = c("ratio", "diff"), #, "p.value"), # p.value is slow, do later
                        gamma_name = c("auto", 1),
                        stringsAsFactors = FALSE)

  # Create marker lists and save them. We don't need to re-calculate markers
  # every time we change other parameters so this is calculated outside the
  # analysis loop.
  # Note: "regression" is also an option for marker method but it needs more
  # than 64 GB of memory so I haven't run it.
  for (R in 1:nrow(params)) {
    normtype <- params$normtype[R]
    marker_meth <- params$marker_meth[R]
    gamma_name <- params$gamma_name[R]

    if (gamma_name == "auto") {
      gamma = NULL
    }
    else {
      gamma <- as.numeric(gamma_name)
    }

    print(paste0("Finding markers for method = ", marker_meth,
                 ", gamma = ", gamma_name, " ..."))

    markers <- find_markers(Y = t(sc_mats[[normtype]]),
                            pure_samples = pure_samples,
                            marker_method = marker_meth,
                            data_type = "rna-seq",
                            gamma = gamma)

    saveRDS(markers, file = file.path(dir_markers,
                                      str_glue(marker_file_format)))
    rm(markers)
    gc()
    print("Done.")
  }

  for (datatype in datatypes) {
    if (cellclasstype == "fine") {
      load(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_finecelltypes_30percentlimit.rda")))
    }
    if (cellclasstype=="broad") {
      pseudobulk <- readRDS(file.path(dir_pseudobulk,
                                      paste0("pseudobulk_", sndata, "_",
                                             datatype, "_broadcelltypes.rds")))
    }

    pseudobulk <- assays(pseudobulk)[["counts"]]

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(sc_mats[[1]]), rownames(pseudobulk))

    pseudobulk_cpm <- calculateCPM(pseudobulk)

    pb_mats <- list()

    pb_mats[["logcpm"]] <- log2(pseudobulk_cpm + 1)
    #pb_mats[["logcpm_scuttle"]] <- normalizeCounts(pseudobulk, log = TRUE, pseudo.count = 1)
    #pb_mats[["logcounts"]] <- log2(pseudobulk + 1)

    # Clear up as much memory as possible
    rm(pseudobulk, pseudobulk_cpm)
    gc()

    ###test dtangle parameters on pseudo-bulk combinations

    params_run <- expand.grid(marker_meth = unique(params$marker_meth),
                              gamma_name = unique(params$gamma_name),
                              sum_fn_type = c("mean", "median"),
                              stringsAsFactors = FALSE)

    dtangle_list <- list()

    for (normtype in unique(params$normtype)) {
      # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
      # Single cell data must be first so indices in pure_samples are correct
      # TODO consider quantile normalization of the matrix
      Y <- t(cbind(sc_mats[[normtype]][keepgene,],
                   pb_mats[[normtype]][keepgene,]))
      gc()

      for (R in 1:nrow(params_run)) {
        marker_meth <- params_run$marker_meth[R]
        gamma_name <- params_run$gamma_name[R]
        sum_fn_type <- params_run$sum_fn_type[R]

        if (gamma_name == "auto") {
          gamma <- NULL
        }
        else {
          gamma <- as.numeric(gamma_name)
        }

        sum_fn <- mean
        if (sum_fn_type == "median") {
          sum_fn <- median
        }

        markers <- readRDS(file.path(dir_markers, str_glue(marker_file_format)))
        markers_use <- list(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, # Percent of markers in each cell type
                            lengths(markers$L)) # All markers for each cell type

        for (n_markers in markers_use) {
          n_markers_name <- n_markers
          if (is(n_markers, "integer")) { # the list of lengths
            n_markers_name <- "all"
          }

          name <- str_glue(list_name_format)

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

        # Periodically save the list in case of crash
        print("Saving list checkpoint...")
        saveRDS(dtangle_list, file = file.path(dir_output,
                                               str_glue(list_file_format)))
      } # End params_run loop

      # Next iteration will start with new data, remove the old data
      rm(Y)
      gc()
    } # End normtype loop

    # Save the completed list
    print("Saving final list...")
    saveRDS(dtangle_list, file = file.path(dir_output,
                                           str_glue(list_file_format)))
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
