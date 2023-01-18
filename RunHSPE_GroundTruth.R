library(hspeSparse) # HSPE with my mods to make it sparse matrix-friendly
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(stringr)

source("Filenames.R")

marker_file_format <- paste0("hspe_markers_{sndata}_{cellclasstype}",
                             "_normalization_{normtype}_method_{marker_meth}.rds")

list_name_format <- paste0("{sndata}_{cellclasstype}_method_{marker_meth}",
                           "_lossfn_{loss_fn}_nmarkers_{n_markers_name}",
                           "_normalization_{normtype}")

list_file_format <- "hspe_list_{sndata}_{datatype}_{cellclasstype}.rds"

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito") #,
              #"seaRef", "seaAD")

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
                        stringsAsFactors = FALSE)

  # Create marker lists and save them. We don't need to re-calculate markers
  # every time we change other parameters so this is calculated outside the
  # analysis loop.
  # Note: "regression" is also an option for marker method but it needs more
  # than 64 GB of memory so I haven't run it.
  for (R in 1:nrow(params)) {
    normtype <- params$normtype[R]
    marker_meth <- params$marker_meth[R]

    print(paste0("Finding markers for method = ", marker_meth, " ..."))
    markers <- find_markers(Y = t(sc_mats[[normtype]]),
                            pure_samples = pure_samples,
                            marker_method = marker_meth)

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

    ###test hspe parameters on pseudo-bulk combinations
    params_run <- expand.grid(marker_meth = unique(params$marker_meth),
                              loss_fn = c("var", "L2"),
                              stringsAsFactors = FALSE)

    hspe_list <- list()

    for (normtype in unique(params$normtype)) {
      # Pre-combine matrices so this isn't repeatedly done on every hspe call.
      # Single cell data must be first so indices in pure_samples are correct
      # TODO consider quantile normalization of the matrix
      Y <- t(cbind(sc_mats[[normtype]][keepgene,],
                   pb_mats[[normtype]][keepgene,]))
      gc()

      for (R in 1:nrow(params_run)) {
        marker_meth <- params_run$marker_meth[R]
        loss_fn <- params_run$loss_fn[R]

        markers <- readRDS(file.path(dir_markers, str_glue(marker_file_format)))
        markers_use <- list(0.01, 0.02, 0.05, 0.1, 0.2, 0.5) #, 0.75, # Percent of markers in each cell type
                            #lengths(markers$L)) # All markers for each cell type

        for (n_markers in markers_use) {
          n_markers_name <- n_markers
          if (is(n_markers, "integer")) { # the list of lengths
            n_markers_name <- "all"
          }

          name <- str_glue(list_name_format)

          result <- hspe(Y = Y,
                         pure_samples = pure_samples,
                         n_markers = n_markers,
                         markers = markers,
                         loss_fn = loss_fn,
                         seed = 12345)

          # Only keep results for pseudobulk samples, and get rid of "diag",
          # which is huge and unneeded
          result$estimates <- result$estimates[colnames(pb_mats[[normtype]]), ]
          result <- result[1:4] # "diag" is #5

          hspe_list[[name]] <- result

          gc()
          print(name)
        } # End n_markers loop

        # Periodically save the list in case of crash
        print("Saving list checkpoint...")
        saveRDS(hspe_list, file = file.path(dir_output,
                                            str_glue(list_file_format)))
      } # End marker_meth loop

      # Next iteration will start with new data, remove the old data
      rm(Y)
      gc()
    } # End normtype loop

    print("Saving final list...")
    saveRDS(hspe_list, file = file.path(dir_output,
                                        str_glue(list_file_format)))
  } # End datatype loop
}
