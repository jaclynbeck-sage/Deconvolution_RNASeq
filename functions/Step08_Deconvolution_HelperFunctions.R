Modify_Music_Input <- function(data, params) {
  data$reference <- as(data$reference, "SingleCellExperiment")

  sc_basis <- Load_MusicBasis(params$reference_data_name, params$granularity)

  if (is.null(sc_basis)) {
    # Pre-compute sc.basis to save time. Use all genes in the single cell
    # data, not just the ones that exist in both bulk and sc data sets
    sce_unfiltered <- Load_SingleCell(dataset = params$reference_data_name,
                                      granularity = params$granularity,
                                      normalization = params$normalization)

    sc_basis <- MuSiC::music_basis(sce_unfiltered,
                                   non.zero = TRUE,
                                   markers = rownames(sce_unfiltered),
                                   clusters = "celltype",
                                   samples = "sample",
                                   select.ct = NULL,
                                   ct.cov = FALSE,
                                   verbose = TRUE)

    # Save for later use
    Save_MusicBasis(sc_basis,
                    params$reference_data_name,
                    params$granularity)

    rm(sce_unfiltered)
    gc()
  } else {
    message("Using pre-computed sc_basis")
  }

  data$sc_basis <- sc_basis
  return(data)
}


# Params is unused here but passed in just in case we need to change that
Modify_DtangleHSPE_Input <- function(data, params) {
  metadata <- colData(data$reference)

  celltypes <- levels(metadata$celltype)
  pure_samples <- lapply(celltypes, function(ct) {
    which(metadata$celltype == ct)
  })
  names(pure_samples) <- celltypes

  data$pure_samples <- pure_samples

  data$reference <- assay(data$reference, "counts")

  # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
  # Input data must be first so indices in pure_samples are correct.
  data$Y <- t(cbind(data$reference, data$test))

  data$reference <- NULL
  data$test <- NULL
  gc()

  return(data)
}
