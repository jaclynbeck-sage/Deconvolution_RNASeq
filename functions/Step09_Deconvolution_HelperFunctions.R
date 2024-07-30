Modify_Music_Input <- function(data, params) {
  data$reference <- as(data$reference, "SingleCellExperiment")

  sc_basis <- Load_MusicBasis(params$reference_data_name, params$granularity)

  if (is.null(sc_basis)) {
    # Pre-compute sc.basis to save time. Use all genes in the single cell
    # data, not just the ones that exist in both bulk and sc data sets
    sce_unfiltered <- Load_SingleCell(dataset = params$reference_data_name,
                                      granularity = params$granularity,
                                      output_type = params$normalization)

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


# If there is a pre-computed corrected signature, copy it to the directory that
# the Cibersort docker container can read. If not, create a CibersortX-formatted
# data file for single cell reference input so it can compute the adjusted
# signature, then replace the reference object with the filename. Cell type
# names also can't have any "." in them, so this does a string replace to change
# them to "_". To avoid differences in how R and C sort strings with multiple
# capital letters (i.e. "Oligodendrocyte" vs "OPC" get sorted differently
# between the two languages), we also change all cell type names to lower case.
Modify_CibersortX_Input <- function(data, params) {
  # Copy adjusted signature file to cibersort directory if the file exists
  sig_params <- params %>% select(-reference_input_type)
  sig_filename <- list.files(dir_cibersort_corrected_signatures,
                             pattern = paste(sig_params, collapse = "_"),
                             full.names = TRUE)

  if (length(sig_filename) == 1) {
    file.copy(from = sig_filename,
              to = file.path(dir_cibersort, basename(sig_filename)),
              overwrite = TRUE)
    data$singlecell_filename <- ""
  } else {
    # Otherwise CibersortX needs the original single cell data to compute an
    # adjusted signature
    message("Copying single cell data to CibersortX directory...")
    sce <- Load_SingleCell(dataset = params$reference_data_name,
                           granularity = params$granularity,
                           output_type = params$normalization)
    sce$celltype <- str_replace(as.character(sce$celltype), "\\.", "_")
    sce$celltype <- str_to_lower(sce$celltype)

    f_name <- Save_SingleCellToCibersort(sce = sce,
                                         dataset_name = params$reference_data_name,
                                         granularity = params$granularity)
    data$singlecell_filename <- f_name
    rm(sce)
    gc()
  }

  # Lower-casing the reference signature has to be handled in the inner loop
  # so the inner loop can convert back to the original names
  colnames(data$reference) <- str_replace(colnames(data$reference), "\\.", "_")

  return(data)
}
