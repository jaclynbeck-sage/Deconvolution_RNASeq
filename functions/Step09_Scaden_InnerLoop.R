# Uses the "omnideconv" library to run Scaden, which has a 'train' step and
# a 'predict' step. Scaden trains the bulk data against the single cell data,
# so to improve results we input one tissue at a time instead of the whole
# bulk data set at once.
#
# Arguments:
#   data = a named list that must contain the following items:
#     reference = a SingleCellExperiment object containing the reference data.
#                 Must have a colData column named "celltype", denoting the cell
#                 type assignment of each cell.
#     test = a matrix (genes x samples) of bulk data. Must be a dense matrix.
#     metadata = a data frame with columns for "sample" and "tissue"
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            marker_subtype, marker_input_type. This variable is unused except
#            to get added to the results object.
#   algorithm = the name of the algorithm being run. This value is unused.
#
# Returns:
#   a list containing entries for the celltype percentage estimates ("estimates"),
#   "params", which is the parameter set used for this run, and "markers", which
#   is the list of genes shared by the single cell and bulk data sets (since
#   Scaden doesn't actually use markers).
Scaden_InnerLoop <- function(data, params, algorithm = NULL) {
  # Run each tissue separately
  res_pcts_tissue <- lapply(levels(data$bulk_metadata$tissue), function(tiss) {
    meta_sub <- subset(data$bulk_metadata, tissue == tiss)
    data_sub <- data$test[, as.character(meta_sub$sample)]

    # Check that there are no NA samples
    stopifnot(all(colnames(data_sub) == meta_sub$sample))

    # directory names for the model folders need to be passed in as full paths
    full_filepath <- file.path(getwd(), dir_scaden_models)
    folder_name <- paste0(paste(params, collapse = "_"),
                          "_", tiss)
    mod_path <- file.path(full_filepath, folder_name)
    temp_path <- file.path(full_filepath, "tmp", folder_name)

    dir.create(mod_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(temp_path, recursive = TRUE, showWarnings = FALSE)

    # Check if simulated data already exists -- if so, it can be re-used.
    dir_simulated <- file.path(full_filepath, "simulated_data")
    file_simulated <- file.path(dir_simulated,
                                paste0(paste(params$reference_data_name,
                                             params$granularity,
                                             params$normalization,
                                             params$n_cells,
                                             sep = "_"),
                                       ".h5ad"))

    # TPM and CPM normalization are the same for single cell
    file_simulated <- str_replace(file_simulated, "_tpm_", "_cpm_")

    # If it doesn't exist, make simulated data
    if (!file.exists(file_simulated)) {
      # Needs to be a dense matrix
      dense_sce <- t(as.matrix(counts(data$reference)))
      gc()

      simulated_h5 <- omnideconv::scaden_simulate(
        cell_type_annotations = as.character(data$reference$celltype),
        gene_labels = rownames(data$reference),
        single_cell_object = dense_sce,
        temp_dir = temp_path,
        cells = params$n_cells,
        verbose = TRUE
      )

      rm(dense_sce)
      gc()

      # Move the simulated data file out of the temp directory to its expected
      # location for use next time
      dir.create(dir_simulated, recursive = TRUE, showWarnings = FALSE)
      file.copy(from = file.path(temp_path, "data.h5ad"),
                to = file_simulated)
    } else {
      message(str_glue("Using existing simulated data to train model: {basename(file_simulated)}"))
      simulated_h5 <- anndata::read_h5ad(file_simulated)
    }

    # Train the model with this tissue's data set
    processed <- omnideconv::scaden_process(h5ad = simulated_h5,
                                            temp_dir = temp_path,
                                            bulk_gene_expression = data_sub,
                                            verbose = TRUE)

    model_path <- omnideconv::scaden_train(h5ad_processed = processed,
                                           temp_dir = temp_path,
                                           model_path = mod_path,
                                           verbose = TRUE)

    # Clear memory
    gc()

    # Predict step
    res_pcts <- omnideconv::deconvolute_scaden(signature = model_path,
                                               bulk_gene_expression = data_sub,
                                               temp_dir = temp_path,
                                               verbose = TRUE)

    # Remove temp directory with dense matrix files
    unlink(temp_path, recursive = TRUE)

    return(res_pcts)
  })

  # Combine results from each tissue. Params and markers are not tissue-specific.
  res <- list("estimates" = do.call(rbind, res_pcts_tissue),
              "params" = params,
              "markers" = intersect(rownames(data$reference), rownames(data$test)))

  # Put the samples back in their original order
  res$estimates <- res$estimates[colnames(data$test), levels(data$reference$celltype)]

  return(res)
}
