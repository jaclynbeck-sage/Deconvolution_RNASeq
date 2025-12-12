# Uses the "omnideconv" library to run Scaden, which has a 'train' step and
# a 'predict' step.
#
# Arguments:
#   data = a named list that must contain the following items:
#     reference = a Summarized object containing the reference data, which is
#                 a simulated pseudobulk data set generated in Step 05.
#     test = a matrix (genes x samples) of bulk data. Must be a dense matrix.
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: algorithm, reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            and marker_subtype. This variable is unused except to get added to
#            the results object.
#
# Returns:
#   a list containing entries for the celltype percentage estimates ("estimates"),
#   "params", which is the parameter set used for this run, and "markers", which
#   is the list of genes shared by the single cell and bulk data sets (since
#   Scaden doesn't actually use markers).
Scaden_InnerLoop <- function(data, params) {
  reticulate::use_virtualenv("r-omnideconv")

  # directory names for the model folders need to be passed in as full paths
  full_filepath <- file.path(getwd(), dir_scaden_models)
  folder_name <- paste(params, collapse = "_")
  mod_path <- file.path(full_filepath, folder_name)
  temp_path <- file.path(full_filepath, "tmp", folder_name)

  dir.create(mod_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_path, recursive = TRUE, showWarnings = FALSE)

  # Put the simulated data in the AnnData format expected by Scaden
  adata <- anndata::AnnData(X = t(assay(data$reference, "counts")),
                            obs = as.data.frame(colData(data$reference)),
                            var = list(symbol = rownames(data$reference)))
  adata$uns <- list(cell_types = colnames(adata$obs))

  # Train the model with this tissue's data set
  processed <- omnideconv::scaden_process(h5ad = adata,
                                          temp_dir = temp_path,
                                          bulk_gene_expression = data$test,
                                          verbose = TRUE)

  model_path <- omnideconv::scaden_train(h5ad_processed = processed,
                                         temp_dir = temp_path,
                                         model_path = mod_path,
                                         verbose = TRUE)

  # Clear memory
  gc()

  # Predict step
  ests <- omnideconv::deconvolute_scaden(signature = model_path,
                                         bulk_gene_expression = data$test,
                                         temp_dir = temp_path,
                                         verbose = TRUE)

  # Cell type names get changed with make.names() when making the reference
  # data, so we undo that change here
  stopifnot("celltype_names" %in% names(metadata(data$reference)))

  orig_celltypes <- metadata(data$reference)$celltype_names
  mn <- make.names(orig_celltypes)

  stopifnot(all(colnames(ests) %in% mn))
  ests <- ests[, mn]
  colnames(ests) <- orig_celltypes

  # Remove temp directory with dense matrix files
  unlink(temp_path, recursive = TRUE)

  res <- list("estimates" = ests,
              "params" = params,
              "markers" = intersect(rownames(data$reference), rownames(data$test)))

  # Put the samples back in their original order
  res$estimates <- res$estimates[colnames(data$test), ]

  return(res)
}
