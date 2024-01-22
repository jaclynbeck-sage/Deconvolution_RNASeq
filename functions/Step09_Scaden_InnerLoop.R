# Uses the "omnideconv" library to run Scaden, which has a 'train' step and
# a 'predict' step.
#
# Arguments:
#   sce = a SingleCellExperiment object containing the reference data. Must have
#         a colData column named "celltype", denoting the cell type assignment
#         of each cell.
#   bulk_mat = a matrix (genes x samples) of bulk data. Must be a dense matrix.
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
Scaden_InnerLoop <- function(sce, bulk_mat, params, algorithm) {
  # directory names need to be passed in as full paths
  full_filepath <- file.path(getwd(), dir_scaden_models)
  folder_name <- paste(params, collapse = "_")
  mod_path <- file.path(full_filepath, folder_name)
  temp_path <- file.path(full_filepath, "tmp", folder_name)

  dir.create(mod_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_path, recursive = TRUE, showWarnings = FALSE)

  # Training step requires a dense matrix of single cell counts
  model_path <- build_model_scaden(single_cell_object = as.matrix(counts(sce)),
                                   cell_type_annotations = as.character(sce$celltype),
                                   bulk_gene_expression = bulk_mat,
                                   model_path = mod_path,
                                   temp_dir = temp_path,
                                   verbose = TRUE)

  # Clear memory
  gc()

  # Predict step
  res_pcts <- deconvolute_scaden(signature = model_path,
                                 bulk_gene_expression = bulk_mat,
                                 temp_dir = temp_path,
                                 verbose = TRUE)

  res <- list("estimates" = res_pcts,
              "params" = params,
              "markers" = intersect(rownames(sce), rownames(bulk_mat)))

  # Remove temp directory with dense matrix files
  unlink(temp_path, recursive = TRUE)

  return(res)
}
