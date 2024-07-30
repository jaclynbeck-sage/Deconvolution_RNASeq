# Runs Dtangle or HSPE on the given bulk data, using either single cell or
# pseudobulk data as input and a set of parameters that defines what markers
# to use.
#
# Arguments:
#   data = a named list that must contain the following items:
#     Y = a matrix of combined singlecell (or pseudobulk) data and bulk data.
#         rows are samples, columns are genes. Single cell data should be first.
#     pure_samples = a named list where each entry is a vector of indices into
#                    "Y" that correspond to that cell type, from the single cell
#                    data.
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            marker_subtype, marker_input_type.
#   algorithm = the name of the algorithm being run (options are "Dtangle"
#               and "HSPE")
#
# Returns:
#   a list containing entries for the celltype percentage estimates (name varies
#   between algorithms), "params", which is the parameter set used for this run,
#   and "markers", which is the list of genes used as markers for this run
DtangleHSPE_InnerLoop <- function(data, params, algorithm) {
  markers <- FilterMarkers_FromParams(available_genes = colnames(data$Y),
                                      params)

  if (Check_MissingMarkers(markers, params) ||
      Check_TooFewMarkers(markers, params, 3) ||
      # Check_TooManyMarkers(markers, params, 5000) || # No longer checking this
      Check_NotEnoughNewMarkers(markers, params)) {
    return(NULL)
  }

  # Dtangle-specific function call ---------------------------------------------
  if (algorithm == "Dtangle") {
    result <- dtangle(Y = data$Y[, unlist(markers)],
                      pure_samples = data$pure_samples,
                      data_type = "rna-seq",
                      n_markers = lengths(markers), # pass the actual number of markers we have after filtering
                      markers = markers)

    # Only keep results for bulk test samples
    test_samples <- setdiff(1:nrow(data$Y), unlist(data$pure_samples))
    result$estimates <- result$estimates[test_samples, ]
  }
  # HSPE-specific function call ------------------------------------------------
  else if (algorithm == "HSPE") {
    result <- hspe(Y = data$Y[, unlist(markers)],
                   pure_samples = data$pure_samples,
                   n_markers = lengths(markers), # pass the actual number of markers we have after filtering
                   markers = markers,
                   loss_fn = "L2", # 'L2' usually converges faster than 'var' with nearly identical results
                   seed = 12345,
                   verbose = TRUE)

    # Get rid of "diag" (index 5), which is huge and unneeded
    result <- result[1:4]
  }

  # Add the params we used to generate this run
  result$params <- params
  result$markers <- unlist(result$markers)
  print(paste(result$params, collapse = "  "))

  return(result)
}
