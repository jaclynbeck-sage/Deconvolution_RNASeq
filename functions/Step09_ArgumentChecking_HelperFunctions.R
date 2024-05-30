# Helper functions for checking arguments during step 9. These exist to avoid
# copy/pasting the same checks and error messages across all the algorithm files.

# For readability
Is_SignatureMatrix <- function(markers_obj) {
  return(is(markers_obj, "matrix") || is(markers_obj, "data.frame"))
}


# For readability
Is_MarkerList <- function(markers_obj) {
  return(is(markers_obj, "list"))
}


# If any cell type has 0 markers, the FilterMarkers and FilterSignature functions
# return NULL, so this function checks for a NULL markers_obj and prints out an
# message.
#
# Arguments:
#   markers_obj = either a signature matrix (or dataframe) or a list of cell type markers
#   params = a named vector or one-row data frame with the parameters being used
#            for this call of the algorithm
#
# Returns:
#   TRUE if a cell type has no markers with this parameter set, FALSE otherwise
Check_MissingMarkers <- function(markers_obj, params) {
  if (is.null(markers_obj)) {
    param_set <- paste(params, collapse = "  ")
    msg <- paste("*** Missing markers for at least one cell type for param set:",
                 param_set, "/  *** skipping ***")
    message(msg)
    return(TRUE)
  }
  return(FALSE)
}


# If there's less than ~3 markers per cell type, this isn't useful
# information. Checking total markers (signature) or mean markers (marker list)
# instead of markers per cell type allows some leeway for rarer cell types to
# have less than 3 markers as long as the other cell types have enough.
#
# Arguments:
#   same as Check_MissingMarkers
#   low_threshold = the lowest acceptable number of markers per cell type
#
# Returns:
#   TRUE if there are fewer than low_threshold markers per cell type, FALSE otherwise
Check_TooFewMarkers <- function(markers_obj, params, low_threshold = 3) {
  if (Is_SignatureMatrix(markers_obj)) {
    if (params$marker_type != "None") {
      markers <- Load_Markers(params$reference_data_name, params$granularity,
                              params$marker_type, params$marker_subtype,
                              input_type = params$marker_input_type)
      totals <- sapply(markers, function(X) {
        sum(X %in% rownames(markers_obj))
      })
      not_enough <- (mean(totals) < low_threshold)
    } else {
      not_enough <- (nrow(markers_obj) < (low_threshold * ncol(markers_obj)))
    }
  }
  else if (Is_MarkerList(markers_obj)) {
    not_enough <- (mean(lengths(markers_obj)) < low_threshold)
  }
  else {
    stop("Invalid input for marker_obj.")
  }

  if (not_enough) {
    param_set <- paste(params, collapse = "  ")
    msg <- paste("*** Too few markers for param set:", param_set,
                 "/  *** skipping ***")
    message(msg)
    return(TRUE)
  }
  return(FALSE)
}


# Early testing showed that for Dtangle/HSPE, it's not useful to use a large
# number of markers, so if the total markers used is > 5000, don't test. 5000 is
# pretty generous, and allows for 500 markers for 10 cell types. Algorithms
# that are designed to use all/most genes (DeconRNASeq, MuSiC) do not call this
# function.
# Arguments:
#   same as Check_MissingMarkers
#   high_threshold = the highest acceptable number of total markers
#
# Returns:
#   TRUE if there are more than high_threshold total markers, FALSE otherwise
Check_TooManyMarkers <- function(markers_obj, params, high_threshold = 5000) {
  if (Is_SignatureMatrix(markers_obj)) {
    too_many <- (nrow(markers_obj) > high_threshold)
  }
  else if (Is_MarkerList(markers_obj)) {
    too_many <- (sum(lengths(markers_obj)) > high_threshold)
  }
  else {
    stop("Invalid input for marker_obj.")
  }

  if (too_many) {
    param_set <- paste(params, collapse = "  ")
    msg <- paste("*** Marker set too large for param set:",
                 param_set, "/  *** skipping ***")
    message(msg)
    return(TRUE)
  }
  return(FALSE)
}


# For whole-number n_markers parameters, the n_markers argument (mostly)
# doubles each time. If there aren't enough markers in each cell type to do
# anything new with this n_markers value, skip testing.
# Arguments:
#   same as Check_MissingMarkers
#
# Returns:
#   TRUE if there are not enough new markers to be useful, FALSE otherwise
Check_NotEnoughNewMarkers <- function(markers_obj, params) {
  # Not applicable for fractional n_marker values
  if (params$n_markers <= 1) {
    return(FALSE)
  }

  low_threshold <- params$n_markers / 2

  if (Is_SignatureMatrix(markers_obj)) {
    not_enough <- nrow(markers_obj) <= (low_threshold * ncol(markers_obj))
  }
  else if (Is_MarkerList(markers_obj)) {
    not_enough <- all(lengths(markers_obj) <= low_threshold)
  }
  else {
    stop("Invalid input for marker_obj.")
  }

  if (not_enough) {
    param_set <- paste(params, collapse = "  ")
    msg <- paste("*** Less than half of the requested", params$n_markers,
                 "markers left after filtering for param set:",
                 param_set, "/  *** skipping ***")
    message(msg)
    return(TRUE)
  }
  return(FALSE)
}
