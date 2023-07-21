# Helper functions for running Dtangle and HSPE. The functions in this file
# help clear up clutter in the main algorithm script and make it clearer
# what is happening there.

source(file.path("functions", "General_HelperFunctions.R"))

# Get_DtangleHSPEInput: loads in single cell or pseudobulk pure sample data as
# the reference set, loads in the test pseudobulk data, and combines the
# reference and test data in the format dtangle/hspe expects.
#
# Arguments:
#   reference_data_name = the name of the data set to load
#   test_data_name = either "donors" or "training", for testing on pseudobulk
#                    data, or one of "Mayo", "MSBB", or "ROSMAP" to test on
#                    bulk data
#   granularity = either "broad" or "fine", for which level of cell types to
#                 use as the reference
#   reference_input_type = either "singlecell" or "pseudobulk", for whether the
#                          reference data is single cell data or pseudobulked
#                          pure samples
#
# Returns:
#   a list with entries for "Y", which is the combined reference + test data,
#   and "pure_samples", which is a list of indices into Y that correspond to
#   samples from each cell type.
Get_DtangleHSPEInput <- function(reference_data_name, test_data_name,
                                 granularity, reference_input_type, normalization) {

  data <- Load_AlgorithmInputData(reference_data_name, test_data_name,
                                  granularity, reference_input_type,
                                  output_type = normalization)

  ##### Get indices of pure samples #####
  metadata <- colData(data$reference)

  celltypes <- levels(metadata$celltype)
  pure_samples <- lapply(celltypes, function(ct) {
    which(metadata$celltype == ct)
  })
  names(pure_samples) <- celltypes

  data$reference <- assay(data$reference, "counts")
  data$test <- as(assay(data$test, "counts"), "matrix")

  # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
  # Input data must be first so indices in pure_samples are correct.
  Y <- t(cbind(data$reference, data$test))
  return(list("Y" = Y, "pure_samples" = pure_samples))
}


# Get_Gamma: helper function to turn gamma_name into a value that dtangle/hspe
# can use.
#
# Arguments:
#   gamma_name = either "auto" or a numerical value between 0 and 1
#
# Returns:
#   NULL if gamma_name is "auto", otherwise the numerical value
Get_Gamma <- function(gamma_name) {
  gamma <- NULL
  if (gamma_name != "auto") {
    gamma <- as.numeric(gamma_name)
  }

  return(gamma)
}


# Get_SumFn: helper function to turn sum_fun_type into the actual function name
#
# Arguments:
#   sum_fn_type = either "mean" or "median", as a character string
#
# Returns:
#   reference to the function mean or the function median
Get_SumFn <- function(sum_fn_type) {
  sum_fn <- mean
  if (sum_fn_type == "median") {
    sum_fn <- median
  }

  return(sum_fn)
}


# Get_NMarkers: helper function to get the actual value for n_markers to pass
# to dtangle/HSPE. If n_markers is < 1, it's a fraction and can be passed in
# as-is. If n_markers = 1, that signifies to use all markers, but we need to
# pass in a vector corresponding to the full number of genes in each cell type.
# If n_markers > 1, that means use the same number of markers for all cell types
# rather than a relative amount, so a vector with n_markers repeated needs to
# be passed in, but accounting for cases where a cell type has less than
# n_markers genes in its marker list.
#
# Arguments:
#   markers = a named list, where the names are cell types and each item in
#             the list is a vector of genes that dtangle/hspe have determined
#             to be markers for that cell type
#   n_markers = a number that can either be a decimal that is > 0 and <= 1,
#               to signify to use that fraction of markers for each cell type,
#               or a whole number > 1 to signify using that same number of
#               markers for each cell type (or less, if a cell type has fewer)
#
# Returns:
#   a vector of whole numbers with the length of each marker list. If n_markers
#   was a fraction, this fraction is converted to an integer number of markers
#   for each cell type. This is done because the dtangle/HSPE code rounds DOWN
#   for fractional numbers of markers, which can sometimes result in 0 markers
#   for a cell type.
Get_NMarkers <- function(markers, n_markers) {
  n_markers_new <- n_markers

  # dtangle/hspe rounds down if length * n_markers is a decimal, so we need to
  # input a list of the length of each marker set instead, rounded UP to ensure
  # at least one marker per cell type
  if (n_markers <= 1) {
    n_markers_new <- ceiling(lengths(markers) * n_markers)
  }
  else if (n_markers > 1) {
    n_markers_new <- sapply(lengths(markers), min, n_markers)
  }

  return(n_markers_new)
}


# Check_NMarkers - checks that after filtering for genes that exist in both
# data sets and subsetting to n_markers, there are still enough markers to be
# useful and not too many markers to confuse the algorithm.
#
# Arguments:
#   n_markers_orig - the original n_markers argument from the parameter set,
#                    which may be a whole number or a fraction
#   n_markers - a vector of the number of markers per cell type (integers),
#               which was calculated using n_markers_orig
#   err_string - The header for the error to print out if a check fails
#
# Returns:
#   TRUE if n_markers passed all checks, FALSE otherwise
Check_NMarkers <- function(n_markers_orig, n_markers, err_string) {
  # For whole-number n_markers_orig arguments, the n_markers argument (mostly)
  # doubles each time. If there aren't enough markers in each cell type to do
  # anything new with this n_markers value, skip testing.
  if (n_markers_orig > 1 & all(n_markers <= (n_markers_orig/2))) {
    print(paste0(err_string, "Not enough total markers. ***"))
    return(FALSE)
  }

  # If there's less than ~3 markers per cell type, this isn't useful
  # information. Taking the mean still allows some leeway for rarer cell types
  # to have less than 3 markers as long as the other cell types have enough.
  if (mean(n_markers) < 3) {
    print(paste0(err_string, "Not enough markers per cell type. ***"))
    return(FALSE)
  }

  # Early testing showed it's not useful to use a large number of markers, so
  # if the total markers used is > 5000, don't test. 5000 is pretty generous,
  # and allows for 500 markers for 10 cell types.
  if (sum(n_markers) > 5000) {
    print(paste0(err_string, "Marker set too large. ***"))
    return(FALSE)
  }

  # All checks passed
  return(TRUE)
}
