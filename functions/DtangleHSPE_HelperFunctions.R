# Helper functions for running Dtangle and HSPE. The functions in this file
# help clear up clutter in the main algorithm script and make it clearer
# what is happening there.

source(file.path("functions", "FileIO_HelperFunctions.R"))

# Get_DtangleHSPEInput: loads in single cell or pseudobulk pure sample data as
# the reference set, loads in the test pseudobulk data, and combines the
# reference and test data in the format dtangle/hspe expects.
#
# Arguments:
#   dataset = the name of the data set to load
#   datatype = either "donors" or "training", for which type of test pseudobulk
#              data to load
#   granularity = either "broad" or "fine", for which level of cell types to
#                 use as the reference
#   input_type = either "singlecell" or "pseudobulk", for whether the reference
#                data is single cell data or pseudobulk pure samples
#
# Returns:
#   a list with entries for "Y", which is the combined reference + test data,
#   and "pure_samples", which is a list of indices into Y that correspond to
#   samples from each cell type.
Get_DtangleHSPEInput <- function(dataset, datatype, granularity, input_type) {

  ##### Reference data set #####
  if (input_type == "singlecell") {
    input_obj <- Load_SingleCell(dataset, granularity, output_type = "logcpm")
  }
  else { # Input is pseudobulk pure samples
    input_obj <- Load_PseudobulkPureSamples(dataset, granularity,
                                            output_type = "logcpm")
  }

  input_mat <- assay(input_obj, "counts")
  metadata <- colData(input_obj)

  ##### Get indices of pure samples #####
  celltypes <- levels(metadata$celltype)
  pure_samples <- lapply(celltypes, function(ct) {
    which(metadata$celltype == ct)
  })
  names(pure_samples) <- celltypes

  ##### Test data #####

  # Ground truth pseudobulk sets
  if (datatype == "donors" | datatype == "training") {
    pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, "logcpm")
    bulk_mat <- assay(pseudobulk, "counts")
  }
  # ROSMAP, Mayo, or MSBB
  else {
    genes <- Load_GeneConversion(dataset)
    bulk_mat <- Load_BulkData(datatype, genes, output_type = "logcpm")
  }

  # These SHOULD have the same rownames, but just in case.
  keepgene <- intersect(rownames(input_mat), rownames(bulk_mat))

  # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
  # Input data must be first so indices in pure_samples are correct.
  Y <- t(cbind(input_mat[keepgene,],
               bulk_mat[keepgene,]))

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
#   a single number (for fraction) or a vector of whole numbers (for "all
#   markers" or "same number of markers" conditions)
Get_NMarkers <- function(markers, n_markers) {
  n_markers_new <- n_markers

  # dtangle/hspe don't interpret "1" as 100%, so we need to input a list of
  # the length of each marker set instead
  if (n_markers == 1) {
    n_markers_new <- lengths(markers)
  }
  else if (n_markers > 1) {
    n_markers_new <- sapply(lengths(markers), min, n_markers)
  }

  return(n_markers_new)
}
