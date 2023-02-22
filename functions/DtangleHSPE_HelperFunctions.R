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

  # TEMPORARY: dtangle code will not work with a DelayedArray.
  # The seaRef dataset will fit in memory all at once, so this converts it
  # to a sparse matrix. The seaAD data set will NOT fit so this won't work on it.
  if (is(input_mat, "DelayedArray")) {
    input_mat <- as(input_mat, "CsparseMatrix")
  }

  ##### Get indices of pure samples #####
  celltypes <- levels(metadata$celltype)
  pure_samples <- lapply(celltypes, function(ct) {
    which(metadata$celltype == ct)
  })
  names(pure_samples) <- celltypes

  ##### Test data #####

  pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, "logcpm")
  bulk_mat <- assay(pseudobulk, "counts")

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
