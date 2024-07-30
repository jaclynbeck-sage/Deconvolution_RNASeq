# Configuration for running MuSiC in Step09
# Variables that need to be defined are:
#   normalizations - vector with any combination of 'counts', 'cpm', 'tmm',
#                    'tpm', 'log_cpm', 'log_tmm', or 'log_tpm', depending on
#                    which scale the algorithm expects values to be on and what
#                    it makes sense to run with the algorithm
#   reference_input_types - vector with any of 'signature', 'singlecell', or
#                    'pseudobulk' (for inputting a signature matrix, single cell
#                    data, or pseudobulk of the single cell data, respectively).
#                    DWLS and DeconRNASeq use 'signature' only, the other
#                    algorithms can use 'singlecell' or 'pseudobulk', or both of
#                    those.
#   params_markers - a matrix of combinations of marker-related arguments as
#                    output by 'CreateParams_FilterableSignature' or
#                    'CreateParams_MarkerTypes'. See those functions for input
#                    options. Any arguments not defined in the function call in
#                    this config file will be set to default values (see the
#                    functions for defaults).
#   additional_args - a matrix of combinations of algorithm-related arguments
#                    that are for the algorithm's function call. If the function
#                    doesn't need other arguments or you don't want to test
#                    them, this value should be set to NULL.
#   inner_loop_file - which file defines the algorithm-specific inner loop function
#   inner_loop_func - a string, the name of the inner loop function to call
#   required_libraries - a vector of libraries needed for execution of the
#                    inner loop function, which is usually just the algorithm's
#                    library
#   cores - an integer value of how many CPUs to use when running in parallel.
#           For MuSiC, use 1/4 to 1/2 the available cores, as this algorithm
#           multi-threads within each process
alg_config <- list(
  # Used to create params_loop1 (data-specific arguments) in main function
  normalizations = c("counts"), # MuSiC works on linear-scale data
  reference_input_types = c("singlecell"),

  # For params_loop2 (algorithm-specific arguments) in main function.
  # MuSiC was designed to use the whole set of genes, but I am testing whether
  # narrowing the gene set helps. Although MuSiC inputs singlecell/pseudobulk
  # data directly and not a signature, calling this function is the easiest
  # way to get the genes to use both for the whole gene set and for a subset
  # of markers.
  params_markers = CreateParams_FilterableSignature(
                      filter_level = c(0, 1, 2, 3),
                      n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                                    3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
                    ), # all other args are default
  additional_args = expand_grid(ct.cov = c(FALSE), # Their ct.cov methods don't work as intended and/or take forever to converge
                                centered = c(TRUE, FALSE),
                                normalize = c(TRUE, FALSE)),

  # Define the function that runs each param set
  inner_loop_file = file.path("functions", "Step09_Music_InnerLoop.R"),
  inner_loop_func = "Music_InnerLoop", # must be a string
  required_libraries = c("MuSiC", "SingleCellExperiment"),

  # How many cores to use in parallel
  cores = 6 # Machine this was set up for has 16 cores
)
