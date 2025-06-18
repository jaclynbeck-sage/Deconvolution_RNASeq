# Configuration for running Dtangle in Step09
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
#           For Dtangle, as many cores as will fit in RAM. Assume between 5-20 GB
#           of RAM needed per core, depending on the dataset.
alg_config <- list(
  # Used to create params_loop1 (data-specific arguments) in main function
  normalizations = c("log_cpm", "log_tmm", "log_tpm"), # Dtangle works on log-scale data
  reference_input_types = c("pseudobulk"),

  # For params_loop2 (algorithm-specific arguments) in main function
  params_markers = CreateParams_MarkerTypes(), # all args are default
  additional_args = NULL, # dtangle has additional arguments but we aren't testing them

  # Define the function that runs each param set
  inner_loop_file = file.path("functions", "Step09_DtangleHSPE_InnerLoop.R"),
  inner_loop_func = "DtangleHSPE_InnerLoop", # must be a string
  required_libraries = c("dtangle"),

  # How many cores to use in parallel
  cores = parallel::detectCores() - 1 # Machine this was set up for has 16 cores
)
