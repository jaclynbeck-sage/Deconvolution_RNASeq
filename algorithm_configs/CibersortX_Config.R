# Configuration for running CibersortX in Step09
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
#           For CibersortX, this should be 1 since it can't run in parallel. It
#           uses too much memory and disk space to be able to stack processes.
alg_config <- list(
  # Used to create params_loop1 (data-specific arguments) in main function
  normalizations = c("cpm", "tpm"), # CibersortX works on linear-scale data
  reference_input_types = c("cibersortx", "signature"), # CibersortX can use a signature matrix it creates,
                                                        # or the signature created by this pipeline

  # For params_loop2 (algorithm-specific arguments) in main function
  # Note: If the input type is "cibersortx", only filter_level = 0 will be used.
  params_markers = CreateParams_FilterableSignature(filter_levels = c(0, 3)), # All other args are default

  additional_args = expand_grid(batch_correct = c(TRUE, FALSE)), # Use rmbatchSmode to correct signature or not

  # Define the function that runs each param set
  inner_loop_file = file.path("functions", "Step09_CibersortX_InnerLoop.R"),
  inner_loop_func = "CibersortX_InnerLoop", # must be a string
  required_libraries = c("omnideconv"),

  # How many cores to use in parallel
  cores = 1 # Only use more than one core if batch-adjusted signature matrices are pre-computed
)
