library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step09_Deconvolution_HelperFunctions.R"))

# Parameter setup --------------------------------------------------------------

# options: "CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "Music", "Scaden"
algorithm <- "CibersortX"

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

ct_ad_only <- TRUE

# Force these three algorithms to use all data all the time. For the paper, I
# also ran DeconRNASeq and Dtangle with ct_ad_only = FALSE.
if ((algorithm %in% c("CibersortX", "HSPE", "Scaden"))) {
  ct_ad_only <- FALSE
}

source(file.path("algorithm_configs", str_glue("{algorithm}_Config.R")))

# datasets and normalization parameters
params_loop1 <- tidyr::expand_grid(
  algorithm = algorithm,
  reference_data_name = datasets,
  test_data_name = c("Mayo", "MSBB", "ROSMAP"),
  granularity = c("broad_class", "sub_class"),
  reference_input_type = alg_config$reference_input_types,
  normalization = alg_config$normalizations,
  regression_method = c("none", "edger", "lme", "dream")
) %>%
  arrange(normalization)

# Algorithm-specific parameters -- marker types, number of markers, plus any
# additional arguments to the algorithm itself. alg_config$additional_args can
# be NULL if no other additional args.
params_loop2 <- tidyr::expand_grid(alg_config$params_markers,
                                   alg_config$additional_args)


# Parallel execution setup -----------------------------------------------------

# NOTE: Recommendations for number of cores, per algorithm:
#   DWLS: 1/2 the available cores, as this doesn't need much RAM but multi-threads a little
#   DeconRNASeq and Music: 1/4 to 1/2 the available cores, as these algorithms multi-thread
#   Dtangle and HSPE: as many cores as will fit in RAM. Assume between 5-20 GB
#                     of RAM needed per core, depending on the dataset.
#   CibersortX and Scaden: only 1 core because of the memory usage
# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- alg_config$cores
cl <- makeCluster(cores, type = "FORK", outfile = str_glue("{algorithm}_output.txt"))
registerDoParallel(cl)

required_libraries <- alg_config$required_libraries


# Iterate through parameters ---------------------------------------------------

# Outer loop - each row of params_loop1 represents a single/unique call to
# Load_AlgorithmInputData. The inner loop then runs through all parameter sets
# on that data.
# NOTE: the helper functions have to be sourced inside the foreach loop
#       so they exist in each newly-created parallel environment

for (P in 1:nrow(params_loop1)) {
  data <- Load_AlgorithmInputData_FromParams(params_loop1[P, ])

  if (ct_ad_only) {
    message("Running CT and AD cases only.")
    data$test <- data$test[, data$test$diagnosis %in% c("CT", "AD")]
  } else {
    message("Running with all samples.")
  }

  bulk_metadata <- colData(data$test)
  data$test <- as.matrix(assay(data$test, "counts"))

  # Extra processing for CibersortX: Some re-formatting of the input, plus
  # re-use of batch-corrected signature if it exists
  if (algorithm == "CibersortX") {
    data <- Modify_CibersortX_Input(data, params_loop1[P, ])

  } else if (algorithm %in% c("Dtangle", "HSPE")) {
    # Extra pre-processing needed for Dtangle/HSPE -- reformat the input
    data <- Modify_DtangleHSPE_Input(data, params_loop1[P, ])

  } else if (algorithm == "Music") {
    # Extra pre-processing needed for MuSiC -- calculate or load sc_basis
    data <- Modify_Music_Input(data, params_loop1[P,])

  } else if (algorithm == "Scaden") {
    # Scaden needs to run each tissue separately so we need to bring that info
    # into the inner loop
    data$bulk_metadata <- bulk_metadata
  }

  ## Loop through algorithm-specific arguments ---------------------------------
  # Inner loop - each row of params_loop2 represents a single/unique call to
  # the algorithm with specific parameters like which markers to use, how many
  # from each cell type, and any changes to arguments in the function call.

  results_list <- foreach(R = 1:nrow(params_loop2), .packages = required_libraries) %dopar% {
    source(file.path("functions", "General_HelperFunctions.R"))
    source(file.path("functions", "Step09_ArgumentChecking_HelperFunctions.R"))
    source(alg_config$inner_loop_file) # defined in the config

    params <- cbind(params_loop1[P, ], params_loop2[R, ])

    # There are some invalid parameter combinations for CibersortX: if the
    # reference_input_type = cibersortx, only use filter_level = 0 because
    # CibersortX already filtered its signature. If reference_input_type =
    # signature, only use filter_level = 3 because CibersortX expects a filtered
    # signature.
    if (algorithm == "CibersortX") {
      if ((params$reference_input_type == "cibersortx") && (params$filter_level != 0)) {
        return(NULL)
      } else if ((params$reference_input_type == "signature") && (params$filter_level != 3)) {
        return(NULL)
      }
    }

    # If we are picking up from a failed/crashed run, and we've already run
    # this parameter set, load the result instead of re-running the algorithm
    prev_res <- Load_AlgorithmIntermediate(params)
    if (!is.null(prev_res)) {
      message(paste0("Using previously-run result for ",
                     paste(params, collapse = "_")))
      prev_res$param_id <- paste(paste(params_loop1[P, ], collapse = "_"),
                                 R, sep = "_")
      return(prev_res)
    }

    # Otherwise, call the algorithm-specific function to run it with this
    # set of parameters
    set.seed(12345)
    inner_loop_func <- match.fun(alg_config$inner_loop_func)

    res <- inner_loop_func(data, params)

    if (!is.null(res)) {
      res$param_id <- paste(paste(params_loop1[P, ], collapse = "_"),
                            R, sep = "_")
    }

    # Save each result in case of crashing
    Save_AlgorithmIntermediate(res)
    return(res)
  } # end foreach loop

  # It's possible for some items in results_list to be null if there was an error.
  # Filter them out.
  results_list <- results_list[lengths(results_list) > 0]

  # Give every result in the list a unique name
  names(results_list) <- sapply(results_list, "[[", "param_id")
  results_list <- results_list[sort(names(results_list))]

  # Save the completed list
  name_base <- paste(params_loop1[P, ], collapse = "_")
  print(str_glue("Saving final list for {name_base}..."))
  Save_AlgorithmOutputList(results_list, algorithm,
                           test_dataset = params_loop1$test_data_name[P],
                           name_base = name_base)

  rm(results_list, data)
  gc()
}

stopCluster(cl)
