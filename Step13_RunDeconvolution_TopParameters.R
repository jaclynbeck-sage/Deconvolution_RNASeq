library(foreach)
library(doParallel)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step09_Deconvolution_HelperFunctions.R"))

# Parameter setup --------------------------------------------------------------

# options: "CibersortX", "DWLS", "DeconRNASeq", "Dtangle", "HSPE", "Music", "Scaden"
# Note that if an algorithm was run in Step09 with "ct_ad_only" = FALSE, all this
# script will do is subset the existing output to the results from the top
# parameter sets and save those as a new file.
algorithm <- "CibersortX"

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

source(file.path("algorithm_configs", str_glue("{algorithm}_Config.R")))

# datasets and normalization parameters
params_loop1 <- expand_grid(reference_data_name = datasets,
                            test_data_name = c("Mayo", "MSBB", "ROSMAP"),
                            granularity = c("broad_class", "sub_class"),
                            reference_input_type = alg_config$reference_input_types,
                            normalization = alg_config$normalizations,
                            regression_method = c("none", "edger", "lme", "dream"),
                            samples = "all_samples") %>%
  arrange(normalization)


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
  # Check if there is a file of top parameters -- if there isn't, either errors
  # haven't been calculated for this set of data or there were no valid estimates
  # for this set of data, so we can skip running this one.
  file_params <- select(params_loop1[P, ], -samples)
  file_id <- paste(c(algorithm, file_params), collapse = "_")
  top_param_file <- file.path(dir_top_parameters,
                              params_loop1$test_data_name[P],
                              algorithm,
                              str_glue("top_parameters_{file_id}.rds"))

  if (!file.exists(top_param_file)) {
    warning(paste(top_param_file, "doesn't exist! Skipping..."))
    next
  }

  top_params <- readRDS(top_param_file)

  # Use all the parameter sets in the file. We remove columns that also exist
  # in params_loop1 and also remove unused columns "total_markers_used" and
  # "algorithm".
  params_loop2 <- top_params$params %>%
    select(-reference_data_name, -test_data_name, -granularity,
           -reference_input_type, -normalization, -regression_method,
           -total_markers_used, -algorithm)

  # Load input data
  data <- Load_AlgorithmInputData_FromParams(params_loop1[P, ])

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
    data <- Modify_Music_Input(data, params_loop1[P, ])
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

    # Not very memory efficient but we don't want to alter "data" and affect
    # other threads using it
    data_filt <- data

    # If we've already run this script on this parameter set, the intermediate
    # file will have "all_samples" at the end of the filename to differentiate
    # it from Step09 results that don't have all samples.
    prev_res <- Load_AlgorithmIntermediate(algorithm, params)
    if (!is.null(prev_res)) {
      message(paste0("Using previously-run result for ",
                     paste(params, collapse = " ")))

      prev_res$param_id <- paste(algorithm,
                                 paste(params_loop1[P, ], collapse = "_"),
                                 R, sep = "_")

      return(prev_res)
    }

    # If we haven't already run this script, in most cases we have already run
    # Step09 and there may be a saved estimates file containing results for
    # all samples or CT/AD samples only. If this is the case, we only need to
    # re-run the algorithm on the samples not in the file, and concatenate the
    # results together, instead of re-running all samples.
    prev_res <- Load_AlgorithmOutputList(algorithm,
                                         params$reference_data_name,
                                         params$test_data_name,
                                         params$granularity,
                                         params$reference_input_type,
                                         params$normalization,
                                         params$regression_method)
    if (!is.null(prev_res)) {
      message(paste0("Found Step09 result for ",
                     paste(params, collapse = " ")))

      # Find which result matches these parameters
      prev_params <- do.call(rbind, lapply(prev_res, "[[", "params"))
      prev_params <- prev_params[, colnames(params_loop2)]
      match_params <- plyr::match_df(prev_params, params_loop2[R, ])

      if (nrow(match_params) != 1) {
        stop("Step09 results do not contain this parameter set.")
      }

      prev_res <- prev_res[[rownames(match_params)]]
      gc()

      message(paste("Found Step09 result containing", nrow(prev_res$estimates),
                    "samples."))

      if ("test" %in% names(data)) {
        missing_samples <- setdiff(colnames(data$test),
                                   rownames(prev_res$estimates))
        data_filt$test <- data$test[, missing_samples]
      } else if ("Y" %in% names(data)) {
        missing_samples <- setdiff(rownames(data$Y),
                                   rownames(prev_res$estimates))
        data_filt$Y <- data$Y[missing_samples, ]

        # Since "Y" has both pseudobulk and bulk samples in it, missing_samples
        # will include all pseudobulk samples as well as any missing bulk
        # samples. This is needed for subsetting data_filt, but when looking at
        # the length of missing_samples below it will look like more samples
        # need to be run. This piece of code gets how many bulk samples are
        # actually missing.
        missing_samples <- missing_samples[-unlist(data$pure_samples)]
      }

      # If this file already has all samples in it, no need to re-run
      if (length(missing_samples) == 0) {
        message("No additional samples need to be run.")
        prev_res$param_id <- paste(algorithm,
                                   paste(params_loop1[P, ], collapse = "_"),
                                   R, sep = "_")
        return(prev_res)
      }

      message(paste("Running", length(missing_samples), "additional samples."))
    }

    # Call the algorithm-specific function to run it with this set of parameters
    set.seed(12345)
    inner_loop_func <- match.fun(alg_config$inner_loop_func)

    res <- inner_loop_func(data_filt, params, algorithm)

    if (!is.null(res)) {
      res$param_id <- paste(algorithm,
                            paste(params_loop1[P, ], collapse = "_"),
                            R, sep = "_")

      # Add CT/AD samples back to the result, if applicable
      if (!is.null(prev_res)) {
        res$estimates <- rbind(res$estimates[, colnames(prev_res$estimates)],
                               prev_res$estimates)
      }
    }

    # Save each result in case of crashing
    Save_AlgorithmIntermediate(res, algorithm)
    gc()
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
                           name_base = name_base,
                           top_params = TRUE)

  rm(results_list, data)
  gc()
}

stopCluster(cl)
