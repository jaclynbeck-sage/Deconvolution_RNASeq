library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))

#### Parameter setup #####

# options: "DWLS", "DeconRNASeq", "Dtangle", "HSPE", "Music"
algorithm <- "DWLS"

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

source(file.path("algorithm_configs", str_glue("{algorithm}_Config.R")))

# datasets and normalization parameters
params_loop1 <- expand_grid(reference_data_name = datasets,
                            test_data_name = c("Mayo", "MSBB", "ROSMAP"),
                            granularity = c("broad_class"),
                            reference_input_type = alg_config$reference_input_types,
                            normalization = alg_config$normalizations,
                            regression_method = c("none")) %>% #, "edger", "deseq2", "dream")) %>%
                    arrange(normalization)

# Algorithm-specific parameters -- marker types, number of markers, plus any
# additional arguments to the algorithm itself
params_loop2 <- expand_grid(alg_config$params_markers,
                            alg_config$additional_args) # this value can be NULL if no other additional args


##### Parallel execution setup #####

# NOTE: Recommendations for number of cores, per algorithm:
#   DWLS: 1/2 the available cores, as this doesn't need much RAM but multi-threads a little
#   deconRNASeq and music: 1/4 to 1/2 the available cores, as these algorithms multi-thread
#   dtangle and hspe: as many cores as will fit in RAM. Assume between 5-20 GB
#                     of RAM needed per core, depending on the dataset.
# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- alg_config$cores
cl <- makeCluster(cores, type = "FORK", outfile = str_glue("{algorithm}_output.txt"))
registerDoParallel(cl)

required_libraries <- alg_config$required_libraries


#### Iterate through parameters ####

# Outer loop - each row of params_loop1 represents a single/unique call to
# Load_AlgorithmInputData in its own thread. The inner loop then runs through
# all parameter sets on that data.
# NOTE: the helper functions have to be sourced inside the foreach loop
#       so they exist in each newly-created parallel environment

for (P in 1:nrow(params_loop1)) {
  data <- Load_AlgorithmInputData_FromParams(params_loop1[P,])

  data$test <- as.matrix(assay(data$test, "counts"))

  # Some extra pre-processing needed for MuSiC -- calculate sc_basis
  if (algorithm == "Music") {
    data$reference <- as(data$reference, "SingleCellExperiment")

    # Pre-compute sc.basis to save time
    sc_basis <- MuSiC::music_basis(data$reference, non.zero = TRUE,
                                   markers = rownames(data$reference),
                                   clusters = "celltype", samples = "sample",
                                   select.ct = NULL,
                                   ct.cov = FALSE,
                                   verbose = TRUE)
  }

  # Some extra pre-processing needed for Dtangle/HSPE -- reformat the input
  if (algorithm %in% c("Dtangle", "HSPE")) {
    metadata <- colData(data$reference)

    celltypes <- levels(metadata$celltype)
    pure_samples <- lapply(celltypes, function(ct) {
      which(metadata$celltype == ct)
    })
    names(pure_samples) <- celltypes

    data$reference <- assay(data$reference, "counts")

    # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
    # Input data must be first so indices in pure_samples are correct.
    Y <- t(cbind(data$reference, data$test))
  }


  ##### Loop through algorithm-specific arguments #####
  # Inner loop - each row of params_loop2 represents a single/unique call to
  # the algorithm with specific parameters like which markers to use, how many
  # from each cell type, and any changes to arguments in the function call.

  results_list <- foreach(R = 1:nrow(params_loop2), .packages = required_libraries) %dopar% { #lapply(1:nrow(params_loop2), function(R) {
    source(file.path("functions", "General_HelperFunctions.R"))
    source(file.path("functions", "Step09_ArgumentChecking_HelperFunctions.R"))
    source(alg_config$inner_loop_file) # defined in the config

    params <- cbind(params_loop1[P,], params_loop2[R,])

    # For backward compatibility -- for algorithms that only input one kind of
    # data (DeconRNASeq and DWLS are signature-only), remove the
    # reference_input_type variable since it wasn't there originally for these
    # algorithms
    if ("signature" %in% params_loop1$reference_input_type) {
      params <- params %>% select(-reference_input_type)
    }

    # If we are picking up from a failed/crashed run, and we've already run
    # this parameter set, load the result instead of re-running the algorithm
    prev_res <- Load_AlgorithmIntermediate(algorithm, params)
    if (!is.null(prev_res)) {
      message(paste0("Using previously-run result for ", paste(params, collapse = " ")))
      return(prev_res)
    }

    # Otherwise, call the algorithm-specific function to run it with this
    # set of parameters
    set.seed(12345)
    inner_loop_func <- match.fun(alg_config$inner_loop_func)

    # Music needs sc_basis passed in
    if (algorithm == "Music") {
      res <- inner_loop_func(data$reference, data$test, sc_basis, params,
                             verbose = FALSE)
    }
    # Dtangle/HSPE need "Y" and "pure_samples" passed in
    else if (algorithm %in% c("Dtangle", "HSPE")) {
      res <- inner_loop_func(Y, pure_samples, params, algorithm)
    }
    else {
      res <- inner_loop_func(data$reference, data$test, params, algorithm)
    }

    # Save each result in case of crashing
    Save_AlgorithmIntermediate(res, algorithm)
    return(res)
  } # end foreach loop

  # It's possible for some items in results_list to be null if there was an error.
  # Filter them out.
  results_list <- results_list[lengths(results_list) > 0]

  # Give every result in the list a unique name
  name_base <- paste(params_loop1[P,], collapse = "_")
  names(results_list) <- paste(algorithm, name_base,
                               1:length(results_list), sep = "_")

  # Save the completed list
  print(str_glue("Saving final list for {name_base}..."))
  Save_AlgorithmOutputList(results_list, algorithm,
                           test_dataset = params_loop1$test_data_name[P],
                           name_base = name_base)

  rm(results_list, data)
  gc()
}

stopCluster(cl)
