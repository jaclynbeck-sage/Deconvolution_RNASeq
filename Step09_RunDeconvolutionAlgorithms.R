library(foreach)
library(doParallel)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))

# Parameter setup --------------------------------------------------------------

# options: "CibersortX", "DWLS", "DeconRNASeq", "Dtangle", "HSPE", "Music", "Scaden"
algorithm <- "CibersortX"

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

use_top_params <- FALSE
ct_ad_only <- TRUE

# If we're using the top parameter sets, we have to use all data.
# Force these three algorithms to use all data all the time. For the paper, I
# also ran DeconRNASeq and Dtangle with ct_ad_only = FALSE.
if (use_top_params || (algorithm %in% c("CibersortX", "HSPE", "Scaden"))) {
  ct_ad_only <- FALSE
}

source(file.path("algorithm_configs", str_glue("{algorithm}_Config.R")))

# datasets and normalization parameters
params_loop1 <- expand_grid(reference_data_name = datasets,
                            test_data_name = c("Mayo", "MSBB", "ROSMAP"),
                            granularity = c("broad_class", "sub_class"),
                            reference_input_type = alg_config$reference_input_types,
                            normalization = alg_config$normalizations,
                            regression_method = c("none", "edger", "lme", "dream")) %>%
  arrange(normalization)

# Algorithm-specific parameters -- marker types, number of markers, plus any
# additional arguments to the algorithm itself
params_loop2 <- expand_grid(alg_config$params_markers,
                            alg_config$additional_args) # this value can be NULL if no other additional args


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
  }

  data$test <- as.matrix(assay(data$test, "counts"))

  # Some extra pre-processing needed for MuSiC -- calculate sc_basis
  if (algorithm == "Music") {
    data$reference <- as(data$reference, "SingleCellExperiment")

    sc_basis <- Load_MusicBasis(params_loop1$reference_data_name[P],
                                params_loop1$granularity[P])

    if (is.null(sc_basis)) {
      # Pre-compute sc.basis to save time. Use all genes in the single cell
      # data, not just the ones that exist in both bulk and sc data sets
      sce_unfiltered <- Load_SingleCell(dataset = params_loop1$reference_data_name[P],
                                        granularity = params_loop1$granularity[P],
                                        output_type = params_loop1$normalization[P])

      sc_basis <- MuSiC::music_basis(sce_unfiltered,
                                     non.zero = TRUE,
                                     markers = rownames(sce_unfiltered),
                                     clusters = "celltype",
                                     samples = "sample",
                                     select.ct = NULL,
                                     ct.cov = FALSE,
                                     verbose = TRUE)

      # Save for later use
      Save_MusicBasis(sc_basis,
                      params_loop1$reference_data_name[P],
                      params_loop1$granularity[P])

      rm(sce_unfiltered)
      gc()
    } else {
      message("Using pre-computed sc_basis")
    }
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

    rm(data)
    gc()

    # Load the subset of parameters that were determined by Dtangle
    if (algorithm == "HSPE") {
      params_filename <- paste0("hspe_params_",
                                paste(params_loop1[P, ], collapse = "_"),
                                ".rds")
      params_filename <- file.path(dir_hspe_params, params_filename)

      if (file.exists(params_filename)) {
        params_subset <- readRDS(params_filename)
        params_subset$string_id <- apply(params_subset, 1, paste, collapse = "_")
      } else {
        message(str_glue("WARNING: {params_filename} does not exist. Skipping..."))
        next # We have to skip this file
      }
    }
  }

  # Extra processing for CibersortX: If there is a pre-computed corrected
  # signature, copy it to the directory that the Cibersort docker container can
  # read. If not, create a CibersortX-formatted data file for single cell
  # reference input so it can compute the adjusted signature, then replace the
  # reference object with the filename. Cell type names also can't have any "."
  # in them, so this does a string replace to change them to "_". To avoid
  # differences in how R and C sort strings with multiple capital letters (i.e.
  # "Oligodendrocyte" vs "OPC" get sorted differently between the two
  # languages), we also change all cell type names to lower case.
  if (algorithm == "CibersortX") {
    # Copy adjusted signature file to cibersort directory if the file exists
    sig_params <- params_loop1[P, ] %>% select(-reference_input_type)
    sig_filename <- list.files(dir_cibersort_corrected_signatures,
                               pattern = paste(sig_params, collapse = "_"),
                               full.names = TRUE)

    if (length(sig_filename) == 1) {
      file.copy(from = sig_filename,
                to = file.path(dir_cibersort, basename(sig_filename)),
                overwrite = TRUE)
      data$singlecell_filename <- ""
    } else {
      # Otherwise CibersortX needs the original single cell data to compute an
      # adjusted signature
      message("Copying single cell data to CibersortX directory...")
      sce <- Load_SingleCell(dataset = params_loop1$reference_data_name[P],
                             granularity = params_loop1$granularity[P],
                             output_type = params_loop1$normalization[P])
      sce$celltype <- str_replace(as.character(sce$celltype), "\\.", "_")
      sce$celltype <- str_to_lower(sce$celltype)

      f_name <- Save_SingleCellToCibersort(sce = sce,
                                           dataset_name = params_loop1$reference_data_name[P],
                                           granularity = params_loop1$granularity[P])
      data$singlecell_filename <- f_name
      rm(sce)
      gc()
    }

    # Lower-casing the reference signature has to be handled in the inner loop
    # so the inner loop can convert back to the original names
    colnames(data$reference) <- str_replace(colnames(data$reference), "\\.", "_")
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

    # For HSPE, skip this parameter set if it isn't in the specified subset of
    # parameters to run
    if (algorithm == "HSPE") {
      if (nrow(plyr::match_df(params, params_subset, on = colnames(params))) == 0) {
        return(NULL)
      }
    }

    # There are some invalid parameter combinations for CibersortX: if the
    # reference_input_type = cibersortx, only use filter_level = 0 because
    # CibersortX already filtered its signature. If reference_input_type =
    # signature, only use filter_level = 3 because CibersortX expects a filtered
    # signature.
    if (algorithm == "CibersortX") {
      if ((params$reference_input_type == "cibersortx") & (params$filter_level != 0)) {
        return(NULL)
      } else if ((params$reference_input_type == "signature") & (params$filter_level != 3)) {
        return(NULL)
      }
    }

    # If we are picking up from a failed/crashed run, and we've already run
    # this parameter set, load the result instead of re-running the algorithm
    prev_res <- Load_AlgorithmIntermediate(algorithm, params)
    if (!is.null(prev_res)) {
      message(paste0("Using previously-run result for ",
                     paste(params, collapse = " ")))
      prev_res$param_id <- paste(algorithm,
                                 paste(params_loop1[P, ], collapse = "_"),
                                 R, sep = "_")
      return(prev_res)
    }

    # Otherwise, call the algorithm-specific function to run it with this
    # set of parameters
    set.seed(12345)
    inner_loop_func <- match.fun(alg_config$inner_loop_func)

    # Music needs sc_basis passed in
    if (algorithm == "Music") {
      res <- inner_loop_func(data$reference, data$test,
                             sc_basis, params, verbose = FALSE)
    }
    # Dtangle/HSPE need "Y" and "pure_samples" passed in
    else if (algorithm %in% c("Dtangle", "HSPE")) {
      res <- inner_loop_func(Y, pure_samples, params, algorithm)
    }
    # CibersortX needs the signature and single cell filename passed in
    else if (algorithm == "CibersortX") {
      res <- inner_loop_func(data$reference, data$test,
                             data$singlecell_filename, params)
    } else {
      res <- inner_loop_func(data$reference, data$test, params, algorithm)
    }

    if (!is.null(res)) {
      res$param_id <- paste(algorithm,
                            paste(params_loop1[P, ], collapse = "_"),
                            R, sep = "_")
    }

    # Save each result in case of crashing
    Save_AlgorithmIntermediate(res, algorithm)
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
