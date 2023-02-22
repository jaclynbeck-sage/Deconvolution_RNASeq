# Runs Dtangle on a variety of parameters and saves the list of results to a
# file.
#
# NOTE: This script relies on Dtangle markers, so these must be calculated
# before-hand.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(foreach)
library(doParallel)

##### Parallel execution setup #####

cores <- 4
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("dtangleSparse", "Matrix", "SummarizedExperiment",
                        "SingleCellExperiment", "stringr", "dplyr")

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito")#,
              "seaRef") #, "seaAD")

input_types = c("singlecell", "pseudobulk")

params_loop1 <- expand.grid(dataset = datasets,
                            granularity = c("broad"), #, "fine"),
                            datatype = c("donors", "training"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

params_loop2 <- expand.grid(marker_method = c("ratio", "diff", "p.value", "regression"),
                            gamma_name = c("auto", 0.5),
                            sum_fn_type = c("mean", "median"),
                            n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2,
                                          0.5, 0.75, 1.0),
                          stringsAsFactors = FALSE)

#### Iterate through parameters in parallel ####

foreach (P = 1:nrow(params_loop1), .packages = required_libraries) %dopar% {
  # This needs to be sourced inside the loop for parallel processing
  source(file.path("functions",  "FileIO_HelperFunctions.R"))

  dataset <- params_loop1$dataset[P]
  granularity <- params_loop1$granularity[P]
  datatype <- params_loop1$datatype[P]

  # Each dataset / datatype / granularity combo gets its own list
  dtangle_list <- list()

  ##### Run on both single cell and pseudobulk input matrices #####
  # NOTE: input_types and all parameters from params_loop2 are executed in the
  # same thread because they all use the same input data, which can be large
  # and it isn't feasible to have multiple copies of it between multiple threads.

  for (input_type in input_types) {

    ##### Prepare input matrix #####

    if (input_type == "singlecell") {
      input_obj <- Load_SingleCell(dataset, granularity, output_type = "logcpm")
    }
    else { # Input is pseudobulk pure samples
      input_obj <- Load_PseudobulkPureSamples(dataset, granularity,
                                              output_type = "logcpm")
    }

    metadata <- colData(input_obj)
    input_mat <- assay(input_obj, "counts")

    # Clear up as much memory as possible
    rm(input_obj)
    gc()

    celltypes <- levels(metadata$celltype)
    pure_samples <- lapply(celltypes, function(ct) {
      which(metadata$celltype == ct)
    })
    names(pure_samples) <- celltypes

    # TEMPORARY: dtangle code below will not work with a DelayedArray.
    # The seaRef dataset will fit in memory all at once, so this converts it
    # to a sparse matrix. The seaAD data set will NOT fit so this won't work on it.
    if (is(input_mat, "DelayedArray")) {
      input_mat <- as(input_mat, "CsparseMatrix")
    }

    ##### Prepare test data matrix #####

    pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, "logcpm")
    bulk_mat <- assays(pseudobulk)[["counts"]]

    # Clear up as much memory as possible
    rm(pseudobulk)
    gc()

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(input_mat), rownames(bulk_mat))

    # Filter params_loop2 down: Early testing showed we don't need to test
    # higher percentages for method = "ratio", since this method returns the
    # full set of genes and smaller lists do better.
    params_run <- subset(params_loop2, !(marker_method == "ratio" & n_markers > 0.2))

    # "p.value" and "regression" aren't feasible for single cell input
    if (input_type == "singlecell") {
      params_run <- subset(params_run, marker_method == "ratio" | marker_method == "diff")
    }

    # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
    # Input data must be first so indices in pure_samples are correct.
    Y <- t(cbind(input_mat[keepgene,],
                 bulk_mat[keepgene,]))
    gc()

    ##### Iterate through Dtangle parameters #####

    for (R in 1:nrow(params_run)) {
      marker_method <- params_run$marker_method[R]
      gamma_name <- params_run$gamma_name[R]
      sum_fn_type <- params_run$sum_fn_type[R]
      n_markers <- params_run$n_markers[R]

      gamma <- NULL
      if (gamma_name != "auto") {
        gamma <- as.numeric(gamma_name)
      }

      sum_fn <- mean
      if (sum_fn_type == "median") {
        sum_fn <- median
      }

      markers <- Load_DtangleMarkers(dataset, granularity, input_type,
                                     marker_method)
      n_markers_name <- n_markers

      # dtangle doesn't interpret "1" as 100%, so we need to input a list of
      # the length of each marker set instead
      if (n_markers == 1) {
        n_markers_name <- "all"
        n_markers <- lengths(markers$L)
      }

      name <- str_glue(paste0("{dataset}_{granularity}_{datatype}_",
                              "{input_type}_{R}"))

      result <- dtangle(Y = Y,
                        pure_samples = pure_samples,
                        data_type = "rna-seq",
                        gamma = gamma, # If gamma is not NULL, it will override data_type argument
                        n_markers = n_markers,
                        markers = markers,
                        summary_fn = sum_fn)

      # Only keep results for pseudobulk samples
      result$estimates <- result$estimates[colnames(bulk_mat), ]

      # Add the params we used to generate this run
      result$params <- cbind(params_loop1[P,], "input_type" = input_type,
                             params_run[R,])

      dtangle_list[[name]] <- result

      gc()
      print(paste(result$params, collapse = "  "))
    } # End params_run loop

    # Next iteration will start with new data, remove the old data
    rm(Y)
    gc()

    # Save the completed list
    print("Saving final list...")
    Save_AlgorithmOutputList(dtangle_list, dataset, datatype, granularity)
  } # end input_types loop

  return(NULL)
} # end params_loop1 parallel loop

stopCluster(cl)
