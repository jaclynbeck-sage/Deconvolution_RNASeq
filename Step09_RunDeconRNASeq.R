# Runs DeconRNASeq on a variety of parameters and saves the list of results to a
# file.
#
# NOTE: This script relies on Dtangle markers, so these must be calculated
# before-hand.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the inner foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(SummarizedExperiment)
library(stringr)
library(reshape2)

source(file.path("functions", "General_HelperFunctions.R"))

##### Parallel execution setup #####

# NOTE: Assume about 5-20 GB RAM needed per core. DeconRNASeq multi-threads,
#       so only use ~half the cores available on the machine.
# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- 6
cl <- makeCluster(cores, type = "FORK", outfile = "deconRNAseq_output.txt")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("DeconRNASeq", "scuttle")

#### Parameter setup #####

datasets <- c("cain", "lau", "leng", "mathys", "seaRef") #, "morabito", "seaAD")

params_loop1 <- expand_grid(reference_data_name = datasets,
                            test_data_name = c("Mayo", "MSBB", "ROSMAP"),
                            granularity = c("broad_class"),
                            normalization = c("cpm", "tmm", "tpm"),
                            regression_method = c("none", "edger", "deseq2", "dream")) %>%
                  arrange(normalization)

marker_types <- list("dtangle" = c("ratio", "diff", "p.value", "regression"),
                     "autogenes" = c("correlation", "distance", "combined"),
                     "seurat" = c("None"),
                     "deseq2" = c("DESeq2"))

params_tmp <- CreateParams_FilterableSignature(
                filter_level = c(1, 2, 3),
                n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                              3, 5, 10, 20, 50, 100, 200),
                marker_types,
                marker_input_type = c("singlecell", "pseudobulk"),
                marker_order = c("distance", "correlation")
)

params_loop2 <- expand_grid(params_tmp,
                            use_scale = c(TRUE, FALSE),
                            recalc_cpm = c(FALSE)) #c(TRUE, FALSE))

#### Iterate through parameters ####

for (P in 1:nrow(params_loop1)) {
  reference_data_name <- params_loop1$reference_data_name[P]
  test_data_name <- params_loop1$test_data_name[P]
  granularity <- params_loop1$granularity[P]
  normalization <- params_loop1$normalization[P]
  regression_method <- params_loop1$regression_method[P]

  data <- Load_AlgorithmInputData(reference_data_name, test_data_name,
                                  granularity,
                                  reference_input_type = "signature",
                                  output_type = normalization,
                                  regression_method = regression_method)

  data$reference <- as.data.frame(data$reference)
  data$test <- as.data.frame(assay(data$test, "counts"))

  ##### Filter levels, number of markers, and DeconRNASeq arguments #####
  # NOTE: the helper functions have to be sourced inside the foreach loop
  #       so they exist in each newly-created parallel environment

  decon_list <- foreach (R = 1:nrow(params_loop2),
                         .packages = required_libraries) %dopar% {
    source(file.path("functions", "DeconRNASeq_InnerLoop.R"))
    source(file.path("functions", "General_HelperFunctions.R"))

    params = cbind(params_loop1[P,], params_loop2[R,])

    # If we are picking up from a failed/crashed run, and we've already run
    # this parameter set, load the result instead of re-running the algorithm
    prev_res <- Load_AlgorithmIntermediate("deconRNASeq", params)
    if (!is.null(prev_res)) {
      print(paste0("Using previously-run result for ",
                   paste(params, collapse = " ")))
      return(prev_res)
    }

    set.seed(12345)
    res <- DeconRNASeq_InnerLoop(data$reference, data$test,
                                 cbind(params_loop1[P,], params_loop2[R,]))

    Save_AlgorithmIntermediate(res, "deconRNASeq")
    return(res)
  }

  # It's possible for some items in decon_list to be null if there was an error.
  # Filter them out.
  decon_list <- decon_list[lengths(decon_list) > 0]

  name_base <- str_glue(paste0("{reference_data_name}_{test_data_name}_",
                               "{granularity}_{normalization}_{regression_method}"))
  names(decon_list) <- paste("deconRNASeq",
                              name_base,
                              1:length(decon_list), sep = "_")

  # Save the completed list
  print(str_glue("Saving final list for {name_base}..."))
  Save_AlgorithmOutputList(decon_list, "deconRNASeq", reference_data_name,
                           test_data_name, granularity, normalization,
                           regression_method)

  rm(decon_list, data)
  gc()
}

stopCluster(cl)
