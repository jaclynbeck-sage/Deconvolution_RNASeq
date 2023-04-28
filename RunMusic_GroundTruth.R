# Runs MuSiC on a variety of parameters and saves the list of results to a
# file.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))

##### Parallel execution setup #####

# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- 8
# Music is very verbose so we make an output file instead of putting on the
# console, because RStudio won't display everything in the console otherwise.
cl <- makeCluster(cores, type = "FORK", outfile = "music_output.txt")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("MuSiC", "SingleCellExperiment")

#### Parameter setup #####

datasets <- c("cain", "lau", "leng", "mathys", "morabito", "seaRef") #, "seaAD")

params_loop1 <- expand_grid(reference_data_name = datasets,
                            test_data_name = c("Mayo", "MSBB", "ROSMAP"), #c("donors", "training"),
                            granularity = c("broad")) %>% arrange(test_data_name)

marker_types <- list("dtangle" = c("diff"), #c("ratio", "diff", "p.value", "regression"),
                     "autogenes" = c("correlation", "distance", "combined"),
                     "seurat" = c("None"))

params_tmp <- CreateParams_FilterableSignature(
                  filter_level = c(0, 1, 2, 3),
                  n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                                5, 10, 20, 50, 100, 200, 500),
                  marker_types,
                  marker_input_type = c("singlecell", "pseudobulk")
)

params_loop2 <- expand_grid(params_tmp,
                            ct.cov = c(FALSE), # Their code for ct.cov=TRUE doesn't work as intended
                            centered = c(TRUE, FALSE),
                            normalize = c(TRUE, FALSE))

#### Iterate through parameters in parallel ####

for (P in 1:nrow(params_loop1)) {
  reference_data_name <- params_loop1$reference_data_name[P]
  test_data_name <- params_loop1$test_data_name[P]
  granularity <- params_loop1$granularity[P]

  data <- Load_AlgorithmInputData(reference_data_name, test_data_name,
                                  granularity,
                                  reference_input_type = "singlecell",
                                  output_type = "counts")

  data$test <- as.matrix(assay(data$test, "counts"))

  A <- Load_AvgLibSize(reference_data_name, granularity)

  ##### Iterate through combinations of MuSiC arguments in parallel #####
  # NOTE: the helper functions have to be sourced inside the foreach loop
  #       so they exist in each newly-created parallel environment

  music_list <- foreach (R = 1:nrow(params_loop2),
                         .packages = required_libraries) %dopar% {
    source(file.path("functions", "General_HelperFunctions.R"))
    source(file.path("functions", "Music_InnerLoop.R"))

    set.seed(12345)
    res <- Music_InnerLoop(data$reference, data$test, A,
                           cbind(params_loop1[P,], params_loop2[R,]))

    Save_AlgorithmIntermediate(res, "music")
    return(res)
  }

  # It's possible for some items in music_list to be null if there was an error.
  # Filter them out.
  music_list <- music_list[lengths(music_list) > 0]

  names(music_list) <- paste0("music_",
                              str_glue("{reference_data_name}_{test_data_name}_{granularity}_"),
                              1:length(music_list))

  print(str_glue("Saving final list for {reference_data_name} {test_data_name} {granularity}..."))
  Save_AlgorithmOutputList(music_list, "music", reference_data_name,
                           test_data_name, granularity)

  rm(music_list, data)
  gc()
}

stopCluster(cl)
