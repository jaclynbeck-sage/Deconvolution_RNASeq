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
source("Music_Edits.R")

##### Parallel execution setup #####

# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- 6
# Music is very verbose so we make an output file instead of putting on the
# console, because RStudio won't display everything in the console otherwise.
cl <- makeCluster(cores, type = "FORK", outfile = "music_output.txt")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("MuSiC", "SingleCellExperiment")

#### Parameter setup #####

datasets <- c("cain", "lau", "leng", "mathys", "morabito", "seaRef") #, "seaAD")

reference_input_types = c("singlecell", "pseudobulk")

params_loop1 <- expand_grid(reference_data_name = datasets,
                            test_data_name = c("Mayo", "MSBB", "ROSMAP"),
                            granularity = c("broad_class"),
                            normalization = c("counts", "cpm", "tmm", "tpm")) %>%
                  arrange(test_data_name)

marker_types <- list("dtangle" = c("ratio", "diff", "p.value", "regression"),
                     "autogenes" = c("correlation", "distance", "combined"),
                     "seurat" = c("None"),
                     "deseq2" = c("DESeq2"))

params_tmp <- CreateParams_FilterableSignature(
                  filter_level = c(0, 1, 2, 3),
                  n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                                3, 5, 10, 20, 50, 100, 200, 500),
                  marker_types,
                  marker_input_types = c("singlecell", "pseudobulk"),
                  marker_order = c("distance", "correlation")
)

params_loop2 <- expand_grid(params_tmp,
                            ct.cov = c(FALSE), # Their ct.cov methods don't work as intended and/or take forever to converge
                            centered = c(TRUE, FALSE),
                            normalize = c(TRUE, FALSE))


#### Iterate through parameters in parallel ####

for (P in 1:nrow(params_loop1)) {
  reference_data_name <- params_loop1$reference_data_name[P]
  test_data_name <- params_loop1$test_data_name[P]
  granularity <- params_loop1$granularity[P]
  normalization <- params_loop1$normalization[P]

  music_list <- list()

  for (input_type in reference_input_types) {

    data <- Load_AlgorithmInputData(reference_data_name, test_data_name,
                                    granularity,
                                    reference_input_type = input_type,
                                    output_type = normalization)

    if (input_type == "pseudobulk") {
      data$reference$sample <- str_replace(data$reference$sample, ".*_", "")
    }

    data$reference <- as(data$reference, "SingleCellExperiment")
    data$test <- as.matrix(assay(data$test, "counts"))

    # Pre-compute sc.basis to save time
    # Hack to get the function redirect to work: The "<<-" operator creates a
    # global variable that can be used inside the modified music_basis function.
    sc_basis_precomputed <<- NULL

    sc_basis <- music_basis(data$reference, non.zero = TRUE,
                            markers = rownames(data$reference),
                            clusters = "celltype", samples = "sample",
                            select.ct = NULL,
                            ct.cov = FALSE,
                            verbose = TRUE)

    # Add the covariance variable (equivalent to re-calling music_basis with ct.cov = TRUE)
    Sigma.ct <- music_Sigma.ct(data$reference,
                               non.zero = TRUE,
                               markers = NULL,
                               clusters = "celltype",
                               samples = "sample",
                               select.ct = NULL)
    colnames(Sigma.ct) <- rownames(data$reference)

    sc_basis$Sigma.ct <- Sigma.ct

    # Save in case of crashing
    saveRDS(sc_basis, file = file.path(dir_tmp,
                                       str_glue("sc_basis_{reference_data_name}_{granularity}_{input_type}_{normalization}.rds")))

    ##### Iterate through combinations of MuSiC arguments in parallel #####
    # NOTE: the helper functions have to be sourced inside the foreach loop
    #       so they exist in each newly-created parallel environment

    music_list_tmp <- foreach (R = 1:nrow(params_loop2),
                               .packages = required_libraries) %dopar% {
      source(file.path("functions", "General_HelperFunctions.R"))
      source(file.path("functions", "Music_InnerLoop.R"))

      set.seed(12345)
      res <- Music_InnerLoop(data$reference, data$test,
                             sc_basis,
                             cbind(params_loop1[P,],
                                   "reference_input_type" = input_type,
                                   params_loop2[R,]),
                             verbose = FALSE)

      Save_AlgorithmIntermediate(res, "music")
      return(res)
    }

    # It's possible for some items in music_list to be null if there was an error.
    # Filter them out.
    music_list_tmp <- music_list_tmp[lengths(music_list_tmp) > 0]

    name_base <- str_glue("{reference_data_name}_{test_data_name}_{granularity}_{input_type}_{normalization}")
    names(music_list_tmp) <- paste("music",
                                   name_base,
                                   1:length(music_list_tmp), sep = "_")

    music_list <- append(music_list, music_list_tmp)

    rm(data)
    gc()
  }

  print(str_glue("Saving final list for {reference_data_name} {test_data_name} {granularity} {normalization}..."))
  Save_AlgorithmOutputList(music_list, "music", reference_data_name,
                           test_data_name, granularity, normalization)

  rm(music_list)
  gc()
}

stopCluster(cl)
