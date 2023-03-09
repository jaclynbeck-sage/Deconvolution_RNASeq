# Runs MuSiC on a variety of parameters and saves the list of results to a
# file.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(foreach)
library(doParallel)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

##### Parallel execution setup #####

# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- 8
cl <- makeCluster(cores, type = "FORK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("MuSiC", "SingleCellExperiment")

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

params_loop1 <- expand.grid(dataset = datasets,
                            datatype = c("donors", "training"),
                            granularity = c("broad", "fine"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

params_loop2 <- expand.grid(ct.cov = c(TRUE, FALSE),
                            centered = c(TRUE, FALSE),
                            normalize = c(TRUE, FALSE))

#### Iterate through parameters in parallel ####

for (P in 1:nrow(params_loop1)) {
  dataset <- params_loop1$dataset[P]
  datatype <- params_loop1$datatype[P]
  granularity <- params_loop1$granularity[P]

  sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

  pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity,
                                output_type = "counts")
  pseudobulk <- assay(pseudobulk, "counts")

  A <- Load_AvgLibSize(dataset, granularity)

  # These SHOULD have the same rownames, but just in case.
  keepgene <- intersect(rownames(sce), rownames(pseudobulk))
  pseudobulk <- as.matrix(pseudobulk[keepgene, ])
  sce <- sce[keepgene, ]

  ##### Iterate through combinations of MuSiC arguments in parallel #####
  # NOTE: the helper functions have to be sourced inside the foreach loop
  #       so they exist in each newly-created parallel environment

  music_list <- foreach (R = 1:nrow(params_loop2),
                         .packages = required_libraries) %dopar% {
    source(file.path("functions", "General_HelperFunctions.R"))
    source(file.path("functions", "Music_InnerLoop.R"))
    set.seed(12345)

    res <- Music_InnerLoop(sce, pseudobulk, A,
                           cbind(params_loop1[P,], params_loop2[R,]))

    Save_AlgorithmIntermediate(res, "music")
    return(res)
  }

  # It's possible for some items in music_list to be null if there was an error.
  # Filter them out.
  music_list <- music_list[lengths(music_list) > 0]

  names(music_list) <- paste0("music_",
                              str_glue("{dataset}_{datatype}_{granularity}_"),
                              1:length(music_list))

  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(music_list, "music", dataset, datatype, granularity)

  rm(music_list, pseudobulk, sce)
  gc()
}

stopCluster(cl)
