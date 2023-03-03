# Runs MuSiC on a variety of parameters and saves the list of results to a
# file.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(foreach)
library(doParallel)

##### Parallel execution setup #####

cores <- 12
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("MuSiC", "SummarizedExperiment", "stringr", "dplyr")

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

foreach (P = 1:nrow(params_loop1), .packages = required_libraries) %dopar% {
  # These need to be sourced inside the loop for parallel processing
  source(file.path("functions", "General_HelperFunctions.R"))
  source(file.path("functions", "FileIO_HelperFunctions.R"))
  source(file.path("functions", "Music_InnerLoop.R"))

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

  # Each dataset / datatype / granularity combo gets its own list
  music_list <- list()

  ##### Iterate through combinations of MuSiC arguments #####
  # NOTE: This set of parameters (params_loop2) are all executed in the same
  # thread because they use the same single cell and pseudobulk data

  music_list <- foreach (R = 1:nrow(params_loop2)) %do% {
    res <- Music_InnerLoop(sce, pseudobulk, A,
                           cbind(params_loop1[P,], params_loop2[R,]))
    return(res)
  }

  # It's possible for some items in music_list to be null if there was an error.
  # Filter them out.
  music_list <- music_list[!is.null(music_list)]

  names(music_list) <- paste0("music_",
                              str_glue("{dataset}_{granularity}_{datatype}_"),
                              1:nrow(params_loop2))

  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(music_list, "music", dataset, datatype, granularity)

  rm(music_list, pseudobulk, sce)
  gc()
}

stopCluster(cl)
