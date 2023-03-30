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
library(foreach)
library(doParallel)
library(SummarizedExperiment)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

##### Parallel execution setup #####

# NOTE: Assume about 5-20 GB RAM needed per core. DeconRNASeq multi-threads,
#       so only use ~half the cores available on the machine.
# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- 2
cl <- makeCluster(cores, type = "FORK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("DeconRNASeq")

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

params_loop1 <- expand.grid(dataset = datasets,
                            datatype = c("donors", "training"),
                            granularity = c("broad", "fine"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

params_loop2 <- expand.grid(filter_level = c(0, 1, 2, 3),
                            n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                                          10, 20, 50, 100, 200, 500),
                            marker_type = c("dtangle", "autogenes", "seurat"),
                            use_scale = c(TRUE, FALSE),
                            stringsAsFactors = FALSE)

marker_types <- list("dtangle" = c("ratio", "diff", "p.value", "regression"),
                     "autogenes" = c("correlation", "distance", "combined"),
                     "seurat" = c("None"))
marker_types <- do.call(rbind, lapply(names(marker_types), function(N) {
  data.frame(marker_type = N, marker_subtype = marker_types[[N]])
}))

params_loop2 <- merge(params_loop2, marker_types, by = "marker_type")

# Some filter_type / n_markers combos are not valid, get rid of them
# (filter levels 1 & 2 don't use n_markers or marker_type arguments)
low_filt <- params_loop2$filter_level < 3
params_loop2$marker_type[low_filt] <- "None"
params_loop2$marker_subtype[low_filt] <- "None"
params_loop2$n_markers[low_filt] <- -1

params_loop2 <- params_loop2 %>% distinct() %>% arrange(filter_level)


#### Iterate through parameters ####

for (P in 1:nrow(params_loop1)) {
  dataset <- params_loop1$dataset[P]
  datatype <- params_loop1$datatype[P]
  granularity <- params_loop1$granularity[P]

  # Input data
  signature <- Load_SignatureMatrix(dataset, granularity)
  signature <- as.data.frame(signature)

  # Test data
  pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, output_type = "cpm")
  pseudobulk <- assay(pseudobulk, "counts")
  pseudobulk <- as.data.frame(as.matrix(pseudobulk))

  ##### Filter levels, number of markers, and DeconRNASeq arguments #####
  # NOTE: the helper functions have to be sourced inside the foreach loop
  #       so they exist in each newly-created parallel environment

  decon_list <- foreach (R = 1:nrow(params_loop2),
                         .packages = required_libraries) %dopar% {
    source(file.path("functions", "DeconRNASeq_InnerLoop.R"))
    source(file.path("functions", "General_HelperFunctions.R"))

    set.seed(12345)
    res <- DeconRNASeq_InnerLoop(signature, pseudobulk,
                                 cbind(params_loop1[P,], params_loop2[R,]))

    Save_AlgorithmIntermediate(res, "deconRNASeq")
    return(res)
  }

  # It's possible for some items in decon_list to be null if there was an error.
  # Filter them out.
  decon_list <- decon_list[lengths(decon_list) > 0]

  names(decon_list) <- paste0("deconRNASeq_",
                              str_glue("{dataset}_{datatype}_{granularity}_"),
                              1:length(decon_list))

  # Save the completed list
  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(decon_list, "deconRNASeq", dataset, datatype, granularity)

  rm(decon_list)
  gc()
}

stopCluster(cl)
