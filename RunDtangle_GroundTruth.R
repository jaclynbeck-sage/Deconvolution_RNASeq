# Runs Dtangle or HSPE on a variety of parameters and saves the list of results
# to a file. The code to run these two algorithms is nearly identical so it has
# been combined into one file.
#
# NOTE: This script relies on Dtangle markers, so these must be calculated
# before-hand.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "DtangleHSPE_HelperFunctions.R"))

##### Edit this variable to run either dtangle or hspe #####

algorithm <- "dtangle" # "dtangle" or "hspe", all lower case

##### Parallel execution setup #####

# NOTE: Dtangle and HSPE are very memory-intensive. The amount of memory used
#       per core will vary by data set, usually about 2-3x the size of the data
#       set in memory. For most datasets, it's between 5-20 GB per core for
#       single cell input, much less for pseudobulk input. Adjust accordingly.
# NOTE: HSPE multi-threads and will use all available CPU, so only use ~2 cores
#       regardless of memory constraints.
# NOTE: "FORK" is more memory-efficient but only works on Unix systems. For
#       other systems, use "PSOCK" and reduce the number of cores.
cores <- 8
cl <- makeCluster(cores, type = "FORK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c(str_glue("{algorithm}Sparse"))

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

input_types = c("singlecell", "pseudobulk")

params_loop1 <- expand.grid(dataset = datasets,
                            granularity = c("broad", "fine"),
                            datatype = c("donors", "training"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

if (algorithm == "dtangle") {
  params_loop2 <- expand.grid(marker_type = c("autogenes", "dtangle", "seurat"),
                              marker_input_type = c("singlecell", "pseudobulk"),
                              gamma_name = c("auto"),
                              sum_fn_type = c("mean", "median"),
                              n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2,
                                            0.5, 0.75, 1.0,
                                            5, 10, 20, 50, 100, 200, 500),
                              stringsAsFactors = FALSE)
  marker_types <- list("dtangle" = c("ratio", "diff", "p.value", "regression"),
                       "autogenes" = c("correlation", "distance", "combined"),
                       "seurat" = c("None"))
  marker_types <- do.call(rbind, lapply(names(marker_types), function(N) {
    data.frame(marker_type = N, marker_subtype = marker_types[[N]])
  }))

  params_loop2 <- merge(params_loop2, marker_types, by = "marker_type")

} else if (algorithm == "hspe") {
  # Needs updating for marker_subtype
  params_loop2 <- expand.grid(marker_method = c("ratio", "diff", "p.value", "regression"),
                              loss_fn = c("var", "L2"),
                              n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2,
                                            0.5, 0.75, 1.0,
                                            5, 10, 20, 50, 100, 200, 500),
                              stringsAsFactors = FALSE)
}

#### Iterate through parameters ####

for (P in 1:nrow(params_loop1)) {
  dataset <- params_loop1$dataset[P]
  granularity <- params_loop1$granularity[P]
  datatype <- params_loop1$datatype[P]

  # Each dataset / datatype / granularity combo gets its own list
  dtangle_list <- list()

  ##### Run on both single cell and pseudobulk input matrices #####

  for (input_type in input_types) {

    ##### Prepare reference and test matrices #####

    input_list <- Get_DtangleHSPEInput(dataset, datatype, granularity, input_type)

    # Free up unused memory
    gc()

    ##### Iterate through Dtangle/HSPE parameters in parallel #####
    # NOTE: the helper functions have to be sourced inside the foreach loop
    #       so they exist in each newly-created parallel environment

    dtangle_list_tmp <- foreach (R = 1:nrow(params_loop2),
                                 .packages = required_libraries) %dopar% {
      source(file.path("functions", "DtangleHSPE_HelperFunctions.R"))
      source(file.path("functions", "DtangleHSPE_InnerLoop.R"))

      set.seed(12345)
      res <- DtangleHSPE_InnerLoop(Y = input_list[["Y"]],
                                   pure_samples = input_list[["pure_samples"]],
                                   params = cbind(params_loop1[P,],
                                                  "input_type" = input_type,
                                                  params_loop2[R,]),
                                   algorithm = algorithm,
                                   limit_n_markers = TRUE)

      Save_AlgorithmIntermediate(res, algorithm)
      return(res)
    }

    # It's possible for some items in dtangle_list to be null if a parameter
    # set was skipped. Filter them out.
    dtangle_list_tmp <- dtangle_list_tmp[lengths(dtangle_list_tmp) > 0]

    name <- str_glue("{algorithm}_{dataset}_{datatype}_{granularity}_{input_type}_")
    names(dtangle_list_tmp) <- paste0(name, 1:length(dtangle_list_tmp))

    dtangle_list <- append(dtangle_list, dtangle_list_tmp)

    # Next iteration will start with new data, remove the old data
    rm(input_list)
    gc()
  } # end input_types loop

  # Save the completed list
  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(dtangle_list, algorithm, dataset, datatype, granularity)
} # end params_loop1 parallel loop

stopCluster(cl)
