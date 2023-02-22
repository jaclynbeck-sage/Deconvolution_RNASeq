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

##### Edit this variable to run either dtangle or hspe #####

algorithm <- "dtangle" # "dtangle" or "hspe", all lower case

##### Parallel execution setup #####

cores <- 12
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("Matrix", "SummarizedExperiment", "stringr", "dplyr",
                        "SingleCellExperiment",
                        str_glue("{algorithm}Sparse"))

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito")#,
              #"seaRef") #, "seaAD")

input_types = c("singlecell", "pseudobulk")

params_loop1 <- expand.grid(dataset = datasets,
                            granularity = c("broad"), #, "fine"),
                            datatype = c("donors", "training"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

if (algorithm == "dtangle") {
  params_loop2 <- expand.grid(marker_method = c("ratio", "diff", "p.value", "regression"),
                              gamma_name = c("auto", 0.5),
                              sum_fn_type = c("mean", "median"),
                              n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2,
                                            0.5, 0.75, 1.0),
                              stringsAsFactors = FALSE)
} else if (algorithm == "hspe") {
  params_loop2 <- expand.grid(marker_method = c("ratio", "diff", "p.value", "regression"),
                              loss_fn = c("var", "L2"),
                              n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2,
                                            0.5, 0.75, 1.0),
                              stringsAsFactors = FALSE)
}

#### Iterate through parameters in parallel ####

foreach (P = 1:nrow(params_loop1), .packages = required_libraries) %dopar% {
  # This needs to be sourced inside the loop for parallel processing
  source(file.path("functions", "FileIO_HelperFunctions.R"))
  source(file.path("functions", "DtangleHSPE_HelperFunctions.R"))

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

    ##### Prepare reference and test matrices #####

    input_list <- Get_DtangleHSPEInput(dataset, datatype, granularity, input_type)
    pure_samples <- input_list[["pure_samples"]]

    # Free up unused memory
    gc()

    # Filter params_loop2 down: Early testing showed we don't need to test
    # higher percentages for method = "ratio", since this method returns the
    # full set of genes and smaller lists do better.
    params_run <- subset(params_loop2, !(marker_method == "ratio" & n_markers > 0.2))

    # "p.value" and "regression" aren't feasible for single cell input
    if (input_type == "singlecell") {
      params_run <- subset(params_run, marker_method == "ratio" | marker_method == "diff")
    }

    ##### Iterate through Dtangle/HSPE parameters #####

    for (R in 1:nrow(params_run)) {
      marker_method <- params_run$marker_method[R]
      n_markers <- params_run$n_markers[R]

      markers <- Load_DtangleMarkers(dataset, granularity, input_type,
                                     marker_method)

      # dtangle/hspe don't interpret "1" as 100%, so we need to input a list of
      # the length of each marker set instead
      if (n_markers == 1) {
        n_markers <- lengths(markers$L)
      }

      name <- str_glue(paste0("{dataset}_{granularity}_{datatype}_",
                              "{input_type}_{R}"))

      ##### Dtangle-specific function call #####
      if (algorithm == "dtangle") {
        gamma <- Get_Gamma(params_run$gamma_name[R])
        sum_fn <- Get_SumFn(params_run$sum_fn_type[R])

        result <- dtangle(Y = input_list[["Y"]],
                          pure_samples = pure_samples,
                          data_type = "rna-seq",
                          gamma = gamma, # If gamma is not NULL, it will override data_type argument
                          n_markers = n_markers,
                          markers = markers,
                          summary_fn = sum_fn)

        # Only keep results for pseudobulk test samples
        test_samples <- setdiff(1:nrow(input_list[["Y"]]), unlist(pure_samples))
        result$estimates <- result$estimates[test_samples, ]
      }
      ##### HSPE-specific function call #####
      else if (algorithm == "hspe") {
        loss_fn <- params_run$loss_fn[R]

        result <- hspe(Y = input_list[["Y"]],
                       pure_samples = pure_samples,
                       n_markers = n_markers,
                       markers = markers,
                       loss_fn = loss_fn,
                       seed = 12345)

        # Get rid of "diag" (index 5), which is huge and unneeded
        result <- result[1:4]
      }

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
  } # end input_types loop

  # Save the completed list
  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(dtangle_list, algorithm, dataset, datatype, granularity)

  return(NULL)
} # end params_loop1 parallel loop

stopCluster(cl)
