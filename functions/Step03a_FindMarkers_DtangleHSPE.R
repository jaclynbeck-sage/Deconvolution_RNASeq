# Dtangle and HSPE calculate their own markers based on the reference data.
# Both algorithms use the same method/code to find markers, so only one run
# per parameter set is needed for use with both algorithms.
#
# This script runs through multiple parameters and data input types and saves
# the resulting markers to files for later usage. This is done outside the
# main algorithm loops because both algorithms can use the same markers file for
# the same parameters, so this removes redundancy. Additionally, we use Dtangle
# markers for some of the cell type signature options for DeconRNASeq, so these
# need to be available before running that algorithm.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(foreach)
library(doParallel)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("hspeSparse", "SingleCellExperiment",
                        "SummarizedExperiment", "stringr", "dplyr")

# Assume this needs up to 20 GB of RAM per core depending on dataset, and adjust
# accordingly.
cores <- 4
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# The set of parameters to iterate over
params_data <- expand.grid(dataset = c("cain", "lau", "lengEC", "lengSFG",
                                       "mathys", "morabito", "seaRef"), #, "seaAD"))
                           granularity = c("broad", "fine"),
                           input_type = c("singlecell", "pseudobulk"),
                           stringsAsFactors = FALSE) %>% arrange(dataset)

# Each param set is completely independent of others, so we run it in parallel.
foreach (R = 1:nrow(params_data), .packages = required_libraries) %dopar% {
  # This needs to be sourced inside the loop for parallel processing
  source(file.path("functions", "FileIO_HelperFunctions.R"))

  dataset <- params_data$dataset[R]
  granularity <- params_data$granularity[R]
  input_type <- params_data$input_type[R]

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

  marker_methods <- c("ratio", "diff")

  # We can use p.value and regression with pseudobulk but not single cell due to
  # memory and/or time constraints
  if (input_type == "pseudobulk") {
    marker_methods <- c("ratio", "diff", "p.value", "regression")
  }

  # Create marker lists and save them.
  for (marker_meth in marker_methods) {
    print(str_glue(paste0("Finding markers for {dataset} {input_type} ",
                          "({granularity}), method = {marker_meth} ...")))

    markers <- find_markers(Y = t(input_mat),
                            pure_samples = pure_samples,
                            marker_method = marker_meth)

    Save_DtangleMarkers(markers, dataset, granularity, input_type, marker_meth)

    rm(markers)
    gc()
  }

  return(NULL)
}

stopCluster(cl)
