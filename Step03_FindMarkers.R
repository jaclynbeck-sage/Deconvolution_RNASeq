library(dplyr)
library(foreach)
library(doParallel)
library(reticulate)
library(stringr)
source("Filenames.R")

##### Settings #####

# Which algorithms to run for marker finding
marker_types_run <- c("dtangle", "seurat", "autogenes")

# Run multiple parameter sets in parallel for dtangle/HSPE marker finding?
# Assume most data sets use <20 GB of RAM per core, but seaRef will use ~50.
dtangle_do_parallel <- TRUE
dtangle_n_cores <- 4
dtangle_clust_type <- "FORK" # Use PSOCK for non-Unix systems

# Which datasets to run on
datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD"))

# What granularities?
granularities <- c("broad", "fine")


# Dtangle/HSPE only: input types
dtangle_input_types <- c("singlecell", "pseudobulk")


##### Dtangle/HSPE markers #####

if ("dtangle" %in% marker_types_run) {

  if (dtangle_do_parallel) {
    cl <- makeCluster(dtangle_n_cores, type = dtangle_clust_type, outfile = "")
    registerDoParallel(cl)
  }

  source(file.path("functions", "Step03a_FindMarkers_DtangleHSPE.R"))
  FindMarkers_DtangleHSPE(datasets, granularities, dtangle_input_types)

  if (dtangle_do_parallel) {
    stopCluster(cl)
  }
}


##### Seurat / MAST markers #####

# Needs all available CPUs, so this isn't run in parallel.

if ("seurat" %in% marker_types_run) {
  source(file.path("functions", "Step03b_FindMarkers_Seurat.R"))
  FindMarkers_Seurat(datasets, granularities)
}


##### AutogeneS markers (requires python/reticulate) #####

# This section is a little janky to make R and Python work together without
# crashing. We need to use reticulate to establish the conda environment but
# then need to call the python script as a system command, because otherwise
# reticulate/RStudio crashes when trying to call ag.optimize directly from R.

if ("autogenes" %in% marker_types_run) {
  # This activates the conda environment we need, and running the useless command
  # somehow forces things to load properly for when we call the system command
  use_condaenv("autogenes_env")
  py_run_string("import os")

  for (dataset in datasets) {
    for (granularity in granularities) {
      print(str_glue("Finding markers for {dataset} / {granularity}..."))
      output_file_prefix <- file.path(dir_markers,
                                      str_glue("autogenes_markers_{dataset}_{granularity}"))

      python_file <- file.path("functions", "Step03c_FindMarkers_AutogeneS.py")

      cmd <- str_glue("python {python_file} {dataset} {granularity} {output_file_prefix}")
      system(cmd)
    }
  }
}
