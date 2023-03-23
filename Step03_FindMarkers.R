library(dplyr)
library(foreach)
library(doParallel)

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


if ("seurat" %in% marker_types_run) {
  source(file.path("functions", "Step03b_FindMarkers_Seurat.R"))
  FindMarkers_Seurat(datasets, granularities)
}

if ("autogenes" %in% marker_types_run) {

}
