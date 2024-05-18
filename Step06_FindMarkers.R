# Finds cell type markers in single cell data using 4 different algorithms:
# Dtangle, Seurat/MAST, AutogeneS, and DESeq2

library(dplyr)
library(foreach)
library(doParallel)
library(stringr)

# Settings ---------------------------------------------------------------------

# Which algorithms to run for marker finding
marker_types_run <- c("dtangle", "seurat", "autogenes", "deseq2", "excluded_genes")

# Run multiple parameter sets in parallel for dtangle/HSPE marker finding if
# dtangle_do_parallel is TRUE. Other algorithms can't be run in parallel.
# Assume most data sets use <20 GB of RAM per core, but seaRef will use > 100 GB
# for single cell.
dtangle_do_parallel <- TRUE
dtangle_n_cores <- 4
dtangle_clust_type <- "FORK" # Use PSOCK for non-Unix systems

# Which datasets to run on
datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

# What granularities?
granularities <- c("broad_class", "sub_class")

# Dtangle/HSPE only: input types
dtangle_input_types <- c("singlecell", "pseudobulk")


# Dtangle/HSPE markers ---------------------------------------------------------

if ("dtangle" %in% marker_types_run) {
  if (dtangle_do_parallel) {
    cl <- makeCluster(dtangle_n_cores, type = dtangle_clust_type, outfile = "")
    registerDoParallel(cl)
  }

  source(file.path("functions", "Step06a_FindMarkers_DtangleHSPE.R"))
  FindMarkers_DtangleHSPE(datasets, granularities, dtangle_input_types)

  if (dtangle_do_parallel) {
    stopCluster(cl)
  }
}


# Seurat / MAST markers --------------------------------------------------------

# Uses all available CPUs already, so this isn't run in parallel.

if ("seurat" %in% marker_types_run) {
  source(file.path("functions", "Step06b_FindMarkers_Seurat.R"))
  FindMarkers_Seurat(datasets, granularities)
}


# AutogeneS markers (requires python/reticulate) -------------------------------

# This section uses reticulate, so it's not run in parallel in case that causes
# issues with python environments.

if ("autogenes" %in% marker_types_run) {
  source(file.path("functions", "Step06c_FindMarkers_AutogeneS.R"))
  FindMarkers_AutogeneS(datasets, granularities)
}


# DESeq2 markers from pseudobulk -----------------------------------------------

if ("deseq2" %in% marker_types_run) {
  source(file.path("functions", "Step06d_FindMarkers_DESeq2.R"))
  FindMarkers_DESeq2(datasets, granularities)
}


# Find marker genes that change with diagnosis ---------------------------------

if ("excluded_genes" %in% marker_types_run) {
  source(file.path("functions", "Step06e_FindMarkerExclusions.R"))
  FindMarkerExclusions(datasets, granularities)
}
