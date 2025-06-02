# Finds cell type markers in single cell data using 4 different algorithms:
# Dtangle, Seurat/MAST, AutogeneS, and DESeq2

library(dplyr)
library(foreach)
library(doParallel)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

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


# Filter markers down to exclude diagnosis-influenced genes --------------------

for (dataset in datasets) {
  exclusions_file <- file.path(dir_metadata,
                               str_glue("{dataset}_excluded_genes.rds"))
  if (file.exists(exclusions_file)) {
    exclusions <- readRDS(exclusions_file)
    exclusions <- exclusions$genes
  } else { # seaRef won't have exclusions based on diagnosis
    exclusions <- c()
  }

  marker_files <- list.files(path = dir_markers,
                             pattern = dataset,
                             full.names = TRUE)

  for (file in marker_files) {
    markers <- readRDS(file)
    markers$filtered_for_dx <- lapply(markers$filtered, setdiff, exclusions)
    saveRDS(markers, file)

    loss <- (1 - sum(lengths(markers$filtered_for_dx)) / sum(lengths(markers$filtered))) * 100
    print(str_glue("{basename(file)}: removed {round(loss)}% of total markers"))
  }
}


# Calculate ordering by correlation --------------------------------------------

# For every combination of reference data set, bulk data set, normalization,
# regression, and granularity

params_loop <- expand_grid(test_data_name = c("Mayo_TCX", "Mayo_CBE", "ROSMAP_ACC", "ROSMAP_DLPFC", "ROSMAP_PCC"), # TODO MSBB
                           normalization = c("cpm", "tmm", "tpm"),
                           regression_method = c("none", "edger", "lme", "dream"))

for (granularity in granularities) {
  for (dataset in datasets) {
    marker_files <- list.files(path = dir_markers,
                               pattern = str_glue("{dataset}_{granularity}"),
                               full.names = TRUE)

    for (file in marker_files) {
      markers <- readRDS(file)
      markers$ordered_by_correlation <- list()

      for (P in 1:nrow(params_loop)) {
        params <- params_loop[P, ]
        data <- Load_BulkData(params$test_data_name, params$normalization,
                              params$regression_method)
        markers_ordered <- OrderMarkers_ByCorrelation(markers$filtered,
                                                      as.matrix(assay(data, "counts")))

        markers$ordered_by_correlation[[paste(params, collapse = "_")]] <- markers_ordered
      }

      print(file)
      saveRDS(markers, file)
    }
  }
}
