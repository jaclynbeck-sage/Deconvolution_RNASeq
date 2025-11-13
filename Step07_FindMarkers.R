# Finds cell type markers in single cell data using 4 different algorithms:
# Dtangle, Seurat, AutogeneS, and DESeq2

library(dplyr)
library(stringr)
library(parallel)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

# Settings ---------------------------------------------------------------------

# Which algorithms to run for marker finding
marker_types_run <- c("dtangle", "seurat", "autogenes", "deseq2", "excluded_genes")

# Run multiple parameter sets in parallel where possible. Dtangle and Deseq2
# marker finding can be run in parallel.
do_parallel <- TRUE
clust_type <- "FORK" # Use PSOCK for non-Unix systems
n_cores <- detectCores() - 1

# Which datasets to run on
datasets <- all_singlecell_datasets()

# What granularities?
granularities <- c("broad_class", "sub_class")


# Dtangle/HSPE markers ---------------------------------------------------------

if ("dtangle" %in% marker_types_run) {
  if (do_parallel) {
    cl <- makeCluster(n_cores, type = clust_type, outfile = "")
  } else {
    cl <- NULL
  }

  source(file.path("functions", "Step07a_FindMarkers_DtangleHSPE.R"))
  FindMarkers_DtangleHSPE(datasets, granularities, cl)

  if (do_parallel) {
    stopCluster(cl)
  }
}


# Seurat markers ---------------------------------------------------------------

# This is pretty fast already and needs a lot of memory, so this isn't run in parallel.

if ("seurat" %in% marker_types_run) {
  source(file.path("functions", "Step07b_FindMarkers_Seurat.R"))
  FindMarkers_Seurat(datasets, granularities)
}


# AutogeneS markers (requires python/reticulate) -------------------------------

# This section uses reticulate, so it's not run in parallel in case that causes
# issues with python environments.

if ("autogenes" %in% marker_types_run) {
  source(file.path("functions", "Step07c_FindMarkers_AutogeneS.R"))
  FindMarkers_AutogeneS(datasets, granularities)
}


# DESeq2 markers from pseudobulk -----------------------------------------------

if ("deseq2" %in% marker_types_run) {
  source(file.path("functions", "Step07d_FindMarkers_DESeq2.R"))
  if (!do_parallel) {
    n_cores <- 1
  }
  FindMarkers_DESeq2(datasets, granularities, n_cores)
}


# Find marker genes that change with diagnosis ---------------------------------

if ("excluded_genes" %in% marker_types_run) {
  source(file.path("functions", "Step07e_FindMarkerExclusions.R"))
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
    markers$ad_gene_exclusions <- exclusions

    saveRDS(markers, file)

    loss1 <- (1 - length(setdiff(unlist(markers$all), markers$ad_gene_exclusions)) /
                sum(lengths(markers$all))) * 100
    loss2 <- (1 - length(setdiff(unlist(markers$filtered), markers$ad_gene_exclusions)) /
                sum(lengths(markers$filtered))) * 100
    print(str_glue("{basename(file)}: removed {round(loss1)}% of ",
                   "non-logfc-filtered and {round(loss2)}% of logfc-filtered ",
                   "markers."))
  }
}


# Calculate ordering by correlation --------------------------------------------

# For every combination of reference data set, bulk data set, normalization,
# regression, and granularity

params_loop <- expand_grid(test_data_name = all_bulk_datasets(),
                           normalization = c("cpm", "tmm", "tpm"),
                           regression_method = c("none", "edger", "lme", "combat"))

for (dataset in datasets) {
  marker_files <- list.files(path = dir_markers,
                             pattern = dataset,
                             full.names = TRUE)

  if (do_parallel) {
    cl <- makeCluster(n_cores, type = clust_type, outfile = "")
  } else {
    cl <- NULL
  }

  parallel::parLapply(cl, marker_files, function(file) {
    markers <- readRDS(file)
    markers$ordered_by_correlation <- list(all = list(),
                                           filtered = list())

    for (P in 1:nrow(params_loop)) {
      params <- params_loop[P, ]
      data <- Load_BulkData(params$test_data_name, params$normalization,
                            params$regression_method)
      data <- as.matrix(assay(data, "counts"))
      marker_name <- paste(params, collapse = "_")
      cat(marker_name, "\n")

      markers_all_ordered <- OrderMarkers_ByCorrelation(markers$all, data)
      markers_filt_ordered <- OrderMarkers_ByCorrelation(markers$filtered, data)

      markers$ordered_by_correlation$all[[marker_name]] <- markers_all_ordered
      markers$ordered_by_correlation$filtered[[marker_name]] <- markers_filt_ordered
    }

    print(file)
    saveRDS(markers, file)
  })

  if (do_parallel) {
    stopCluster(cl)
  }
}
