# Dtangle and HSPE calculate their own markers based on the reference data.
# Both algorithms use the same method/code to find markers, so only one run
# per parameter set is needed for use with both algorithms.
#
# This function runs through multiple parameters and data input types and saves
# the resulting markers to files for later usage with all of the deconvolution
# algorithms.

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("hspe", "SingleCellExperiment", "Matrix",
                        "SummarizedExperiment", "stringr", "dplyr")

FindMarkers_DtangleHSPE <- function(datasets, granularities, input_types) {
  params_data <- expand.grid(dataset = datasets,
                             granularity = granularities,
                             input_type = input_types,
                             stringsAsFactors = FALSE) %>% dplyr::arrange(dataset)

  # Each param set is completely independent of others, so we run it in parallel.
  # If no parallel cluster is set up before running this function, the foreach
  # loop will run one by one instead of in parallel.
  foreach(R = 1:nrow(params_data), .packages = required_libraries) %dopar% {
    # This needs to be sourced inside the loop for parallel processing
    source(file.path("functions", "FileIO_HelperFunctions.R"))

    dataset <- params_data$dataset[R]
    granularity <- params_data$granularity[R]
    input_type <- params_data$input_type[R]

    if (input_type == "singlecell") {
      input_obj <- Load_SingleCell(dataset, granularity, output_type = "log_cpm")
    } else { # pseudobulk pure samples
      input_obj <- Load_PseudobulkPureSamples(dataset, granularity,
                                              output_type = "log_cpm")
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
    for (method in marker_methods) {
      message(str_glue(paste0("Finding markers for {dataset} {input_type} ",
                              "({granularity}), method = {method} ...")))

      markers <- hspe::find_markers(Y = t(input_mat),
                                    pure_samples = pure_samples,
                                    marker_method = method)

      # Filter the marker list to genes that meet certain criteria:
      #   ratio: all Vs higher than the mean (there isn't an intuitive criteria for this one)
      #   diff: all Vs with >= 1 log-fold change
      #   p.value: all Vs >= 0.95 (signifying p <= 0.05)
      #   regression: all Vs with >= 1 log-fold change
      markers[["filtered"]] <- lapply(markers$V, function(vals) {
        new_vals <- switch(method,
                           "ratio" = vals[vals >= mean(vals)],
                           "diff" = vals[vals >= 1],
                           "p.value" = vals[vals >= 0.95],
                           "regression" = vals[vals >= 1])
        return(names(new_vals))
      })

      # Format like the other marker finding algorithms
      markers$L <- lapply(markers$L, names)
      markers_final <- list("all" = markers$L,
                            "filtered" = markers$filtered)

      Save_Markers(markers_final, dataset, granularity,
                   marker_type = "dtangle",
                   input_type = input_type,
                   marker_subtype = method)

      rm(markers)
      gc()
    }

    return(NULL)
  }
}
