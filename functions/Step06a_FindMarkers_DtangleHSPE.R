# Dtangle and HSPE calculate their own markers based on the reference data.
# Both algorithms use the same method/code to find markers, so only one run
# per parameter set is needed for use with both algorithms.
#
# This function runs through multiple parameters and data input types and saves
# the resulting markers to files for later usage with all of the deconvolution
# algorithms.

FindMarkers_DtangleHSPE <- function(datasets, granularities, cl = NULL) {
  params_data <- expand.grid(dataset = datasets,
                             granularity = granularities,
                             stringsAsFactors = FALSE) |>
    dplyr::arrange(dataset)

  # Each param set is completely independent of others, so we run it in
  # parallel. If no parallel cluster is set up before running this function, the
  # loop will run one by one instead of in parallel.
  parallel::parApply(cl, params_data, 1, function(params) {
    # This needs to be sourced inside the loop for parallel processing
    source(file.path("functions", "FileIO_HelperFunctions.R"))

    dataset <- params["dataset"]
    granularity <- params["granularity"]

    pb <- Load_PseudobulkPureSamples(dataset, granularity,
                                     output_type = "log_cpm")
    metadata <- colData(pb)
    input_mat <- assay(pb, "counts")

    celltypes <- levels(metadata$celltype)
    pure_samples <- lapply(celltypes, function(ct) {
      which(metadata$celltype == ct)
    })

    names(pure_samples) <- celltypes
    marker_methods <- c("ratio", "diff", "p.value", "regression")

    # Create marker lists and save them.
    for (method in marker_methods) {
      message(str_glue("Finding markers for {dataset} ({granularity}), ",
                       "method = {method} ..."))

      markers <- dtangle::find_markers(Y = t(input_mat),
                                       pure_samples = pure_samples,
                                       marker_method = method)

      # Filter the marker list to genes that meet certain criteria:
      #   ratio: all Vs higher than the mean (there isn't an intuitive criteria for this one)
      #   diff: all Vs with >= 1 log-fold change
      #   p.value: all Vs >= 0.95 (signifying p <= 0.05)
      #   regression: all Vs with >= 1 log-fold change

      # LFC needs a smaller threshold for sub class
      lfc_thresh <- ifelse(granularity == "broad_class", 1, 0.5)

      markers[["filtered"]] <- lapply(markers$V, function(vals) {
        new_vals <- switch(method,
                           "ratio" = vals[vals >= mean(vals)],
                           "diff" = vals[vals >= lfc_thresh],
                           "p.value" = vals[vals >= 0.95],
                           "regression" = vals[vals >= lfc_thresh])
        return(names(new_vals))
      })

      # Format like the other marker finding algorithms
      markers$L <- lapply(markers$L, names)
      markers_final <- list("all" = markers$L,
                            "filtered" = markers$filtered,
                            "V" = markers$V)

      Save_Markers(markers_final, dataset, granularity,
                   marker_type = "dtangle",
                   marker_subtype = method)

      rm(markers)
      gc()
    }

    return(NULL)
  })
}
