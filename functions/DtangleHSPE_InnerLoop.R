DtangleHSPE_InnerLoop <- function(Y, pure_samples, params, algorithm, limit_n_markers = FALSE) {
  dataset <- params$dataset
  granularity <- params$granularity
  input_type <- params$input_type
  marker_input_type <- params$marker_input_type

  marker_type <- params$marker_type
  marker_subtype <- params$marker_subtype
  n_markers <- params$n_markers

  # "p.value" and "regression" aren't feasible for single cell input
  if (input_type == "singlecell" & marker_type == "dtangle") {
    if (marker_subtype == "p.value" | marker_subtype == "regression") {
      print(str_glue("*** Skipping combination singlecell / dtangle / {marker_subtype}"))
      return(NULL)
    }
  }

  markers <- Load_Markers(dataset, granularity, marker_type, marker_subtype,
                          marker_input_type)

  markers <- lapply(markers, function(X) {intersect(X, colnames(Y))})

  if (any(lengths(markers) == 0)) {
    print(str_glue(paste0("*** Skipping marker method {marker_type} / ",
                          "{marker_subtype}: Missing markers for some cell ",
                          "types. ***")))
    return(NULL)
  }

  n_markers_orig <- n_markers
  n_markers <- Get_NMarkers(markers, n_markers_orig)

  if (limit_n_markers) {
    # For whole-number n_markers arguments, the n_markers argument (mostly)
    # doubles each time. If there aren't enough markers in each cell type to do
    # anything new with this n_markers value, skip testing.
    if (is.vector(n_markers) & all(n_markers <= (n_markers_orig/2))) {
      print(str_glue("*** Skipping n_markers {n_markers_orig}: Not enough total markers. ***"))
      return(NULL)
    }
    # Early testing showed it's not useful to use a large number of markers, so
    # if the total markers used is > 5000, don't test.
    if ( (is.vector(n_markers) & (sum(n_markers) > 5000)) |
         (!is.vector(n_markers) & (sum(lengths(markers) * n_markers) > 5000)) ) {
      print(str_glue("*** Skipping n_markers {n_markers_orig}: Marker set too large. ***"))
      return(NULL)
    }
  }

  ##### Dtangle-specific function call #####
  if (algorithm == "dtangle") {
    gamma <- Get_Gamma(params$gamma_name)
    sum_fn <- Get_SumFn(params$sum_fn_type)

    result <- dtangle(Y = Y,
                      pure_samples = pure_samples,
                      data_type = "rna-seq",
                      gamma = gamma, # If gamma is not NULL, it will override data_type argument
                      n_markers = n_markers,
                      markers = markers,
                      summary_fn = sum_fn)

    # Only keep results for bulk test samples
    test_samples <- setdiff(1:nrow(Y), unlist(pure_samples))
    result$estimates <- result$estimates[test_samples, ]
  }
  ##### HSPE-specific function call #####
  else if (algorithm == "hspe") {
    loss_fn <- params$loss_fn

    result <- hspe(Y = Y,
                   pure_samples = pure_samples,
                   n_markers = n_markers,
                   markers = markers,
                   loss_fn = loss_fn,
                   seed = 12345)

    # Get rid of "diag" (index 5), which is huge and unneeded
    result <- result[1:4]
  }

  # Add the params we used to generate this run
  result$params <- params
  print(paste(result$params, collapse = "  "))

  return(result)
}
