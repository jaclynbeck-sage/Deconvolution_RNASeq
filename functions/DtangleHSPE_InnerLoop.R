DtangleHSPE_InnerLoop <- function(Y, pure_samples, params, algorithm,
                                  limit_n_markers = FALSE) {
  reference_data_name <- params$reference_data_name
  granularity <- params$granularity
  reference_input_type <- params$reference_input_type
  marker_input_type <- params$marker_input_type

  marker_type <- params$marker_type
  marker_subtype <- params$marker_subtype
  n_markers <- as.numeric( params$n_markers )
  marker_order <- params$marker_order

  # "p.value" and "regression" aren't feasible for single cell input
  if (marker_input_type == "singlecell" & marker_type == "dtangle") {
    if (marker_subtype == "p.value" | marker_subtype == "regression") {
      print(str_glue("*** Skipping combination singlecell / dtangle / {marker_subtype}"))
      return(NULL)
    }
  }

  markers <- Load_Markers(reference_data_name, granularity, marker_type,
                          marker_subtype, marker_input_type)

  markers <- lapply(markers, function(X) {intersect(X, colnames(Y))})

  if (any(lengths(markers) == 0)) {
    print(str_glue(paste0("*** Skipping marker method {marker_type} / ",
                          "{marker_subtype}: Missing markers for some cell ",
                          "types. ***")))
    return(NULL)
  }

  if (marker_order == "correlation") {
    markers <- OrderMarkers_ByCorrelation(markers, t(Y[-unlist(pure_samples),]))
  }

  n_markers_orig <- n_markers
  n_markers <- Get_NMarkers(markers, n_markers_orig)

  if (limit_n_markers) {
    err_string <- str_glue("*** Skipping {marker_type} / {marker_subtype} / {n_markers_orig}: ")
    ok <- Check_NMarkers(n_markers_orig, n_markers, err_string)

    if (ok == FALSE) {
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
    result <- result %>% select(-diag)
  }

  # Add the params we used to generate this run
  result$params <- params
  result$markers <- unlist(result$markers)
  print(paste(result$params, collapse = "  "))

  return(result)
}
