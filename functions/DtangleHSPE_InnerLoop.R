DtangleHSPE_InnerLoop <- function(Y, pure_samples, params, algorithm, limit_n_markers = FALSE) {
  dataset <- params$dataset
  granularity <- params$granularity
  input_type <- params$input_type

  marker_method <- params$marker_method
  n_markers <- params$n_markers

  markers <- Load_DtangleMarkers(dataset, granularity, input_type,
                                 marker_method)

  markers <- lapply(markers, function(X) {intersect(X, colnames(Y))})

  if (any(lengths(markers) == 0)) {
    print(str_glue("*** Skipping marker method {marker_method}: Missing markers for some cell types. ***"))
    return(NULL)
  }

  n_markers_orig <- n_markers
  n_markers <- Get_NMarkers(markers, n_markers_orig)

  # For whole-number n_markers arguments, the n_markers argument (mostly)
  # doubles each time. If there aren't enough markers in each cell type to do
  # anything new with this n_markers value, skip testing.
  if (limit_n_markers & is.vector(n_markers) & all(n_markers <= (n_markers_orig/2))) {
    print(str_glue("*** Skipping n_markers {n_markers_orig}: Not enough total markers. ***"))
    return(NULL)
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
