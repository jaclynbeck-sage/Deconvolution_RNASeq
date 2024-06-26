# Runs a signature-based deconvolution algorithm on the given bulk data, using a
# signature matrix and a set of parameters. Current algorithms that use this
# method are DeconRNASeq and DWLS.
#
# Arguments:
#   signature = a data.frame describing the average expression of each gene (in
#               cpm) for each cell type (rows = genes, columns = cell types).
#               Must be a data.frame, not a matrix.
#   bulk_mat = a matrix of bulk expression (rows = genes, columns = samples).
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            marker_subtype, marker_input_type, and any algorithm-specific
#            variables.
#   algorithm = the name of the algorithm being run (options are "DeconRNASeq"
#               and "DWLS")
#
# Returns:
#   a list containing entries for the celltype percentage estimates (name varies
#   between algorithms), "params", which is the parameter set used for this run,
#   and "markers", which is the list of genes used as markers for this run
SignatureBased_InnerLoop <- function(signature, bulk_mat, params, algorithm) {
  # Ensure signature is a matrix
  signature <- as.matrix(signature)

  # Filter signature matrix according to parameters
  signature_filt <- FilterSignature_FromParams(signature, params, bulk_mat)

  if (Check_MissingMarkers(signature, params) ||
    Check_TooFewMarkers(signature, params, 3) ||
    Check_NotEnoughNewMarkers(signature, params)) {
    return(NULL)
  }

  # Signature and bulk data should have the same rows post-filtering.
  keepgene <- intersect(rownames(signature_filt), rownames(bulk_mat))
  bulk_mat_filt <- bulk_mat[keepgene, ]
  signature_filt <- signature_filt[keepgene, ]

  tryCatch(
    {
      # Run DeconRNASeq --------------------------------------------------------
      if (algorithm == "DeconRNASeq") {
        # DeconRNASeq-specific argument
        use_scale <- as.logical(params$use_scale)

        # These matrices need to be data.frames for this algorithm
        bulk_mat_filt <- as.data.frame(bulk_mat_filt)
        signature_filt <- as.data.frame(signature_filt)

        res <- DeconRNASeq(bulk_mat_filt, signature_filt,
                           proportions = NULL,
                           known.prop = FALSE,
                           use.scale = use_scale,
                           fig = FALSE)

        rownames(res$out.all) <- colnames(bulk_mat_filt)
      } # end DeconRNASeq

      # Run DWLS ---------------------------------------------------------------
      else if (algorithm == "DWLS") {
        # DWLS-specific argument
        solver_type <- params$solver_type

        # There is a lot more setup for this algorithm so for readability, it's
        # moved to a function defined below
        res <- RunDWLS(signature_filt, bulk_mat_filt, solver_type, params)
      } # End DWLS
      else {
        stop(str_glue("Unsupported algorithm name '{algorithm}'."))
      }

      res$params <- params
      res$markers <- rownames(signature_filt)

      # End of loop
      print(paste(res$params, collapse = "  "))
      return(res)
    },
    ##### Error handling #####
    error = function(err) {
      param_set <- paste(params, collapse = "  ")
      msg <- paste("*** Error running param set: ", param_set,
                   "/  *** skipping ***")
      print(msg)
      print(err)

      return(NULL)
    }
  )
}

# Helper function for above
RunDWLS <- function(signature_filt, bulk_mat_filt, solver_type, params) {
  if (solver_type == "DWLS") {
    solve_fn <- solveDampenedWLS
  } else if (solver_type == "SVR") {
    solve_fn <- solveSVR
  } else {
    stop(str_glue("Unsupported solver_type '{solver_type}'"))
  }

  # Calls solveDampenedWLS or solveSVR. These functions can only be called on
  # one sample at a time so this loops over all samples. DWLS prints the output
  # of every call so we redirect this to a file to remove it from the main
  # status file from threading
  sink(file.path(dir_estimates_tmp,
                 paste0("DWLS_tmp_", paste(params, collapse = "_"), ".txt")))
  res_pcts <- apply(bulk_mat_filt, 2, function(B) {
    return(solve_fn(signature_filt, B))
  })
  sink()

  res_pcts <- t(res_pcts) # puts cell types as columns, samples as rows

  # DWLS can produce negative numbers but these are usually really close to
  # zero. Include a check for larger negative numbers just in case. DWLS can
  # also produce NAs, and these should get thrown out as well.
  if (any(is.na(res_pcts)) | any(res_pcts < -1e-3)) {
    param_set <- paste(params, collapse = "  ")
    msg <- paste("*** Negative numbers or NA in result for for param set:",
                 param_set, "/  *** skipping ***")
    return(NULL)
  }

  # Fix any numbers really close to 0 that are negative
  res_pcts[res_pcts < 0] <- 0
  res_pcts <- sweep(res_pcts, 1, rowSums(res_pcts), "/")

  res <- list("estimates" = res_pcts)

  return(res)
}
