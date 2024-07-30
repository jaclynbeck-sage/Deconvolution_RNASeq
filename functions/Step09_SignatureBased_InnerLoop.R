# Runs a signature-based deconvolution algorithm on the given bulk data, using a
# signature matrix and a set of parameters. Current algorithms that use this
# method are DeconRNASeq and DWLS.
#
# Arguments:
#   data = a named list that must contain the following items:
#     reference = a signature matrix describing the average expression of each
#                 gene (in cpm) for each cell type (rows = genes, columns = cell types).
#     test = a matrix of bulk expression (rows = genes, columns = samples).
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
SignatureBased_InnerLoop <- function(data, params, algorithm) {
  # Ensure signature is a matrix
  signature <- as.matrix(data$reference)

  # Filter signature matrix according to parameters
  signature_filt <- FilterSignature_FromParams(signature, params)

  if (Check_MissingMarkers(signature_filt, params) ||
      Check_TooFewMarkers(signature_filt, params, 3) ||
      Check_NotEnoughNewMarkers(signature_filt, params)) {
    return(NULL)
  }

  # Signature and bulk data should have the same rows post-filtering.
  keepgene <- intersect(rownames(signature_filt), rownames(data$test))
  bulk_mat_filt <- data$test[keepgene, ]
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

        res_orig <- DeconRNASeq(bulk_mat_filt,
                                signature_filt,
                                proportions = NULL,
                                known.prop = FALSE,
                                use.scale = use_scale,
                                fig = FALSE)

        rownames(res_orig$out.all) <- colnames(bulk_mat_filt)
        res <- list("estimates" = res_orig$out.all)
      } # end DeconRNASeq

      # Run DWLS ---------------------------------------------------------------
      else if (algorithm == "DWLS") {
        # DWLS-specific argument
        solver_type <- params$solver_type

        res_pcts <- omnideconv::deconvolute_dwls(bulk_mat_filt,
                                                 signature_filt,
                                                 dwls_submethod = solver_type,
                                                 verbose = FALSE)

        if (any(is.na(res_pcts)) || any(res_pcts < 0)) {
          message(paste("NA or negative numbers in",
                        paste(params, collapse = "_")))
        }
        res_pcts[res_pcts < 0] <- 0
        res_pcts <- sweep(res_pcts, 1, rowSums(res_pcts), "/")

        res <- list("estimates" = res_pcts)
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
    ## Error handling ----------------------------------------------------------
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
