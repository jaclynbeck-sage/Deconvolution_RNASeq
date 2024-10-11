# Runs MuSiC on the given bulk data, using a reference single cell data set and
# a set of parameters.
#
# Arguments:
#   data = a named list that must contain the following elements:
#     reference = a SingleCellExperiment object containing the reference data.
#                 Must have a colData column named "celltype", denoting the cell
#                 type assignment of each cell.
#     test = a matrix of bulk expression (rows = genes, columns = samples).
#            Must be a matrix, not a data.frame.
#     sc_basis = pre-computed list output from music_basis()
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: algorithm, reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            marker_subtype, marker_input_type, and algorithm-specific
#            variables ct_cov, centered, normalize
#   verbose = whether Music should have verbose output
#
# Returns:
#   a list containing entries for "Est.prop.weighted" and "Est.prop.allgene",
#   which are output by MuSiC, "Est.pctRNA.weighted" and "Est.pctRNA.allgene",
#   which are calculated based on the "M.S" matrix from sc_basis, "params",
#   which is the parameter set used for this run, and "markers", which is the
#   list of genes used as markers for this run
Music_InnerLoop <- function(data, params, verbose = FALSE) {
  # Unpack variables for readability, enforce they are the correct types
  ct_cov <- as.logical(params$ct.cov)
  centered <- as.logical(params$centered)
  normalize <- as.logical(params$normalize)

  celltypes <- levels(data$reference$celltype)
  n_celltypes <- length(celltypes) # Needed for the Check_ functions to work

  # We can use the FilterSignature function to get the list of genes to use
  tmp <- FilterSignature_FromParams(data$test[, 1:n_celltypes], params)

  if (Check_MissingMarkers(tmp, params) ||
      Check_TooFewMarkers(tmp, params, 3) ||
      Check_NotEnoughNewMarkers(tmp, params)) {
    return(NULL)
  }

  markers_use <- rownames(tmp)

  # Modify sc_basis for this set of markers.
  sc_basis_precomputed <- data$sc_basis
  sc_basis_precomputed$Sigma <- sc_basis_precomputed$Sigma[markers_use, ]
  sc_basis_precomputed$Sigma.ct <- sc_basis_precomputed$Sigma.ct[, markers_use]
  sc_basis_precomputed$Disgn.mtx <- sc_basis_precomputed$Disgn.mtx[markers_use, ]
  sc_basis_precomputed$M.theta <- sc_basis_precomputed$M.theta[markers_use, ]

  # Run MuSiC ------------------------------------------------------------------
  # Sometimes MuSiC will produce too many NAs if too many cell types are missing
  # from too many samples, and this eventually throws errors. The best thing to
  # do is just ignore the error and continue the loop
  tryCatch(
    {
      cs <- data.frame(cells = names(data$sc_basis$M.S),
                       sizes = data$sc_basis$M.S)

      result <- music_prop(bulk.mtx = data$test,
                           sc.sce = data$reference,
                           markers = markers_use,
                           clusters = "celltype",
                           samples = "sample",
                           verbose = verbose,
                           cell_size = cs,
                           ct.cov = ct_cov,
                           centered = centered,
                           normalize = normalize,
                           sc.basis = sc_basis_precomputed)

      # Remove "Weight.gene", "r.squared.full", and "Var.prop". "Weight.gene"
      # especially is a very large array and is unneeded, so this reduces
      # output size.
      result <- result[c("Est.prop.weighted", "Est.prop.allgene")]
      result$Est.prop.weighted <- result$Est.prop.weighted[, celltypes]
      result$Est.prop.allgene <- result$Est.prop.allgene[, celltypes]

      # Convert proportion of cells to percent RNA, and assign this as the
      # estimates output. We don't use Est.pctRNA.allgene downstream but want to
      # keep it for reference
      M.S <- data$sc_basis$M.S[celltypes]
      result$estimates <- ConvertPropCellsToPctRNA(
        result$Est.prop.weighted, M.S
      )
      result$Est.pctRNA.allgene <- ConvertPropCellsToPctRNA(
        result$Est.prop.allgene, M.S
      )

      # Remove Est.prop.weighted and Est.prop.allgene, we only need pctRNA
      # from this point forward.
      result$Est.prop.weighted <- NULL
      result$Est.prop.allgene <- NULL

      result$params <- params
      result$markers <- markers_use
      print(paste(result$params, collapse = "  "))

      return(result)
    },
    # Error handling -----------------------------------------------------------
    error = function(err) {
      param_set <- paste(params, collapse = "  ")
      print(paste("*** Error running param set", param_set))
      print(err)
      print("*** skipping ***")

      return(NULL)
    }
  )
}
