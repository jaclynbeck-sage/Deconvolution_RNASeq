# Runs DeconRNASeq on the given bulk data, using a signature matrix and a set
# of parameters. This is used for both training / ground truth data and testing
# data, which each have unique setups but use this same piece of code for
# execution.
#
# Arguments:
#   signature = a data.frame describing the average expression of each gene (in
#               cpm) for each cell type (rows = genes, columns = cell types).
#               Must be a data.frame, not a matrix.
#   bulk_df = a data.frame of bulk expression (rows = genes, columns = samples).
#             Must be a data.frame, not a matrix.
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            marker_subtype, marker_input_type, and use_scale.
#
# Returns:
#   a list containing entries for "out.all" and "out.pca", which are output by
#   DeconRNASeq, and "params", which is the parameter set used for this run
DeconRNASeq_InnerLoop <- function(signature, bulk_df, params) {

  # Unpack variables for readability, enforce they are the correct types
  reference_data_name <- params$reference_data_name
  granularity <- params$granularity

  filter_level <- as.numeric( params$filter_level )
  n_markers <- as.numeric( params$n_markers )
  marker_type <- params$marker_type
  marker_subtype <- params$marker_subtype
  marker_input_type <- params$marker_input_type
  marker_order <- params$marker_order
  use_scale <- as.logical( params$use_scale )
  recalc_cpm <- as.logical( params$recalc_cpm )

  if (recalc_cpm) {
    signature <- as.data.frame(scuttle::calculateCPM(signature))
    bulk_df <- as.data.frame(scuttle::calculateCPM(bulk_df))
  }

  # Filter signature matrix according to parameters
  signature_filt <- FilterSignature(signature, filter_level, reference_data_name,
                                    granularity, n_markers, marker_type,
                                    marker_subtype, marker_input_type, marker_order,
                                    bulk_df)

  if (is.null(signature_filt)) {
    param_set <- paste(params, collapse = "  ")
    print(c("*** Missing markers for at least one cell type for param set", param_set))
    print("*** skipping ***")
    return(NULL)
  }

  # We want at least ~3 markers per cell type or there isn't enough information
  # to work from
  if (nrow(signature_filt) < 3*ncol(signature_filt)) {
    param_set <- paste(params, collapse = "  ")
    print(c("*** Too few markers for param set", param_set))
    print("*** skipping ***")
    return(NULL)
  }

  # Signature and bulk data should have the same rows post-filtering
  keepgene <- intersect(rownames(signature_filt), rownames(bulk_df))
  bulk_df_filt <- bulk_df[keepgene,]
  signature_filt <- signature_filt[keepgene,]

  # Run DeconRNASeq
  tryCatch({
    res <- DeconRNASeq(bulk_df_filt, signature_filt, proportions = NULL,
                       known.prop = FALSE, use.scale = use_scale, fig = FALSE)

    rownames(res$out.all) <- colnames(bulk_df_filt)
    res$params <- params
    res$markers <- rownames(signature_filt)

    print(paste(res$params, collapse = "  "))
    return(res)
  },
  error = function(err) {
    param_set <- paste(params, collapse = "  ")
    print(c("*** Error running param set", param_set))
    print(err)
    print("*** skipping ***")

    return(NULL)
  })
}
