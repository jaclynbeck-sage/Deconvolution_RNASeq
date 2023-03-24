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
#            containing the following variables: dataset, granularity,
#            filter_level, n_markers, marker_type, marker_subtype, and use_scale.
#
# Returns:
#   a list containing entries for "out.all" and "out.pca", which are output by
#   DeconRNASeq, and "params", which is the parameter set used for this run
DeconRNASeq_InnerLoop <- function(signature, bulk_df, params) {

  # Unpack variables for readability, enforce they are the correct types
  dataset <- params$dataset
  granularity <- params$granularity

  filter_level <- as.numeric( params$filter_level )
  n_markers <- as.numeric( params$n_markers )
  marker_type <- params$marker_type
  marker_subtype <- params$marker_subtype
  use_scale <- as.logical( params$use_scale )

  # Filter signature matrix according to parameters
  signature_filt <- FilterSignature(signature, filter_level, dataset,
                                    granularity, n_markers, marker_type,
                                    marker_subtype)

  # Signature and bulk data should have the same rows post-filtering
  keepgene <- intersect(rownames(signature_filt), rownames(bulk_df))
  bulk_df_filt <- bulk_df[keepgene,]
  signature_filt <- signature_filt[keepgene,]

  # Run DeconRNASeq
  res <- DeconRNASeq(bulk_df_filt, signature_filt, proportions = NULL,
                     known.prop = FALSE, use.scale = use_scale, fig = FALSE)

  rownames(res$out.all) <- colnames(bulk_df_filt)
  res$params <- params

  print(paste(res$params, collapse = "  "))
  return(res)
}
