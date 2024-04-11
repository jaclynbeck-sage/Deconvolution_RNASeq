# Merges a list of Seurat objects back together, keeping the mapped dimension
# reductions but removing the integrated counts, which aren't needed and take up
# needed space in memory.
#
# Arguments:
#   query_list = a list of Seurat objects
#
# Returns:
#   A single Seurat object
Merge_Queries <- function(query_list) {
  merged <- merge(query_list[[1]], query_list[2:length(query_list)],
                  merge.data = FALSE,
                  merge.dr = c("ref.pca", "ref.umap"))
  return(merged)
}


# Wrapper around FindTransferAnchors and MapQuery, because this is used twice
# in the script with slightly different arguments and this helps avoid making
# copy/paste errors.
#
# Arguments:
#   query = the Seurat object to be mapped
#   refernce = the Seurat object used as a reference
#   dims_use = the dimensions to use for FindTransferAnchors
#   map_cols = a string or a vector of strings indicating which metadata fields
#              to transfer labels for (i.e. "broad_class" or "sub_class")
#   recompute_residuals = if TRUE, FindTransferAnchors will compute residuals
#              for the query dataset using the reference dataset's SCT model.
#              If FALSE, the query dataset must be SCTransformed and have its
#              pearson residuals already calculated.
#
# Returns:
#   a Seurat object that has been mapped onto the reference, with predicted
#   labels set in the metadata and mapped dimension reductions "ref.pca" and
#   "ref.umap".
Map_Cells <- function(query, reference, dims_use = 1:30, map_cols = NULL,
                      recompute_residuals = TRUE) {
  anchors <- FindTransferAnchors(reference = reference,
                                 query = query,
                                 dims = dims_use,
                                 normalization.method = "SCT",
                                 reference.reduction = "pca",
                                 recompute.residuals = recompute_residuals)

  # Doing a full MapQuery instead of just TransferData isn't strictly necessary
  # but is useful for examining how well the mapping went
  mapped <- MapQuery(anchorset = anchors,
                     reference = reference,
                     query = query,
                     refdata = map_cols,
                     reference.reduction = "pca",
                     reduction.model = "umap")
  return(mapped)
}
