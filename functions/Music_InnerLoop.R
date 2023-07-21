# Runs MuSiC on the given bulk data, using a reference single cell data set and
# a set of parameters. This is used for both training / ground truth data and
# testing data, which each have unique setups but use this same piece of code
# for execution.
#
# Arguments:
#   sce = a SingleCellExperiment object containing the reference data. Must have
#         a colData column named "celltype", denoting the cell type assignment
#         of each cell.
#   bulk_df = a matrix of bulk expression (rows = genes, columns = samples).
#             Must be a matrix, not a data.frame.
#   A = a named vector of the average expected library size of each cell type,
#       normalized so sum(A) = 1.
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: ct_cov, centered, normalize
#
# Returns:
#   a list containing entries for "Est.prop.weighted" and "Est.prop.allgene",
#   which are output by MuSiC, "Est.pctRNA.weighted" and "Est.pctRNA.allgene",
#   which are calculated based on the "A" matrix, and "params", which is the
#   parameter set used for this run
Music_InnerLoop <- function(sce, bulk_mtx, A, params) {
  # Unpack variables for readability, enforce they are the correct types
  reference_data_name <- params$reference_data_name
  granularity <- params$granularity

  filter_level <- as.numeric( params$filter_level )
  n_markers <- as.numeric( params$n_markers )
  marker_type <- params$marker_type
  marker_subtype <- params$marker_subtype
  marker_input_type <- params$marker_input_type
  marker_order <- params$marker_order

  ct_cov <- as.logical( params$ct.cov )
  centered <- as.logical( params$centered )
  normalize <- as.logical( params$normalize )

  # We can use the FilterSignature function to get the list of genes to use
  tmp <- FilterSignature(bulk_mtx[,1:2], filter_level, reference_data_name,
                         granularity, n_markers, marker_type, marker_subtype,
                         marker_input_type, marker_order, bulk_mtx)

  if (is.null(tmp)) {
    param_set <- paste(params, collapse = "  ")
    print(c("*** Missing markers for at least one cell type for param set", param_set))
    print("*** skipping ***")
    return(NULL)
  }

  markers_use <- rownames(tmp)

  # We want at least ~3 markers per cell type or there isn't enough information
  # to work from
  if (length(markers_use) < 3*length(A)) {
    param_set <- paste(params, collapse = "  ")
    print(c("*** Too few markers for param set", param_set))
    print("*** skipping ***")
    return(NULL)
  }

  # Sometimes ct.cov = TRUE will produce too many NAs if too many cell types
  # are missing from too many donors, and this eventually throws errors. The
  # best thing to do is just ignore the error and continue the loop
  tryCatch({
    result <- music_prop(bulk.mtx = bulk_mtx, sc.sce = sce,
                         markers = markers_use,
                         clusters = "celltype",
                         samples = "donor", verbose = FALSE,
                         ct.cov = ct_cov, centered = centered,
                         normalize = normalize)

    # Remove "Weight.gene", "r.squared.full", and "Var.prop". "Weight.gene"
    # especially is a very large array and is unneeded, so this reduces
    # output size.
    result <- result[c("Est.prop.weighted", "Est.prop.allgene")]
    result$Est.prop.weighted <- result$Est.prop.weighted[,names(A)]
    result$Est.prop.allgene <- result$Est.prop.allgene[,names(A)]

    # Convert proportion of cells to percent RNA
    result$Est.pctRNA.weighted <- ConvertPropCellsToPctRNA(result$Est.prop.weighted, A)
    result$Est.pctRNA.allgene <- ConvertPropCellsToPctRNA(result$Est.prop.allgene, A)

    result$params <- params
    result$markers <- markers_use
    print(paste(result$params, collapse = "  "))

    return(result)
  },
  error = function(err) {
    param_set <- paste(params, collapse = "  ")
    print(c("*** Error running param set", param_set))
    print(err)
    print("*** skipping ***")

    return(NULL)
  })
}
