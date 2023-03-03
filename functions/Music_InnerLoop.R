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
  ct_cov <- params$ct.cov
  centered <- params$centered
  normalize <- params$normalize

  # Sometimes ct.cov = TRUE will produce too many NAs if too many cell types
  # are missing from too many donors, and this eventually throws errors. The
  # best thing to do is just ignore the error and continue the loop
  tryCatch({
    result <- music_prop(bulk.mtx = bulk_mtx, sc.sce = sce,
                         clusters = "celltype",
                         samples = "donor", verbose = TRUE,
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
