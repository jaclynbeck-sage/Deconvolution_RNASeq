# Helper functions that are used in multiple scripts.
# NOTE: For the calculations below that use model.matrix and matrix
# multiplication, using this method to sum counts is way faster than iterating
# over cell type/donor combos and using sum(), colSums(), or rowSums(). We do
# it this way even though it's less readable, because for data sets this large
# the decrease in run-time is significant.

library(Matrix)
library(stringr)
library(scuttle)

source(file.path("functions", "FileIO_HelperFunctions.R"))

# CalculatePercentRNA: Calculates the percentage of RNA contributed by each cell
# type to each sample as:
#   (sum of RNA counts over all cells of type T in sample S) /
#     (sum of RNA counts over all cells of all types in sample S)
#
# Raw, un-normalized RNA counts are added together rather than adding normalized
# CPMs, in order to be consistent with what we would get from a bulk sample of
# purified cells.
#
# Arguments:
#   singlecell_counts = matrix of counts (rows = genes, cols = cells)
#   samples = vector of sample IDs that is the same length as
#             ncol(singlecell_counts). Must either be a factor (contains
#             multiple sample IDs) or a vector with the same sample ID repeated
#             (for all cells belonging to a single sample).
#   celltypes = vector of cell type assignments that is the same length as
#               ncol(singlecell_counts). Must be a factor.
#
# Returns:
#   matrix of percentages (rows = samples, cols = cell types), whose rows sum
#   to 1.
CalculatePercentRNA <- function(singlecell_counts, samples, celltypes) {
  # Sum all counts per gene for each donor/celltype combination
  y <- model.matrix(~0 + samples:celltypes)
  count_sums <- singlecell_counts %*% y

  # Sum over all genes for each donor/celltype combination
  count_sums <- colSums(count_sums)

  # names are of the format "samples<sample>:celltypes<celltype>", split
  # them apart by the ":" into a data frame w/ 2 columns
  col_info <- str_split(names(count_sums), ":", simplify = TRUE)
  sample_list <- unique(col_info[,1])

  pctRNA <- sapply(sample_list, function (samp) {
    # All entries for this donor = 1 entry per cell type
    cols <- which(col_info[,1] == samp)
    pct <- count_sums[cols] / sum(count_sums[cols])

    # Remove extra labels added by model.matrix
    names(pct) <- str_replace(names(pct), "samples.*celltypes", "")
    return(pct)
  })

  colnames(pctRNA) <- str_replace(colnames(pctRNA), "samples", "")
  return(t(pctRNA))
}


# CalculateA: Calculates the A-matrix (average library size of each cell type)
# to provide a way to convert from proportion of cells to percent RNA in
# samples where ground truth is not known. The calculation is as follows:
#   For each sample S:
#     A_s_T = average total RNA count of all cells of type T in sample S
#     A_s = (vector of A_s_T for all cell types) normalized to sum to 1
#   A_T = average normalized A_s_T over all samples
#   A = (vector of A_T for all cell types), which sums to 1
#
# Because A is an estimate of average RNA count per cell type, raw RNA counts
# are averaged rather than using normalized CPMs. Normalizing A *after*
# averaging counts removes scaling issues with test samples and results in a
# vector describing each cell type's average total count relative to other types.
#
# Arguments:
#   sce = SingleCellExperiment containing raw RNA counts
#   samples = vector of sample IDs that is the same length as
#             ncol(sce). Must be a factor.
#   celltypes = vector of cell type assignments that is the same length as
#               ncol(sce). Must be a factor.
#
# Returns:
#   vector of average library size of each cell type, normalized to sum to 1.
CalculateA <- function(sce, samples, celltypes) {
  # Sum all counts per gene for each donor/celltype combination
  y <- model.matrix(~0 + samples:celltypes)
  count_sums <- counts(sce) %*% y

  # Sum over all genes for each donor/celltype combination
  count_sums <- colSums(count_sums)

  # For each sample
  A_s <- sapply(levels(samples), function(samp) {
    # For each cell of type <ct> in the sample
    sapply(levels(celltypes), function(ct) {
      name <- str_glue("samples{samp}:celltypes{ct}")
      cells <- celltypes == ct & samples == samp

      # Average library size = sum over all genes for this donor/celltype
      # divided by number of cells for this donor/celltype
      return(as.numeric(count_sums[name] / sum(cells)))
    })
  })

  # Normalize each donor row to sum to 1
  A_s <- t(sweep(A_s, 2, colSums(A_s, na.rm = TRUE), "/"))

  # Take the means but exclude entries for donors who don't have a certain cell type
  A_s[A_s == 0] <- NA
  A <- colMeans(A_s, na.rm = TRUE)
  A <- A / sum(A) # Enforce sum to 1
}


# CalculateSignature: Creates a cell-type-specific "signature" matrix that
# describes the expected count of each gene for each cell type:
#   Sig_s_T = (sum of RNA counts per gene, for all cells of type T in sample S),
#             converted to CPM after summation
#   Sig_T = average Sig_s_T over all samples
#
# Raw, un-normalized RNA counts are added together rather than adding normalized
# CPMs, in order to be consistent with what we would get from a bulk sample of
# purified cells. The sum is then converted to CPM to remove scaling issues
# between samples.
#
# Arguments:
#   sce = SingleCellExperiment containing raw RNA counts
#   samples = vector of sample IDs that is the same length as
#             ncol(sce). Must be a factor.
#   celltypes = vector of cell type assignments that is the same length as
#               ncol(sce). Must be a factor.
#
# Returns:
#   matrix of average CPMs for each gene, for each cell type (rows = genes,
#   columns = cell types)
CalculateSignature <- function(sce, samples, celltypes) {
  # Sum all counts per gene for each donor/celltype combination
  y <- model.matrix(~0 + samples:celltypes)
  count_sums <- counts(sce) %*% y

  # Remove columns where a donor doesn't have the specified cell type
  count_sums <- count_sums[,colSums(count_sums) > 0]

  # Convert to CPM -- creates an "average" profile in CPM for each cell type
  count_sums <- calculateCPM(count_sums)

  # colnames are of the format samples<sample>:celltypes<celltype>, extract
  # "<celltype>" as a vector
  cts <- str_replace(colnames(count_sums), ".*:celltypes", "")

  # Get the mean over all donors, for each gene and cell type
  sig <- sapply(levels(celltypes), function(ct) {
    cols <- which(cts == ct)
    rowMeans(count_sums[,cols], na.rm = TRUE)
  })

  return(sig)
}


# Filter Signature: filters the signature matrix, which includes all genes
# present in the single cell dataset, down to a set of genes which are more
# informative. Because we don't know exactly what will be the most informative
# for each algorithm, this function has multiple filtering settings:
#   Filter level 0 = don't filter at all
#                1 = keep genes where at least one cell type has > 1 cpm
#                2 = keep genes where at least one cell type has > 10 cpm
#                3 = use Dtangle-determined markers, and filter the markers
#                    further to use the top <filt_percent> genes of each cell
#                    type (Dtangle lists markers in descending order of variance)
#
# Arguments:
#   signature = signature matrix of cpm values for each gene for each cell type
#               (rows = genes, columns = cell types)
#   filter_level = the filter setting to use (range 0-3)
#   dataset = the data set used to get Dtangle markers. Only applicable if
#             filter_level = 3.
#   granularity = the level of cell types used ("broad" or "fine") to get
#                 Dtangle markers. Only applicable if filter_level = 3.
#   filt_percent = percent of Dtangle markers to keep (range 0-1). Only
#                  applicable if filter_level = 3.
#
# Returns:
#   signature matrix containing only the genes that pass the specified filters.
#   rows = genes, columns = cell types.
FilterSignature <- function(signature, filter_level = 1, dataset = NULL, granularity = NULL, filt_percent = 1.0) {
  # Filter for genes where at least one cell type expresses at > 1 cpm
  if (filter_level == 1) {
    ok <- which(rowSums(signature >= 1) > 0)
    signature <- signature[ok, ]
  }

  # Filter for genes where at least one cell type expresses at > 10 cpm
  else if (filter_level == 2) {
    ok <- which(rowSums(signature >= 10) > 0)
    signature <- signature[ok, ]
  }

  else if (filter_level == 3 & !is.null(dataset) & !is.null(granularity)) {
    # TODO is this the best one?
    markers <- Load_DtangleMarkers(dataset, granularity, "pseudobulk", "diff")

    n_markers <- ceiling(lengths(markers$L) * filt_percent)
    markers <- lapply(names(markers$L), function(N) {
      names(markers$L[[N]][1:n_markers[N]])
    })

    signature <- signature[unique(unlist(markers)),]
  }

  return(signature)
}


# ConvertPropCellsToPctRNA: for algorithms that output proportion of cells
# as (number of cells of type T) / (number of total cells), convert these
# values to percent of RNA contributed by each cell type:
#   percent RNA = (proportion of cells * average library size of the cell type),
#                 normalized so all percentages sum to 1
#
# Arguments:
#   propCells: a matrix of calculated cell-type proportions (rows = samples,
#              columns = cell types)
#   A: the A-matrix (average library size per cell type) calculated from the
#      reference data set
#
# Returns:
#   a matrix of percent RNA contributed by each cell type (rows = samples,
#   columns = cell types)
ConvertPropCellsToPctRNA <- function(propCells, A) {
  pct <- sweep(propCells, 2, A, "*")
  pct <- sweep(pct, 1, rowSums(pct), "/")
  return(pct)
}
