# Helper functions that are used in multiple scripts.
# NOTE: For the calculations below that use model.matrix and matrix
# multiplication, using this method to sum counts is way faster than iterating
# over cell type/donor combos and using sum(), colSums(), or rowSums(). We do
# it this way even though it's less readable, because for data sets this large
# the decrease in run-time is significant.

library(Matrix)
library(stringr)
library(scuttle)
library(dplyr)
library(tidyr)
library(reshape2)

source(file.path("functions", "FileIO_HelperFunctions.R"))

# Load_AlgorithmInputData: Ease-of-use function that gets both the reference
# data and the test data from a set of different formats, and makes sure
# both data sets have the same genes in the same order.
#
# Arguments:
#   reference_data_name = the name of the reference singlecell data set
#   test_data_name = the name of the test data set to deconvolve
#   granularity = either "broad" or "fine", for the level of cell types used
#   reference_input_type = the type of reference data to load. Options are:
#                           "singlecell": the full single cell data set as a
#                                         SingleCellExperiment object
#                           "pseudobulk": pseudobulked "pure" samples as a
#                                         SummarizedExperiment object
#                           "signature": the signature matrix
#   output_type = if reference_input_type is "singlecell" or "pseudobulk",
#                 specifies how the counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm) of non-zero cpm entries
#                   "log1p_cpm" will take the log2(cpm+1) of cpms
#
# Returns:
#   a list with entries "reference" and "test", containing the reference data
#   and the test data, respectively. Because different algorithms need the
#   input data in different formats, the objects are returned in the format they
#   originated in (a SingleCellExperiment, SummarizedExperiment, or matrix)
#   and it is assumed that the individual algorithms will manipulate the formats
#   afterward.
Load_AlgorithmInputData <- function(reference_data_name, test_data_name,
                                    granularity = "broad",
                                    reference_input_type = "singlecell",
                                    output_type = "counts") {
  # Reference input
  if (reference_input_type == "singlecell") {
    reference_obj <- Load_SingleCell(reference_data_name, granularity,
                                     output_type)
  }
  else if (reference_input_type == "pseudobulk") {
    reference_obj <- Load_PseudobulkPureSamples(reference_data_name,
                                                granularity, output_type)
  }
  else if (reference_input_type == "signature") {
    reference_obj <- Load_SignatureMatrix(reference_data_name, granularity)
  }
  else {
    print("*** Error: Invalid reference_input_type specified! ***")
    return(NULL)
  }

  # Test data
  if (test_data_name == "donors" | test_data_name == "training") {
    test_obj <- Load_Pseudobulk(reference_data_name, test_data_name,
                                granularity, output_type)
  }
  # ROSMAP, Mayo, or MSBB
  else {
    test_obj <- Load_BulkData(test_data_name, output_type)
  }

  genes <- intersect(rownames(reference_obj), rownames(test_obj))
  reference_obj <- reference_obj[genes,]
  test_obj <- test_obj[genes,]

  return(list("reference" = reference_obj, "test" = test_obj))
}


# CreateParams_MarkerTypes - Creates a parameter matrix using tidyr::expand_grid
# (which can take data frames in the input) that has all variables required for
# loading different combinations of markers: n_markers, marker_type,
# marker_subtype, and marker_input_type. Invalid combinations of these variables
# are removed from the parameter set before returning.
#
# Arguments:
#   n_markers = a vector containing one or more percentages (range 0-1.0) and/or
#               one or more integers (range 2-Inf) specifying how many markers
#               per cell type to use.
#   marker_types = a list where the names of the entries are one of "autogenes",
#                  "dtangle", or "seurat", and the items in each entry are a
#                  list of marker subtypes to use for that algorithm. See
#                  FilterSignature for more detail. Must be a list and not a
#                  vector.
#   marker_input_types = dtangle-specific: a vector of one or all of
#                        c("singlecell", "pseudobulk") designating whether to
#                        test dtangle marker sets from singlecell input,
#                        pseudobulk input, or both
#   marker_order = "distance" or "correlation", whether markers should be
#                  ordered by largest expression difference between cell types
#                  or by correlation with other markers for the same cell type
#
# Returns:
#   a tibble containing all possible valid combinations of the arguments
CreateParams_MarkerTypes <- function(n_markers, marker_types,
                                     marker_input_types, marker_order) {
  marker_types <- melt(marker_types) %>% dplyr::rename(marker_type = "L1",
                                                       marker_subtype = "value")

  params <- tidyr::expand_grid(n_markers = n_markers,
                               marker_types,
                               marker_input_type = marker_input_types,
                               marker_order = marker_order)

  # marker_input_type only applies to dtangle markers
  params$marker_input_type[params$marker_type != "dtangle"] <- "None"

  # We don't use dtangle 'p.value' and 'regression' markers for `singlecell`
  # input because the compute time for these is too high
  params <- subset(params, !(marker_input_type == "singlecell" &
                               marker_subtype %in% c("regression", "p.value")))

  # We don't need to re-order markers when we're using the whole marker set
  params$marker_order[params$n_markers == 1] <- "None"

  params <- params %>% distinct()
  return(params)
}


# CreateParams_FilterableSignature - Creates a parameter matrix using
# tidyr::expand_grid (which can take data frames in the input) that has all
# variables required for filtering a signature matrix: filter_level, n_markers,
# marker_type, marker_subtype, and marker_input_type. Invalid combinations of
# these variables are removed from the parameter set before returning.
#
# Arguments:
#   filter_levels = a vector containing one or more filter levels, as used in
#                   the function FilterSignature. Options are c(0, 1, 2, 3).
#   The rest of these arguments are only relevant for filter_level = 3:
#   n_markers = a vector containing one or more percentages (range 0-1.0) and/or
#               one or more integers (range 2-Inf) specifying how many markers
#               per cell type to use.
#   marker_types = a list where the names of the entries are one of "autogenes",
#                  "dtangle", or "seurat", and the items in each entry are a
#                  list of marker subtypes to use for that algorithm. See
#                  FilterSignature for more detail. Must be a list and not a
#                  vector.
#   marker_input_types = dtangle-specific: a vector of one or all of
#                        c("singlecell", "pseudobulk") designating whether to
#                        test dtangle marker sets from singlecell input,
#                        pseudobulk input, or both
#   marker_order = "distance" or "correlation", whether markers should be
#                  ordered by largest expression difference between cell types
#                  or by correlation with other markers for the same cell type
#
# Returns:
#   a tibble containing all possible valid combinations of the arguments
CreateParams_FilterableSignature <- function(filter_levels, n_markers,
                                             marker_types, marker_input_types,
                                             marker_order) {
  params_tmp <- CreateParams_MarkerTypes(n_markers, marker_types,
                                         marker_input_types, marker_order)

  params <- tidyr::expand_grid(filter_level = filter_levels,
                               params_tmp)

  # Some filter_type / n_markers combos are not valid, get rid of them
  # (filter levels 1 & 2 don't use n_markers or marker_type arguments)
  low_filt <- params$filter_level < 3
  params$marker_type[low_filt] <- "None"
  params$marker_input_type[low_filt] <- "None"
  params$marker_subtype[low_filt] <- "None"
  params$n_markers[low_filt] <- -1
  params$marker_order[low_filt] <- "None"

  params <- params %>% distinct()
  return(params)
}


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
#                3 = use algorithm-determined markers, and filter the markers
#                    further to use the top <filt_percent> genes of each cell
#                    type (markers are in descending order of variance or logFC)
#
# Arguments:
#   signature = signature matrix of cpm values for each gene for each cell type
#               (rows = genes, columns = cell types)
#   filter_level = the filter setting to use (range 0-3)
#   The rest of the arguments are only used if filter_level == 3:
#   reference_data_name = the name of the data set used to get markers.
#   granularity = the level of cell types used ("broad" or "fine") to get
#                 markers.
#   filt_percent = percent of markers to keep (range 0-1) OR an integer
#                  describing the number of markers to use per cell type
#                  (range 2-Inf). For example, passing in 10 will use 10
#                  markers for each cell type, while passing in 0.5 will use
#                  50% of each cell type's markers (different number of markers
#                  per cell type).
#   marker_type = one of "autogenes", "dtangle", or "seurat", indicating which
#                 algorithm's marker set to use.
#   marker_subtype = the subtype of markers to load, specific to the marker_type.
#                    For dtangle, one of "ratio", "diff", "p.value", or
#                    "regression", which was used as the input to
#                    dtangle::find_markers.
#                    For autogenes, one of "correlation", "distance", or
#                    "combined", specifying which weighting scheme was used to
#                    pick markers.
#                    Seurat doesn't have a subtype so this can be left as NULL.
#   marker_input_type = for marker_type == "dtangle" only, either "singlecell"
#                       or "pseudobulk", for which type of input was used to
#                       generate the markers. For other marker_types, this
#                       argument is ignored.
#   marker_order = "distance" or "correlation". Whether the markers should be
#                  ordered by distance (prioritize markers with the highest
#                  expression difference between the cell type and all other
#                  cell types) or correlation (prioritize markers with the
#                  highest correlation to other markers for that cell type)
#   test_data = a data.frame or matrix of expression data. Only used if
#               marker_order is "correlation".
#
# Returns:
#   signature matrix containing only the genes that pass the specified filters.
#   rows = genes, columns = cell types.
FilterSignature <- function(signature, filter_level = 1, reference_data_name = NULL,
                            granularity = NULL, filt_percent = 1.0,
                            marker_type = "dtangle", marker_subtype = "diff",
                            marker_input_type = "pseudobulk",
                            marker_order = "distance", test_data = NULL) {
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

  # Filter for genes specific to cell-type marker sets
  else if (filter_level == 3 & !is.null(reference_data_name) & !is.null(granularity)) {
    markers <- Load_Markers(reference_data_name, granularity, marker_type,
                            marker_subtype, input_type = marker_input_type)

    if (is.null(markers)) {
      return(NULL)
    }

    # Remove genes that don't exist in the signature matrix
    markers <- lapply(markers, function(X) {intersect(X, rownames(signature))})

    if (marker_order == "correlation") {
      if (is.null(test_data)) {
        print(paste0("Error! No data provided for calculating marker ",
                     "correlation. Markers will be ordered by distance ",
                     "instead of correlation."))
      }
      else {
        markers <- OrderMarkers_ByCorrelation(markers, test_data)
      }
    }

    # Percentage of each cell type's markers
    if (filt_percent <= 1) {
      n_markers <- ceiling(lengths(markers) * filt_percent)
    }
    # Fixed number of markers for each cell type
    else {
      n_markers <- sapply(lengths(markers), min, filt_percent)
    }

    if (any(n_markers == 0)) {
      return(NULL)
    }

    markers <- lapply(names(markers), function(N) {
      markers[[N]][1:n_markers[N]]
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


GetNMarkers_Optimal <- function(marker_list, signature, score = "correlation", filter = "percent") {
  if (score == "correlation") {
    score_fun <- function(sig_filt) {
      tmp <- cor(sig_filt)
      return(mean(abs(tmp[upper.tri(tmp)])))
    }
  }
  else if (score == "condition" ) {
    score_fun <- kappa
  }

  if (filter == "percent") {
    incs <- seq(0.001, 1, 0.001)
    filter_fun <- function(mkrs, inc) { mkrs[1:ceiling(length(mkrs)*inc)] }
  }
  else { # Fixed
    incs <- seq(1, 500, 1)
    filter_fun <- function(mkrs, inc) { mkrs[1:min(length(mkrs), inc)] }
  }

  vals <- sapply(incs, function(X) {
    markers_filt <- lapply(marker_list, filter_fun, X)
    sig_filt <- signature[unlist(markers_filt),]
    return(score_fun(sig_filt))
  })

  ind <- which.min(vals)

  return(incs[ind])
}


ConvertToROSMAPCelltypes <- function(df, remove_unused = TRUE) {
  orig_cols <- colnames(df)

  # Some datasets are missing one or more of these vascular types
  for (col in c("Endo", "Peri", "VLMC")) {
    if (!(col %in% colnames(df))) {
      df <- cbind(df, rep(0, nrow(df)))
      colnames(df)[ncol(df)] <- col
    }
  }

  df <- as.data.frame(df) %>%
    mutate(Neuro = Exc + Inh,
           Oligo = Oligo + OPC,
           Endo = Endo + Peri + VLMC)

  cols <- c("Astro", "Endo", "Micro", "Neuro", "Oligo")

  if (remove_unused == FALSE) {
    cols <- sort(unique(c(orig_cols, cols)))
  }

  return(df %>% select(all_of(cols)))
}


# OrderMarkers_ByCorrelation: By default markers are ordered from highest to
# lowest difference in expression between the target cell type and other cell
# types. However this doesn't guarantee that the top markers are at all
# correlated, or that the markers are the most informative for the data set
# being tested. This function filters the marker list and re-orders it such
# that for each cell type, the marker list is now the largest possible set of
# marker genes that are all positively correlated with each other in the test
# data set. The remaining markers are then ordered from highest to lowest
# average correlation with each other.
# NOTE: This problem can be solved exactly by treating the correlation matrix
#       as an adjacency graph and using igraph::largest_cliques(), however the
#       run-time increases exponentially with number of genes and isn't
#       realistic for our needs. The code below is a greedy approximation.
#
# Arguments:
#   marker_list = a list of marker gene names, one list entry per cell type
#   data = a gene x sample expression matrix (or data.frame) for calculating
#          correlation. This is usually the test data set.
#
# Returns:
#   an updated list of marker genes, one list entry per cell type
OrderMarkers_ByCorrelation <- function(marker_list, data) {
  # Assumption that if the values in data aren't very large, this is log-scale
  # data that needs to be put into linear scale
  if (max(data) < 100) {
    data <- 2^data-1
  }

  # For each cell type, update the marker list
  new_list <- sapply(names(marker_list), function(N) {
    markers <- marker_list[[N]]
    markers <- markers[markers %in% rownames(data)]
    corr_mat <- cor(as.matrix(t(data[markers,])))

    # Re-order genes by average correlation, highest first
    corr_means <- rowMeans(corr_mat)
    corr_means <- sort(corr_means, decreasing = TRUE)

    # Start at the most positively correlated genes, iteratively add genes that
    # have only positive correlations with the existing set of genes
    new_markers <- c(names(corr_means)[1])
    for (m in 2:length(corr_means)) {
      marker <- names(corr_means)[m]
      tmp <- corr_mat[marker, new_markers]
      if (all(tmp >= 0)) {
        new_markers <- c(new_markers, marker)
      }
    }

    # Sort by average correlation within this group of markers
    new_means <- rowMeans(corr_mat[new_markers, new_markers])
    return(new_markers[order(new_means, decreasing = TRUE)])
  })

  return(new_list)
}
