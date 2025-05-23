# Helper functions that are used in multiple scripts.
# NOTE: For the calculations below that use model.matrix and matrix
# multiplication, using this method to sum counts is way faster than iterating
# over cell type/sample combos and using sum(), colSums(), or rowSums(). We do
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
#   granularity = either "broad_class" or "sub_class", for the level of cell types used
#   reference_input_type = the type of reference data to load. Options are:
#                           "singlecell": the full single cell data set as a
#                                         SingleCellExperiment object
#                           "pseudobulk": pseudobulked "pure" samples as a
#                                         SummarizedExperiment object
#                           "signature": the signature matrix
#                           "cibersortx": the signature matrix created by
#                                         CibersortX PLUS the full single cell
#                                         data set as a SingleCellExperiment
#                                         object
#   output_type = one of "counts", "cpm", "tmm", "tpm", "log_cpm", "log_tmm", or
#                 "log_tpm". See Load_CountsFile for description.
#   regression_method = "none", if raw uncorrected counts should be used for
#                       bulk data, or one of "edger", "deseq2", or "dream", to
#                       use batch-corrected counts from one of those methods.
#                       Applies to bulk data only.
#
# Returns:
#   a list with entries "reference" and "test", containing the reference data
#   and the test data, respectively. Because different algorithms need the
#   input data in different formats, the objects are returned in the format they
#   originated in (a SingleCellExperiment, SummarizedExperiment, or matrix)
#   and it is assumed that the individual algorithms will manipulate the formats
#   afterward. If the reference_input_type is "cibersortx", the "reference" slot
#   will be the single cell data and the list will have an additional field
#   called "cibersortx_signature" with CibersortX's calculated signature matrix.
Load_AlgorithmInputData <- function(reference_data_name, test_data_name,
                                    granularity = "broad_class",
                                    reference_input_type = "singlecell",
                                    output_type = "counts",
                                    regression_method = "none") {
  # Reference input
  if (reference_input_type == "singlecell") {
    reference_obj <- Load_SingleCell(reference_data_name, granularity,
                                     output_type)
  } else if (reference_input_type == "pseudobulk") {
    reference_obj <- Load_PseudobulkPureSamples(reference_data_name,
                                                granularity, output_type)
  } else if (reference_input_type == "signature") {
    reference_obj <- Load_SignatureMatrix(reference_data_name, granularity,
                                          output_type)
  } else if (reference_input_type == "cibersortx") {
    reference_obj <- Load_SignatureMatrix(reference_data_name, granularity,
                                          output_type = "cibersortx")
  } else {
    stop("Invalid reference_input_type specified!")
  }

  # Test data
  if (test_data_name == "sc_samples" || test_data_name == "training") {
    test_obj <- Load_Pseudobulk(reference_data_name, test_data_name,
                                granularity, output_type)
  }
  # ROSMAP, Mayo, or MSBB
  else {
    test_obj <- Load_BulkData(test_data_name, output_type, regression_method)
  }

  # CibersortX signature doesn't need gene filtering, everything else does
  if (reference_input_type != "cibersortx") {
    genes <- intersect(rownames(reference_obj), rownames(test_obj))
    reference_obj <- reference_obj[genes, ]
    test_obj <- test_obj[genes, ]
  }

  return(list("reference" = reference_obj, "test" = test_obj))
}


# Load_AlgorithmInputData_FromParams: Ease-of-use shortcut function that takes in
# a list of parameters and calls Load_AlgorithmInputData from the broken-out list.
# See 'Load_AlgorithmInputData' for arguments and return.
Load_AlgorithmInputData_FromParams <- function(params) {
  return(Load_AlgorithmInputData(params$reference_data_name,
                                 params$test_data_name,
                                 params$granularity,
                                 params$reference_input_type,
                                 params$normalization, # 'output_type' arg in function
                                 params$regression_method))
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
#               per cell type to use. If NULL, default values of c(0.01, 0.02,
#               0.05, 0.1, 0.2, 0.5, 0.75, 1.0, 3, 5, 10, 20, 50, 100, 200, 500)
#               will be used.
#   marker_types = a list where the names of the entries are one of "autogenes",
#                  "dtangle", "seurat", or "deseq2" and the items in each entry
#                  are a list of marker subtypes to use for that algorithm. See
#                  FilterSignature for more detail. Must be a list and not a
#                  vector. If NULL, all 4 marker types with all possible
#                  subtypes will be used.
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
CreateParams_MarkerTypes <- function(n_markers = NULL, marker_types = NULL,
                                     marker_input_types = c("singlecell", "pseudobulk"),
                                     marker_order = c("distance", "correlation")) {
  # Default values for n_markers and marker_types if they are not defined
  if (is.null(n_markers)) {
    n_markers <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                   3, 5, 10, 20, 50, 100, 200, 500)
  }
  if (is.null(marker_types)) {
    marker_types <- list("dtangle" = c("ratio", "diff", "p.value", "regression"),
                         "autogenes" = c("correlation", "distance", "combined"),
                         "seurat" = c("None"),
                         "deseq2" = c("None"))
  }

  marker_types <- melt(marker_types) %>%
    dplyr::rename(marker_type = "L1", marker_subtype = "value")

  params <- tidyr::expand_grid(n_markers = n_markers,
                               marker_types, # this is data frame
                               marker_input_type = marker_input_types,
                               marker_order = marker_order)

  # marker_input_type only applies to dtangle markers
  params$marker_input_type[params$marker_type != "dtangle"] <- "None"

  # We don't use dtangle 'p.value' and 'regression' markers for 'singlecell'
  # input because the compute time for these is too high
  params <- subset(params, !(marker_input_type == "singlecell" &
    marker_subtype %in% c("regression", "p.value")))

  # We don't need to re-order markers when we're using the whole marker set
  params$marker_order[params$n_markers == 1] <- "distance"

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
#   The rest of the arguments are only relevant for filter_level = 3. See
#   their descriptions in CreateParams_MarkerTypes.
#
# Returns:
#   a tibble containing all possible valid combinations of the arguments
CreateParams_FilterableSignature <- function(filter_levels = c(1, 2, 3),
                                             n_markers = NULL,
                                             marker_types = NULL,
                                             marker_input_types = c("singlecell", "pseudobulk"),
                                             marker_order = c("distance", "correlation")) {
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
  params$marker_order[low_filt] <- "distance"

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
  # Sum all counts per gene for each sample/celltype combination
  y <- model.matrix(~ 0 + samples:celltypes)
  count_sums <- singlecell_counts %*% y

  # Sum over all genes for each sample/celltype combination
  count_sums <- colSums(count_sums)

  # names are of the format "samples<sample>:celltypes<celltype>", split
  # them apart by the ":" into a data frame w/ 2 columns
  col_info <- str_split(names(count_sums), ":", simplify = TRUE)
  sample_list <- unique(col_info[, 1])

  pctRNA <- sapply(sample_list, function(samp) {
    # All entries for this sample = 1 entry per cell type
    cols <- which(col_info[, 1] == samp)
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
#   dataset = name of the dataset to load
#   granularity = either "broad_class" or "sub_class"
#
# Returns:
#   vector of average library size of each cell type, normalized to sum to 1.
CalculateA <- function(dataset, granularity) {
  pb <- Load_PseudobulkPureSamples(dataset, granularity, output_type = "counts")

  pb_counts <- assay(pb, "counts")
  count_sums <- colSums(pb_counts)

  samples <- str_replace(names(count_sums), ".*_", "")

  # For each sample
  A_s <- sapply(unique(samples), function(samp) {
    cols <- grepl(samp, names(count_sums))

    # Average library size = sum over all genes for each sample/celltype
    # divided by number of cells for that sample/celltype
    a_s <- count_sums[cols] / pb$n_cells[cols]
    names(a_s) <- pb$celltype[cols]

    # Some samples don't have all cell types
    missing <- setdiff(levels(pb$celltype), names(a_s))
    a_s[missing] <- NA

    a_s <- a_s[levels(pb$celltype)]

    return(a_s)
  })

  # Normalize each sample row to sum to 1
  A_s <- t(sweep(A_s, 2, colSums(A_s, na.rm = TRUE), "/"))

  # Take the means but exclude entries for samples who don't have a certain cell type
  A_s[A_s == 0] <- NA
  A <- colMeans(A_s, na.rm = TRUE)
  A <- A / sum(A) # Enforce sum to 1
}


# CalculateSignature: Creates a cell-type-specific "signature" matrix that
# describes the expected count of each gene for each cell type.
#
# Uses pseudobulk pure samples, where each sample is the sum of all raw counts
# of a specific cell type from a specific donor. The pseudobulk counts are
# normalized to CPM, and the CPM values of each gene for each cell type are
# averaged together.
#
# Arguments:
#   dataset = the name of the dataset to load
#   granularity = either "broad_class" or "sub_class"
#   output_type = either "cpm" or "tmm", for whether to use pure CPM or normalize
#                 using tmm factors
#   geom_mean = whether to take the geometric mean instead of the arithmetic mean
#
# Returns:
#   matrix of average CPMs for each gene, for each cell type (rows = genes,
#   columns = cell types)
CalculateSignature <- function(dataset, granularity, output_type, geom_mean = FALSE) {
  pb <- Load_PseudobulkPureSamples(dataset, granularity, output_type)

  pb_cpm <- assay(pb, "counts")

  # Get the mean over all samples, for each gene and cell type
  sig <- sapply(levels(pb$celltype), function(ct) {
    cols <- which(pb$celltype == ct)
    if (geom_mean) {
      log_means <- rowMeans(log(pb_cpm[, cols] + 1), na.rm = TRUE)
      cpm_means <- exp(log_means) - 1
      cpm_means[cpm_means < 0] <- 0

    } else {
      cpm_means <- rowMeans(pb_cpm[, cols], na.rm = TRUE)
    }
    return(cpm_means)
  })

  return(sig)
}


# FilterSignature: filters the signature matrix, which includes all genes
# present in the single cell dataset, down to a set of genes which are more
# informative. Because we don't know exactly what will be the most informative
# for each algorithm, this function has multiple filtering settings:
#   Filter level 0 = don't filter at all
#                1 = keep genes where at least one cell type has > 1 cpm
#                2 = keep genes where at least one cell type has > 10 cpm
#                3 = use algorithm-determined markers, and filter the markers
#                    further to use the top <n_markers> genes of each cell
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
#   n_markers = percent of markers to keep (range 0-1) OR an integer describing
#               the number of markers to use per cell type (range 2-Inf). For
#               example, passing in 10 will use 10 markers for each cell type,
#               while passing in 0.5 will use 50% of each cell type's markers
#               (different number of markers per cell type).
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
#   test_data_name = the name of the test data set. Only used if marker_order
#                    is "correlation".
#   normalization = the normalization strategy. Only used if marker_order is
#                   "correlation".
#   regression_method = the regression method for the bulk data. Only used if
#                   marker_order is "correlation".
#
# Returns:
#   signature matrix containing only the genes that pass the specified filters.
#   rows = genes, columns = cell types.
FilterSignature <- function(signature, filter_level = 1, reference_data_name = NULL,
                            granularity = NULL, n_markers = 1.0,
                            marker_type = "dtangle", marker_subtype = "diff",
                            marker_input_type = "pseudobulk",
                            marker_order = "distance", test_data_name = NULL,
                            normalization = NULL, regression_method = NULL) {
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
    markers <- FilterMarkers(reference_data_name, granularity, n_markers,
                             marker_type, marker_subtype, marker_input_type,
                             marker_order,
                             available_genes = rownames(signature),
                             test_data_name = test_data_name,
                             normalization = normalization,
                             regression_method = regression_method)
    if (is.null(markers)) {
      return(NULL)
    }

    signature <- signature[unique(unlist(markers)), ]
  }

  return(signature)
}


# Shortcut function for FilterSignature where the 'params' object can get passed
# in instead of unpacking all the variables inside it.
#
# Arguments:
#   signature = signature matrix of cpm values for each gene for each cell type
#               (rows = genes, columns = cell types)
#   params = a named vector or one-row data frame with the parameters to use.
#            Must contain variables with the same names as the arguments to
#            Filter_Signature (except for "signature").
#
# Returns:
#   signature matrix containing only the genes that pass the specified filters.
#   rows = genes, columns = cell types.
FilterSignature_FromParams <- function(signature, params) {
  return(FilterSignature(signature,
                         filter_level = params$filter_level,
                         reference_data_name = params$reference_data_name,
                         granularity = params$granularity,
                         n_markers = params$n_markers,
                         marker_type = params$marker_type,
                         marker_subtype = params$marker_subtype,
                         marker_input_type = params$marker_input_type,
                         marker_order = params$marker_order,
                         test_data_name = params$test_data_name,
                         normalization = params$normalization,
                         regression_method = params$regression_method))
}


# FilterMarkers: filters the list of cell type markers, which includes all
# markers found for each cell type, down to a set of genes which are more
# informative.
#
# Arguments:
#   See descriptions for FilterSignature arguments
#   available_genes = a vector of gene names that are valid for the current data
#                     set. This is used to filter out marker genes that don't
#                     exist in the data.
#
# Returns:
#   a list where each item is a vector of marker genes for a given cell type
FilterMarkers <- function(reference_data_name, granularity, n_markers,
                          marker_type, marker_subtype, marker_input_type,
                          marker_order, available_genes, test_data_name = NULL,
                          normalization = NULL, regression_method = NULL) {
  markers <- Load_Markers(reference_data_name, granularity, marker_type,
                          marker_subtype,
                          input_type = marker_input_type)

  if (is.null(markers)) {
    return(NULL)
  }

  if (marker_order == "correlation") {
    data_name <- paste(test_data_name, normalization, regression_method,
                       sep = "_")
    data_name <- str_replace(data_name, "log_", "")
    data_name <- str_replace(data_name, "counts", "cpm")

    markers <- markers$ordered_by_correlation[[data_name]]
  } else {
    markers <- markers$filtered
  }

  # Remove genes that don't exist in the data
  markers <- lapply(markers, function(X) {
    intersect(X, available_genes)
  })

  # We can't use these markers if any cell type has 1 or 0 markers after filtering
  if (any(lengths(markers) < 2)) {
    return(NULL)
  }

  # Percentage of each cell type's markers
  if (n_markers <= 1) {
    n_markers <- ceiling(lengths(markers) * n_markers)
  }
  # Fixed number of markers for each cell type
  else {
    n_markers <- sapply(lengths(markers), min, n_markers)
  }

  if (any(n_markers == 0)) {
    return(NULL)
  }

  markers <- lapply(names(markers), function(N) {
    markers[[N]][1:n_markers[N]]
  })

  return(markers)
}


# Shortcut function for FilterMarkers where the 'params' object can get passed
# in instead of unpacking all the variables inside it.
#
# Arguments:
#   available_genes = a vector of genes that are valid for this data. Usually
#                     it's a vector of genes that exist in both the single cell
#                     and bulk data set.
#   params = a named vector or one-row data frame with the parameters to use.
#            Must contain variables with the same names as the arguments to
#            Filter_Markers (except for "available_genes").
#
# Returns:
#   a list where each item is a vector of marker genes for a given cell type
FilterMarkers_FromParams <- function(available_genes, params) {
  return(FilterMarkers(reference_data_name = params$reference_data_name,
                       granularity = params$granularity,
                       n_markers = params$n_markers,
                       marker_type = params$marker_type,
                       marker_subtype = params$marker_subtype,
                       marker_input_type = params$marker_input_type,
                       marker_order = params$marker_order,
                       available_genes = available_genes,
                       test_data_name = params$test_data_name,
                       normalization = params$normalization,
                       regression_method = params$regression_method))
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


# Converts the cell type names used in this pipeline to the cell type names used
# in the ROSMAP IHC data. This function only works for broad_class cell types,
# as the IHC data has no subclasses.
#
# Arguments:
#   df = a data frame where rows are samples and columns are cell types
#   remove_unused = whether to remove this pipeline's cell type names from the
#                   result when they differ from ROSMAP's names, or whether to
#                   leave them in. (For example, whether to remove "Excitatory"
#                   and "Inhibitory" columns since they are replaced by "Neuro")
#
# Returns:
#   a data frame where rows are samples and columns are cell types with names
#   that match the ROSMAP IHC data (if remove_unused is TRUE) or may
#   additionally contain cell type names from this pipeline that do not overlap
#   with ROSMAP names (if remove_unused is FALSE).
ConvertToROSMAPCelltypes <- function(df, remove_unused = TRUE) {
  orig_cols <- colnames(df)

  # ROSMAP IHC doesn't distinguish between neuronal types or between
  # oligodendrocytes and OPCs. We assume here that "Endo" is approximately
  # equivalent to our pipeline's "Vascular" cell type, although Vascular
  # contains endothelial cells, pericytes, and VLMCs.
  df <- as.data.frame(df) %>%
    mutate(Neuro = Excitatory + Inhibitory,
           Oligo = Oligodedrocyte + OPC,
           Endo = Vascular)

  cols <- c("Astro", "Endo", "Micro", "Neuro", "Oligo")

  if (remove_unused == FALSE) {
    cols <- sort(unique(c(orig_cols, cols)))
  }

  return(df %>% select(all_of(cols)))
}


# OrderMarkers_ByCorrelation: By default, markers are ordered from highest to
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
    data <- 2^data - 1
  }

  # For each cell type, update the marker list
  new_list <- sapply(names(marker_list), function(N) {
    markers <- marker_list[[N]]
    markers <- markers[markers %in% rownames(data)]

    # Markers don't need to be ordered if there are 3 or less after filtering
    if (length(markers) <= 3) {
      return(markers)
    }

    markers <- markers[rowSums(data[markers, ] > 0) >= 3]

    # Check again for 3 or less markers after filtering by gene expression
    if (length(markers) <= 3) {
      return(markers)
    }

    corr_mat <- cor(as.matrix(t(data[markers, ])))

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


# Clean_BulkCovariates - takes a covariates dataframe for one of the bulk
# datasets, makes sure that categorical variables are factors, scales numerical
# variables, and merges the cleaned covariates dataframe with the metadata
# dataframe.
#
# Arguments:
#   metadata - the metadata dataframe (colData()) from a SummarizedExperiment
#   covariates - a dataframe of covariates, where rows are samples and columns
#                are the covariates
#   scale_numerical - TRUE or FALSE, whether to scale numeric columns
#
# Returns:
#   a dataframe with merged metadata and cleaned covariates
Clean_BulkCovariates <- function(metadata, covariates, scale_numerical = TRUE) {
  covariates <- subset(covariates, specimenID %in% rownames(metadata))

  covariates <- covariates |>
    mutate(
      ageDeath_num = suppressWarnings(as.numeric(ageDeath)),
      ageDeath = case_when(
        !is.na(ageDeath_num) ~ cut(
          ageDeath_num,
          breaks = c(0, 65, 70, 75, 80, 85, 90),
          labels = c("Under 65", "65 - 69", "70 - 74", "75 - 79", "80 - 84", "85 - 89"),
          right = FALSE,
          include.lowest = TRUE
        ),
        ageDeath == "90_or_over" ~ "90+", # Fix for Mayo
        .default = ageDeath # 90+ or NA
      )
    ) |>
    select(-ageDeath_num)

  covariates$ageDeath <- factor(covariates$ageDeath,
                                 levels = c("Under 65", "65 - 69", "70 - 74",
                                            "75 - 79", "80 - 84", "85 - 89", "90+"))

  for (col in c("diagnosis", "sex", "race", "spanish", "ethnicity",
                "individualID", "apoeGenotype", "projid", "batch", "Braak",
                "Thal", "CERAD")) {
    if (col %in% colnames(covariates)) {
      covariates[, col] <- factor(covariates[, col])
    }
  }

  # Remove duplicate columns that already exist in colData
  covariates <- covariates %>% select(-diagnosis, -tissue)

  # Merge covariates into the metadata
  sample_order <- rownames(metadata)
  metadata <- merge(metadata, covariates,
                    by.x = "sample", by.y = "specimenID",
                    sort = FALSE)
  rownames(metadata) <- metadata$sample

  # Put the data frame back in the original order, as merge might change it
  metadata <- data.frame(metadata[sample_order, ])

  # Scale numerical covariates if scale == TRUE
  if (scale_numerical) {
    for (colname in colnames(metadata)) {
      if (is.numeric(metadata[, colname])) {
        metadata[, colname] <- as.numeric(scale(metadata[, colname]))
      }
    }
  }

  # Remove columns that are all NA or all the same value
  all_na <- sapply(colnames(metadata), function(col_name) {
    all(is.na(metadata[, col_name]))
  })

  all_same <- sapply(colnames(metadata), function(col_name) {
    all(metadata[, col_name] == metadata[1, col_name])
  })
  all_same[is.na(all_same)] <- FALSE

  cols_keep <- !all_na & !all_same

  return(metadata[, cols_keep])
}


# FileParams_FromParams - takes a parameters data frame and extracts the
# columns that are used to name output/error files. If params_df is from a
# single output or error file, all values for those columns are the same for
# each row so this function will return a one-row data frame.
#
# Arguments:
#   params_df - a data frame of any number of rows that must contain the columns
#               in the select statement below
#
# Returns:
#   a data frame
FileParams_FromParams <- function(params_df) {
  params_df %>%
    dplyr::select_at(Get_ParameterColumnNames()) %>%
    dplyr::distinct() %>%
    as.data.frame()
}


# Get_ParameterColumns - Helper function to get the names of the parameters that
# exist across all deconvolution trials. These are used all over the code so
# this avoids hard-coding the list in multiple places.
#
# Arguments: none
# Returns: a vector of column names
Get_ParameterColumnNames <- function() {
  return(c("algorithm", "reference_data_name", "test_data_name", "granularity",
           "reference_input_type", "normalization", "regression_method"))
}


# List_to_DF - a helper function for turning a list of data frames into a single
# data frame via rbind. Providing a sublist_name will extract all items with that
# name from each item in the list and rbind those together instead.
#
# Arguments:
#   input_list - a list of data frames or a list of lists of data frames
#   sublist_name - if NULL, input_list must be a list of data frames and this
#                  function will rbind all items in input_list together.
#                  If a value is provided, input_list must be a list where each
#                  item contains a named list of data frames whose name matches
#                  sublist_name. In that case, sublists with that name are
#                  extracted from each item in input_list and rbind-ed together.
#
# Returns:
#   a data frame
List_to_DF <- function(input_list, sublist_name = NULL) {
  if (is.null(sublist_name)) {
    return(do.call(rbind, input_list))
  }

  return(do.call(rbind, lapply(input_list, "[[", sublist_name)))
}
