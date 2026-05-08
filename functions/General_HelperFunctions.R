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

# Shortcut functions for listing out all single cell and bulk data sets
all_singlecell_datasets <- function() {
  c("cain", "lau", "mathys", "seaRef")
}

all_bulk_datasets <- function() {
  c("Mayo_CBE", "Mayo_TCX",
    "MSBB_FP", "MSBB_IFG", "MSBB_PHG", "MSBB_STG",
    "ROSMAP_ACC", "ROSMAP_DLPFC_1", "ROSMAP_DLPFC_2", "ROSMAP_PCC")
}

is_bulk <- function(dataset) {
  return(dataset %in% c("Mayo", "MSBB", "ROSMAP", all_bulk_datasets()))
}

is_singlecell <- function(dataset) {
  return(dataset %in% all_singlecell_datasets())
}


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
#   normalization = one of "counts", "cpm", "tmm", "tpm", "log_cpm", "log_tmm",
#                   or "log_tpm". See Load_CountsFile for description.
#   regression_method = "none", if raw uncorrected counts should be used for
#                       bulk data, or one of "edger", "deseq2", or "combat", to
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
                                    normalization = "counts",
                                    regression_method = "none") {
  # Reference input
  if (reference_input_type == "singlecell") {
    reference_obj <- Load_SingleCell(reference_data_name, granularity,
                                     normalization)
  } else if (reference_input_type == "pseudobulk") {
    reference_obj <- Load_PseudobulkPureSamples(reference_data_name,
                                                granularity, normalization)
  } else if (reference_input_type == "signature") {
    reference_obj <- Load_SignatureMatrix(reference_data_name, granularity,
                                          normalization)
  } else if (reference_input_type == "cibersortx") {
    reference_obj <- Load_SignatureMatrix(reference_data_name, granularity,
                                          normalization = "cibersortx")
  } else if (reference_input_type == "scaden") {
    reference_obj <- Load_SimulatedScadenData(reference_data_name, granularity,
                                              normalization)
  } else {
    stop("Invalid reference_input_type specified!")
  }

  # Test data
  if (test_data_name == "sc_samples" || test_data_name == "training") {
    test_obj <- Load_Pseudobulk(reference_data_name, test_data_name,
                                granularity, normalization)
  }
  # ROSMAP, Mayo, or MSBB
  else {
    test_obj <- Load_BulkData(test_data_name, normalization, regression_method)
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
                                 params$normalization,
                                 params$regression_method))
}


# CreateParams_MarkerTypes - Creates a parameter matrix using tidyr::expand_grid
# (which can take data frames in the input) that has all variables required for
# loading different combinations of markers: n_markers, marker_type, and
# marker_subtype. Invalid combinations of these variables are removed from the
# parameter set before returning.
#
# Arguments:
#   n_markers = a vector containing one or more integers (range 3-Inf)
#               specifying how many markers per cell type to use. If NULL,
#               default values of c(3, 5, 10, 20, 50, 100, 200, 500) will be
#               used.
#   marker_types = a list where the names of the entries are one of "autogenes",
#                  "dtangle", "seurat", or "deseq2" and the items in each entry
#                  are a list of marker subtypes to use for that algorithm. See
#                  FilterSignature for more detail. Must be a list and not a
#                  vector. If NULL, all 4 marker types with all possible
#                  subtypes will be used.
#   marker_order = "distance" or "correlation", whether markers should be
#                  ordered by largest expression difference between cell types
#                  or by correlation with other markers for the same cell type
#   filter_ad_genes = whether to filter genes that change with AD from the set
#                  of markers
#
# Returns:
#   a tibble containing all possible valid combinations of the arguments
CreateParams_MarkerTypes <- function(n_markers = NULL, marker_types = NULL,
                                     marker_order = c("distance", "correlation"),
                                     filter_logfc_genes = TRUE,
                                     filter_ad_genes = FALSE) {
  # Default values for n_markers and marker_types if they are not defined
  if (is.null(n_markers)) {
    n_markers <- c(3, 5, 10, 20, 50, 100, 200, 500)
  }
  if (is.null(marker_types)) {
    marker_types <- list("dtangle" = c("ratio", "diff"),
                         "autogenes" = c("correlation", "distance", "combined"),
                         "seurat" = c("None"),
                         "deseq2" = c("None"))
  }

  marker_types <- melt(marker_types) %>%
    dplyr::rename(marker_type = "L1", marker_subtype = "value")

  params <- tidyr::expand_grid(n_markers = n_markers,
                               marker_types, # this is a data frame
                               marker_order = marker_order,
                               filter_logfc_genes = filter_logfc_genes,
                               filter_ad_genes = filter_ad_genes)

  params <- params %>% distinct()
  return(params)
}


# CreateParams_FilterableSignature - Creates a parameter matrix using
# tidyr::expand_grid (which can take data frames in the input) that has all
# variables required for filtering a signature matrix: filter_level, n_markers,
# marker_type, and marker_subtype. Invalid combinations of these variables are
# removed from the parameter set before returning.
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
                                             marker_order = c("distance", "correlation"),
                                             filter_logfc_genes = TRUE,
                                             filter_ad_genes = FALSE) {
  params_tmp <- CreateParams_MarkerTypes(n_markers, marker_types, marker_order,
                                         filter_logfc_genes, filter_ad_genes)

  params <- tidyr::expand_grid(filter_level = filter_levels,
                               params_tmp)

  # Some filter_type / n_markers combos are not valid, get rid of them
  # (filter levels 0, 1 & 2 don't use n_markers or marker_type arguments)
  low_filt <- params$filter_level < 3
  params$marker_type[low_filt] <- "None"
  params$marker_subtype[low_filt] <- "None"
  params$n_markers[low_filt] <- -1
  params$marker_order[low_filt] <- "distance"
  params$filter_logfc_genes[low_filt] <- FALSE
  params$filter_ad_genes[low_filt] <- FALSE


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
#   sce = a SingleCellExperiment, which must have "sample", "broad_class", and
#         "sub_class" columns in its metadata.
#   granularity = either "broad_class" or "sub_class"
#
# Returns:
#   matrix of percentages (rows = samples, cols = cell types), whose rows sum
#   to 1.
CalculatePercentRNA <- function(sce, granularity) {
  sce$celltype <- colData(sce)[, granularity]
  aggr <- scuttle::aggregateAcrossCells(sce,
                                        ids = paste(sce$sample, sce$celltype),
                                        statistics = "sum")

  count_sums <- colSums(counts(aggr))
  count_sums <- data.frame(ids = names(count_sums), counts = as.numeric(count_sums))

  pcts <- merge(as.data.frame(colData(aggr)), count_sums) |>
    group_by(sample) |>
    mutate(pct = counts / sum(counts)) |>
    select(sample, celltype, pct) |>
    tidyr::pivot_wider(id_cols = "sample",
                       names_from = "celltype",
                       values_from = "pct") |>
    tibble::column_to_rownames("sample") |>
    as.matrix()

  # Fill in 0's for any samples that didn't have a particular cell type
  pcts[is.na(pcts)] <- 0

  return(pcts)
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
  pb <- Load_PseudobulkPureSamples(dataset, granularity, normalization = "counts")

  count_sums <- colSums(assay(pb, "counts"))
  count_sums <- data.frame(sample = names(count_sums), count_sums = count_sums)

  A <- merge(as.data.frame(colData(pb)), count_sums) |>
    group_by(sample_orig) |>
    mutate(
      # sum of counts for a specific cell type / number of cells of that cell type
      # gives average library size
      avg_lib_size = count_sums / ncells,
      # Normalize average library size across each sample so
      # sum(lib_size of every cell type) for each sample = 1
      norm_lib_size = avg_lib_size / sum(avg_lib_size)
    ) |>
    select(sample_orig, celltype, norm_lib_size) |>
    # Make a sample x celltype matrix
    tidyr::pivot_wider(id_cols = sample_orig,
                       names_from = "celltype",
                       values_from = "norm_lib_size") |>
    tibble::column_to_rownames("sample_orig") |>
    # Mean library size for each cell type across all samples, ignoring NAs
    # where a sample doesn't have that cell type.
    colMeans(na.rm = TRUE)

  A <- A / sum(A) # Enforce sum to 1
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
#   marker_order = "distance" or "correlation". Whether the markers should be
#                  ordered by distance (prioritize markers with the highest
#                  expression difference between the cell type and all other
#                  cell types) or correlation (prioritize markers with the
#                  highest correlation to other markers for that cell type)
#   filter_logfc_genes = if TRUE, use a filtered set of markers where each gene
#                  has a >1 logfc (>0.25 for subclass) between the target cell
#                  type and any other cell type. If FALSE, use all discovered
#                  markers.
#   filter_ad_genes = remove genes that change with AD in the reference data set
#                  from the list of markers, if markers are being used.
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
                            granularity = NULL, n_markers = 3,
                            marker_type = "dtangle", marker_subtype = "diff",
                            marker_order = "distance", filter_logfc_genes = TRUE,
                            filter_ad_genes = FALSE, test_data_name = NULL,
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
                             marker_type, marker_subtype, marker_order,
                             available_genes = rownames(signature),
                             filter_logfc_genes = filter_logfc_genes,
                             filter_ad_genes = filter_ad_genes,
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
                         marker_order = params$marker_order,
                         filter_logfc_genes = params$filter_logfc_genes,
                         filter_ad_genes = params$filter_ad_genes,
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
                          marker_type, marker_subtype, marker_order,
                          available_genes, filter_logfc_genes = TRUE,
                          filter_ad_genes = FALSE, test_data_name = NULL,
                          normalization = NULL, regression_method = NULL) {
  markers <- Load_Markers(reference_data_name, granularity, marker_type,
                          marker_subtype)

  if (is.null(markers)) {
    return(NULL)
  }

  marker_category <- ifelse(filter_logfc_genes, "filtered", "all")

  if (marker_order == "correlation") {
    data_name <- paste(test_data_name, normalization, regression_method,
                       sep = "_")
    data_name <- str_replace(data_name, "log_", "") |>
      str_replace("counts_tpm", "tpm") |>
      str_replace("counts", "cpm")

    markers_out <- markers$ordered_by_correlation[[marker_category]][[data_name]]
  } else {
    markers_out <- markers[[marker_category]]
  }

  markers_out <- lapply(markers_out, function(X) {
    markers$genes[X]
  })

  # Remove genes that change in AD for this data set if applicable
  if (filter_ad_genes) {
    markers_out <- lapply(markers_out, function(mkrs) {
      setdiff(mkrs, markers$ad_gene_exclusions)
    })
  }

  # Remove genes that don't exist in the data
  markers_out <- lapply(markers_out, function(X) {
    intersect(X, available_genes)
  })

  # Get the minimum of n_markers and the actual number of markers for each cell type
  n_markers <- sapply(lengths(markers_out), min, n_markers)

  # We can't use these markers if any cell type has < 3 markers after filtering
  if (any(n_markers < 3)) {
    return(NULL)
  }

  markers_filt <- lapply(names(markers_out), function(N) {
    markers_out[[N]][1:n_markers[N]]
  })

  names(markers_filt) <- names(markers_out)

  return(markers_filt)
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
#            FilterMarkers (except for "available_genes").
#
# Returns:
#   a list where each item is a vector of marker genes for a given cell type
FilterMarkers_FromParams <- function(available_genes, params) {
  return(FilterMarkers(reference_data_name = params$reference_data_name,
                       granularity = params$granularity,
                       n_markers = params$n_markers,
                       marker_type = params$marker_type,
                       marker_subtype = params$marker_subtype,
                       marker_order = params$marker_order,
                       available_genes = available_genes,
                       filter_logfc_genes = params$filter_logfc_genes,
                       filter_ad_genes = params$filter_ad_genes,
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
# that for each cell type, the marker list is now the set of marker genes with
# a positive average correlation with all other marker genes. These markers are
# then ordered from highest to lowest average correlation with each other.
#
# Arguments:
#   marker_list = a list of marker gene names, one list entry per cell type
#   data = a gene x sample expression matrix (or data.frame) for calculating
#          correlation. This is usually the test data set. The data should be
#          log-normalized.
#
# Returns:
#   an updated list of marker genes, one list entry per cell type
OrderMarkers_ByCorrelation <- function(marker_list, data) {
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
    diag(corr_mat) <- NA

    # Re-order genes by average correlation, highest first
    corr_means <- rowMeans(corr_mat, na.rm = TRUE)
    corr_means <- sort(corr_means, decreasing = TRUE)

    new_markers <- names(corr_means)[corr_means >= 0]

    # If there aren't enough markers that are positively correlated with each
    # other, just return all markers sorted by average correlation
    if (length(new_markers) < 3) {
      return(names(corr_means))
    }

    # Sort by average correlation within this group of markers now that negative
    # markers have been removed
    new_means <- rowMeans(corr_mat[new_markers, new_markers], na.rm = TRUE) |>
      sort(decreasing = TRUE)
    return(names(new_means))
  })

  return(new_list)
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


# Helper function for marker finding. Given a data frame of average expression
# for each gene/cell type, filter to genes where the log2-FC between the
# target celltype and the highest-expressing non-target celltype is above a
# certain threshold. The threshold is 1 for broad class and 0.25 for sub class.
# A gene can only be a marker for exactly one cell type for broad class, but can
# be a marker for up to 2 cell types for sub class. This function handles both
# cases: for broad class, the log2-FC is expr[celltype] - max(expr[!celltype]),
# while for sub class it is similar except that that instead of looking for the
# max among all non-target celltypes, it looks for the max among all non-target
# celltypes for which this gene is not a marker. This allows for two cell types
# with the same marker gene to have similar expression as long as the rest of
# the cell types have lower expression.
Get_QualityMarkers <- function(expr_df, markers, granularity) {
  # For broad class, get the top-expressing cell type. For sub
  # class, get the top two expressing cell types.
  #cell_thresh <- 1 #ifelse(granularity == "broad_class", 1, 2)
  celltypes <- colnames(expr_df) |> sort()

  # Turn expression df into long format
  expr_df <- expr_df[unique(markers), ] |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(cols = -gene,
                        names_to = "celltype",
                        values_to = "expr")

  # Find the top expressor of each marker gene
  expr_mod <- expr_df |>
    group_by(gene) |>
    slice_max(order_by = expr, n = 1) # cell_thresh)

  # Assign marker status to top expressors only
  expr_df <- expr_df |>
    group_by(celltype) |>
    mutate(
      is_marker = gene %in% expr_mod$gene[expr_mod$celltype == unique(celltype)]
    ) |>
    ungroup()

  # Helper function
  getLog2FC <- function(celltype, expr) { #, is_marker) {
    sapply(celltype, function(ct) {
      expr[celltype == ct] - max(expr[celltype != ct]) # & !is_marker])
    })
  }

  # Filtered marker list sorted by log2FC
  sorted_logfc <- expr_df |>
    group_by(gene) |>
    mutate(log2FC = getLog2FC(celltype, expr)) |> #, is_marker)) |>
    subset(is_marker == TRUE) |>
    dplyr::arrange(desc(log2FC))

  # Create one full list containing all genes that were identified as
  # markers, and one list filtered to (log2FC between highest and
  # second-highest expression) > 1 (or 0.25 for sub_class)
  markers_all <- sapply(celltypes, function(ct) {
    return(sorted_logfc$gene[sorted_logfc$celltype == ct])
  })

  thresh <- config::get("step07_find_markers")$lfc_threshold[[granularity]]
  markers_filt <- sapply(celltypes, function(ct) {
    return(sorted_logfc$gene[sorted_logfc$celltype == ct &
                               sorted_logfc$log2FC >= thresh])
  })

  return(list("all" = markers_all,
              "filtered" = markers_filt,
              "sorted_logfc" = sorted_logfc))
}
