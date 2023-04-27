# Helper functions for loading or saving data sets and markers.

library(SummarizedExperiment)
library(SingleCellExperiment)
library(scuttle)
library(stringr)
source("Filenames.R")

##### Single cell / Pseudobulk objects #####

# Load_SingleCell: reads a SingleCellExperiment data set from a file and
# transforms the counts according to output_type. Single-cell-specific wrapper
# function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad" or "fine", for which level of cell types to
#                 use in the metadata
#   output_type = either "counts", "cpm", "logcpm", or "log1p_cpm" to determine
#                  how the counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm) of non-zero cpm entries
#                   "log1p_cpm" will take the log2(cpm+1) of cpms
#
# Returns:
#   a SingleCellExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable. The colData DataFrame also has a column added called
#   "celltype", which is populated with either the broad or fine cell type
#   assignments based on granularity.
Load_SingleCell <- function(dataset, granularity, output_type = "counts") {
  sc_file <- file.path(dir_input, str_glue("{dataset}_sce.rds"))

  singlecell <- Load_CountsFile(sc_file, output_type)
  metadata <- colData(singlecell)

  # Assign the column "celltype" to be either the broad or fine cell types
  if (granularity == "broad") {
    metadata$celltype <- metadata$broadcelltype
  }
  else {
    metadata$celltype <- metadata$subcluster
  }

  colData(singlecell) <- metadata
  return(singlecell)
}

# Load_PseudobulkPureSamples: reads a SummarizedExperiment data set from a file
# and transforms the counts according to output_type. Pure sample pseudobulk-
# specific wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load in.
#   output_type = either "counts", "cpm", "logcpm", or "log1p_cpm" to determine
#                 how the counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm) of non-zero cpm entries
#                   "log1p_cpm" will take the log2(cpm+1) of cpms
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable. It also gets colData set to a DataFrame containing one column
#   for "celltype", which is populated with cell type assignments, which can be
#   extracted from the sample names in the data set.
Load_PseudobulkPureSamples <- function(dataset, granularity, output_type = "counts") {
  pb_file <- str_glue(paste0("pseudobulk_{dataset}_puresamples_",
                             "{granularity}.rds"))
  pb_file <- file.path(dir_pseudobulk, pb_file)

  pseudobulk <- Load_CountsFile(pb_file, output_type)
  return(pseudobulk)
}


# Save_PseudobulkPureSamples: saves a SummarizedExperiment object to a file with
# a specific filename format (which should match the format used in
# Load_PseudobulkPureSamples).
#
# Arguments:
#   se = a SummarizedExperiment object
#   dataset = the name of the data set to load in
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load in.
#
# Returns:
#   nothing
Save_PseudobulkPureSamples <- function(se, dataset, granularity) {
  filename <- str_glue("pseudobulk_{dataset}_puresamples_{granularity}.rds")
  saveRDS(se, file = file.path(dir_pseudobulk, filename))
}


# Load_Pseudobulk: reads a SummarizedExperiment data set from a file and
# transforms the counts according to output_type. Donor or training pseudobulk-
# specific wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   data_type = either "donors" or "training", for which type of pseudobulk to
#               load in.
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load in.
#   output_type = either "counts", "cpm", "logcpm", or "log1p_cpm to determine
#                 how the counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm) of non-zero cpm entries
#                   "log1p_cpm" will take the log2(cpm+1) of cpms
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable.
Load_Pseudobulk <- function(dataset, data_type, granularity, output_type = "counts") {
  pb_file <- str_glue(paste0("pseudobulk_{dataset}_{data_type}_",
                             "{granularity}.rds"))
  pb_file <- file.path(dir_pseudobulk, pb_file)

  pseudobulk <- Load_CountsFile(pb_file, output_type)
  return(pseudobulk)
}


# Save_Pseudobulk: saves a SummarizedExperiment object to a file with a specific
# filename format (which should match the format used in Load_Pseudobulk).
#
# Arguments:
#   se = a SummarizedExperiment object
#   dataset = the name of the data set to load in
#   data_type = either "donors" or "training", for which type of pseudobulk to
#               load in.
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load in.
#
# Returns:
#   nothing
Save_Pseudobulk <- function(se, dataset, data_type, granularity) {
  filename <- str_glue("pseudobulk_{dataset}_{data_type}_{granularity}.rds")
  saveRDS(se, file = file.path(dir_pseudobulk, filename))
}


# Load_CountsFile: reads a SummarizedExperiment or SingleCellExperiment from
# a file and transforms the counts according to output_type.
#
# Arguments:
#   filename = the name of of the file to read in, including file path
#   output_type = either "counts", "cpm", "logcpm", or "log1p_cpm" to determine
#                 how the counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm) of non-zero cpm entries
#                   "log1p_cpm" will take the log2(cpm+1) of cpms
#
# Returns:
#   a SummarizedExperiment or SingleCellExperiment object that is the exact same
#   as what was read from the file, except that the "counts" slot is set with
#   the transformed values if applicable, or NULL if there is an error.
#
# NOTE: SingleCellExperiment is a subclass of SummarizedExperiment, so this
#       method works for both single cell and pseudobulk data.
# NOTE: Both objects allow for putting cpm/log2 values in other named slots to
#       preserve the counts slot, but some of these single cell datasets are so
#       large that it is impractical or impossible to have a raw counts matrix
#       *and* the transformed values held in memory at the same time. Instead
#       we overwrite the counts slot. To be consistent with single cell data and
#       allow for inter-operability, pseudobulk data is treated the same way.
Load_CountsFile <- function(filename, output_type) {
  if (!file.exists(filename)) {
    print(str_glue("Error: {filename} doesn't exist!"))
    return(NULL)
  }

  if (!(output_type %in% c("counts", "cpm", "logcpm", "log1p_cpm"))) {
    print(paste0("Error! output_type should be either \"counts\", \"cpm\",",
                 "\"logcpm\", or \"log1p_cpm\"."))
    return(NULL)
  }

  se_obj <- readRDS(filename)

  # TODO seaRef only, this loads it all into memory
  if (is(assay(se_obj, "counts"), "DelayedArray")) {
    if (!file.exists(path(assay(se_obj, "counts")))) {
      path(assay(se_obj, "counts")) <- file_searef_h5
    }
    assay(se_obj, "counts") <- as(assay(se_obj, "counts"), "CsparseMatrix")
  }

  if (output_type == "counts") {
    return(se_obj)
  }

  # The other three output_types all need CPM calculation
  cpms <- calculateCPM(se_obj)

  # Log2 of non-zero entries, no pseudocount
  if (output_type == "logcpm") {
    if (is(cpms, "matrix")) {
      cpms[cpms != 0] <- log2(cpms[cpms != 0])
    }
    else { # Sparse matrix
      cpms@x <- log2(cpms@x)
    }
  }

  # log2(x+1)
  if (output_type == "log1p_cpm") {
    if (is(cpms, "matrix")) {
      cpms <- log2(cpms + 1)
    }
    else { # Sparse matrix
      cpms@x <- log2(cpms@x+1)
    }
  }

  # Replace raw counts with 'cpms', which is either CPMs or log2 values
  assays(se_obj)[["counts"]] <- cpms
  return(se_obj)
}


##### Bulk (test) data #####

# Arguments:
#   dataset = the name of the dataset ("ROSMAP", "Mayo", or "MSBB")
#   output_type = either "counts", "cpm", or "logcpm", to determine how the
#                 counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm)
Load_BulkData <- function(dataset, output_type = "counts") {
  bulk_file <- str_glue(paste0("{dataset}_se.rds"))
  bulk_file <- file.path(dir_input, bulk_file)

  bulk <- Load_CountsFile(bulk_file, output_type)
  return(bulk)
}


##### Helper matrix load functions #####

# Load_AvgLibSize: Loads the average library size per cell type (called "A")
# that was calculated from single cell data.
#
# Arguments:
#   dataset = the name of the data set to load
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load the A vector for.
#
# Returns:
#   a named vector where the names are cell types and the values are the
#   normalized average library size of each cell type. Sums to 1.
Load_AvgLibSize <- function(dataset, granularity) {
  # This is a list with entries for "A_broad" and "A_fine"
  A <- readRDS(file.path(dir_input, str_glue("{dataset}_A_matrix.rds")))
  return(A[[str_glue("A_{granularity}")]])
}


# Load_SignatureMatrix: loads the cell-type signature matrix calculated from
# single cell data.
#
# Arguments:
#   dataset = the name of the data set to load
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load the signature matrix for.
#
# Returns:
#   a matrix with rows = genes and columns = cell types, where the values are
#   in counts per million, describing the expected value of each gene for each
#   cell type.
Load_SignatureMatrix <- function(dataset, granularity) {
  sig_matrix <- readRDS(file.path(dir_input, str_glue("{dataset}_signature.rds")))
  return(sig_matrix[[str_glue("sig_{granularity}")]])
}


# Load_GeneConversion: loads the data frame that contains mappings between
# Ensembl ID and gene symbol.
#
# Arguments:
#   dataset = the name of the data set to load the gene info from
#
# Returns:
#   a DataFrame with columns "Ensembl.ID" and "Symbol", used to convert
#   between gene naming systems
Load_GeneConversion <- function(dataset) {
  return(readRDS(file.path(dir_input, str_glue("{dataset}_gene_names.rds"))))
}


##### General algorithm I/O #####

# Save_AlgorithmIntermediate: saves a single output from an algorithm to
# a temporary folder, to safeguard against having to re-start processing on the
# full parameter set list in the event of a crash. The file is named with the
# parameters used to generate the result.
#
# Arguments:
#   result = a named list output by one of the deconvolution algorithms. Can
#            also be NULL instead of a list. If it's a list, it must contain an
#            entry called "params", which is a single-row dataframe containing
#            the parameters used to generate this result.
#   algorithm = the name of the algorithm to pre-pend to the file name
#
# Returns:
#   nothing
Save_AlgorithmIntermediate <- function(result, algorithm) {
  if (!is.null(result)) {
    filename <- paste(result$params, collapse = "_")
    saveRDS(result, file.path(dir_params_lists_tmp,
                              paste0(algorithm, "_", filename, ".rds")))
  }
}

# Save_AlgorithmOutputList: saves a list of output from one of the algorithms
# to an RDS file, named with a specific format.
#
# Arguments:
#   output_list = a list of outputs from one of the deconvolution algorithms,
#                 which contains output run under different parameter sets
#   algorithm = the name of the algorithm
#   dataset = the name of the data set
#   datatype = either "donors" or "training", if the algorithm was run on donor
#              or training pseudobulk, OR one of the bulk data sets ("Mayo",
#              "MSBB", "ROSMAP")
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#
# Returns:
#   Nothing
Save_AlgorithmOutputList <- function(output_list, algorithm, dataset, datatype, granularity) {
  list_file_format <- "{algorithm}_list_{dataset}_{datatype}_{granularity}.rds"

  out_directory <- switch(datatype,
                          "Mayo" = dir_mayo_output,
                          "MSBB" = dir_msbb_output,
                          "ROSMAP" = dir_rosmap_output,
                          dir_params_lists
  )

  saveRDS(output_list, file = file.path(out_directory,
                                        str_glue(list_file_format)))
}


# Load_AlgorithmOutputList: loads a list of output from one of the algorithms,
# named with a specific format.
#
# Arguments:
#   algorithm = the name of the algorithm
#   dataset = the name of the data set
#   datatype = either "donors" or "training", if the algorithm was run on donor
#              or training pseudobulk, OR one of the bulk data sets ("Mayo",
#              "MSBB", "ROSMAP")
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#
# Returns:
#   a list of outputs from one of the deconvolution algorithms, which contains
#   output run under different parameter sets
Load_AlgorithmOutputList <- function(algorithm, dataset, datatype, granularity) {
  list_file_format <- "{algorithm}_list_{dataset}_{datatype}_{granularity}.rds"

  out_directory <- switch(datatype,
                          "Mayo" = dir_mayo_output,
                          "MSBB" = dir_msbb_output,
                          "ROSMAP" = dir_rosmap_output,
                          dir_params_lists
  )

  params_file <- file.path(out_directory, str_glue(list_file_format))

  if (!file.exists(params_file)) {
    print(paste(params_file, "doesn't exist!"))
    return(list())
  }

  return(readRDS(params_file))
}


# Save_ErrorList: saves a list of calculated error metrics for a single
# algorithm for a single data set, named with a specific format.
#
# Arguments:
#   error_list = a list of errors containing entries for mean errors, errors
#                by celltype, errors by subject, goodness-of-fit, and parameters
#                for an algorithm / dataset / datatype / granularity combo
#   algorithm = the name of the algorithm
#   dataset = the name of the data set
#   datatype = either "donors" or "training", to signify if the algorithm was
#              run on donor or training pseudobulk
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#
# Returns:
#   Nothing
Save_ErrorList <- function(error_list, algorithm, dataset, datatype, granularity) {
  error_file_format <- "errors_{algorithm}_{dataset}_{datatype}_{granularity}.rds"
  saveRDS(error_list, file = file.path(dir_errors, str_glue(error_file_format)))
}


# Load_ErrorList: loads a list of calculated error metrics for a single
# algorithm for a single data set.
#
# Arguments:
#   algorithm = the name of the algorithm
#   dataset = the name of the data set
#   datatype = either "donors" or "training", to signify if the algorithm was
#              run on donor or training pseudobulk
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#
# Returns:
#   a list of errors containing entries for mean errors, errors by celltype,
#   errors by subject, goodness-of-fit, and parameters for an algorithm /
#   dataset / datatype / granularity combo
Load_ErrorList <- function(algorithm, dataset, datatype, granularity) {
  error_file_format <- "errors_{algorithm}_{dataset}_{datatype}_{granularity}.rds"
  error_file <- file.path(dir_errors, str_glue(error_file_format))

  if (!file.exists(error_file)) {
    print(paste(error_file, "doesn't exist!"))
    return(list())
  }

  return(readRDS(error_file))
}


##### Dtangle/HSPE #####

# Save_DtangleMarkers: saves a set of Dtangle/HSPE markers to an RDS file, named
# with a specific format.
#
# Arguments:
#   markers = the output of hspe::find_markers, which is a list.
#   dataset = the name of the data set
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used to generate the markers.
#   input_type = either "singlecell" or "pseudobulk", for which type of input
#                was used to generate the markers.
#   marker_method = one of "ratio", "diff", "p.value", or "regression", which
#                   was used as the input to find_markers.
#
# Returns:
#   Nothing
Save_DtangleMarkers <- function(markers, dataset, granularity, input_type, marker_method) {
  marker_file_format <- paste0("dtangle_markers_{dataset}_{granularity}_",
                               "input_{input_type}_method_{marker_method}.rds")
  saveRDS(markers, file = file.path(dir_markers, str_glue(marker_file_format)))
}


# Load_Markers: a generic function that reads a set of markers from an RDS file,
# named with a specific format. Markers can be from dtangle/HSPE, AutogeneS, or
# Seurat.
#
# Arguments:
#   dataset = the name of the data set
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used to generate the markers.
#   marker_type = one of "dtangle", "autogenes", or "seurat" to indicate which
#                 markers to load
#   marker_subtype = the subtype of markers to load, specific to the marker_type.
#                    For dtangle, one of "ratio", "diff", "p.value", or
#                    "regression", which was used as the input to
#                    dtangle::find_markers.
#                    For autogenes, one of "correlation", "distance", or
#                    "combined", specifying which weighting scheme was used to
#                    pick markers.
#                    Seurat doesn't have a subtype so this can be left as NULL.
#   input_type = for marker_type == "dtangle" only, either "singlecell" or
#                "pseudobulk", for which type of input was used to generate the
#                markers. For other marker_types, leave as NULL.
#
# Returns:
#   a named list where each entry is a list of marker gene names for a given
#   cell type
Load_Markers <- function(dataset, granularity, marker_type, marker_subtype = NULL,
                         input_type = NULL) {
  if (!(marker_type %in% c("dtangle", "autogenes", "seurat"))) {
    print("marker_type must be one of 'dtangle', 'autogenes', or 'seurat')")
    return(NULL)
  }

  marker_file_format <- NULL
  if (marker_type == "dtangle") {
    marker_file_format <- paste0("dtangle_markers_{dataset}_{granularity}_",
                                 "input_{input_type}_method_{marker_subtype}.rds")
  }
  else if (marker_type == "autogenes") {
    marker_file_format <- "autogenes_markers_{dataset}_{granularity}_{marker_subtype}.rds"
  }
  else if (marker_type == "seurat") {
    marker_file_format <- "seurat_markers_{dataset}_{granularity}.rds"
  }

  marker_file <- file.path(dir_markers, str_glue(marker_file_format))
  if (!file.exists(marker_file)) {
    print(str_glue("Marker file {marker_file} doesn't exist!"))
    return(NULL)
  }

  markers <- readRDS(file = marker_file)

  if (marker_type == "dtangle") {
    markers <- markers$filtered #lapply(markers$L, names)
  }

  return(markers)
}
