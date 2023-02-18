# Helper functions for loading or saving data sets and markers.

library(SummarizedExperiment)
library(SingleCellExperiment)
library(scuttle)
library(stringr)

##### Single cell / Pseudobulk objects #####

# Load_SingleCell: reads a SingleCellExperiment data set from a file and
# transforms the counts according to output_type. Single-cell-specific wrapper
# function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad" or "fine", for which level of cell types to
#                 use in the metadata
#   output_type = either "counts", "cpm", or "logcpm", to determine how the
#                 counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm + 1)
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
#   output_type = either "counts", "cpm", or "logcpm", to determine how the
#                 counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm + 1)
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable. It also gets colData set to a DataFrame containing one column
#   for "celltype", which is populated with cell type assignments, which can be
#   extracted from the sample names in the data set.
Load_PseudobulkPureSamples <- function(dataset, granularity, output_type = "counts") {
  pb_file <- str_glue(paste0("pseudobulk_{dataset}_puresamplesbydonor_",
                             "{granularity}celltypes.rds"))
  pb_file <- file.path(dir_pseudobulk, pb_file)

  pseudobulk <- Load_CountsFile(pb_file, output_type)

  # TODO make metadata for these files in generate_pseudobulk
  celltypes <- str_replace(colnames(pseudobulk), "puresample_", "")
  celltypes <- str_replace(celltypes, "_.*", "")
  celltypes <- factor(celltypes)

  metadata <- DataFrame(celltype = celltypes)
  rownames(metadata) <- colnames(pseudobulk)

  colData(pseudobulk) <- metadata
  return(pseudobulk)
}


# Load_Pseudobulk: reads a SummarizedExperiment data set from a file and
# transforms the counts according to output_type. Donor or training pseudobulk-
# specific wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad" or "fine", for which level of cell types to
#                 load in.
#   output_type = either "counts", "cpm", or "logcpm", to determine how the
#                 counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm + 1)
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable.
Load_Pseudobulk <- function(dataset, data_type, granularity, output_type = "counts") {
  pb_file <- str_glue(paste0("pseudobulk_{dataset}_{data_type}_",
                             "{granularity}celltypes.rds"))
  pb_file <- file.path(dir_pseudobulk, pb_file)

  pseudobulk <- Load_CountsFile(pb_file, output_type)
  return(pseudobulk)
}


# Load_CountsFile: reads a SummarizedExperiment or SingleCellExperiment from
# a file and transforms the counts according to output_type.
#
# Arguments:
#   filename = the name of of the file to read in, including file path
#   output_type = either "counts", "cpm", or "logcpm", to determine how the
#                 counts are transformed:
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "logcpm" will take the log2(cpm + 1)
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

  se_obj <- readRDS(filename)

  if (output_type == "counts") {
    return(se_obj)
  }
  else if (output_type == "cpm" | output_type == "logcpm") {
    cpms <- calculateCPM(se_obj)

    if (output_type == "logcpm") {
      # TODO consider using scuttle normalizeCounts()
      if (is(cpms, "DelayedArray") | is(cpms, "matrix")) {
        cpms <- log2(cpms + 1)
      }
      else { # Sparse matrix
        cpms@x <- log2(cpms@x + 1)
      }
    }

    # Replace raw counts with 'cpms', which is either CPMs or log2(CPMs)
    assays(se_obj)[["counts"]] <- cpms
    return(se_obj)
  }
  else {
    print("Error! output_type should be either \"counts\", \"cpm\", or \"logcpm\".")
    return(NULL)
  }
}


##### Dtangle/HSPE Markers #####

marker_file_format <- paste0("dtangle_markers_{dataset}_{granularity}_",
                             "input_{input_type}_method_{marker_method}.rds")

# Save_DtangleMarkers: saves a set of Dtangle/HSPE markers to an RDS file, named
# with a specific format. This is an ease-of-use function to prevent having to
# edit multiple files whenver the format of the filename changes.
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
  saveRDS(markers, file = file.path(dir_markers, str_glue(marker_file_format)))
}


# Load_DtangleMarkers: saves a set of Dtangle/HSPE markers to an RDS file, named
# with a specific format. This is an ease-of-use function to prevent having to
# edit multiple files whenver the format of the filename changes.
#
# Arguments:
#   dataset = the name of the data set
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used to generate the markers.
#   input_type = either "singlecell" or "pseudobulk", for which type of input
#                was used to generate the markers.
#   marker_method = one of "ratio", "diff", "p.value", or "regression", which
#                   was used as the input to find_markers.
#
# Returns:
#   a list that was the output of hspe::find_markers for this parameter set
Load_DtangleMarkers <- function(dataset, granularity, input_type, marker_method) {
  markers <- readRDS(file = file.path(dir_markers, str_glue(marker_file_format)))
  return(markers)
}
