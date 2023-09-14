# Helper functions for loading or saving data sets and markers.

library(SummarizedExperiment)
library(SingleCellExperiment)
library(scuttle)
library(stringr)
library(DESeq2)
library(preprocessCore)
source("Filenames.R")

##### Single cell / Pseudobulk objects #####

# Load_SingleCell: reads a SingleCellExperiment data set from a file and
# transforms the counts according to output_type. Single-cell-specific wrapper
# function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to use in the metadata
#   output_type = one of "counts", "vst", "cpm", "tmm", "log_cpm", "log_tmm",
#                 "qn_cpm", "qn_tmm", "qn_log_cpm", or "qn_log_tmm". See
#                 Load_CountsFile for description.
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

  if (!(granularity %in% c("broad_class", "sub_class"))) {
    stop("Error! 'granularity' should be either 'broad_class' or 'sub_class'.")
  }

  # Assign the column "celltype" to be either the broad or fine cell types
  metadata$celltype <- metadata[,granularity]

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
#   output_type = one of "counts", "vst", "cpm", "tmm", "log_cpm", "log_tmm",
#                 "qn_cpm", "qn_tmm", "qn_log_cpm", or "qn_log_tmm". See
#                 Load_CountsFile for description.
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable.
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
#   output_type = one of "counts", "vst", "cpm", "tmm", "log_cpm", "log_tmm",
#                 "qn_cpm", "qn_tmm", "qn_log_cpm", or "qn_log_tmm". See
#                 Load_CountsFile for description.
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
#   output_type = one of "counts", "vst", "cpm", "tmm", "log_cpm", "log_tmm",
#                 "qn_cpm", "qn_tmm", "qn_log_cpm", or "qn_log_tmm":
#                   "counts" will return raw, unaltered counts
#                   "vst" will return the variance stabilized transform of the
#                         counts, using DESeq2::vst()
#                   "cpm" will normalize the counts to counts per million
#                   "tmm" will normalize using TMM factors from edgeR
#                   "log_cpm" will take the log2(cpm+1)
#                   "log_tmm" will take the log2(cpm+1)
#                   "qn_cpm" will quantile normalize cpms
#                   "qn_tmm" will quantile normalize tmms
#                   "qn_log_cpm" will quantile normalize log_cpm values
#                   "qn_log_tmm" will quantile normalize log_tmm values
#                 Quantile normalization is done by diagnosis (bulk) or
#                 diagnosis + celltype (single cell), where the expression data
#                 is split by those categories, quantile normalized within each
#                 category, and re-combined.
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

  output_opts <- expand.grid(norm = c("cpm", "tmm"),
                             log = c("", "log_"),
                             qn = c("", "qn_"))
  output_opts <- c("counts", "vst",
                   paste0(output_opts$qn, output_opts$log, output_opts$norm))

  if (!(output_type %in% output_opts)) {
    print(paste0("Error! output_type should be one of 'counts', 'vst', 'cpm', ",
                 "'tmm', 'log_cpm', 'log_tmm', 'qn_cpm', 'qn_tmm', ",
                 "'qn_log_cpm', or 'qn_log_tmm'."))
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

  # VST returns log2-scale values
  # TODO this doesn't work on single cell data
  if (output_type == "vst") {
    norm_counts <- DESeqDataSetFromMatrix(assay(se_obj, "counts"),
                                          colData(se_obj),
                                          design = ~1)
    norm_counts <- vst(norm_counts)
    assays(se_obj)[["counts"]] <- assay(norm_counts)
    return(se_obj)
  }

  # The remaining output_types all need CPM or TMM
  norm_counts <- calculateCPM(se_obj)

  if (grepl("tmm", output_type)) {
    # This is equivalent to counts * 1e6 / (libSize * tmm).
    # Since norm_counts is already (counts * 1e6 / libSize), calling this
    # function just divides by tmm, preserving sparse matrices.
    norm_counts <- normalizeCounts(norm_counts, se_obj$tmm_factors,
                                   log = FALSE, center.size.factors = FALSE)
  }

  # log2(x+1)
  if (grepl("log", output_type)) {
    if (is(norm_counts, "matrix")) {
      norm_counts <- log2(norm_counts + 1)
    }
    else { # Sparse matrix
      norm_counts@x <- log2(norm_counts@x+1)
    }
  }

  # Quantile normalize by class. For bulk data this is by diagnosis. For
  # single cell data this is by cell type and diagnosis.
  if (grepl("qn", output_type)) {
    metadata <- colData(se_obj)
    if ("broadcelltype" %in% colnames(colData(se_obj))) {
      metadata$group <- paste(metadata$broadcelltype, metadata$diagnosis)
    }
    else {
      metadata$group <- metadata$diagnosis
    }

    for (G in unique(metadata$group)) {
      sub <- subset(metadata, group == G)
      q <- normalize.quantiles(as.matrix(norm_counts[,rownames(sub)]),
                               copy = FALSE,
                               keep.names = TRUE)
      #q <- limma::normalizeQuantiles(norm_counts[,rownames(sub)])
      if (is(norm_counts, "matrix")) {
        norm_counts[,rownames(sub)] <- q
      }
      else {
        norm_counts[,rownames(sub)] <- as(q, "CsparseMatrix")
      }
    }
  }

  # Replace raw counts with the normalized counts
  assays(se_obj)[["counts"]] <- norm_counts
  return(se_obj)
}


##### Bulk (test) data #####

# Arguments:
#   dataset = the name of the dataset ("ROSMAP", "Mayo", or "MSBB")
#   output_type = one of "counts", "vst", "cpm", "tmm", "log_cpm", "log_tmm",
#                 "qn_cpm", "qn_tmm", "qn_log_cpm", or "qn_log_tmm". See
#                 Load_CountsFile for description.
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


# Save_Covariates: saves a dataframe of covariates to a file with a specific
# filename format.
#
# Arguments:
#   dataset_name = the name of the data set
#   covariates = a data frame of covariates
#
# Returns:
#   nothing
Save_Covariates <- function(dataset_name, covariates) {
  saveRDS(covariates, file.path(dir_covariates,
                                str_glue("{dataset_name}_covariates.rds")))
}

# Load_Covariates: loads a dataframe of covariates from a file with a specific
# filename format.
#
# Arguments:
#   dataset_name = the name of the data set
#
# Returns:
#   a data frame of covariates
Load_Covariates <- function(dataset_name) {
  return(readRDS(file.path(dir_covariates,
                           str_glue("{dataset_name}_covariates.rds"))))
}


# Save_MapReference: saves a reference data set for cell type mapping to a
# file with a specific filename format.
#
# Arguments:
#   dataset_name = the name of the single cell dataset
#   seurat_object = Either a single Seurat object (broad_class) or a list of
#                   Seurat objects (sub_class). Seurat object(s) should be
#                   SCTransform-ed, and have PCA and UMAP dimension reductions.
#                   Run_UMAP needs to be run with 'return.model = TRUE' in order
#                   for this object to work as a reference.
#   granularity = either 'broad_class' or 'sub_class', for whether the
#                 Seurat object(s) have been prepared for broad class mapping or
#                 sub class mapping.
#
# Returns:
#   Nothing
Save_MapReference <- function(dataset_name, seurat_object, granularity) {
  saveRDS(seurat_object,
          file.path(dir_map_reference,
                    str_glue("reference_{dataset_name}_{granularity}.rds")))
}


# Load_MapReference: loads a reference data set for cell type mapping from a
# file with a specific filename format.
#
# Arguments:
#   dataset_name = the name of the single cell dataset
#   granularity = either 'broad_class' or 'sub_class', for whether the
#                 Seurat object(s) have been prepared for broad class mapping or
#                 sub class mapping.
#
# Returns:
#   a single Seurat object (if granularity = "broad_class") or a list of
#   Seurat objects (if granularity = "sub_class")
Load_MapReference <- function(dataset_name, granularity) {
  readRDS(file.path(dir_map_reference,
                    str_glue("reference_{dataset_name}_{granularity}.rds")))
}


# Save_PreprocessedData: Saves single cell or bulk RNA seq data that has been
# pre-processed (run through Step02) but has not had cell types mapped or
# normalization applied.
#
# Arguments:
#   dataset_name = the name of the dataset to load
#   se_object = a SummarizedExperiment or SingleCellExperiment object
#
# Returns:
#   Nothing
Save_PreprocessedData <- function(dataset_name, se_object) {
  saveRDS(se_object,
          file = file.path(dir_preprocessed,
                           str_glue("{dataset_name}_preprocessed.rds")))
}


# Load_PreprocessedData: Loads single cell or bulk RNA seq data that has been
# pre-processed (run through Step02) but has not had cell types mapped or
# normalization applied.
#
# Arguments:
#   dataset_name = the name of the dataset to load
#
# Returns:
#   A SummarizedExperiment (or SingleCellExperiment)
Load_PreprocessedData <- function(dataset_name) {
  data <- readRDS(file.path(dir_preprocessed,
                            str_glue("{dataset_name}_preprocessed.rds")))
  if (is(assay(data, "counts"), "DelayedMatrix")) {
    assay(data, "counts") <- as(assay(data, "counts"), "CsparseMatrix")
  }
  return(data)
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
#   reference_dataset = the name of the reference data set
#   test_dataset = either "donors" or "training", if the algorithm was run on
#                  donor or training pseudobulk, OR one of the bulk data sets
#                  ("Mayo", "MSBB", "ROSMAP")
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#   normalization = the normalization strategy used. Same as the 'output_type'
#                   argument to Load_CountsFile.
#
# Returns:
#   Nothing
Save_AlgorithmOutputList <- function(output_list, algorithm, reference_dataset,
                                     test_dataset, granularity, normalization) {
  list_file_format <- paste0("estimates_{algorithm}_{reference_dataset}_",
                             "{test_dataset}_{granularity}_{normalization}.rds")

  out_directory <- switch(test_dataset,
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
#   reference_dataset = the name of the reference data set
#   test_dataset = either "donors" or "training", if the algorithm was run on
#                  donor or training pseudobulk, OR one of the bulk data sets
#                  ("Mayo", "MSBB", "ROSMAP")
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#   normalization = the normalization strategy used. Same as the 'output_type'
#                   argument to Load_CountsFile.
#
# Returns:
#   a list of outputs from one of the deconvolution algorithms, which contains
#   output run under different parameter sets
Load_AlgorithmOutputList <- function(algorithm, reference_dataset, test_dataset,
                                     granularity, normalization) {
  list_file_format <- paste0("estimates_{algorithm}_{reference_dataset}_",
                             "{test_dataset}_{granularity}_{normalization}.rds")

  out_directory <- switch(test_dataset,
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
#   reference_dataset = the name of the reference data set
#   test_dataset = either "donors" or "training", to signify if the algorithm was
#                  run on donor or training pseudobulk, or one of "Mayo",
#                  "MSBB", or "ROSMAP" if the algorithm was run on bulk data
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#   normalization = the type of normalization used. See Load_CountsFile
#                   output_type parameter explanation for valid values.
#
# Returns:
#   Nothing
Save_ErrorList <- function(error_list, algorithm, reference_dataset, test_dataset,
                           granularity, normalization) {
  error_file_format <- paste0("errors_{algorithm}_{reference_dataset}_",
                              "{test_dataset}_{granularity}_{normalization}.rds")
  saveRDS(error_list, file = file.path(dir_errors, str_glue(error_file_format)))
}


# Load_ErrorList: loads a list of calculated error metrics for a single
# algorithm for a single data set.
#
# Arguments:
#   algorithm = the name of the algorithm
#   reference_dataset = the name of the reference data set
#   test_dataset = either "donors" or "training", to signify if the algorithm was
#                  run on donor or training pseudobulk, or one of "Mayo",
#                  "MSBB", or "ROSMAP" if the algorithm was run on bulk data
#   granularity = either "broad" or "fine", for which level of cell types was
#                 used for markers and pseudobulk creation.
#   normalization = the type of normalization used. See Load_CountsFile
#                   output_type parameter explanation for valid values.
#
# Returns:
#   a list of errors containing entries for mean errors, errors by celltype,
#   errors by subject, goodness-of-fit, and parameters for an algorithm /
#   dataset / datatype / granularity / normalization combo
Load_ErrorList <- function(algorithm, reference_dataset, test_dataset,
                           granularity, normalization) {
  error_file_format <- paste0("errors_{algorithm}_{reference_dataset}_",
                              "{test_dataset}_{granularity}_{normalization}.rds")
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
  markers <- markers$filtered # All 3 marker types have a 'filtered' list

  return(markers)
}
