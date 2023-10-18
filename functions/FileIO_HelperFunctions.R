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
#   "celltype", which is populated with either the broad or subclass cell type
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

# Save_SingleCell: Saves a SingleCellExperiment object to the input folder. This
# object should be finalized data that has passed QC and has all cell types
# mapped to a reference.
#
# Arguments:
#   dataset = the name of the dataset
#   sce = a SingleCellExperiment object
#
# Returns:
#   nothing
Save_SingleCell <- function(dataset, sce) {
  sc_file <- file.path(dir_input, str_glue("{dataset}_sce.rds"))
  saveRDS(sce, sc_file)
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
# transforms the counts according to output_type. Sample or training pseudobulk-
# specific wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   data_type = either "sc_samples" or "training", for which type of pseudobulk
#               to load in.
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load in.
#   output_type = one of "counts", "cpm", "tpm", "tmm",
#                 "log_cpm", "log_tpm", "log_tmm",
#                 "qn_cpm", "qn_tpm", "qn_tmm",
#                 "qn_log_cpm", "qn_log_tpm", or "qn_log_tmm". See
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
#   data_type = either "sc_samples" or "training", for which type of pseudobulk
#               to load in.
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load in.
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
#   output_type = one of "counts", "cpm", "tpm", "tmm",
#                 "log_cpm", "log_tpm", "log_tmm",
#                 "qn_cpm", "qn_tpm", "qn_tmm",
#                 "qn_log_cpm", "qn_log_tpm", or "qn_log_tmm":
#                   "counts" will return raw, unaltered counts
#                   "cpm" will normalize the counts to counts per million
#                   "tpm" will normalize the counts to transcripts per million
#                         for bulk data only. Single cell data will default to
#                         cpm since UMI data is similar to tpm already.
#                   "tmm" will normalize using TMM factors from edgeR
#                   "log_cpm" will take the log2(cpm+1)
#                   "log_tpm" will take the log2(tpm+1)
#                   "log_tmm" will take the log2(tmm+1)
#                   "qn_cpm" will quantile normalize cpms
#                   "qn_tpm" will quantile normalize tpms
#                   "qn_tmm" will quantile normalize tmms
#                   "qn_log_cpm" will quantile normalize log_cpm values
#                   "qn_log_tpm" will quantile normalize log_tpm values
#                   "qn_log_tmm" will quantile normalize log_tmm values
#                 Quantile normalization is done by diagnosis + tissue (bulk) or
#                 diagnosis + celltype (single cell), where the expression data
#                 is split by those categories, quantile normalized within each
#                 category, and re-combined.
#   regression_method = which type of regressed/corrected counts to use, or "none".
#                 The corrected counts (or raw counts if "None") will then be
#                 transformed as usual according to `output_type`. Valid values:
#                   "edger" - regression via edgeR::glmQLFit (bulk only)
#                   "deseq2" - regression via DESeq2::DESeq (bulk only)
#                   "dream" - regression via voom/dream (bulk only)
#                   "none" - use raw counts
#                 NOTE: I don't know how well tmm or tpm will work with regressed
#                 data.
#
# Returns:
#   a SummarizedExperiment or SingleCellExperiment object that is the exact same
#   as what was read from the file, except that the "counts" slot is set with
#   the transformed values if applicable, or NULL if there is an error.
#
# NOTE: SingleCellExperiment is a subclass of SummarizedExperiment, so this
#       method works for both single cell and bulk data.
# NOTE: Both objects allow for putting cpm/log2 values in other named slots to
#       preserve the counts slot, but some of these single cell datasets are so
#       large that it is impractical or impossible to have a raw counts matrix
#       *and* the transformed values held in memory at the same time. Instead
#       we overwrite the counts slot. To be consistent with single cell data and
#       allow for inter-operability, bulk data is treated the same way.
Load_CountsFile <- function(filename, output_type, regression_method = "none") {
  if (!file.exists(filename)) {
    print(str_glue("Error: {filename} doesn't exist!"))
    return(NULL)
  }

  output_opts <- expand.grid(norm = c("cpm", "tpm", "tmm"),
                             log = c("", "log_"),
                             qn = c("", "qn_"))
  output_opts <- c("counts",
                   paste0(output_opts$qn, output_opts$log, output_opts$norm))

  if (!(output_type %in% output_opts)) {
    stop(paste0("Error! 'output_type' should be one of ", output_opts, "."))
  }

  if (!(regression_method %in% c("none", "edger", "deseq2", "dream"))) {
    stop("Error! 'regression_method' should be one of 'none', 'edger', 'deseq2', or 'dream'.")
  }

  se_obj <- readRDS(filename)

  # Applies to seaRef only, this loads it all into memory.
  if (is(assay(se_obj, "counts"), "DelayedArray")) {
    if (!file.exists(path(assay(se_obj, "counts")))) {
      path(assay(se_obj, "counts")) <- file_searef_h5
    }
    assay(se_obj, "counts") <- as(assay(se_obj, "counts"), "CsparseMatrix")
  }

  # If using regressed counts, replace the 'counts' matrix with those counts,
  # and replace the tmm_factors field with the TMM factors for the regressed
  # data
  if (regression_method != "none") {
    assay(se_obj, "counts") <- assay(se_obj, paste0("corrected_", regression_method))
    se_obj$tmm_factors <- colData(se_obj)[, paste0("tmm_factors_", regression_method)]
  }

  if (output_type == "counts") {
    return(se_obj)
  }

  # The remaining output_types all need CPM, TPM, or TMM

  # TPM for bulk data only
  if (grepl("tpm", output_type) & "exon_length" %in% colnames(rowData(se_obj))) {
    norm_counts <- calculateTPM(se_obj, lengths = rowData(se_obj)$exon_length)
  }
  # CPM for single cell, and for bulk if TPM is not specified
  else {
    norm_counts <- calculateCPM(se_obj)
  }

  if (grepl("tmm", output_type)) {
    # This is equivalent to counts * 1e6 / (libSize * tmm).
    # Since norm_counts is already (counts * 1e6 / libSize), calling this
    # function just divides by tmm, while preserving sparse matrices.
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

  # Quantile normalize by class. For bulk data this is by diagnosis and tissue.
  # For single cell data this is by cell type and diagnosis.
  if (grepl("qn", output_type)) {
    metadata <- colData(se_obj)
    if ("broadcelltype" %in% colnames(colData(se_obj))) {
      metadata$group <- paste(metadata$broadcelltype, metadata$diagnosis)
    }
    else {
      metadata$group <- paste(metadata$diagnosis, metadata$tissue)
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
#   regression_method = "none", if raw uncorrected counts should be used for
#                       bulk data, or one of "edger", "deseq2", or "dream", to
#                       use batch-corrected counts from one of those methods.
#
# Returns:
#   a SummarizedExperiment object with (possibly normalized) data in the
#   "counts" slot
Load_BulkData <- function(dataset, output_type = "counts", regression_method = "none") {
  bulk_file <- str_glue(paste0("{dataset}_se.rds"))
  bulk_file <- file.path(dir_input, bulk_file)

  bulk <- Load_CountsFile(bulk_file, output_type, regression_method)
  return(bulk)
}

# Arguments:
#   dataset = the name of the dataset ("ROSMAP", "Mayo", or "MSBB")
#   se = a SummarizedExperiment object
#
# Returns:
#   nothing
Save_BulkData <- function(dataset, se) {
  bulk_file <- str_glue(paste0("{dataset}_se.rds"))
  bulk_file <- file.path(dir_input, bulk_file)

  saveRDS(se, file = bulk_file)
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
#   output_type = any output type from the 'output_type' argument of
#                 Load_CountsFile. If output_type contains 'tmm', the signature
#                 matrix will be normalized with TMM factors. Otherwise it is
#                 in CPM. If output_type contains 'log', the signature matrix
#                 will be transformed as log2(x+1).
#
# Returns:
#   a matrix with rows = genes and columns = cell types, where the values are
#   in counts per million (or TMM-normalized), describing the expected value of
#   each gene for eachcell type.
Load_SignatureMatrix <- function(dataset, granularity, output_type) {
  sig_matrix <- readRDS(file.path(dir_input, str_glue("{dataset}_signature.rds")))

  if (grepl("tmm", output_type)) {
    sig <- sig_matrix$tmm[[granularity]]
  }
  else {
    sig <- sig_matrix$cpm[[granularity]]
  }

  if (grepl("log", output_type)) {
    sig <- log2(sig + 1)
  }
  return(sig)
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


# Save_ModelFormulas: saves a list of formulas for linear modeling to a file
# with a specific filename format.
#
# Arguments:
#   dataset_name = the name of the data set
#   formula_list = a list of formulas ('formula_fixed' and 'formula_mixed'). The
#                  formulas should be strings, not the 'formula' object.
#
# Returns:
#   nothing
Save_ModelFormulas <- function(dataset_name, formula_list) {
  saveRDS(formula_list, file.path(dir_metadata,
                                  str_glue("model_formulas_{dataset_name}.rds")))
}

# Load_ModelFormulas: loads a list of formulas for linear modeling from a file
# with a specific filename format.
#
# Arguments:
#   dataset_name = the name of the data set
#
# Returns:
#   a list of formulas ('formula_fixed' and 'formula_mixed'). The formulas are
#   strings, not the 'formula' object.
Load_ModelFormulas <- function(dataset_name) {
  readRDS(file.path(dir_metadata, str_glue("model_formulas_{dataset_name}.rds")))
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
Load_PreprocessedData <- function(dataset_name, remove_excluded = TRUE) {
  data <- readRDS(file.path(dir_preprocessed,
                            str_glue("{dataset_name}_preprocessed.rds")))

  if (remove_excluded) {
    # Remove mitochondrial and non-coding genes, which are marked as 'exclude'
    genes <- rowData(data)
    data <- data[!genes$exclude, ]
  }

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

# Load_AlgorithmIntermediate: loads a single output from an algorithm from the
# temporary folder, in order to re-load previous results when restarting due to
# crash or shutdown. The file is named with the parameters used to generate the
# result.
#
# Arguments:
#   algorithm = the name of the algorithm to pre-pend to the file name
#   params = a named vector or single-row dataframe with the parameters used to
#            generate the result
#
# Returns:
#   a result object from one of the algorithms if a file matching the input
#   parameters exists, or NULL if not
Load_AlgorithmIntermediate <- function(algorithm, params) {
  filename <- paste(params, collapse = "_")
  filename <- file.path(dir_params_lists_tmp,
                        paste0(paste0(algorithm, "_", filename, ".rds")))
  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  return(NULL)
}

# Save_AlgorithmOutputList: saves a list of output from one of the algorithms
# to an RDS file, named with a specific format.
#
# Arguments:
#   output_list = a list of outputs from one of the deconvolution algorithms,
#                 which contains output run under different parameter sets
#   algorithm = the name of the algorithm
#   test_dataset = the name of the test dataset
#   name_base = a string that uniquely identifies this list of outputs, which
#               currently is of the format:
#               '<reference_data_name>_<test_data_name>_<granularity>_<reference_input_type>_<normalization>_<regression_method>'
#
# Returns:
#   Nothing
Save_AlgorithmOutputList <- function(output_list, algorithm, test_dataset, name_base) {
  list_file_format <- paste0("estimates_{algorithm}_{name_base}.rds")

  out_directory <- switch(test_dataset,
                          "Mayo" = dir_mayo_output,
                          "MSBB" = dir_msbb_output,
                          "ROSMAP" = dir_rosmap_output,
                          dir_params_lists)

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
#   regression_method = the type of regression used for bulk counts. Same as
#                       the 'regression_method' argument to Load_CountsFile.
#
# Returns:
#   a list of outputs from one of the deconvolution algorithms, which contains
#   output run under different parameter sets
Load_AlgorithmOutputList <- function(algorithm, reference_dataset, test_dataset,
                                     granularity, normalization,
                                     regression_method = "none") {
  list_file_format <- paste0("estimates_{algorithm}_{reference_dataset}_",
                             "{test_dataset}_{granularity}_{normalization}_",
                             "{regression_method}.rds")

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
#   dataset_name = the name of the bulk dataset for these errors
#   error_list = a list of errors containing entries for mean errors, errors
#                by celltype, errors by subject, goodness-of-fit, and parameters
#                for an algorithm / dataset / datatype / granularity combo
#   algorithm = the name of the algorithm
#   params_data = a one-row dataframe or named list of parameters that were
#                 used to generate the estimates for this list
#
# Returns:
#   Nothing
Save_ErrorList <- function(dataset_name, error_list, algorithm, params_data) {
  name_base <- paste(params_data, collapse = "_")
  error_file_format <- paste0("errors_{algorithm}_{name_base}.rds")

  dir_errors_alg <- file.path(dir_errors, dataset_name, algorithm)
  dir.create(dir_errors_alg, recursive = TRUE, showWarnings = FALSE)

  saveRDS(error_list, file = file.path(dir_errors_alg, str_glue(error_file_format)))
}


# Load_ErrorList: loads a list of calculated error metrics for a single
# algorithm for a single data set.
#
# Arguments:
#   algorithm = the name of the algorithm
#   params = a one-row dataframe or named list of parameters used to generate
#            the file
#
# Returns:
#   a list of errors containing entries for mean errors, errors by celltype,
#   errors by subject, goodness-of-fit, and parameters for an algorithm /
#   dataset / datatype / granularity / normalization combo
Load_ErrorList <- function(algorithm, params) {
  name_base <- paste(params, collapse = "_")
  error_file_format <- paste0("errors_{algorithm}_{name_base}.rds")
  error_file <- file.path(dir_errors, params$test_data_name, algorithm,
                          str_glue(error_file_format))

  if (!file.exists(error_file)) {
    #print(paste(error_file, "doesn't exist!"))
    return(NULL)
  }

  return(readRDS(error_file))
}


# Save_ErrorIntermediate: saves the error calculation from a single param set to
# a temporary folder, to safeguard against having to re-start processing on the
# full parameter set list in the event of a crash. The file is named with the
# parameters used to generate the result.
#
# Arguments:
#   error_obj = a named list containing both the output of one of the
#            deconvolution algorithms plus the error calculations. It must
#            contain an entry called "params", which is a single-row dataframe
#            containing the parameters used to generate this result.
#   algorithm = the name of the algorithm to pre-pend to the file name
#
# Returns:
#   nothing
Save_ErrorIntermediate <- function(error_obj, algorithm) {
  params <- error_obj$params %>% select(-total_markers_used)
  filename <- paste(params, collapse = "_")
  saveRDS(error_obj, file.path(dir_errors_tmp, str_glue("{algorithm}_{filename}.rds")))
}


# Load_ErrorIntermediate: loads an error calculation from a single param set
# from a temporary folder, in order to re-load previous calculations when
# restarting due to crash or shutdown. The file is named with the parameters
# used to generate the result.
#
# Arguments:
#   algorithm = the name of the algorithm to pre-pend to the file name
#   params = a named vector or single-row dataframe with the parameters used to
#            generate the result
#
# Returns:
#   a named list containing both the output of one of the deconvolution
#   algorithms plus the error calculations if a file matching the input
#   parameters exists, or NULL if not
Load_ErrorIntermediate <- function(algorithm, params) {
  filename <- paste(params, collapse = "_")
  filename <- file.path(dir_errors_tmp, str_glue("{algorithm}_{filename}.rds"))

  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  return(NULL)
}


##### Dtangle/HSPE #####

# Save_DtangleMarkers: saves a set of Dtangle/HSPE markers to an RDS file, named
# with a specific format.
#
# Arguments:
#   markers = the output of hspe::find_markers, which is a list.
#   dataset = the name of the data set
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types was used to generate the markers.
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
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types was used to generate the markers.
#   marker_type = one of "dtangle", "autogenes", "seurat", or "deseq2", to
#                 indicate which markers to load
#   marker_subtype = the subtype of markers to load, specific to the marker_type.
#                    For dtangle, one of "ratio", "diff", "p.value", or
#                    "regression", which was used as the input to
#                    dtangle::find_markers.
#                    For autogenes, one of "correlation", "distance", or
#                    "combined", specifying which weighting scheme was used to
#                    pick markers.
#                    Seurat doesn't have a subtype so this can be left as NULL.
#                    For deseq2, one of "DESeq2" or "glmGamPoi".
#   input_type = for marker_type == "dtangle" only, either "singlecell" or
#                "pseudobulk", for which type of input was used to generate the
#                markers. For other marker_types, leave as NULL.
#
# Returns:
#   a named list where each entry is a list of marker gene names for a given
#   cell type
Load_Markers <- function(dataset, granularity, marker_type, marker_subtype = NULL,
                         input_type = NULL) {
  if (!(marker_type %in% c("dtangle", "autogenes", "seurat", "deseq2"))) {
    print("marker_type must be one of 'dtangle', 'autogenes', 'seurat', or 'deseq2')")
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
  else if (marker_type == "deseq2") {
    marker_file_format <- "deseq2_markers_{dataset}_{granularity}_{marker_subtype}.rds"
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
