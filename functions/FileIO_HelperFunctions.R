# Helper functions for loading or saving data sets and markers.

library(SummarizedExperiment)
library(SingleCellExperiment)
library(scuttle)
library(stringr)
library(DESeq2)
#library(preprocessCore)
library(vroom)
source("Filenames.R")


# Generic read functions -------------------------------------------------------

# Load_CountsFile: reads a SummarizedExperiment or SingleCellExperiment from
# a file and transforms the counts according to 'normalization'.
#
# Arguments:
#   filename = the name of of the file to read in, including file path
#   normalization = one of "counts", "counts_tpm", "cpm", "tpm", "tmm",
#                   "log_cpm", "log_tpm", or "log_tmm":
#                     "counts" will return raw, unaltered counts
#                     "counts_tpm" will return the TPM matrix for bulk data,
#                         scaled to the same order of magnitude as the counts
#                         matrix. Single cell data will default to "counts"
#                         since UMI data is similar to TPM already.
#                     "cpm" will normalize the counts to counts per million
#                     "tpm" will use the transcripts per million matrix for bulk
#                         data only, normalized to CPM. The extra CPM
#                         normalization is done because the TPM matrix values no
#                         longer add up to 1e6 due to the large number of genes
#                         removed from the data. Single cell data will default
#                         to "cpm".
#                     "tmm" will normalize using TMM factors from edgeR
#                     "log_cpm" will take the log2(cpm+1)
#                     "log_tpm" will take the log2(tpm+1)
#                     "log_tmm" will take the log2(tmm+1)
#   regression_method = which type of regressed/corrected counts to use, or "none".
#                   The corrected counts (or raw counts if "none") will then be
#                   transformed as usual according to "normalization". Valid values:
#                     "edger" - regression via edgeR::glmQLFit (bulk only)
#                     "lme" - regression via lme4::lmer or lm (bulk only)
#                     "combat" - regression via ComBat (bulk only)
#                     "none" - use raw counts
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
Load_CountsFile <- function(filename, normalization, regression_method = "none") {
  if (!file.exists(filename)) {
    message(str_glue("Error: {filename} doesn't exist!"))
    return(NULL)
  }

  norm_opts <- expand.grid(norm = c("cpm", "tpm", "tmm"),
                           log = c("", "log_"))
  norm_opts <- c("counts", "counts_tpm",
                 paste0(norm_opts$log, norm_opts$norm))

  if (!(normalization %in% norm_opts)) {
    stop(paste0("Error! 'normalization' should be one of ",
                paste(norm_opts, collapse = ", "),
                "."))
  }

  if (!(regression_method %in% c("none", "edger", "lme", "combat"))) {
    stop("Error! 'regression_method' should be one of 'none', 'edger', 'lme', or 'combat'.")
  }

  se_obj <- readRDS(filename)

  # Bulk data only -- if TPM is part of the normalization, replace the "counts"
  # assay with the TPM assay and remove tmm factors.
  if (grepl("tpm", normalization) && "tpm" %in% names(assays(se_obj))) {
    # The TPM matrices go to two decimal places, so we multiply by 100 to get
    # count-like data. This has the advantage of putting the TPM data on
    # roughly the same scale as the actual count data as well.
    assay(se_obj, "counts") <- assay(se_obj, "tpm") * 100
    se_obj$tmm_factors <- 1
  }

  # If using regressed counts, replace the 'counts' matrix with those counts,
  # and replace the tmm_factors field with the TMM factors for the regressed
  # data
  if (regression_method != "none") {
    field_name <- regression_method
    tmm_factors <- colData(se_obj)[, paste0("tmm_factors_", regression_method)]

    if (grepl("tpm", normalization)) {
      field_name <- str_glue("{regression_method}_tpm")
      tmm_factors <- 1 # TMM factors ignored for TPM data
    }

    assay(se_obj, "counts") <- assay(se_obj, paste0("corrected_", field_name))
    se_obj$tmm_factors <- tmm_factors
  }

  if (normalization == "counts" || normalization == "counts_tpm") {
    return(se_obj)
  }

  # The remaining normalizations all need normalization by CPM, TPM, or TMM

  # If normalization is TPM, the appropriate TPM data for the regression method
  # will be in the "counts" assay and get normalized to sum to 1e6. Otherwise,
  # this acts on raw counts.
  norm_counts <- scuttle::calculateCPM(se_obj)

  if (grepl("tmm", normalization)) {
    # This is equivalent to counts * 1e6 / (libSize * tmm).
    # Since norm_counts is already (counts * 1e6 / libSize), calling this
    # function just divides by tmm, while preserving sparse matrices.
    norm_counts <- scuttle::normalizeCounts(norm_counts, se_obj$tmm_factors,
                                            log = FALSE,
                                            center.size.factors = FALSE)
  }

  # log2(x+1)
  if (grepl("log", normalization)) {
    if (is(norm_counts, "matrix")) {
      norm_counts <- log2(norm_counts + 1)
    }
    else { # Sparse matrix
      norm_counts@x <- log2(norm_counts@x + 1)
    }
  }

  # Replace raw counts with the normalized counts
  assays(se_obj)[["counts"]] <- norm_counts
  return(se_obj)
}


# Single cell / Pseudobulk wrapper functions -----------------------------------

# Load_SingleCell: reads a SingleCellExperiment data set from a file and
# transforms the counts according to 'normalization'. Single-cell-specific
# wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to use in the metadata
#   normalization = one of "counts", "counts_tpm", "cpm", "tmm", "tpm",
#                   "log_cpm", "log_tmm", or "log_tpm". See Load_CountsFile for
#                   description.
#
# Returns:
#   a SingleCellExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable. The colData DataFrame also has a column added called
#   "celltype", which is populated with either the broad or subclass cell type
#   assignments based on granularity.
Load_SingleCell <- function(dataset, granularity, normalization = "counts") {
  sc_file <- file.path(dir_singlecell, str_glue("{dataset}_sce.rds"))

  singlecell <- Load_CountsFile(sc_file, normalization, regression_method = "none")
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
  sc_file <- file.path(dir_singlecell, str_glue("{dataset}_sce.rds"))
  saveRDS(sce, sc_file)
}


# Load_PseudobulkPureSamples: reads a SummarizedExperiment data set from a file
# and transforms the counts according to 'normalization'. Pure sample
# pseudobulk-specific wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load in.
#   normalization = one of "counts", "counts_tpm", "cpm", "tmm", "tpm",
#                   "log_cpm", "log_tmm", or "log_tpm". See Load_CountsFile for
#                   description.
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable.
Load_PseudobulkPureSamples <- function(dataset, granularity, normalization = "counts") {
  pb_file <- str_glue("pseudobulk_{dataset}_puresamples_{granularity}.rds")
  pb_file <- file.path(dir_pseudobulk, pb_file)

  pseudobulk <- Load_CountsFile(pb_file, normalization, regression_method = "none")
  return(pseudobulk)
}


# Save_PseudobulkPureSamples: saves a SummarizedExperiment object to a file with
# a specific filename format (which should match the format used in
# Load_PseudobulkPureSamples).
#
# Arguments:
#   se = a SummarizedExperiment object
#   dataset = the name of the data set to load in
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load in.
#
# Returns:
#   nothing
Save_PseudobulkPureSamples <- function(se, dataset, granularity) {
  filename <- str_glue("pseudobulk_{dataset}_puresamples_{granularity}.rds")
  saveRDS(se, file = file.path(dir_pseudobulk, filename))
}


# Load_Pseudobulk: reads a SummarizedExperiment data set from a file and
# transforms the counts according to 'normalization'. Sample or training
# pseudobulk-specific wrapper function for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the data set to load in
#   data_type = either "sc_samples" or "training", for which type of pseudobulk
#               to load in.
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load in.
#   normalization = one of "counts", "counts_tpm", "cpm", "tmm", "tpm",
#                   "log_cpm", "log_tmm", or "log_tpm". See Load_CountsFile for
#                   description.
#
# Returns:
#   a SummarizedExperiment object that is the exact same as what was read from
#   the file, except that the "counts" slot is set with the transformed values
#   if applicable.
Load_Pseudobulk <- function(dataset, data_type, granularity, normalization = "counts") {
  pb_file <- str_glue("pseudobulk_{dataset}_{data_type}_{granularity}.rds")
  pb_file <- file.path(dir_pseudobulk, pb_file)

  pseudobulk <- Load_CountsFile(pb_file, normalization, regression_method = "none")
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


Load_SimulatedScadenData <- function(dataset, granularity, normalization) {
  sim_file <- str_glue("simulated_{dataset}_{granularity}_{normalization}.rds")
  sim_file <- file.path(dir_scaden_pseudobulk, sim_file)

  # Always use the CPM file if normalization is "tpm", since TPM doesn't apply
  # to single cell data.
  sim_file <- str_replace(sim_file, "_tpm", "_cpm")

  # Simulated scaden data is always normalized to CPM. If the normalization was
  # TMM, the adjustment was already accounted for in the simulated data and it
  # just needs to be normalized to library size.
  sim_data <- Load_CountsFile(sim_file, normalization = "cpm", regression_method = "none")
  return(sim_data)
}

save_SimulatedScadenData <- function(se, dataset, granularity, normalization) {
  filename <- str_glue("simulated_{dataset}_{granularity}_{normalization}.rds")
  saveRDS(se, file = file.path(dir_scaden_pseudobulk, filename))
}


# Bulk wrapper functions -------------------------------------------------------

# Load_BulkData: reads a SummarizedExperiment object from a file and transforms
# it according to the output type and regression method. Bulk-specific wrapper
# for Load_CountsFile.
#
# Arguments:
#   dataset = the name of the dataset ("ROSMAP", "Mayo", or "MSBB")
#   normalization = one of "counts", "counts_tpm", "cpm", "tmm", "tpm",
#                   "log_cpm", "log_tmm", or "log_tpm". See Load_CountsFile for
#                   description.
#   regression_method = "none", if raw uncorrected counts should be used for
#                       bulk data, or one of "edger", "lme", or "combat", to
#                       use batch-corrected counts from one of those methods.
#
# Returns:
#   a SummarizedExperiment object with (possibly normalized) data in the
#   "counts" slot
Load_BulkData <- function(dataset, normalization = "counts", regression_method = "none") {
  bulk_file <- str_glue(paste0("{dataset}_se.rds"))
  bulk_file <- file.path(dir_bulk, bulk_file)

  bulk <- Load_CountsFile(bulk_file, normalization, regression_method)
  return(bulk)
}

# Save_BulkData: Save a SummarizedExperiment object to the input folder. This
# object should be finalized data that has passed QC and has slots for normalized
# data from edgeR, lme4, and combat
#
# Arguments:
#   dataset = the name of the dataset ("ROSMAP", "Mayo", or "MSBB")
#   se = a SummarizedExperiment object
#
# Returns:
#   nothing
Save_BulkData <- function(dataset, se) {
  bulk_file <- str_glue(paste0("{dataset}_se.rds"))
  bulk_file <- file.path(dir_bulk, bulk_file)

  saveRDS(se, file = bulk_file)
}


# Helper matrix load/save functions --------------------------------------------

# Load_AvgLibSize: Loads the average library size per cell type (called "A")
# that was calculated from single cell data.
#
# Arguments:
#   dataset = the name of the data set to load
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load the A vector for.
#
# Returns:
#   a named vector where the names are cell types and the values are the
#   normalized average library size of each cell type. Sums to 1.
Load_AvgLibSize <- function(dataset, granularity) {
  # This is a list with entries for "A_broad_class" and "A_sub_class"
  A <- readRDS(file.path(dir_singlecell, str_glue("{dataset}_A_matrix.rds")))
  return(A[[str_glue("A_{granularity}")]])
}


# Load_SignatureMatrix: loads the cell-type signature matrix calculated from
# single cell data.
#
# Arguments:
#   dataset = the name of the data set to load
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types to load the signature matrix for.
#   normalization = any value from the 'normalization' argument of
#                 Load_CountsFile, or "cibersortx". If 'normalization' contains
#                 'tmm', the signature matrix will be normalized with TMM
#                 factors. Otherwise it is in CPM. If 'normalization' contains
#                 'log', the signature matrix will be transformed as log2(x+1).
#                 If 'normalization' is "cibersortx", the signature matrix
#                 generated by CibersortX will be loaded instead of the
#                 signature generated by this pipeline.
#
# Returns:
#   a matrix with rows = genes and columns = cell types, where the values are
#   in counts per million (or TMM-normalized), describing the expected value of
#   each gene for eachcell type.
Load_SignatureMatrix <- function(dataset, granularity, normalization) {
  sig_matrix <- readRDS(file.path(dir_signatures, str_glue("{dataset}_signature.rds")))

  if (grepl("tmm", normalization)) {
    sig <- sig_matrix$tmm[[granularity]]
  }
  else if (normalization == "cibersortx") {
    sig <- sig_matrix$cibersortx[[granularity]]
  }
  else {
    sig <- sig_matrix$cpm[[granularity]]
  }

  if (grepl("log", normalization)) {
    sig <- log2(sig + 1)
  }
  return(sig)
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
  saveRDS(covariates,
          file.path(dir_covariates, str_glue("{dataset_name}_covariates.rds")))
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


# Pre-processing load/save functions -------------------------------------------

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
  saveRDS(formula_list,
          file.path(dir_metadata, str_glue("model_formulas_{dataset_name}.rds")))
}


# Load_ModelFormulas: loads a list of formulas for linear modeling from a file
# with a specific filename format.
#
# Arguments:
#   dataset_name = the name of the data set
#
# Returns:
#   a list of formulas ('formula_fixed' and 'formula_mixed'). The formulas are
#   strings, not the 'formula' object. If the formula file doesn't exist, the
#   function returns NULL.
Load_ModelFormulas <- function(dataset_name) {
  model_filename <- file.path(dir_metadata,
                              str_glue("model_formulas_{dataset_name}.rds"))
  if (file.exists(model_filename)) {
    return(readRDS(model_filename))
  }
  return(NULL)
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
  # Use base::saveRDS instead of saveRDS (which defaults to BiocGenerics::saveRDS)
  # so we can save seaRef data, which is not loaded into memory
  base::saveRDS(se_object,
                file = file.path(dir_preprocessed,
                                 str_glue("{dataset_name}_preprocessed.rds")))
}


# Load_PreprocessedData: Loads single cell or bulk RNA seq data that has been
# pre-processed (run through Step02) but has not had cell types mapped or
# normalization applied.
#
# Arguments:
#   dataset_name = the name of the dataset to load
#   remove_excluded = whether to remove genes marked as "exclude" (which include
#                     mitochondrial and non-coding genes). Default is TRUE.
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

  # seaRef only -- load the counts matrix into memory
  if (is(assay(data, "counts"), "DelayedMatrix")) {
    assay(data, "counts") <- as(assay(data, "counts"), "CsparseMatrix")
  }

  return(data)
}


# Algorithm and error calculation I/O ------------------------------------------

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
#
# Returns:
#   nothing
Save_AlgorithmIntermediate <- function(result) {
  if (!is.null(result)) {
    filename <- paste(result$params, collapse = "_")
    saveRDS(result,
            file.path(dir_estimates_tmp,
                      paste0(filename, ".rds")))
  }
}


# Load_AlgorithmIntermediate: loads a single output from an algorithm from the
# temporary folder, in order to re-load previous results when restarting due to
# crash or shutdown. The file is named with the parameters used to generate the
# result.
#
# Arguments:
#   params = a named vector or single-row dataframe with the parameters used to
#            generate the result
#
# Returns:
#   a result object from one of the algorithms if a file matching the input
#   parameters exists, or NULL if not
Load_AlgorithmIntermediate <- function(params) {
  filename <- paste(params, collapse = "_")
  filename <- file.path(dir_estimates_tmp,
                        paste0(filename, ".rds"))
  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  return(NULL)
}


# Save_AlgorithmOutputList: saves a list of final output from one of the
# algorithms to an RDS file, named with a specific format.
#
# Arguments:
#   output_list = a list of outputs from one of the deconvolution algorithms,
#                 which contains output run under different parameter sets
#   algorithm = the name of the algorithm
#   test_dataset = the name of the test dataset
#   name_base = a string that uniquely identifies this list of outputs, which
#               currently is of the format:
#               '<reference_data_name>_<test_data_name>_<granularity>_<reference_input_type>_<normalization>_<regression_method>'
#   top_params = TRUE or FALSE. If TRUE, files will be saved to 'dir_top_estimates'
#                as defined in Filenames.R. If FALSE, files will be saved to
#                'dir_estimates'.
#
# Returns:
#   Nothing
Save_AlgorithmOutputList <- function(output_list, algorithm, test_dataset,
                                     name_base, top_params = FALSE) {
  list_file_name <- str_glue("estimates_{name_base}.rds")

  if (top_params) {
    dir_alg <- file.path(dir_top_estimates, test_dataset, algorithm)
  } else {
    dir_alg <- file.path(dir_estimates, test_dataset, algorithm)
  }

  if (!dir.exists(dir_alg)) {
    dir.create(dir_alg, recursive = TRUE)
  }

  saveRDS(output_list,
          file = file.path(dir_alg, list_file_name))
}


# Load_AlgorithmOutputList: loads a list of output from one of the algorithms,
# named with a specific format.
#
# Arguments:
#   algorithm = the name of the algorithm
#   reference_dataset = the name of the reference data set
#   test_dataset = either "sc_samples" or "training", if the algorithm was run
#                  on single cell samples or training pseudobulk, OR one of the
#                  bulk data sets ("Mayo", "MSBB", "ROSMAP")
#   granularity = either "broad_class" or "sub_class", for which level of cell
#                 types was used for markers and pseudobulk creation.
#   normalization = the normalization strategy used. Same as the 'normalization'
#                   argument to Load_CountsFile.
#   regression_method = the type of regression used for bulk counts. Same as
#                       the 'regression_method' argument to Load_CountsFile.
#
# Returns:
#   a list of outputs from one of the deconvolution algorithms, which contains
#   output run under different parameter sets. If the file doesn't exist, the
#   function returns an empty list.
Load_AlgorithmOutputList <- function(algorithm, reference_dataset, test_dataset,
                                     granularity, reference_input_type,
                                     normalization, regression_method = "none") {
  list_file_format <- paste0("estimates_{algorithm}_{reference_dataset}_",
                             "{test_dataset}_{granularity}_",
                             "{reference_input_type}_{normalization}_",
                             "{regression_method}.rds")

  params_file <- file.path(dir_estimates, test_dataset, algorithm,
                           str_glue(list_file_format))

  if (!file.exists(params_file)) {
    message(paste(params_file, "doesn't exist!"))
    return(list())
  }

  return(readRDS(params_file))
}


# Save_ErrorList: saves a list of calculated error metrics for a single
# algorithm for a single data set, named with a specific format.
#
# Arguments:
#   dataset_name = the name of the bulk dataset for these errors
#   error_list = a list of errors containing entries for mean errors, errors by
#                sample, parameters, and some statistics about the errors
#   algorithm = the name of the algorithm
#   params_data = a one-row dataframe or named list of parameters that were
#                 used to generate the estimates for this list
#   top_params = if TRUE, errors will be saved to the "best errors" directory,
#                otherwise if FALSE they will be saved to the normal errors
#                directory.
#
# Returns:
#   Nothing
Save_ErrorList <- function(dataset_name, error_list, algorithm, params_data,
                           top_params = FALSE) {
  name_base <- paste(params_data, collapse = "_")
  error_file_name <- str_glue("errors_{name_base}.rds")

  if (top_params) {
    dir_errors_alg <- file.path(dir_best_errors, dataset_name, algorithm)
  } else {
    dir_errors_alg <- file.path(dir_errors, dataset_name, algorithm)
  }
  dir.create(dir_errors_alg, recursive = TRUE, showWarnings = FALSE)

  saveRDS(error_list,
          file = file.path(dir_errors_alg, error_file_name))
}


# Load_ErrorList: loads a list of calculated error metrics for a single
# algorithm for a single data set.
#
# Arguments:
#   params = a one-row dataframe or named list of parameters used to generate
#            the file
#   top_params = if TRUE, errors will be loaded from the "best errors" directory,
#                otherwise if FALSE they will be loaded from the normal errors
#                directory.
#
# Returns:
#   a list of errors containing entries for mean errors, errors by sample,
#   parameters, and some statistics about the errors. If the error file doesn't
#   exist, the function returns NULL.
Load_ErrorList <- function(params, top_params = FALSE) {
  name_base <- paste(params, collapse = "_")
  error_file_name <- str_glue("errors_{name_base}.rds")

  if (top_params) {
    error_file <- file.path(dir_best_errors, params$test_data_name,
                            params$algorithm, error_file_name)
  } else {
    error_file <- file.path(dir_errors, params$test_data_name, params$algorithm,
                            error_file_name)
  }

  if (!file.exists(error_file)) {
    #message(paste(error_file, "doesn't exist!"))
    return(NULL)
  }

  return(readRDS(error_file))
}


# Get_ErrorFiles: For a given combination of algorithm / reference dataset /
# bulk dataset / granularity, get all error filenames from the errors directory.
#
# Arguments:
#   bulk_dataset = the name of the bulk dataset
#   algorithm = the name of the algorithm
#   granularity = either "broad_class" or "sub_class"
#   reference_dataset = the name of the single cell dataset used as reference.
#                       If NULL, all error files for the bulk_dataset/algorithm/
#                       granularity are returned.
#
# Returns:
#   a vector of filenames (including full file path)
Get_ErrorFiles <- function(bulk_dataset, algorithm, granularity, reference_dataset = NULL) {
  dir_alg <- file.path(dir_errors, bulk_dataset, algorithm)
  patt <- paste0(reference_dataset, ".*", granularity)
  files <- list.files(dir_alg, pattern = patt, full.names = TRUE)
  return(files)
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
#
# Returns:
#   nothing
Save_ErrorIntermediate <- function(error_obj) {
  params <- error_obj$params %>% select(-total_markers_used)
  filename <- paste(params, collapse = "_") |>
    str_replace_all(" |/", "_") # Get rid of spaces and slashes
  filename <- file.path(dir_errors_tmp, paste0(filename, ".rds"))

  saveRDS(error_obj, filename)
}


# Load_ErrorIntermediate: loads an error calculation from a single param set
# from a temporary folder, in order to re-load previous calculations when
# restarting due to crash or shutdown. The file is named with the parameters
# used to generate the result.
#
# Arguments:
#   params = a named vector or single-row dataframe with the parameters used to
#            generate the result
#
# Returns:
#   a named list containing both the output of one of the deconvolution
#   algorithms plus the error calculations if a file matching the input
#   parameters exists, or NULL if not
Load_ErrorIntermediate <- function(params) {
  filename <- paste(params, collapse = "_") |>
    str_replace_all(" |/", "_") # Get rid of spaces and slashes
  filename <- file.path(dir_errors_tmp, paste0(filename, ".rds"))

  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  return(NULL)
}


# Load_TopParams: loads a top parameters file as calculated in Step 11.
#
# Arguments:
#   params = a named vector or single-row dataframe with the parameters used to
#            generate the top parameters file. Must have the columns named
#            below in the select_cols variable declaration.
#
# Returns:
#   a named list containing the contents of the top parameters file, or NULL
#   if the file doesn't exist.
Load_TopParams <- function(params) {
  select_cols <- c("algorithm", "reference_data_name", "test_data_name",
                   "granularity", "reference_input_type", "normalization",
                   "regression_method")

  file_params <- params[, select_cols]

  file_id <- paste(file_params, collapse = "_")
  top_param_file <- file.path(dir_top_parameters,
                              params$test_data_name,
                              params$algorithm,
                              str_glue("top_parameters_{file_id}.rds"))

  if (!file.exists(top_param_file)) {
    return(NULL)
  }

  top_params <- readRDS(top_param_file)
  return(top_params)
}


# Marker save/load functions ---------------------------------------------------

# Helper function for Save_Markers and Load_Markers
#
# Arguments: see descriptions for Save_Markers arguments
#
# Returns:
#   marker_file_format = a string formatted so that can be used with str_glue()
#                        to insert the proper variable values into the file name
Get_MarkerFileFormat <- function(dataset, granularity, marker_type,
                                 marker_subtype = NULL) {
  if (!(marker_type %in% c("dtangle", "autogenes", "seurat", "deseq2"))) {
    message("marker_type must be one of 'dtangle', 'autogenes', 'seurat', or 'deseq2'")
    return(NULL)
  }

  marker_file_format <- paste0("{marker_type}_markers_{dataset}_{granularity}")

  # Dtangle and AutogeneS need extra information added to the name
  if (marker_type %in% c("dtangle", "autogenes")) {
    marker_file_format <- paste0(marker_file_format, "_{marker_subtype}")
  }

  marker_file_format <- paste0(marker_file_format, ".rds")

  return(marker_file_format)
}


# Save_Markers: a generic function that saves a set of markers to an RDS file,
# named with a specific format. Markers can be from Dtangle/HSPE, AutogeneS,
# Seurat, or DESeq2.
#
# Arguments:
#   markers = a list of markers where each item contains marker genes for one cell type
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
#                    Seurat and DESeq2 don't have a subtype so this can be left
#                    as NULL.
#
# Returns:
#   nothing
Save_Markers <- function(markers, dataset, granularity, marker_type,
                         marker_subtype = NULL) {
  marker_file_format <- Get_MarkerFileFormat(dataset, granularity, marker_type,
                                             marker_subtype)
  if (is.null(marker_file_format)) {
    stop("marker_type must be one of 'autogenes', 'deseq2', 'dtangle', or 'seurat'")
  }

  saveRDS(markers,
          file = file.path(dir_markers, str_glue(marker_file_format)))
}

# Load_Markers: a generic function that reads a set of markers from an RDS file,
# named with a specific format. Markers can be from Dtangle/HSPE, AutogeneS,
# Seurat, or DESeq2.
#
# Arguments: see descriptions for Save_Markers arguments
#
# Returns:
#   a named list where each entry is a list of marker gene names for a given
#   cell type
Load_Markers <- function(dataset, granularity, marker_type,
                         marker_subtype = NULL) {
  marker_file_format <- Get_MarkerFileFormat(dataset, granularity, marker_type,
                                             marker_subtype)

  if (is.null(marker_file_format)) { # This indicates a typo in the input somewhere
    stop("marker_type must be one of 'autogenes', 'deseq2', 'dtangle', or 'seurat'")
  }

  marker_file <- file.path(dir_markers, str_glue(marker_file_format))
  if (!file.exists(marker_file)) {
    message(str_glue("Marker file {marker_file} doesn't exist!"))
    return(NULL) # If a file doesn't exist, the algorithm just moves on
  }

  markers <- readRDS(file = marker_file)

  return(markers)
}


# CibersortX save function -----------------------------------------------------

# Save_SingleCellToCibersort: CibersortX requires passing files to its docker
# container, so the single cell data needs to be written as a dense matrix
# in a specific format. The file is saved to the directory that is shared
# with the CibersortX docker container. omnideconv will do this automatically
# if you pass it a matrix instead of a filename, but it doesn't clear out the
# memory after casting to a dense matrix and doesn't expose its write function,
# so this is a re-implementation.
#
# Arguments:
#   sce - a SingleCellExperiment object, which must have a column for "celltype"
#   output_dir - the path to the directory where the data should be stored
#
# Returns:
#   f_name - the filename of the newly-created file
Save_SingleCellToCibersort <- function(sce, output_dir) {
  # omnideconv requires that the file be named this way
  f_name <- file.path(output_dir, "sample_file_for_cibersort.tsv")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sce <- Dimnames_To_Cibersort(sce)

  # Header (column names) is "GeneSymbol" followed by the cell type labels for
  # each cell
  cts <- as.data.frame(t(c("GeneSymbol", as.character(sce$celltype))))
  vroom_write(cts, file = f_name, delim = "\t", col_names = FALSE)

  # Write in blocks of 500 genes to avoid loading the whole matrix into memory
  # at once, and to speed this process up considerably.
  block_size <- 500
  n_blocks <- ceiling(nrow(sce) / block_size)

  for (block in 1:n_blocks) {
    start <- (block-1) * block_size + 1
    end <- min(block * block_size, nrow(sce))

    df <- as.data.frame(as.matrix(counts(sce[start:end,])))
    df <- cbind(rownames(df), df)
    vroom_write(df, file = f_name, delim = "\t", col_names = FALSE,
                append = TRUE,
                num_threads = parallel::detectCores()-1)
  }

  rm(df)
  gc()

  return(f_name)
}

Genes_To_Cibersort <- function(genes) {
  str_replace_all(as.character(genes), " ENSG", "_ENSG") |> make.names()
}

Cells_To_Cibersort <- function(celltypes, lowercase = FALSE) {
  res <- str_replace_all(as.character(celltypes), "[-/ \\.]", "_") |> make.names()
  if (lowercase) {
    res <- str_to_lower(res)
  }
  return(res)
}

Dimnames_To_Cibersort <- function(obj, col_lower_case = FALSE) {
  # Remove special characters from colnames, remove the space from gene names
  colnames(obj) <- Cells_To_Cibersort(colnames(obj), lowercase = col_lower_case)
  rownames(obj) <- Genes_To_Cibersort(rownames(obj))

  if (is(obj, "SingleCellExperiment")) {
    # Remove special characters from cell type names, and convert them to all
    # lowercase to avoid string sorting issues between R and C.
    obj$celltype <- Cells_To_Cibersort(obj$celltype, lowercase = TRUE)
  }

  return(obj)
}

# Assumes this is a signature matrix where column names are cell types
Cibersort_Celltypes_To_Default <- function(obj, correct_celltype_names) {
  # Cell types come out of CibersortX out of order and lower-case. We put them
  # back in the right order and replace with the correct names.
  cx_names <- Cells_To_Cibersort(correct_celltype_names, lowercase = TRUE)

  stopifnot(all(colnames(obj) %in% cx_names))

  obj <- obj[, cx_names]
  colnames(obj) <- correct_celltype_names
  return(obj)
}

Cibersort_Genes_To_Default <- function(obj, correct_gene_names) {
  cx_genes <- Genes_To_Cibersort(correct_gene_names)

  # There may be fewer than the full number of genes so we have to keep track of
  # the real names of the genes that are kept
  names(cx_genes) <- correct_gene_names
  cx_genes <- cx_genes[cx_genes %in% rownames(obj)]

  stopifnot(length(cx_genes) == nrow(obj))

  obj <- obj[cx_genes, ]
  rownames(obj) <- names(cx_genes)
  return(obj)
}

Cleanup_Cibersort_Docker <- function() {
  find_cmd <- paste0(
    "$(docker ps -a -q --filter ancestor=cibersortx/fractions ",
    "--filter status=exited)"
  )
  system(paste("docker rm", find_cmd))
  gc()
}


# Music save/load functions ----------------------------------------------------

# Save_MusicBasis: Saves the sc_basis object computed by MuSiC::music_basis so
# it doesn't need to be repeatedly calculated.
#
# Arguments:
#   sc_basis = the object returned by MuSiC::music_basis
#   reference_data_name = the name of the single cell data set
#   granularity = either "broad_class" or "sub_class"
#
# Returns:
#   nothing
Save_MusicBasis <- function(sc_basis, reference_data_name, granularity) {
  basis_file <- file.path(dir_music_basis,
                          str_glue("music_basis_{reference_data_name}_{granularity}.rds"))

  saveRDS(sc_basis, basis_file)
}


# Load_MusicBasis: loads the sc_basis object computed by MuSiC::music_basis
#
# Arguments: see descriptions for Save_MusicBasis arguments
#
# Returns:
#   the sc_basis object computed by MuSiC::music_basis. This is a named list of
#   matrices and vectors needed by MuSiC.
Load_MusicBasis <- function(reference_data_name, granularity) {
  basis_file <- file.path(dir_music_basis,
                          str_glue("music_basis_{reference_data_name}_{granularity}.rds"))

  if (!file.exists(basis_file)) {
    return(NULL)
  }

  return(readRDS(basis_file))
}
