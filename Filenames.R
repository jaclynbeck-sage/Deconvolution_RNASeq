# Filenames that are used in multiple files

dir_data <- "data"
dir_figures <- "figures"

dir_input <- file.path(dir_data, "input")
dir_pseudobulk <- file.path(dir_data, "pseudobulk")
dir_output <- file.path(dir_data, "output")
dir_metadata <- file.path(dir_data, "metadata")
dir_tmp <- file.path(dir_data, "tmp")
dir_markers <- file.path(dir_data, "markers")
dir_cibersort <- "/data" # outside this working directory, gets shared with cibersort docker

dir_params_lists <- file.path(dir_output, "params_lists")
dir_errors <- file.path(dir_output, "errors")
dir_best_errors <- file.path(dir_errors, "best_errors")

dir_params_lists_tmp <- file.path(dir_params_lists, "tmp")
dir_errors_tmp <- file.path(dir_errors, "tmp")

dir_raw_data <- file.path(dir_input, "raw_data")
dir_preprocessed <- file.path(dir_input, "preprocessed")
dir_map_reference <- file.path(dir_input, "map_reference")
dir_covariates <- file.path(dir_metadata, "covariates")

dir_cain_raw <- file.path(dir_raw_data, "cain_raw")
dir_lau_raw <- file.path(dir_raw_data, "lau_raw")
dir_leng_raw <- file.path(dir_raw_data, "leng_raw")
dir_mathys_raw <- file.path(dir_raw_data, "mathys_raw")
dir_morabito_raw <- file.path(dir_raw_data, "morabito_raw")
dir_seaad_raw <- file.path(dir_raw_data, "sea-ad_raw")

dir_mayo_raw <- file.path(dir_raw_data, "mayo_raw")
dir_msbb_raw <- file.path(dir_raw_data, "msbb_raw")
dir_rosmap_raw <- file.path(dir_raw_data, "rosmap_raw")

file_searef_h5 <- file.path(dir_seaad_raw, "ref_counts.h5ad")
file_seaad_h5 <- file.path(dir_seaad_raw, "seaad_counts.h5ad")

url_searef_h5 <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
url_seaad_h5 <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad"
url_seaad_merfish <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/MERFISH/SEAAD_MTG_MERFISH_all-nuclei.2023-05-08.h5ad"

dir_mayo_output <- file.path(dir_output, "mayo")
dir_msbb_output <- file.path(dir_output, "msbb")
dir_rosmap_output <- file.path(dir_output, "rosmap")

file_gene_list <- file.path(dir_metadata, "ensembl_gene_list.csv")
file_rosmap_ihc_proportions <- file.path(dir_metadata, "ihc_proportions_normalized.csv")

# Make sure these directories exist
for (D in ls(pattern = "dir_")) {
  dir.create(eval(parse(text = D)), recursive = TRUE, showWarnings = FALSE)
}
