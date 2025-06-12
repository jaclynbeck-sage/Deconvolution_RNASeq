# Filenames that are used in multiple files

dir_data <- "data"
dir_figures <- "figures"

dir_input <- file.path(dir_data, "input")
dir_pseudobulk <- file.path(dir_data, "pseudobulk")
dir_output <- file.path(dir_data, "output")
dir_metadata <- file.path(dir_data, "metadata")
dir_tmp <- file.path(dir_data, "tmp")
dir_markers <- file.path(dir_data, "markers")
dir_signatures <- file.path(dir_data, "signatures")
dir_analysis <- file.path(dir_data, "analysis")
dir_cibersort_corrected_signatures <- file.path(dir_signatures, "cibersortx_batch_corrected")
dir_cibersort <- "/data" # outside this working directory, gets shared with cibersort docker

dir_estimates <- file.path(dir_output, "estimates")
dir_top_estimates <- file.path(dir_output, "top_estimates")
dir_errors <- file.path(dir_output, "errors")
dir_top_parameters <- file.path(dir_output, "top_parameters")
dir_best_errors <- file.path(dir_output, "best_errors")

dir_estimates_tmp <- file.path(dir_estimates, "tmp")
dir_errors_tmp <- file.path(dir_errors, "tmp")

dir_scaden_models <- file.path(dir_output, "scaden_models")
dir_music_basis <- file.path(dir_output, "music_basis")
dir_hspe_params <- file.path(dir_output, "hspe_params")

dir_raw_data <- file.path(dir_input, "raw_data")
dir_preprocessed <- file.path(dir_input, "preprocessed")
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
url_seaad_h5 <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
url_seaad_merfish <- "https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-02-13.h5ad"
url_seaad_donor_metadata <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx"

dir_mayo_output <- file.path(dir_estimates, "Mayo")
dir_msbb_output <- file.path(dir_estimates, "MSBB")
dir_rosmap_output <- file.path(dir_estimates, "ROSMAP")

file_gene_list <- file.path(dir_metadata, "ensembl_gene_list.csv")
file_rosmap_ihc_proportions <- file.path(dir_metadata, "ihc_proportions_normalized.csv")

# Make sure these directories exist
for (D in ls(pattern = "dir_")) {
  dir.create(eval(parse(text = D)), recursive = TRUE, showWarnings = FALSE)
}
