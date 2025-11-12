# Filenames that are used in multiple files. This isn't in the config.yaml file
# because R code in the yaml file can't reference another yaml variable that is
# also defined by R code, so none of the file.path() functions would work.

dir_data <- "data"
dir_figures <- "figures"
dir_downloads <- file.path("data", "downloads")

dir_metadata <- file.path(dir_data, "01_metadata")
dir_preprocessed <- file.path(dir_data, "02_preprocessed")
dir_bulk <- file.path(dir_data, "03_bulk_input")
dir_singlecell <- file.path(dir_data, "04_singlecell_input")
dir_pseudobulk <- file.path(dir_data, "05_pseudobulk")
dir_signatures <- file.path(dir_data, "06_signatures")
dir_markers <- file.path(dir_data, "07_markers")
dir_estimates <- file.path(dir_data, "08_estimates")
dir_algorithm_models <- file.path(dir_data, "08_algorithm_models")
dir_errors <- file.path(dir_data, "10_errors")
dir_top_parameters <- file.path(dir_data, "11_top_parameters")
dir_top_estimates <- file.path(dir_data, "12_top_estimates")
dir_best_errors <- file.path(dir_data, "14_best_error_calculations")

dir_tmp <- file.path(dir_data, "tmp")
dir_estimates_tmp <- file.path(dir_estimates, "tmp")
dir_errors_tmp <- file.path(dir_errors, "tmp")

dir_cibersort_corrected_signatures <- file.path(dir_signatures, "cibersortx_batch_corrected")
dir_cibersort <- "/data" # outside this working directory, gets shared with cibersort docker

dir_scaden_models <- file.path(dir_algorithm_models, "scaden_models")
dir_music_basis <- file.path(dir_algorithm_models, "music_basis")
dir_hspe_params <- file.path(dir_algorithm_models, "hspe_params")

dir_raw_data <- file.path(dir_downloads, "raw_data")
dir_covariates <- file.path(dir_metadata, "covariates")

dir_lau_raw <- file.path(dir_raw_data, "lau_raw")
dir_leng_raw <- file.path(dir_raw_data, "leng_raw")
dir_seaad_raw <- file.path(dir_raw_data, "sea-ad_raw")

file_searef_h5 <- file.path(dir_seaad_raw, "ref_counts.h5ad")
file_seaad_h5 <- file.path(dir_seaad_raw, "seaad_counts.h5ad")

#url_seaad_h5 <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
#url_seaad_merfish <- "https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-02-13.h5ad"
#url_seaad_donor_metadata <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx"

file_gene_list <- file.path(dir_metadata, "ensembl_gene_list.csv")
file_rosmap_ihc_proportions <- file.path(dir_metadata, "ihc_proportions_normalized.csv")

# Make sure these directories exist
for (D in ls(pattern = "dir_")) {
  dir.create(eval(parse(text = D)), recursive = TRUE, showWarnings = FALSE)
}
