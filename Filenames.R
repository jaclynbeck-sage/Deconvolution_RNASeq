# Filenames that are used in multiple files

dir_data <- "DeconvolutionData"
dir_input <- file.path(dir_data, "input")
dir_pseudobulk <- file.path(dir_data, "pseudobulk")
dir_output <- file.path(dir_data, "output")
dir_markers <- file.path(dir_output, "markers")
dir_params_lists <- file.path(dir_output, "params_lists")

dir_cain_raw <- file.path(dir_input, "cain_raw")
dir_lau_raw <- file.path(dir_input, "lau_raw")
dir_leng_raw <- file.path(dir_input, "leng_raw")
dir_mathys_raw <- file.path(dir_input, "mathys_raw")
dir_morabito_raw <- file.path(dir_input, "morabito_raw")
dir_seaad_raw <- file.path(dir_input, "sea-ad_raw")

file_searef_h5 <- file.path(dir_seaad_raw, "ref_counts.h5ad")
file_seaad_h5 <- file.path(dir_seaad_raw, "seaad_counts.h5ad")

url_searef_h5 <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
url_seaad_h5 <- "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2022-08-18.h5ad"

file_rosmap <- file.path(dir_input, "ROSMAP_DLPFC_Counts.tsv")

# Make sure these directories exist
for (D in ls(pat = "dir_")) {
  dir.create(eval(parse(text = D)), recursive = TRUE, showWarnings = FALSE)
}
