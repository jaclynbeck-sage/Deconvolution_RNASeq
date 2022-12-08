# Filenames that are used in multiple files

dir_data <- "DeconvolutionData"
dir_input <- file.path(dir_data, "input")
dir_pseudobulk <- file.path(dir_data, "pseudobulk")
dir_output <- file.path(dir_data, "output")

dir_lau_raw <- file.path(dir_input, "lau_raw")
dir_leng_raw <- file.path(dir_input, "leng_raw")
dir_mathys_raw <- file.path(dir_input, "mathys_raw")
dir_morabito_raw <- file.path(dir_input, "morabito_raw")
dir_seaad_raw <- file.path(dir_input, "sea-ad_raw")

# Make sure these directories exist
for (D in ls(pat = "dir_")) {
  dir.create(eval(parse(text = D)), recursive = TRUE, showWarnings = FALSE)
}
