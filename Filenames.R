# Filenames that are used in multiple files

dir_data <- "DeconvolutionData"
dir_input <- file.path(dir_data, "input")
dir_pseudobulk <- file.path(dir_data, "pseudobulk")
dir_output <- file.path(dir_data, "output")

dir_mathys_raw <- file.path(dir_input, "mathys_raw")

# Make sure these directories exist
for (D in c(dir_pseudobulk, dir_output, dir_mathys_raw)) {
  dir.create(D, recursive = TRUE, showWarnings = FALSE)
}
