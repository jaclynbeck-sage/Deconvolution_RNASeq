# Filenames that are used in multiple files

dir_data <- "DeconvolutionData"
dir_input <- file.path(dir_data, "input")
dir_pseudobulk <- file.path(dir_data, "pseudobulk")
dir_output <- file.path(dir_data, "output")

# Make sure these directories exist
dir.create(dir_pseudobulk, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)
