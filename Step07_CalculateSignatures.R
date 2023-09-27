# This script calculates cell-type gene signatures, using the broad and sub
# class cell type assignments output from mapping in Step 04. It also calculates
# the average cell size (average counts per cell per cell type, normalized) for
# use with MuSiC.
library(SummarizedExperiment)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

for (dataset in datasets) {
  # Calculate the "A" matrix that is needed to convert propCells to pctRNA
  A_broad <- CalculateA(dataset, "broad_class")
  A_sub <- CalculateA(dataset, "sub_class")

  saveRDS(list("A_broad_class" = A_broad, "A_sub_class" = A_sub),
          file = file.path(dir_input, str_glue("{dataset}_A_matrix.rds")))

  # Calculate a signature for each cell type. This matrix includes all genes in
  # the data set and isn't filtered at this point.
  signatures <- lapply(c("cpm", "tmm"), function(output_type) {
    sig_broad <- CalculateSignature(dataset, "broad_class", output_type)
    sig_sub <- CalculateSignature(dataset, "sub_class", output_type)
    return(list("broad_class" = sig_broad, "sub_class" = sig_sub))
  })

  names(signatures) <- c("cpm", "tmm")

  saveRDS(signatures,
          file = file.path(dir_input, str_glue("{dataset}_signature.rds")))
}
