library(dplyr)
library(nnls)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "FileIO_HelperFunctions.R"))

bulk <- Load_BulkData("MSBB", output_type = "cpm")
bulk_cpm <- assay(bulk, "counts")

markers <- Load_Markers(dataset = "cain", granularity = "broad",
                        marker_type = "dtangle", marker_subtype = "diff",
                        input_type = "pseudobulk")

markers <- lapply(markers, function(X) {X[X %in% rownames(bulk_cpm)]})

markers <- OrderMarkers_ByCorrelation(markers, bulk_cpm)

markers_filt <- lapply(markers, function(X) {
  len <- min(length(X), 10)
  X[1:len]
})

R <- sweep(bulk_cpm, 1, bulk_cpm[,1], "/")

R_filt <- R[unlist(markers_filt),]

meds <- apply(R_filt, 2, function(sample) {
  sapply(markers, function(M) {
    median(sample[M], na.rm = TRUE)
  })
})

meds <- t(meds)

res <- nnls(meds, matrix(1, nrow = nrow(meds), ncol = 1))
res$x / sum(res$x)

# I think the samples need to be batch corrected for this to work well
stats::lm.fit()
