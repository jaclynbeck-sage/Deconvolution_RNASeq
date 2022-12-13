library(DeconRNASeq)
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatypes <- list("donors", "training")

###load bulk and snRNA-seq data###
for (sndata in datasets) {
  for (datatype in datatypes) {
    if (cellclasstype == "fine") {
      load(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_finecelltypes_30percentlimit.rda")))
    }
    if (cellclasstype=="broad") {
      pseudobulk <- readRDS(file.path(dir_pseudobulk,
                                      paste0("pseudobulk_", sndata, "_",
                                             datatype, "_broadcelltypes.rds")))
    }

    sce <- readRDS(file.path(dir_input, paste(sndata, "sce.rds", sep = "_")))
    meta <- colData(sce)

    sce_cpm <- calculateCPM(counts(sce))

    signature <- lapply(levels(meta$broadcelltype), function(X) {
      cells <- meta[meta$broadcelltype == X,]
      rowMeans(sce_cpm[,rownames(cells)])
    })
    names(signature) <- levels(meta$broadcelltype)
    signature <- do.call(cbind, signature)

    # Filter for > 1 cpm
    ok <- which(rowSums(signature >= 1) > 0)
    signature <- as.data.frame(signature[ok, ])

    pseudobulk <- assays(pseudobulk)[["counts"]]

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(sce), rownames(pseudobulk))

    pseudobulk_cpm <- calculateCPM(pseudobulk)
    pseudobulk_cpm <- as.data.frame(as.matrix(pseudobulk_cpm))

    # Clear up as much memory as possible
    rm(pseudobulk, sce, sce_cpm)
    gc()

    decon_list <- list()

    for (use.scale in c(TRUE, FALSE)) {
      name <- paste(sndata, cellclasstype,
                    "usescale", use.scale,
                    "normalization", "cpm", sep = "_")

      res <- DeconRNASeq(pseudobulk_cpm[keepgene,], signature, proportions = NULL,
                         known.prop = FALSE, use.scale = use.scale, fig = FALSE)

      res$Est.prop <- res$out.all
      rownames(res$Est.prop) <- colnames(pseudobulk_cpm)
      res <- res[c("Est.prop", "out.pca")]

      decon_list[[name]] <- res

      print(name)
    } # end use.scale loop

    # Save the completed list
    print("Saving final list...")
    saveRDS(decon_list, file = file.path(dir_output,
                                         paste0("deconRNASeq_list_", sndata,
                                                  "_", datatype, "_",
                                                cellclasstype, ".rds")))

    rm(decon_list, meta, pseudobulk_cpm, signature)
    gc()
  } # end datatypes loop
}








