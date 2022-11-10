library(DeconRNASeq)
library(Matrix)
library(SummarizedExperiment)
library(scuttle)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- list("mathys")#,"cain","lau","morabito","lengSFG","lengEC")
datatypes <- list("donors", "training")

###load bulk and snRNA-seq data###
for (sndata in datasets) {
  for (datatype in datatypes) {
    if (cellclasstype == "fine") {
      load(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_finecelltypes_30percentlimit.rda")))
    }
    if (cellclasstype=="broad") {
      if (datatype == "donors") {
        se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", sndata, "_bydonor_broadcelltypes.rds")))
      }
      else {
        se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_",sndata,"_broadcelltypes.rds")))
      }
    }

    load(file.path(dir_input, paste0(sndata,"_counts.rda")))

    # TODO: metadata processing should be done in a function since multiple files do this
    meta <- read.csv(file.path(dir_input, paste0(sndata,"_metadata.csv")), as.is=T)
    rownames(meta) <- meta$cellid

    meta$broadcelltype <- factor(meta$broadcelltype)
    meta$subcluster <- factor(meta$subcluster)

    keep <- intersect(meta$cellid, colnames(fullmat))
    meta <- meta[keep,]
    fullmat <- fullmat[,keep]
    fullmat_cpm <- calculateCPM(fullmat)

    signature <- lapply(levels(meta$broadcelltype), function(X) {
      cells <- meta[meta$broadcelltype == X,]
      rowMeans(fullmat_cpm[,rownames(cells)])
    })
    names(signature) <- levels(meta$broadcelltype)
    signature <- do.call(cbind, signature)

    # Filter for > 1 cpm
    ok <- which(rowSums(signature >= 1) > 0)
    signature <- as.data.frame(signature[ok, ])

    pseudobulk <- assays(se)[["counts"]]

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(fullmat),rownames(pseudobulk))

    pseudobulk_cpm <- calculateCPM(pseudobulk)
    pseudobulk_cpm <- as.data.frame(as.matrix(pseudobulk_cpm))

    # Clear up as much memory as possible
    rm(se, fullmat, pseudobulk)
    gc()

    decon_list <- list()

    for (use.scale in c(TRUE, FALSE)) {
      name <- paste(sndata, cellclasstype,
                    "usescale", use.scale,
                    "normalization", "cpm", sep = "_")

      res <- DeconRNASeq(pseudobulk_cpm, signature, proportions = NULL,
                         known.prop = FALSE, use.scale = use.scale, fig = FALSE)

      res$Est.prop <- res$out.all
      rownames(res$Est.prop) <- colnames(pseudobulk_cpm)
      res <- res[c("Est.prop", "out.pca")]

      decon_list[[name]] <- res

      print(name)
    }

    # Save the completed list
    saveRDS(decon_list, file = file.path(dir_output,
                                         paste0("deconRNASeq_list_", sndata,
                                                  "_", datatype, "_",
                                                cellclasstype, ".rds")))
  }
}








