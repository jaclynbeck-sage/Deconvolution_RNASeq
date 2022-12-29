library(DeconRNASeq)
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(stringr)

source("Filenames.R")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito") #,
              #"seaRef") #, "seaAD")

datatypes <- list("donors", "training")

###load bulk and snRNA-seq data###
for (sndata in datasets) {
  sce <- readRDS(file.path(dir_input, paste(sndata, "sce.rds", sep = "_")))
  meta <- colData(sce)

  # Convert gene names to Ensembl IDs
  genes <- rowData(sce)
  rownames(sce) <- genes[rownames(sce), "Ensembl.ID"]

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

  bulk <- read.table(file_rosmap, header = TRUE, row.names = 1)
  bulk_cpm <- calculateCPM(bulk)
  bulk_cpm <- as.data.frame(as.matrix(bulk_cpm))

  keepgene <- intersect(rownames(sce), rownames(bulk))

  signature <- signature[intersect(rownames(signature), keepgene),]

  # Clear up as much memory as possible
  rm(bulk, sce, sce_cpm)
  gc()

  for (datatype in datatypes) {
    best_params <- readRDS(file.path(dir_output, paste0("best_params_", sndata,
                                                        "_", datatype, "_",
                                                        cellclasstype, ".rds")))
    best_params <- best_params[["deconRNASeq"]]
    best_params <- str_split(best_params, pattern = "_", simplify = TRUE)
    use.scale.params <- as.logical(best_params[,4]) # This is the only variable that changes

    decon_list <- list()

    for (use.scale in use.scale.params) {
      name <- paste(sndata, cellclasstype,
                    "usescale", use.scale,
                    "normalization", "cpm", sep = "_")

      res <- DeconRNASeq(bulk_cpm[keepgene,], signature, proportions = NULL,
                         known.prop = FALSE, use.scale = use.scale, fig = FALSE)

      res$Est.prop <- res$out.all
      rownames(res$Est.prop) <- colnames(bulk_cpm)
      res <- res[c("Est.prop", "out.pca")]

      decon_list[[name]] <- res

      print(name)
    } # end use.scale loop

    # Save the completed list
    print("Saving final list...")
    saveRDS(decon_list, file = file.path(dir_output,
                                         paste0("deconRNASeq_list_", sndata,
                                                  "_", datatype, "_",
                                                cellclasstype, "_ROSMAP.rds")))

    rm(decon_list)
    gc()
  } # end datatypes loop
}








