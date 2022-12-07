library(synapser)
library(biomaRt)
library(Seurat)

# Assumes you have already authenticated a Synapse login or have a Synapse
# config set up with proper credentials
DownloadFromSynapse <- function(synIDs, downloadLocation) {
  synLogin()

  files <- list()

  for (name in names(synIDs)) {
    files[[name]] <- synGet(synIDs[[name]], downloadLocation = downloadLocation)
  }

  return(files)
}

# gene.list is a vector of gene symbols
# TODO: Seurat has a GeneSymbolThesaurus function to get aliases. Might be useful
GeneSymbolToEnsembl_Biomart <- function(gene.list) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ens.genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                     filters = "external_gene_name",
                     values = gene.list,
                     mart = mart)

  dupes <- duplicated(ens.genes$external_gene_name)
  ens.genes <- ens.genes[!dupes,]

  return(ens.genes)
}


##### Custom functions for each data set #####
# Metadata columns must end up in the correct order to be renamed as follows:
#   c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
#
# For gene symbol -> Ensembl ID conversion:
#   Symbols can map to multiple Ensembl IDs. For this application, we use the
#   first Ensembl ID in the list for duplicate symbols, which is what Seurat
#   does. Since most papers used Seurat, we assume this is consistent.

##### Mathys, et al., 2019 #####

ReadMetadata_Mathys <- function(files) {
  clinical <- read.csv(files[["clinical"]]$path)
  cellMeta <- read.table(files[["metadata"]]$path, sep = "\t", header = TRUE)
  metadata <- merge(cellMeta, clinical,
                    by = "projid", all.x = TRUE, all.y = FALSE)

  # See https://www.synapse.org/#!Synapse:syn3191090 for this information (cogdx)
  diagnosis.codes <- list("1" = "Control",
                          "2" = "MCI",
                          "3" = "MCI",
                          "4" = "AD",
                          "5" = "AD",
                          "6" = "Other")
  diagnosis.codes <- melt(diagnosis.codes)
  colnames(diagnosis.codes) <- c("diagnosis", "cogdx")

  metadata <- merge(metadata, diagnosis.codes, by = "cogdx")

  metadata <- metadata[, c("TAG", "projid", "diagnosis", "broad.cell.type", "Subcluster")]

  return(metadata)
}

ReadCounts_Mathys <- function(files) {
  counts <- ReadMtx(mtx = files[["counts"]]$path,
                    cells = files[["metadata"]]$path,
                    features = files[["genes"]]$path,
                    feature.column = 1, skip.cell = 1)
  return(counts)
}

GeneSymbolToEnsembl_Mathys(files) {
  # This file has mappings from Ensembl ID to gene symbol
  genes <- read.table(files[["genes_ensembl"]]$path, sep = "\t")
  colnames(genes) <- c("Ensembl.ID", "Symbol")

  dupes <- duplicated(genes$Symbol)
  genes <- genes[!dupes,]
  rownames(genes) <- genes$Symbol

  return(genes)
}


##### Morabito, et al., 2021 #####

ReadMetadata_Morabito <- function(files) {
  metadata <- read.csv(files[["metadata"]]$path)

  rownames(metadata) <- metadata$X
  metadata <- metadata[, c("X", "Sample.ID", "Diagnosis", "celltype", "cluster")]

  return(metadata)
}

ReadCounts_Morabito <- function(files) {
  counts <- ReadMtx(mtx = files[["counts"]]$path,
                    cells = files[["barcodes"]]$path,
                    features = files[["genes"]]$path,
                    feature.column = 1)
  return(counts)
}

GeneSymbolToEnsembl_Morabito <- function(gene.list) {
  ens.genes <- GeneSymbolToEnsembl_Biomart(gene.list)

  genes <- data.frame(Symbol = gene.list)
  genes <- merge(genes, ens.genes,
                 by.x = "Symbol", by.y = "external_gene_name", all = TRUE)
  colnames(genes) <- c("Symbol", "Ensembl.ID")
  rownames(genes) <- genes$Symbol

  return(genes)
}

