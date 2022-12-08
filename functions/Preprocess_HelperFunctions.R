library(synapser)
library(biomaRt)
library(Seurat)
library(reshape2)
library(readxl)
library(GEOquery)

##### Generic functions #####

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

  genes <- data.frame(Symbol = gene.list)
  genes <- merge(genes, ens.genes,
                 by.x = "Symbol", by.y = "external_gene_name", all = TRUE)
  colnames(genes) <- c("Symbol", "Ensembl.ID")
  rownames(genes) <- genes$Symbol

  return(genes)
}


##### Functions that call dataset-specific functions #####

DownloadData <- function(dataset) {
  files <- switch(dataset,
                  "cain" = DownloadData_Cain(),
                  "lau" = DownloadData_Lau(),
                  "mathys" = DownloadData_Mathys(),
                  "morabito" = DownloadData_Morabito())
  return(files)
}

ReadMetadata <- function(dataset, files) {
  metadata <- switch(dataset,
                     "cain" = ReadMetadata_Cain(files),
                     "lau" = ReadMetadata_Lau(files),
                     "mathys" = ReadMetadata_Mathys(files),
                     "morabito" = ReadMetadata_Morabito(files))
  return(metadata)
}

ReadCounts <- function(dataset, files) {
  counts <- switch(dataset,
                   "cain" = ReadCounts_Cain(files),
                   "lau" = NULL, # Needs more work
                   "mathys" = ReadCounts_Mathys(files),
                   "morabito" = ReadCounts_Morabito(files))
  return(counts)
}

GeneSymbolToEnsembl <- function(dataset, files = NULL, gene.list = NULL) {
  genes <- switch(dataset,
                  "cain" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "lau" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "mathys" = GeneSymbolToEnsembl_Mathys(files),
                  "morabito" = GeneSymbolToEnsembl_Biomart(gene.list))
  return(genes)
}

##### Custom functions for each data set #####
# Metadata columns must end up in the correct order to be renamed as follows:
#   c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
#
# For gene symbol -> Ensembl ID conversion:
#   Symbols can map to multiple Ensembl IDs. For this application, we use the
#   first Ensembl ID in the list for duplicate symbols, which is what Seurat
#   does. Since most papers used Seurat, we assume this is consistent.

##### Cain, et al., 2020 [preprint] #####
# https://doi.org/10.1101/2020.12.22.424084

DownloadData_Cain <- function() {
  # TODO
  files <- list("metadata" = file.path(dir_input, "cain_metadata.csv"),
                "counts" = file.path(dir_input, "cain_counts.rda"))
  return(files)
}

ReadMetadata_Cain <- function(files) {
  metadata <- read.csv(files[["metadata"]])
  metadata <- metadata[, c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster" )]

  return(metadata)
}

ReadCounts_Cain <- function(files) {
  load(files[["counts"]])
  return(fullmat)
}

GeneSymbolToEnsembl_Cain <- function(gene.list) {
  genes <- GeneSymbolToEnsembl_Biomart(gene.list)
  return(genes)
}


##### Lau, et al., 2020 #####
# https://doi.org/10.1073/pnas.2008762117

DownloadData_Lau <- function() {
  gse <- getGEO("GSE157827", destdir = dir_lau_raw)
  geo_metadata <- pData(phenoData(gse[[1]]))

  geo <- getGEOSuppFiles(GEO = "GSE157827", makeDirectory = FALSE,
                         baseDir = dir_lau_raw)

  #untar(rownames(geo)[1], list = TRUE) # List files inside
  untar(rownames(geo)[1], exdir = dir_lau_raw)

  files <- list("metadata" = file.path(dir_input, "lau_metadata.csv"),
                "geo_metadata" = geo_metadata) # TODO
  return(files)
}

ReadMetadata_Lau <- function(files) {
  metadata <- read.csv(files[["metadata"]])
  metadata <- metadata[, c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster" )]

  return(metadata)
}

ReadCounts_Lau <- function(geo_metadata, metadata) {
  samples <- rownames(geo_metadata)
  counts_list <- list()

  for (sample in samples) {
    files <- list.files(path = dir_lau_raw,
                        pattern = paste0(sample, "_.*\\.gz"), full.names = TRUE)
    matrix_file <- grep("matrix", files, value = TRUE)
    barcodes_file <- grep("barcodes", files, value = TRUE)
    features_file <- grep("features", files, value = TRUE)

    counts <- ReadMtx(mtx = matrix_file, cells = barcodes_file, features = features_file)
    colnames(counts) <- paste(colnames(counts), geo_metadata[sample, "title"], sep = "_")

    counts_list[[sample]] <- counts
  }

  counts <- do.call(cbind, counts_list) # Genes are in the same order for each data frame. Should probably write more robust code to explicitly ensure this.
  counts <- counts[, metadata$cellid] # Filter to use only cells in the provided metadata file

  return(counts)
}

GeneSymbolToEnsembl_Lau <- function(gene.list) {
  genes <- GeneSymbolToEnsembl_Biomart(gene.list)
  return(genes)
}


##### Mathys, et al., 2019 #####
# http://dx.doi.org/10.1038/s41586-019-1195-2

DownloadData_Mathys <- function() {
  synIDs <- list("clinical" = "syn3191087",
                 "metadata" = "syn18686372",
                 "genes_ensembl" = "syn18687959",
                 "counts" = "syn18686381",
                 "genes" = "syn18686382")

  files <- DownloadFromSynapse(synIDs, dir_mathys_raw)
  return(files)
}

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

GeneSymbolToEnsembl_Mathys <- function(files) {
  # This file has mappings from Ensembl ID to gene symbol
  genes <- read.table(files[["genes_ensembl"]]$path, sep = "\t")
  colnames(genes) <- c("Ensembl.ID", "Symbol")

  dupes <- duplicated(genes$Symbol)
  genes <- genes[!dupes,]
  rownames(genes) <- genes$Symbol

  return(genes)
}


##### Morabito, et al., 2021 #####
# https://doi.org/10.1038/s41588-021-00894-z

DownloadData_Morabito <- function() {
  synIDs <- list("metadata" = "syn22130806",
                 "barcodes" = "syn24978676",
                 "genes" = "syn24978737",
                 "counts" = "syn22130805")

  files <- DownloadFromSynapse(synIDs, dir_morabito_raw)
  return(files)
}

ReadMetadata_Morabito <- function(files) {
  metadata <- read.csv(files[["metadata"]]$path)
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
  genes <- GeneSymbolToEnsembl_Biomart(gene.list)
  return(genes)
}


##### SEA-AD: Not finished yet #####

ReadMetadata_SEA <- function(files) {
  donor_metadata <- read_excel(path = files[["donor_metadata"]])
  cell_metadata <- read.csv(files[["cell_metadata"]])

  donor_metadata <- subset(donor_metadata, donor_metadata$`Donor ID` %in% cell_metadata$external_donor_name_label)
  # TODO: class_label gives GABA/Glut neuronal split
  metadata <- metadata[, c("sample_name", "external_donor_name_label",
                           "diagnosis", "subclass_label", "cluster_label" )]

  return(metadata)
}

ReadCounts_SEA <- function(files) {
  counts <- readRDS(files[["counts"]]) # Seurat object
  counts <- subset(counts, QCpass == "True")

  return(counts)
}
