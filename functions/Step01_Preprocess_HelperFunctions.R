library(synapser)
library(biomaRt)
library(Seurat)
library(reshape2)
library(stringr)
library(GEOquery)
library(HDF5Array)
library(rhdf5)
library(reticulate)

source("Filenames.R")

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

# This function is unused, saving just in case.
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

# Mayo, MSBB, and ROSMAP all have files in Synapse with the conversion between
# Ensembl ID and gene symbol, in the same format
EnsemblIdToGeneSymbol_BulkData <- function(files) {
  genes <- read.table(files[["genes"]]$path, header = TRUE)
  genes <- genes[,c("hgnc_symbol", "ensembl_gene_id")]
  colnames(genes) <- c("Symbol", "Ensembl.ID")
  rownames(genes) <- genes$Ensembl.ID
  return(genes)
}


##### Functions that call dataset-specific functions #####

DownloadData <- function(dataset) {
  files <- switch(dataset,
                  "cain" = DownloadData_Cain(),
                  "lau" = DownloadData_Lau(),
                  "lengEC" = DownloadData_LengEC(),
                  "lengSFG" = DownloadData_LengSFG(),
                  "mathys" = DownloadData_Mathys(),
                  "morabito" = DownloadData_Morabito(),
                  "seaRef" = DownloadData_SEARef(),
                  "seaAD" = DownloadData_SEAAD(),
                  "Mayo" = DownloadData_Mayo(),
                  "MSBB" = DownloadData_MSBB(),
                  "ROSMAP" = DownloadData_ROSMAP())

  return(files)
}

ReadMetadata <- function(dataset, files) {
  metadata <- switch(dataset,
                     "cain" = ReadMetadata_Cain(files),
                     "lau" = ReadMetadata_Lau(files),
                     "lengEC" = ReadMetadata_Leng(files),
                     "lengSFG" = ReadMetadata_Leng(files),
                     "mathys" = ReadMetadata_Mathys(files),
                     "morabito" = ReadMetadata_Morabito(files),
                     "seaRef" = ReadMetadata_SEARef(files),
                     "seaAD" = ReadMetadata_SEAAD(files),
                     "Mayo" = ReadMetadata_BulkData(files),
                     "MSBB" = ReadMetadata_BulkData(files),
                     "ROSMAP" = ReadMetadata_BulkData(files))
  return(metadata)
}

ReadCounts <- function(dataset, files, metadata) {
  counts <- switch(dataset,
                   "cain" = ReadCounts_Cain(files),
                   "lau" = ReadCounts_Lau(files, metadata),
                   "lengEC" = ReadCounts_Leng(files),
                   "lengSFG" = ReadCounts_Leng(files),
                   "mathys" = ReadCounts_Mathys(files),
                   "morabito" = ReadCounts_Morabito(files),
                   "seaRef" = ReadCounts_SEARef(files),
                   "seaAD" = ReadCounts_SEAAD(files),
                   "Mayo" = ReadCounts_BulkData(files),
                   "MSBB" = ReadCounts_BulkData(files),
                   "ROSMAP" = ReadCounts_BulkData(files))
  return(counts)
}

EnsemblIdToGeneSymbol <- function(dataset, files = NULL, gene.list = NULL) {
  genes <- switch(dataset,
                  "cain" = NULL,
                  "lau" = NULL,
                  "lengEC" = NULL,
                  "lengSFG" = NULL,
                  "mathys" = NULL,
                  "morabito" = NULL,
                  "seaRef" = NULL,
                  "seaAD" = NULL,
                  "Mayo" = EnsemblIdToGeneSymbol_BulkData(files),
                  "MSBB" = EnsemblIdToGeneSymbol_BulkData(files),
                  "ROSMAP" = EnsemblIdToGeneSymbol_BulkData(files))
  return(genes)
}

##### Custom functions for each data set #####


##### Cain, et al., 2020 [preprint] #####
# https://doi.org/10.1101/2020.12.22.424084

DownloadData_Cain <- function() {
  # TODO these are not the original data files
  synIDs <- list("metadata" = "syn38609692",
                 "counts" = "syn38598183")

  files <- DownloadFromSynapse(synIDs, dir_cain_raw)
  return(files)
}

ReadMetadata_Cain <- function(files) {
  metadata <- read.csv(files[["metadata"]]$path)
  metadata <- metadata[, c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster" )]

  noclust <- which(metadata$subcluster == "no.clus")
  metadata$subcluster[noclust] <- metadata$broadcelltype[noclust]

  return(metadata)
}

ReadCounts_Cain <- function(files) {
  load(files[["counts"]]$path)
  return(fullmat)
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

  synIDs <- list("metadata" = "syn38609700")

  files <- DownloadFromSynapse(synIDs, dir_lau_raw)

  files[["geo_metadata"]] = geo_metadata # TODO this is a data frame, not a file
  return(files)
}

ReadMetadata_Lau <- function(files) {
  metadata <- read.csv(files[["metadata"]]$path)
  metadata <- metadata[, c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster" )]

  # TODO Since this isn't the original metadata, I'm not sure about a few things:
  #   1) The number of cells in the unfiltered metadata file is less than what the paper states
  #   2) The paper doesn't classify cells as OPC but this metadata file does
  #   3) The "Exclude" category was added by someone in the deconvolution WG but
  #       I don't know what criteria determined this
  #   4) There are no fine cell types, maybe this data set needs to get mapped
  #       onto another one before use
  metadata <- subset(metadata, broadcelltype != "Exclude")

  return(metadata)
}

# TODO this function is a little hack-y.
ReadCounts_Lau <- function(files, metadata) {
  geo_metadata = files[["geo_metadata"]] # HACK: geo_metadata is a data frame, not a file
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


##### Leng, et al., 2021 #####
# https://doi.org/10.1038/s41593-020-00764-7

DownloadData_LengEC <- function() {
  synIDs <- list("counts" = "syn22722817")
  files <- DownloadFromSynapse(synIDs, dir_leng_raw)
  return(files)
}

DownloadData_LengSFG <- function() {
  synIDs <- list("counts" = "syn22722860")
  files <- DownloadFromSynapse(synIDs, dir_leng_raw)
  return(files)
}

ReadMetadata_Leng <- function(files) {
  sce <- readRDS(files[["counts"]]$path)
  metadata <- colData(sce)

  braak <- list("0" = "Control",
                "2" = "Early Pathology",
                "6" = "AD")

  metadata$cellid <- rownames(metadata)
  metadata$diagnosis <- unlist(braak[metadata$BraakStage])
  metadata$subcluster <- str_replace(metadata$clusterAssignment, "EC:|SFG:", "")
  metadata <- metadata[, c("cellid", "SampleID", "diagnosis", "clusterCellType", "subcluster")]

  return(metadata)
}

ReadCounts_Leng <- function(files) {
  sce <- readRDS(files[["counts"]]$path)
  return(counts(sce))
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

# This function is now unused, saving just in case
GeneSymbolToEnsembl_Mathys <- function(files) {
  # This file has mappings from Ensembl ID to gene symbol
  genes <- read.table(files[["genes_ensembl"]]$path, sep = "\t")
  colnames(genes) <- c("Ensembl.ID", "Symbol")

  dupes <- duplicated(genes$Symbol)
  genes <- genes[!dupes,]
  rownames(genes) <- genes$Symbol

  genes <- genes[,c("Symbol", "Ensembl.ID")]

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


##### SEA-AD #####
# Reference data set (5 donors): https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad
# Full data set (84 donors): https://portal.brain-map.org/explore/seattle-alzheimers-disease/seattle-alzheimers-disease-brain-cell-atlas-download?edit&language=en
# Metadata: https://www.synapse.org/#!Synapse:syn28256462

# These files are h5ad (AnnData) files. There are several R libraries that can
# read this type of file and output the full object with all fields populated
# (e.g. anndata and zellkonverter), however they can't seem to read in these
# particular files correctly, so I use a combination of reticulate and the
# H5ADMatrix library in R to get all the data.

DownloadData_SEARef <- function() {
  synIDs <- list("individual_metadata" = "syn31149116")
  files <- DownloadFromSynapse(synIDs, downloadLocation = dir_seaad_raw)

  files[["counts"]] = file_searef_h5

  if (!file.exists(files[["counts"]])) { # Don't re-download, this file is large
    download.file(url_searef_h5,
                  destfile = files[["counts"]], method = "curl")
  }
  return(files)
}

DownloadData_SEAAD <- function() {
  synIDs <- list("individual_metadata" = "syn31149116")
  files <- DownloadFromSynapse(synIDs, downloadLocation = dir_seaad_raw)

  files[["counts"]] = file_seaad_h5

  if (!file.exists(files[["counts"]])) { # Don't re-download, this file is large
    download.file(url_seaad_h5,
                  destfile = files[["counts"]], method = "curl")
  }
  return(files)
}

ReadMetadata_SEARef <- function(files) {
  donor_metadata <- read.csv(files[["individual_metadata"]]$path)

  ad <- import("anndata")
  adata <- ad$read_h5ad(files[["counts"]], backed = 'r')

  metadata <- adata$obs

  adata$file$close()
  rm(adata)

  metadata$broadcelltype <- as.character(metadata$subclass_label)
  metadata$broadcelltype[metadata$class_label == "Neuronal: GABAergic"] = "GABA"
  metadata$broadcelltype[metadata$class_label == "Neuronal: Glutamatergic"] = "Glut"

  metadata <- merge(metadata, donor_metadata,
                    by.x = "external_donor_name_label", by.y = "individualID")

  metadata <- metadata[, c("sample_name", "external_donor_name_label",
                           "diagnosis", "broadcelltype", "cluster_label" )]

  return(metadata)
}

ReadMetadata_SEAAD <- function(files) {
  donor_metadata <- read.csv(files[["individual_metadata"]]$path)

  ad <- import("anndata")
  adata <- ad$read_h5ad(files[["counts"]], backed = 'r')

  metadata <- adata$obs

  adata$file$close()
  rm(adata)

  metadata$broadcelltype <- as.character(metadata$Subclass)
  metadata$broadcelltype[metadata$Class == "Neuronal: GABAergic"] = "GABA"
  metadata$broadcelltype[metadata$Class == "Neuronal: Glutamatergic"] = "Glut"

  metadata <- merge(metadata, donor_metadata,
                    by.x = "Donor.ID", by.y = "individualID")

  metadata <- metadata[, c("sample_id", "Donor.ID", "diagnosis",
                           "broadcelltype", "Supertype")]
  return(metadata)
}

ReadCounts_SEARef <- function(files) {
  counts <- H5ADMatrix(files[["counts"]])

  col_names <- as.character(HDF5Array(files[["counts"]],
                                      file.path("obs", "sample_name")))
  dimnames(counts)[[2]] <- col_names

  return(counts)
}

ReadCounts_SEAAD <- function(files) {
  counts <- H5ADMatrix(files[["counts"]])

  col_names <- as.character(HDF5Array(files[["counts"]],
                                      file.path("obs", "sample_id")))
  dimnames(counts)[[2]] <- col_names

  return(counts)
}


##### Mayo #####
# Bulk RNA seq data from the Mayo RNA Seq Study:
# https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn5550404
#
# Metadata: https://www.synapse.org/#!Synapse:syn29855549
# Filtered counts: https://www.synapse.org/#!Synapse:syn27024951
# Biomart gene conversion: https://www.synapse.org/#!Synapse:syn27024953

DownloadData_Mayo <- function() {
  synIDs <- list("metadata" = "syn29855549",
                 "counts" = "syn27024951",
                 "genes" = "syn27024953")

  files <- DownloadFromSynapse(synIDs, dir_mayo_raw)
  return(files)
}


##### MSBB #####
# Bulk RNA seq data from the Mount Sinai Brain Bank Study:
# https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn3159438
#
# Metadata: https://www.synapse.org/#!Synapse:syn29855570
# Filtered counts: https://www.synapse.org/#!Synapse:syn27068754
# Biomart gene conversion: https://www.synapse.org/#!Synapse:syn27068755

DownloadData_MSBB <- function() {
  synIDs <- list("metadata" = "syn29855570",
                 "counts" = "syn27068754",
                 "genes" = "syn27068755")

  files <- DownloadFromSynapse(synIDs, dir_msbb_raw)
  return(files)
}


##### ROSMAP #####
# Bulk RNA seq data from the ROSMAP study:
# https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn3219045
#
# Metadata: https://www.synapse.org/#!Synapse:syn29855598
# Filtered counts: https://www.synapse.org/#!Synapse:syn26967451
# Biomart gene conversion: https://www.synapse.org/#!Synapse:syn26967452

DownloadData_ROSMAP <- function() {
  synIDs <- list("metadata" = "syn29855598",
                 "counts" = "syn26967451",
                 "genes" = "syn26967452")

  files <- DownloadFromSynapse(synIDs, dir_rosmap_raw)
  return(files)
}


##### Generic bulk functions #####
# These functions all work on Mayo, MSBB, and ROSMAP since the files all come
# from the harmonization effort

ReadMetadata_BulkData <- function(files) {
  metadata <- read.table(files[["metadata"]]$path, header = T)
  metadata <- metadata[,c("specimenID", "diagnosis", "tissue", "sex")]

  # Necessary because the column names of the counts matrix get converted this
  # way automatically
  metadata$specimenID <- make.names(metadata$specimenID)

  return(metadata)
}

ReadCounts_BulkData <- function(files) {
  counts <- read.table(files[["counts"]]$path, header = TRUE, row.names = 1)
  return(counts)
}
