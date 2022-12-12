library(synapser)
library(biomaRt)
library(Seurat)
library(reshape2)
library(stringr)
library(GEOquery)
library(HDF5Array)
library(rhdf5)


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
                  "lengEC" = DownloadData_LengEC(),
                  "lengSFG" = DownloadData_LengSFG(),
                  "mathys" = DownloadData_Mathys(),
                  "morabito" = DownloadData_Morabito(),
                  "seaRef" = DownloadData_SEARef(),
                  "seaAD" = DownloadData_SEAAD())
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
                     "seaAD" = ReadMetadata_SEAAD(files))
  return(metadata)
}

ReadCounts <- function(dataset, files) {
  counts <- switch(dataset,
                   "cain" = ReadCounts_Cain(files),
                   "lau" = NULL, # Needs more work
                   "lengEC" = ReadCounts_Leng(files),
                   "lengSFG" = ReadCounts_Leng(files),
                   "mathys" = ReadCounts_Mathys(files),
                   "morabito" = ReadCounts_Morabito(files),
                   "seaRef" = ReadCounts_SEARef(files),
                   "seaAD" = ReadCounts_SEAAD(files))
  return(counts)
}

GeneSymbolToEnsembl <- function(dataset, files = NULL, gene.list = NULL) {
  genes <- switch(dataset,
                  "cain" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "lau" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "lengEC" = GeneSymbolToEnsembl_Biomart(gene.list), # TODO they provide the GTF with ID->Symbol mappings
                  "lengSFG" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "mathys" = GeneSymbolToEnsembl_Mathys(files),
                  "morabito" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "seaRef" = GeneSymbolToEnsembl_Biomart(gene.list),
                  "seaAD" = GeneSymbolToEnsembl_Biomart(gene.list))
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
  # TODO these are not the original data files
  synIDs <- list("metadata" = "syn38609692",
                 "counts" = "syn38598183")

  files <- DownloadFromSynapse(synIDs, dir_cain_raw)
  return(files)
}

ReadMetadata_Cain <- function(files) {
  metadata <- read.csv(files[["metadata"]]$path)
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

  files <- list("metadata" = file.path(dir_input, "lau_metadata.csv"),
                "geo_metadata" = geo_metadata) # TODO this is a data frame, not a file
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


##### SEA-AD: Not finished yet #####
# Reference data set (5 donors): https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad
# Full data set (84 donors): https://portal.brain-map.org/explore/seattle-alzheimers-disease/seattle-alzheimers-disease-brain-cell-atlas-download?edit&language=en

# These files are h5ad (AnnData) files. There are several R libraries that can
# read this type of file and output the full object with all fields populated
# (e.g. anndata and zellkonverter), however most of them use reticulate, which
# is incompatible with synapser/PythonEmbedInR, so we can't use them in the
# same R session which makes it difficult to pipeline. So I read pieces by
# hand.

# Note: This function can't read fields that are ENUM type
ReadH5Group <- function(filename, group.name, cols = NULL) {
  if (is.null(cols)) {
    attr <- h5readAttributes(filename, group.name)
    cols <- c(attr[["_index"]], attr[["column-order"]])
  }

  col_list <- list()
  for (column in cols) {
    if (column != "__categories") {
      col_list[[column]] <- HDF5Array(filename, file.path(group.name, column))
    }
  }

  return(col_list)
}

# Note: This function doesn't work with fields that have -1 values
ReplaceColWithCategory <- function(metadata, categories, column.names) {
  if (is.null(column.names)) {
    column.names <- names(categories)
  }

  for (column in column.names) {
    column.safe <- make.names(column) # Replaces invalid characters
    metadata[, column.safe] <- categories[[column]][metadata[, column.safe] + 1]
  }

  return(metadata)
}

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

# Testing calling python from R. Ended up not using this code but saving it
# just in case.
ReadMetadata_SEARef_old <- function(files) {
  pytext <- paste0('import scanpy as sc \n',
                   'data = sc.read_h5ad("', file_searef_h5, '", backed = "r") \n',
                   'metadata = sc.get.obs_df(data, keys = ["sample_name", ',
                   '"external_donor_name_label", "class_label", ',
                   '"subclass_label", "cluster_label"]) \n')
  pyExec(pytext)
  metadata <- pyGet("metadata")

  donor_metadata <- read.csv(files[["individual_metadata"]]$path)

  orig_order <- metadata$sample_name
  metadata <- merge(metadata, donor_metadata,
                    by.x = "external_donor_name_label", by.y = "individualID")

  rownames(metadata) <- metadata$sample_name
  metadata <- metadata[orig_order,]

  metadata$broadcelltype <- metadata$subclass_label
  metadata$broadcelltype[metadata$class_label == "Neuronal: GABAergic"] = "GABA"
  metadata$broadcelltype[metadata$class_label == "Neuronal: Glutamatergic"] = "Glut"

  metadata <- metadata[, c("sample_name", "external_donor_name_label",
                           "diagnosis", "broadcelltype", "cluster_label" )]

  return(metadata)
}

ReadMetadata_SEARef <- function(files) {
  donor_metadata <- read.csv(files[["individual_metadata"]]$path)

  # SEA-Ref specific columns
  cols1 <- c("sample_name", "external_donor_name_label", "class_label",
            "subclass_label", "cluster_label")
  cols2 <- c("class_label", "cluster_label", "external_donor_name_label",
             "subclass_label")

  metadata <- data.frame(ReadH5Group(files[["counts"]], "obs", cols1))
  categories <- ReadH5Group(files[["counts"]], file.path("obs", "__categories"), cols2)

  metadata <- ReplaceColWithCategory(metadata, categories, cols2)

  metadata$broadcelltype <- metadata$subclass_label
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

  # SEA-AD specific columns
  cols1 <- c("sample_id", "Donor ID", "Supertype", "Class", "Subclass")
  cols2 <- c("Donor ID", "Supertype", "Class", "Subclass")

  metadata <- data.frame(ReadH5Group(files[["counts"]], "obs", cols1))
  categories <- ReadH5Group(files[["counts"]], file.path("obs", "__categories"), cols2)

  metadata <- ReplaceColWithCategory(metadata, categories, cols2)

  metadata$broadcelltype <- metadata$Subclass
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
