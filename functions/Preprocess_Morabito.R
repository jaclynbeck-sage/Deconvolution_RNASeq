# Load metadata and counts matrices from Synapse and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
library(Matrix)
library(SingleCellExperiment)
library(edgeR)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))


##### Download files from Synapse #####

synIDs <- list("metadata" = "syn22130806",
               "barcodes" = "syn24978676",
               "genes" = "syn24978737",
               "counts" = "syn22130805")

files <- DownloadFromSynapse(synIDs, dir_morabito_raw)


##### Read in metadata file #####

metadata <- ReadMetadata_Morabito(files)
colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
rownames(metadata) <- metadata$cellid

# Final modifications to metadata
metadata$broadcelltype <- factor(metadata$broadcelltype)
metadata$subcluster <- factor(metadata$subcluster)
metadata$donor <- factor(metadata$donor)
metadata$diagnosis <- factor(metadata$diagnosis)


##### Read in matrix of counts #####

counts <- ReadCounts_Morabito(files)

# Make sure metadata is in the same order as counts
metadata <- metadata[colnames(counts), ]

# Convert gene symbol to Ensembl ID.
genes <- GeneSymbolToEnsembl_Morabito(rownames(counts))
genes <- genes[rownames(counts),]
#rownames(counts) <- genes[rownames(counts), "Ensembl.ID"] # TODO: Can't do this. Not all genes have an Ensembl ID
#rownames(genes) <- genes$Ensembl.ID # See above


##### Create SingleCellExperiment object and save #####

if (!all(metadata$cellid == colnames(counts))) {
  stop("*** Cells in the metadata and counts matrix are in different orders.
       Double-check the data. ***")
}

# edgeR TMM normalization -- Needs a lot of memory
metadata$tmm.size.factors <- calcNormFactors(counts, method = "TMMwsp")

sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata, rowData = genes)

saveRDS(sce, file = file.path(dir_input, "morabito_sce.rds"))

