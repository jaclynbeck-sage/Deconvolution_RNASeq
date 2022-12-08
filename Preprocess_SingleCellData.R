# Load metadata and counts matrices from Synapse and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
library(Matrix)
library(SingleCellExperiment)
library(edgeR)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito")

dataset <- "lengSFG"


##### Download files #####

files <- DownloadData(dataset)


##### Read in metadata file #####

metadata <- ReadMetadata(dataset, files)
colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
rownames(metadata) <- metadata$cellid

# Final modifications to metadata
metadata$broadcelltype <- factor(metadata$broadcelltype)
metadata$subcluster <- factor(metadata$subcluster)
metadata$donor <- factor(metadata$donor)
metadata$diagnosis <- factor(metadata$diagnosis)


##### Read in matrix of counts #####

counts <- ReadCounts(dataset, files)

# Make sure metadata is in the same order as counts
metadata <- metadata[colnames(counts), ]

# Convert gene symbol to Ensembl ID.
genes <- GeneSymbolToEnsembl(dataset, files, rownames(counts))
genes <- genes[rownames(counts),]


##### Create SingleCellExperiment object and save #####

if (!all(metadata$cellid == colnames(counts))) {
  stop("*** Cells in the metadata and counts matrix are in different orders.
       Double-check the data. ***")
}

# edgeR TMM normalization -- Needs a lot of memory
metadata$tmm.size.factors <- calcNormFactors(counts, method = "TMMwsp")

sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata, rowData = genes)

saveRDS(sce, file = file.path(dir_input, paste(dataset, "sce.rds", sep = "_")))

# For reference, using TMM factors:
#tmms <- counts %*% Diagonal(length(norm.factors), norm.factors)
#lib.norm <- colSums(counts) * metadata$tmm.size.factors
#tmm.norm <- counts %*% Diagonal(length(lib.norm), 1e6 / lib.norm)


