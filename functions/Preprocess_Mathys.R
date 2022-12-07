# Load metadata and counts matrices from Synapse and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
library(Matrix)
library(SingleCellExperiment)
library(reshape2)
library(edgeR)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))

##### Download files from Synapse #####

synIDs <- list("clinical" = "syn3191087",
               "metadata" = "syn18686372",
               "genes_ensembl" = "syn18687959",
               "counts" = "syn18686381",
               "genes" = "syn18686382")

files <- DownloadFromSynapse(synIDs, dir_mathys_raw)


##### Read in metadata file #####

metadata <- ReadMetadata_Mathys(files)
colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")

# Final modifications to metadata
metadata$broadcelltype <- factor(metadata$broadcelltype)
metadata$subcluster <- factor(metadata$subcluster)
metadata$donor <- factor(metadata$donor)
metadata$diagnosis <- factor(metadata$diagnosis)


##### Read in matrix of counts #####

counts <- ReadCounts_Mathys(files)

# Make sure metadata is in the same order as counts
metadata <- metadata[colnames(counts), ]

# Convert gene symbol to Ensembl ID.
genes <- GeneSymbolToEnsembl_Mathys(files)
genes <- genes[rownames(counts),]
#rownames(counts) <- genes[rownames(counts), "Ensembl.ID"]
#rownames(genes) <- genes$Ensembl.ID


##### Create SingleCellExperiment object and save #####

if (!all(metadata$cellid == colnames(counts))) {
  stop("*** Cells in the metadata and counts matrix are in different orders.
       Double-check the data. ***")
}

# edgeR TMM normalization
metadata$tmm.size.factors <- calcNormFactors(counts, method = "TMMwsp")
#tmms <- counts %*% Diagonal(length(norm.factors), norm.factors)
#lib.norm <- colSums(counts) * metadata$tmm.size.factors
#tmm.norm <- counts %*% Diagonal(length(lib.norm), 1e6 / lib.norm)

sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata, rowData = genes)

saveRDS(sce, file = file.path(dir_input, "mathys_sce.rds"))

