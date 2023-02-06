# Load metadata and counts matrices from Synapse and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
library(Matrix)
library(SingleCellExperiment)
library(edgeR)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef", "seaAD")

dataset <- "lau"


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

counts <- ReadCounts(dataset, files, metadata)

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

# If the counts matrix is a DelayedArray (i.e. from the SEA-AD files), the
# sce file will contain a pointer to the original data file rather than
# writing the full data to disk again
saveRDS(sce, file = file.path(dir_input, paste(dataset, "sce.rds", sep = "_")))

# Calculate the "A" matrix that is needed to convert propCells to pctRNA
A_broad <- CalculateA(sce, metadata$donor, metadata$broadcelltype)
A_fine <- CalculateA(sce, metadata$donor, metadata$subcluster)

saveRDS(list("A_broad" = A_broad, "A_fine" = A_fine),
        file = file.path(dir_input, paste(dataset, "A_matrix.rds", sep = "_")))

# TODO Calculate a signature for each cell type



# For reference, using TMM factors:
#lib.norm <- colSums(counts) * metadata$tmm.size.factors
#tmm.norm <- counts %*% Diagonal(length(lib.norm), 1e6 / lib.norm)


