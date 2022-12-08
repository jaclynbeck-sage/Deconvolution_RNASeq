# Load metadata and counts matrices from GEO and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.

### Note: this data is unpublished, so this code is temporary until published
# data is on Synapse.

library(Matrix)
library(SingleCellExperiment)
library(edgeR)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))

##### Download files from Synapse #####

# TODO

files <- list("metadata" = file.path(dir_input, "cain_metadata.csv"),
              "counts" = file.path(dir_input, "cain_counts.rda"))


##### Read in metadata file #####

metadata <- ReadMetadata_Cain(files)
colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
rownames(metadata) <- metadata$cellid

# Final modifications to metadata
metadata$broadcelltype <- factor(metadata$broadcelltype)
metadata$subcluster <- factor(metadata$subcluster)
metadata$donor <- factor(metadata$donor)
metadata$diagnosis <- factor(metadata$diagnosis)


##### Read in matrix of counts #####

counts <- ReadCounts_Cain(files)

# Make sure metadata is in the same order as counts
metadata <- metadata[colnames(counts), ]

# Convert gene symbol to Ensembl ID.
genes <- GeneSymbolToEnsembl_Cain(rownames(counts))
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

saveRDS(sce, file = file.path(dir_input, "cain_sce.rds"))

