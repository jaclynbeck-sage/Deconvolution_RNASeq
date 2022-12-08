# Load metadata and counts matrices from GEO and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.

### Note: the metadata for this file was provided by the author, and is not
# publicly available. Some of this code is temporary until a downloadable
# metadata file is available.

library(Matrix)
library(SingleCellExperiment)
library(edgeR)
library(GEOquery)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))

##### Download files from GEO #####

gse <- getGEO("GSE157827", destdir = dir_lau_raw)
geo_metadata <- pData(phenoData(gse[[1]]))

geo <- getGEOSuppFiles(GEO = "GSE157827", makeDirectory = FALSE,
                       baseDir = dir_lau_raw)

#untar(rownames(geo)[1], list = TRUE) # List files inside
untar(rownames(geo)[1], exdir = dir_lau_raw)


##### Read in metadata file #####

# Temporary -- use metadata file provided to me
files <- list("metadata" = file.path(dir_input, "lau_metadata.csv"))
metadata <- ReadMetadata_Lau(files)

colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
rownames(metadata) <- metadata$cellid

# Final modifications to metadata
metadata$broadcelltype <- factor(metadata$broadcelltype)
metadata$subcluster <- factor(metadata$subcluster)
metadata$donor <- factor(metadata$donor)
metadata$diagnosis <- factor(metadata$diagnosis)


##### Read in matrix of counts #####

counts <- ReadCounts_Lau(geo_metadata, metadata)

# Make sure metadata is in the same order as counts
metadata <- metadata[colnames(counts), ]

# Convert gene symbol to Ensembl ID.
genes <- GeneSymbolToEnsembl_Lau(rownames(counts))
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

saveRDS(sce, file = file.path(dir_input, "lau_sce.rds"))

