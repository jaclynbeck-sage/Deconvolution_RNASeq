# Load metadata and counts matrices from Synapse and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
#
# Note: Right now this is just for the Mathys dataset, which needs a lot of
# metadata rearrangement.
library(Matrix)
library(SingleCellExperiment)
library(reshape2)

source("Filenames.R")

##### Read in metadata files #####

# syn3191087
clinical <- read.csv(file.path(dir_mathys_raw, "ROSMAP_clinical.csv"))

# syn18686372
cellMeta <- read.table(file.path(dir_mathys_raw, "filtered_column_metadata.txt"),
                       sep = "\t", header = TRUE)

# syn18687959 -- This file has mappings from Ensembl ID to gene symbol
genes <- read.table(file.path(dir_mathys_raw, "notfiltered_gene_row_names.txt"),
                    sep = "\t")
colnames(genes) <- c("Ensembl.ID", "Symbol")


##### Create metadata data frame #####

metadata <- merge(cellMeta, clinical, by = "projid", all.x = TRUE, all.y = FALSE)

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

# Remove unneeded columns
metadata <- metadata[, c("TAG", "projid", "broad.cell.type", "Subcluster", "diagnosis")]

colnames(metadata) <- c("cell.id", "donor", "broad.cell.type", "fine.cell.type", "diagnosis")

# Make sure metadata is in the same order as the counts matrix will be
rownames(metadata) <- metadata$cell.id
metadata <- metadata[cellMeta$TAG,]

# Obfuscate donor ID for privacy -- only someone with access to the original
# files and this script would be able to recover donor -> projid
projids <- as.character(unique(metadata$donor))
mapping <- paste0("donor", 1:length(projids))
names(mapping) <- projids

metadata$donor <- mapping[as.character(metadata$donor)]

# Final modifications to metadata
metadata$broad.cell.type <- factor(metadata$broad.cell.type)
metadata$fine.cell.type <- factor(metadata$fine.cell.type)
metadata$donor <- factor(metadata$donor)
metadata$diagnosis <- factor(metadata$diagnosis)


##### Read in matrix of counts #####

# syn18686381
counts <- readMM(file.path(dir_mathys_raw, "filtered_count_matrix.mtx"))

# syn18686382
rows <- read.table(file.path(dir_mathys_raw, "filtered_gene_row_names.txt"))
rownames(counts) <- rows$V1
colnames(counts) <- cellMeta$TAG

# Make sure counts is in the same order as metadata
counts <- counts[, metadata$cell.id]

# Convert gene symbol to Ensembl ID. Symbols can map to multiple Ensembl IDs.
# For this application, we use the first Ensembl ID in the list for duplicate
# symbols, which is what Seurat does.
dupes <- duplicated(genes$Symbol)
genes <- genes[!dupes,]
rownames(genes) <- genes$Symbol
genes <- genes[rownames(counts),]
rownames(counts) <- genes[rownames(counts), "Ensembl.ID"]


##### Create SingleCellExperiment object and save #####

rownames(genes) <- genes$Ensembl.ID

sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata, rowData = genes)

saveRDS(sce, file = file.path(dir_input, "mathys_sce.rds"))

