# This script downloads and processes two pieces of information we need in order
# to deconvolve data sets:
#   1. A list of Ensembl ID / gene symbol pairings. The single-cell data sets
#      all use gene symbols while the bulk data sets use Ensembl IDs, so we
#      need this list for conversion.
#   2. Ground-truth cell type proportions for select ROSMAP samples, as
#      determined by IHC on the same tissue as the RNA sequencing samples.
#
# For the gene list, we pull 3 sources of data to make sure we have as much
# coverage as possible: All genes in the current Biomart database, the genes
# used for ROSMAP/Mayo/MSBB in the RNASeq Harmonization Study, and the genes
# used for the Mathys data set. The latter two sources contain genes from
# previous versions of Ensembl, so some symbols are different and some
# Ensembl IDs exist in those lists that are no longer in Biomart. We merge
# the three lists by Ensembl ID and create a list of all possible gene symbols
# associated with each ID.
# NOTE: Mathys is the only single cell data set to provide the mapping between
# their gene symbols and Ensembl IDs, so the other single cell data sets may
# still contain genes not in the created gene list, or have symbols that don't
# map to the correct Ensembl ID used by their data processing pipeline.
#
# For the cell type proportions, the data was determined from separate stains,
# so the proportions do not add up to 1 and some cell types are missing from
# certain donors. We account for this by:
#   a) Removing all donors from the list who are missing at least one cell type
#   b) Normalizing the proportions for each sample to sum to 1.
# NOTE: this is percent of cells, not percent of RNA. There needs to be a
# conversion between percent cells -> percent RNA for deconvolution algorithms
# that output percent RNA. This conversion is handled downstream.

library(stringr)
library(biomaRt)
library(synapser)
library(dplyr)
source("Filenames.R")

synLogin()

##### Gene symbol / Ensembl ID conversions #####

# Biomart query to get all genes in the database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomart_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       mart = mart)
colnames(biomart_genes) <- c("symbol_Biomart", "ensembl_gene_id")

# Gene conversions used for ROSMAP, Mayo, and MSBB. The files from all three
# studies (syn27024953, syn27068755, syn26967452) are identical so we just use
# the one from ROSMAP.
dir.create(file.path(dir_metadata, "gene_files"), showWarnings = FALSE)

filename <- synGet("syn26967452",
                   downloadLocation = file.path(dir_metadata, "gene_files"))

ros_genes <- read.table(filename$path, header = TRUE) %>%
                select(ensembl_gene_id, hgnc_symbol) %>%
                rename(symbol_RNASeq = hgnc_symbol)

# Mathys genes
filename <- synGet("syn18687959",
                   downloadLocation = file.path(dir_metadata, "gene_files"))

mathys_genes <- read.table(filename$path, header = FALSE) %>%
                    rename(ensembl_gene_id = V1, symbol_Mathys = V2)

all_genes <- merge(biomart_genes, ros_genes,
                   by = "ensembl_gene_id", all = TRUE) %>%
             merge(mathys_genes, by = "ensembl_gene_id", all = TRUE)

# Combine the separate symbol fields into a string containing all the symbols
# for each gene.
combine_symbols <- function(...) {
  return(paste(unique( na.omit(c(...)) ), collapse = "|"))
}

all_genes <- all_genes %>% rowwise() %>%
                mutate(hgnc_symbols = combine_symbols(symbol_Biomart, symbol_RNASeq, symbol_Mathys)) %>%
                arrange(ensembl_gene_id) %>% as.data.frame()

write.csv(all_genes, file_gene_list, quote = FALSE, row.names = FALSE)


##### Ground-truth proportions for some ROSMAP samples, generated from IHC #####

dir_ihc_git <- file.path(dir_metadata, "CortexCellDeconv")
system(paste0("git clone https://github.com/ellispatrick/CortexCellDeconv.git ",
              dir_ihc_git))

celltypes <- c("astro", "endo", "microglia", "neuro", "oligo")

# Each cell type has its own .txt file. Read each file in and convert it
# to a data frame.
props <- lapply(celltypes, function(ct) {
  ct_file <- file.path(dir_ihc_git, "CellTypeDeconvAnalysis", "Data",
                       str_glue("IHC.{ct}.txt"))
  ct_props <- read.table(ct_file, header = TRUE, sep = "\t")
  ct_props <- as.data.frame(t(ct_props))
  colnames(ct_props) <- ct
  ct_props$donor <- rownames(ct_props)
  ct_props
})

# Merge the separate cell type data frames into one data frame, by donor
props_df <- props[[1]]
for (i in 2:length(props)) {
  props_df <- merge(props_df, props[[i]], by = "donor", all = TRUE)
}

rownames(props_df) <- str_replace(props_df$donor, "X", "")
props_df <- select(props_df, -donor)

# Save unaltered data to a file for reference
write.csv(props_df, file = file.path(dir_metadata, "ihc_proportions_unnormalized.csv"),
          quote = FALSE, row.names = TRUE)

# Remove donors with missing measurements
good <- rowSums(is.na(props_df)) == 0
props_df <- props_df[good,]

# Adjust all rows to sum to 1
props_df <- sweep(props_df, 1, rowSums(props_df), "/")

colnames(props_df) <- c("Astro", "Endo", "Micro", "Neuro", "Oligo")

# Write normalized proportions file
write.csv(props_df, file = file_rosmap_ihc_proportions,
          quote = FALSE, row.names = TRUE)
