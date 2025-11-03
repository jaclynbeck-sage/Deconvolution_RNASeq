# This script downloads and processes two pieces of information we need in order
# to deconvolve data sets:
#   1. A list of Ensembl ID / gene symbol pairings. The single-cell data sets
#      all use gene symbols while the bulk data sets use Ensembl IDs, so we
#      need this list for conversion. We also need exon lengths of genes from
#      the bulk data sets for calculating TPM normalization.
#   2. Ground-truth cell type proportions for select ROSMAP samples, as
#      determined by IHC on the same tissue as the RNA sequencing samples.
#
# For the gene list, we pull multiple sources of data to make sure we have as
# much coverage as possible:
#   1. The genes used for ROSMAP/Mayo/MSBB in the RNASeq Harmonization Study,
#   2. The genes used for the Seattle Reference Atlas,
#   3. All genes from Ensembl version 98 corresponding to the genome used for
#      the Cain data set,
#   4. The genes used for the Lau, Leng, and Mathys data sets, which are
#      provided by the studies in the raw data
# Some symbols from previous versions of Ensembl are different from the symbols
# in Gencode v43. We merge the multiple lists by Ensembl ID and create a list of
# all possible gene symbols associated with each ID. We then assign a
# "canonical" symbol that will be used in every data set by using the symbol
# from Gencode v43 if it exists, otherwise by finding the most common symbol for
# a given Ensembl ID in the list of symbols from the other data sets.
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
library(synapser)
library(dplyr)
library(purrr)
library(rtracklayer)
library(sageRNAUtils)
library(GEOquery)
library(rhdf5)
source("Filenames.R")

cfg <- config::get("step01_gene_metadata")

synLogin()

dir_gene_files <- file.path(dir_metadata, "gene_files")
dir.create(dir_gene_files, showWarnings = FALSE)

## Bulk RNA Seq genes ----------------------------------------------------------

# Get exon lengths for each gene, for the purpose of calculating TPM on the
# bulk datasets. The bulk datasets were aligned to Gencode release 43.

download.file(cfg$gtf_bulk,
              destfile = file.path(dir_gene_files, basename(cfg$gtf_bulk)),
              method = "curl")
download.file(cfg$fasta_bulk,
              destfile = file.path(dir_gene_files, basename(cfg$fasta_bulk)),
              method = "curl")

# We don't strictly need GC content but this function gets us gene length,
# gene biotype, and symbols
gene_info <- sageRNAUtils::get_gc_content_gtf(
  gtf_file = file.path(dir_gene_files, basename(cfg$gtf_bulk)),
  fasta_file = file.path(dir_gene_files, basename(cfg$fasta_bulk)),
  include_introns = FALSE
)

# Strip version numbers
gene_info$ensembl_gene_id <- str_replace(gene_info$ensembl_gene_id, "\\.[0-9]+", "")

gene_info <- gene_info |>
  dplyr::rename(symbol_bulkRNA = external_gene_name) |>
  # Remove IDs that end in "_PAR_Y"
  subset(!grepl("_PAR_Y", ensembl_gene_id))


## Lau genes -------------------------------------------------------------------

# Download one of the features.tsv files from GEO to get the gene information

res <- getGEOSuppFiles(cfg$geo_lau, makeDirectory = FALSE,
                       baseDir = dir_gene_files,
                       filter_regex = "features")

lau_genes <- read.delim(gzfile(rownames(res)), header = FALSE) |>
  dplyr::select(V1, V2) |>
  dplyr::rename(ensembl_gene_id = V1, symbol_Lau = V2)


## Mathys genes ----------------------------------------------------------------

# They provided their mapping file on Synapse
filename <- synGet(cfg$genes_mathys, version = 1,
                   downloadLocation = dir_gene_files)

mathys_genes <- read.table(filename$path, header = FALSE) |>
  dplyr::rename(ensembl_gene_id = V1, symbol_Mathys = V2)


## Cain and SEA-AD genes -------------------------------------------------------

# GRCh38 Ensembl release 98 (Cain) plus the seaRef GTF file
files <- list("symbol_seaRef" = c(filename = file.path(dir_gene_files, "seaRef_genes.gtf.gz"),
                                  url = cfg$gtf_sea_ad),
              "symbol_v98_cain" = c(filename = file.path(dir_gene_files, "Homo_sapiens.GRCh38.98.gtf.gz"),
                                    url = cfg$gtf_cain))

gtf_genes <- lapply(names(files), function(version) {
  file_info <- files[[version]]
  if (!file.exists(file_info["filename"])) {
    download.file(file_info["url"],
                  destfile = file_info["filename"],
                  method = "curl")
  }

  df <- rtracklayer::import(file_info[["filename"]], format = "gtf") |>
    as.data.frame() |>
    dplyr::select(gene_id, gene_name) |>
    dplyr::distinct()

  colnames(df) <- c("ensembl_gene_id", version)
  return(df)
})

gtf_genes <- purrr::reduce(gtf_genes, dplyr::full_join, by = "ensembl_gene_id")


# Merge all gene sets together -------------------------------------------------

all_genes <- purrr::reduce(list(gene_info, gtf_genes, lau_genes, mathys_genes),
                           dplyr::full_join,
                           by = "ensembl_gene_id") |>
  # Remove some symbols that got set to the Ensembl ID
  mutate(symbol_bulkRNA = ifelse(grepl("ENSG00", symbol_bulkRNA), NA, symbol_bulkRNA))

# For each row/gene that doesn't have a symbol in the bulk RNA seq data, take
# the symbol that appears the most across all data sets. In case of ties,
# which.max returns the first symbol in the tie, which would be the first
# alphabetically.
max_symbol <- all_genes |>
  subset(is.na(symbol_bulkRNA)) |>
  tidyr::pivot_longer(cols = starts_with("symbol"),
                      names_to = "dataset",
                      values_to = "symbol",
                      values_drop_na = TRUE) |>
  group_by(ensembl_gene_id, symbol) |>
  summarize(count = n(), .groups = "drop_last") |>
  summarize(canonical_symbol = symbol[which.max(count)])

# Use bulk RNA seq symbols where possible and concat the most-used symbols from
# the other data sets when bulk RNA symbols are NA
final_symbol <- all_genes |>
  subset(!is.na(symbol_bulkRNA)) |>
  mutate(canonical_symbol = symbol_bulkRNA) |>
  dplyr::select(ensembl_gene_id, canonical_symbol) |>
  rbind(max_symbol)

# Merge back in to the main data frame
all_genes <- merge(all_genes, final_symbol, all = TRUE) |>
  dplyr::arrange(ensembl_gene_id)

all_genes <- subset(all_genes, !is.na(canonical_symbol))
write.csv(all_genes, file_gene_list, quote = FALSE, row.names = FALSE)


# Ground-truth proportions for some ROSMAP samples, generated from IHC ---------

dir_ihc_git <- file.path(dir_metadata, "CortexCellDeconv")
system(paste("git clone", cfg$git_ihc, dir_ihc_git))

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
props_df <- dplyr::select(props_df, -donor)

colnames(props_df) <- c("Astro", "Endo", "Micro", "Neuro", "Oligo")

# Save unaltered data to a file for reference
write.csv(props_df,
          file = file.path(dir_metadata, "ihc_proportions_unnormalized.csv"),
          quote = FALSE,
          row.names = TRUE)

# Remove donors with missing measurements
good <- rowSums(is.na(props_df)) == 0
props_df <- props_df[good, ]

# Adjust all rows to sum to 1
props_df <- sweep(props_df, 1, rowSums(props_df), "/")

# Write normalized proportions file
write.csv(props_df,
          file = file_rosmap_ihc_proportions,
          quote = FALSE,
          row.names = TRUE)
