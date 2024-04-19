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
#   1. All genes in the current Biomart database,
#   2. The genes used for ROSMAP/Mayo/MSBB in the RNASeq Harmonization Study,
#   3. The genes used for the Seattle Reference Atlas,
#   4. All genes from Ensembl versions 98, 93, and 84, corresponding to the
#      versions used in some of the single cell data, and
#   5. The genes used for the Mathys data set
# Some symbols from previous versions of Ensembl are different and some Ensembl
# IDs exist in those lists that are no longer in Biomart. We merge the multiple
# lists by Ensembl ID and create a list of all possible gene symbols associated
# with each ID. We then assign a "canonical" symbol that will be used in every
# data set by finding the first non-null symbol for a given Ensembl ID in the
# list of symbols (ordered by source, as above).
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
library(GenomicFeatures)
library(GenomicRanges)
library(synapser)
library(dplyr)
library(purrr)
library(rtracklayer)
source("Filenames.R")

synLogin()

# URLs for GTF files -----------------------------------------------------------

gtf_seaRef <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/67/39/67390730-a684-47a5-b9f4-89c47cd4e3fc/genesgtf.gz"
gtf_v98 <- "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz"
gtf_v93 <- "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz"
gtf_v84 <- "https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"


# Gene symbol / Ensembl ID conversions -----------------------------------------

# Biomart query to get all genes in the database
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                   version = "110")
biomart_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       mart = mart)

colnames(biomart_genes) <- c("symbol_Biomart", "ensembl_gene_id")
biomart_genes <- subset(biomart_genes, symbol_Biomart != "")

dir_gene_files <- file.path(dir_metadata, "gene_files")
dir.create(dir_gene_files, showWarnings = FALSE)

# Gene conversions used for ROSMAP, Mayo, and MSBB. The files from all three
# studies (syn27024953, syn27068755, syn26967452) are identical so we just use
# the one from ROSMAP.
filename <- synGet("syn26967452", downloadLocation = dir_gene_files)

ros_genes <- read.table(filename$path, header = TRUE) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>%
  dplyr::rename(symbol_RNASeq = hgnc_symbol)

# Mathys genes -- this is the only single cell data set that provides their own
# mapping from gene symbol to Ensembl ID
filename <- synGet("syn18687959", version = 1,
                   downloadLocation = dir_gene_files)

mathys_genes <- read.table(filename$path, header = FALSE) %>%
  dplyr::rename(ensembl_gene_id = V1, symbol_Mathys = V2)

# GRCh38 Ensembl releases 84, 93, and 98, plus the seaRef GTF file
files <- list("symbol_seaRef" = c(filename = file.path(dir_gene_files, "seaRef_genes.gtf.gz"),
                                  url = gtf_seaRef),
              "symbol_v98" = c(filename = file.path(dir_gene_files, "Homo_sapiens.GRCh38.98.gtf.gz"),
                               url = gtf_v98),
              "symbol_v93" = c(filename = file.path(dir_gene_files, "Homo_sapiens.GRCh38.93.gtf.gz"),
                               url = gtf_v93),
              "symbol_v84" = c(filename = file.path(dir_gene_files, "Homo_sapiens.GRCh38.84.gtf.gz"),
                               url = gtf_v84))

gtf_genes <- lapply(names(files), function(version) {
  file_info <- files[[version]]
  if (!file.exists(file_info["filename"])) {
    download.file(file_info["url"],
                  destfile = file_info["filename"],
                  method = "curl")
  }

  df <- rtracklayer::import(gzfile(file_info[["filename"]]), format = "gtf") %>%
    as.data.frame() %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct()

  colnames(df) <- c("ensembl_gene_id", version)
  return(df)
})

gtf_genes <- purrr::reduce(gtf_genes, full_join, by = "ensembl_gene_id")


# Merge all gene sets together -------------------------------------------------

all_genes <- merge(biomart_genes, ros_genes,
                   by = "ensembl_gene_id",
                   all = TRUE) %>%
  merge(gtf_genes, by = "ensembl_gene_id", all = TRUE) %>%
  merge(mathys_genes, by = "ensembl_gene_id", all = TRUE)

# For each row/gene, take the first non-NA symbol. The symbol columns are
# ordered left to right in order of priority (1 - Biomart, 2 - RNASeq,
# 3 - seaRef genes, 4 - GTF genes in descending version order, 5 - Mathys), so
# the first entry from c_across with NAs removed will be the highest-priority
# non-NA gene symbol.
first_symbol <- function(...) {
  vec <- c_across(starts_with("symbol_"))
  return(na.omit(vec)[1])
}

all_genes <- all_genes %>%
  rowwise() %>%
  dplyr::mutate(canonical_symbol = first_symbol()) %>%
  dplyr::arrange(ensembl_gene_id) %>%
  as.data.frame()

# Get exon lengths for each gene, for the purpose of calculating TPM on the
# bulk datasets. The bulk datasets were aligned to Ensembl release 97.
tx <- GenomicFeatures::makeTxDbFromEnsembl(organism = "Homo sapiens",
                                           release = 97)
ex <- GenomicFeatures::exonsBy(tx, by = "gene")
ex <- GenomicRanges::reduce(ex)
exlen <- relist(width(unlist(ex)), ex)
exlens <- sapply(exlen, sum)
exlens <- data.frame(ensembl_gene_id = names(exlens), exon_length = exlens)

all_genes <- merge(all_genes, exlens, by = "ensembl_gene_id", all = TRUE)

all_genes <- subset(all_genes, !is.na(canonical_symbol))
write.csv(all_genes, file_gene_list, quote = FALSE, row.names = FALSE)


# Ground-truth proportions for some ROSMAP samples, generated from IHC ---------

dir_ihc_git <- file.path(dir_metadata, "CortexCellDeconv")
system(paste("git clone https://github.com/ellispatrick/CortexCellDeconv.git",
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
