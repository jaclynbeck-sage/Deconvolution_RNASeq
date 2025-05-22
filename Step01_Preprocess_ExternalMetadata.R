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
#   3. All genes from Ensembl versions 98, 93, and 84, corresponding to the
#      versions used in some of the single cell data, and
#   4. The genes used for the Mathys data set
# Some symbols from previous versions of Ensembl are different and some Ensembl
# IDs exist in those lists that are no longer in Biomart. We merge the multiple
# lists by Ensembl ID and create a list of all possible gene symbols associated
# with each ID. We then assign a "canonical" symbol that will be used in every
# data set by finding the most common symbol for a given Ensembl ID in the list
# of symbols.
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
source("Filenames.R")

synLogin()

# URLs for GTF files -----------------------------------------------------------

gtf_seaRef <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/67/39/67390730-a684-47a5-b9f4-89c47cd4e3fc/genesgtf.gz"
gtf_v98 <- "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz"
gtf_v93 <- "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz"
gtf_v84 <- "https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"


# Gene symbol / Ensembl ID conversions -----------------------------------------

dir_gene_files <- file.path(dir_metadata, "gene_files")
dir.create(dir_gene_files, showWarnings = FALSE)

## Bulk RNA Seq genes ----------------------------------------------------------

# Get exon lengths for each gene, for the purpose of calculating TPM on the
# bulk datasets. The bulk datasets were aligned to Gencode release 31.
gtf_url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
                  "release_31/gencode.v31.primary_assembly.annotation.gtf.gz")
fasta_url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
                    "release_31/GRCh38.primary_assembly.genome.fa.gz")

# We don't strictly need GC content but this function gets us gene length,
# gene biotype, and symbols
gene_info <- sageRNAUtils::get_gc_content_gtf(gtf_url, fasta_url,
                                              include_introns = FALSE)

# Strip version numbers
gene_info$ensembl_gene_id <- str_replace(gene_info$ensembl_gene_id, "\\.[0-9]+", "")

gene_info <- gene_info |>
  dplyr::rename(symbol_RNASeq = external_gene_name) |>
  # Remove IDs that end in "_PAR_Y"
  subset(!grepl("_PAR_Y", ensembl_gene_id))


## Mathys genes ----------------------------------------------------------------

# This is the only single cell data set that provides their own mapping from
# gene symbol to Ensembl ID
filename <- synGet("syn18687959", version = 1,
                   downloadLocation = dir_gene_files)

mathys_genes <- read.table(filename$path, header = FALSE) |>
  dplyr::rename(ensembl_gene_id = V1, symbol_Mathys = V2)


## Other single cell genes -----------------------------------------------------

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

  df <- rtracklayer::import(gzfile(file_info[["filename"]]), format = "gtf") |>
    as.data.frame() |>
    dplyr::select(gene_id, gene_name) |>
    dplyr::distinct()

  colnames(df) <- c("ensembl_gene_id", version)
  return(df)
})

gtf_genes <- purrr::reduce(gtf_genes, dplyr::full_join, by = "ensembl_gene_id")


# Merge all gene sets together -------------------------------------------------

all_genes <- purrr::reduce(list(gene_info, gtf_genes, mathys_genes),
                           dplyr::full_join,
                           by = "ensembl_gene_id")

# For each row/gene, take the symbol that appears the most across all data sets.
# In case of ties, which.max returns the first symbol in the tie, which would
# be the first alphabetically.
max_symbol <- all_genes |>
  tidyr::pivot_longer(cols = starts_with("symbol"),
                      names_to = "dataset",
                      values_to = "symbol",
                      values_drop_na = TRUE) |>
  group_by(ensembl_gene_id, symbol) |>
  summarize(count = n(), .groups = "drop_last") |>
  summarize(canonical_symbol = symbol[which.max(count)])

all_genes <- merge(all_genes, max_symbol, all = TRUE) |>
  dplyr::arrange(ensembl_gene_id)

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
