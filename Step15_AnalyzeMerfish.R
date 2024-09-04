library(anndata)
library(dplyr)
library(readxl)
library(stringr)
library(reshape2)

source(file.path("functions", "FileIO_HelperFunctions.R"))

donor_meta_file <- file.path(dir_seaad_raw, "donor_metadata.xlsx")
merfish_file <- file.path(dir_seaad_raw, "mtg_merfish.h5ad")
cain_mapping_file <- file.path(dir_seaad_raw, "cain_celltype_annotation.csv")

if (!file.exists(merfish_file)) {
  download.file(url_seaad_merfish, merfish_file, method = "curl")
}

download.file(url_seaad_donor_metadata, donor_meta_file, method = "curl")
download.file(url_seaad_cain_map, cain_mapping_file, method = "curl")

donor_metadata <- read_excel(donor_meta_file)
colnames(donor_metadata) <- make.names(colnames(donor_metadata))

# Some cleanup/shortening of consensus diagnosis fields
colnames(donor_metadata) <- str_replace(colnames(donor_metadata),
                                        "Consensus.*choice", "Diagnosis")

diagnosis_fields <- c("Diagnosis.Alzheimers.disease.",
                      "Diagnosis.Alzheimers.Possible..Probable.",
                      "Diagnosis.Control.")

mf_h5 <- anndata::read_h5ad(merfish_file)

metadata <- mf_h5$obs
colnames(metadata) <- make.names(colnames(metadata))

metadata <- subset(metadata, Used.in.analysis == TRUE) %>%
  select(Donor.ID, Cell.ID, Class, Subclass)

donor_metadata <- subset(donor_metadata, Donor.ID %in% metadata$Donor.ID)

for (field in diagnosis_fields) {
  donor_metadata[, field] <- (donor_metadata[, field] == "Checked")
}

donor_metadata$high_thal <- donor_metadata$Thal %in% paste("Thal", 2:5)
donor_metadata$high_braak <- donor_metadata$Braak %in% paste("Braak", c("IV", "V", "VI"))
donor_metadata$high_cerad <- donor_metadata$CERAD.score %in% c("Moderate", "Frequent")

donor_metadata$diagnosis <- ""
donor_metadata$diagnosis[donor_metadata$high_thal & donor_metadata$high_braak] <- "AD"
donor_metadata$diagnosis[!donor_metadata$high_thal & !donor_metadata$high_braak] <- "CT"
donor_metadata$diagnosis[!donor_metadata$high_thal & donor_metadata$high_braak] <- "CT_HIGH_PATH"
donor_metadata$diagnosis[donor_metadata$high_thal & !donor_metadata$high_braak] <- "AD_LOW_PATH"

donor_metadata$diagnosis2 <- ""
donor_metadata$diagnosis2[donor_metadata$high_cerad & donor_metadata$high_braak] <- "AD"
donor_metadata$diagnosis2[!donor_metadata$high_cerad & !donor_metadata$high_braak] <- "CT"
donor_metadata$diagnosis2[!donor_metadata$high_cerad & donor_metadata$high_braak] <- "CT_HIGH_PATH"
donor_metadata$diagnosis2[donor_metadata$high_cerad & !donor_metadata$high_braak] <- "AD_LOW_PATH"

donor_metadata$diagnosis3 <- "OTHER"
donor_metadata$diagnosis3[donor_metadata$Diagnosis.Alzheimers.disease.] <- "AD"
donor_metadata$diagnosis3[donor_metadata$Diagnosis.Alzheimers.Possible..Probable.] <- "AD"
donor_metadata$diagnosis3[donor_metadata$Diagnosis.Control.] <- "CT"

metadata <- merge(metadata, donor_metadata[, c("Donor.ID", "diagnosis", "diagnosis2", "diagnosis3")],
                  by = "Donor.ID")
metadata$Donor.ID <- as.character(metadata$Donor.ID)

cain_map <- read.csv(cain_mapping_file, header = TRUE)
cain_map <- subset(cain_map, Used.in.analysis == "True")

sce <- Load_SingleCell("cain", "broad_class", "counts")

# Some effort required to match cell barcodes
sce$new_barcode <- paste(str_replace(sce$cell_id, ".*_", ""),
                         sce$sample, sep = "-")

sce_meta <- as.data.frame(colData(sce)) %>%
  subset(new_barcode %in% cain_map$X)

cain_map <- subset(cain_map, X %in% sce_meta$new_barcode)

cain_map <- merge(cain_map, sce_meta, by.x = "X", by.y = "new_barcode")

overlap <- table(cain_map$Subclass, cain_map$sub_class)

to_cain <- apply(overlap, 1, function(row) {
  colnames(overlap)[which.max(row)]
})

# One minor fix due to a tie between Pericyte and VLMC here
to_cain["VLMC"] <- "VLMC"

to_seaad <- apply(overlap, 2, function(col) {
  rownames(overlap)[which.max(col)]
})

collapse_c <- list()
collapse_s <- list()

for (N in 1:length(to_cain)) {
  c_name <- names(to_cain)[N]
  s_name <- unique(to_seaad[to_cain[N]])

  if ((length(s_name) > 1) | (s_name != c_name)) {
    print(paste(c_name, "/", s_name))
    collapse_c[[length(collapse_c)+1]] <- c(c_name, s_name)
  }
}

for (N in 1:length(to_seaad)) {
  s_name <- names(to_seaad)[N]
  c_name <- unique(to_cain[to_seaad[N]])

  if ((length(c_name) > 1) | (c_name != s_name)) {
    print(paste(s_name, "/", c_name))
    collapse_s[[length(collapse_s)+1]] <- c(s_name, c_name)
  }
}

to_cain_final <- to_cain
for (item in collapse_s) {
  to_cain_final[to_cain_final %in% item] <- paste(sort(item), collapse = "_")
}

to_seaad_final <- to_seaad
for (item in collapse_c) {
  to_seaad_final[to_seaad_final %in% item] <- paste(sort(item), collapse = "_")
}


saveRDS(list("seaad_to_cain" = to_cain,
             "cain_to_seaad" = to_seaad,
             "seaad_to_cain_collapsed" = to_cain_final,
             "cain_to_seaad_collapsed" = to_seaad_final),
        file.path(dir_metadata, "cain_seaad_celltype_mapping.rds"))


metadata$mapped_sub_class <- to_cain_final[metadata$Subclass]
metadata$mapped_broad_class <- metadata$mapped_sub_class
metadata$mapped_broad_class[grepl("Exc", metadata$mapped_broad_class)] <- "Excitatory"
metadata$mapped_broad_class[grepl("Inh", metadata$mapped_broad_class)] <- "Inhibitory"
metadata$mapped_broad_class[grepl("Endo|VLMC", metadata$mapped_broad_class)] <- "Vascular"

pcts_broad <- table(metadata$Donor.ID, metadata$mapped_broad_class)
pcts_broad <- as.data.frame(sweep(pcts_broad, 1, rowSums(pcts_broad), "/")) %>%
  dplyr::rename(sample = Var1, celltype = Var2, percent = Freq) %>%
  merge(donor_metadata[, c("Donor.ID", "diagnosis", "diagnosis3")],
        by.x = "sample", by.y = "Donor.ID") %>%
  mutate(across(c(sample, celltype, diagnosis), factor))

pcts_sub <- table(metadata$Donor.ID, metadata$mapped_sub_class)
pcts_sub <- as.data.frame(sweep(pcts_sub, 1, rowSums(pcts_sub), "/")) %>%
  dplyr::rename(sample = Var1, celltype = Var2, percent = Freq) %>%
  merge(donor_metadata[, c("Donor.ID", "diagnosis", "diagnosis3")],
        by.x = "sample", by.y = "Donor.ID") %>%
  mutate(across(c(sample, celltype, diagnosis), factor))

broad_stats <- pcts_broad %>%
  group_by(celltype) %>%
  summarize(percent_mean = mean(percent),
            percent_sd = sd(percent),
            percent_rel_sd = percent_sd / percent_mean)

broad_stats_diag <- pcts_broad %>%
  group_by(diagnosis3, celltype) %>%
  summarize(mean_pct = mean(percent),
            sd_pct = sd(percent),
            rel_sd_pct = sd_pct / mean_pct,
            count = n(),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = "diagnosis3",
                     values_from = c("mean_pct", "sd_pct", "rel_sd_pct", "count"))

sub_stats <- pcts_sub %>%
  group_by(celltype) %>%
  summarize(percent_mean = mean(percent),
            percent_sd = sd(percent),
            percent_rel_sd = percent_sd / percent_mean)

sub_stats_diag <- pcts_sub %>%
  group_by(diagnosis3, celltype) %>%
  summarize(mean_pct = mean(percent),
            sd_pct = sd(percent),
            rel_sd_pct = sd_pct / mean_pct,
            count = n(),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = "diagnosis3",
                     values_from = c("mean_pct", "sd_pct", "rel_sd_pct", "count"))

do_t_test <- function(data) {
  tmp = t.test(data$percent[data$diagnosis3 == "CT"], data$percent[data$diagnosis3 == "AD"])
  return(tmp$p.value)
}

tmp = pcts_broad %>%
  subset(diagnosis3 %in% c("AD", "CT")) %>%
  group_by(celltype) %>%
  summarize(p_val = do_t_test(.data))

tmp2 = pcts_sub %>%
  subset(diagnosis3 %in% c("AD", "CT")) %>%
  group_by(celltype) %>%
  summarize(p_val = do_t_test(.data))

