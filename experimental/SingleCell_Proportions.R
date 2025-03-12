library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

### Cell type proportion differences between diagnoses

reference_dataset <- "leng"
granularity <- "broad_class"

sce <- Load_SingleCell(reference_dataset, granularity, "counts")
table(sce$celltype)

if (reference_dataset == "lau") {
  sce$diagnosis <- factor(str_replace(as.character(sce$diagnosis), "healthy c", "C"))
}

props <- table(sce$sample, sce$celltype)
props <- sweep(props, 1, rowSums(props), "/")
props <- as.data.frame(props)
colnames(props) <- c("sample", "celltype", "proportion")

diagnosis <- data.frame(sample = sce$sample, diagnosis = sce$diagnosis) %>%
  distinct()

props <- merge(props, diagnosis, by = "sample")
props_sub <- subset(props, diagnosis %in% c("Control", "AD"))

ann <- aov(proportion ~ celltype * diagnosis, data = props_sub)
summary(ann)

tuk <- TukeyHSD(ann, "celltype:diagnosis")
celltypes <- levels(sce$celltype)
terms <- paste0(celltypes, ":Control-", celltypes, ":AD")
results <- as.data.frame(tuk$`celltype:diagnosis`)[terms,]
print(results)

ggplot(props_sub, aes(x = celltype, y = proportion, color = diagnosis)) +
  geom_boxplot()


### Mismatches between original cell assignments and mapped cell assignments

reference_dataset <- "seaRef"
granularity <- "broad_class"

sce <- Load_SingleCell(reference_dataset, granularity, "counts")

sce$original_broad_class <- str_replace(sce$original_broad_class, "Endothelial", "Vascular")
sce$original_broad_class <- str_replace(sce$original_broad_class, "Pericyte", "Vascular")
sce$original_broad_class <- str_replace(sce$original_broad_class, "VLMC", "Vascular")
sim <- table(sce$broad_class, sce$original_broad_class)

print(sim)

# Things that were originally assigned cell type X and are not mapped to cell type X
false_neg <- sapply(rownames(sim), function(N) {
  (sum(sim[, N]) - sim[N, N]) / sum(sim[, N])
})
print(false_neg)

# Things that were mapped to cell type X that were not originally assigned cell type X
false_pos <- sapply(rownames(sim), function(N) {
  (sum(sim[N, ]) - sim[N, N]) / sum(sim[N, ])
})
print(false_pos)

# Leng doesn't have clear concordance between original sub class and mapped sub class,
# but best guess based on table() for subtypes is:
if (reference_dataset == "leng") {
  mapping <- c(
    "Exc.1" = "Exc.1",
    "Exc.2" = "Exc.1", # some in Exc.2 and Exc.4
    "Exc.3" = "Exc.1", # but a lot in Exc.8/9
    "Exc.4" = "Exc.4", # but a lot in Exc.1, some in Exc.2
    "Exc.5" = "Exc.7", # many in Exc.4 and some in Exc.1 and Exc.2
    "Exc.6" = "Exc.8", # Exc.8, 9, and 6 get merged
    "Exc.7" = "Exc.1",
    "Inh.1" = "Inh.2", # some in Inh.4
    "Inh.2" = "Inh.1",
    "Inh.3" = "Inh.3",
    "Inh.4" = "Inh.5" # Inh.5 and 7 get merged
  )
}
# So it's hard to verify subclasses for Leng. We have to merge a lot of (mapped)
# subclasses to account for the number of mapped subclasses that do not appear
# in the above list because they map to a smaller proportion of the original
# cell type.
# Merge Exc.1 and Exc.10 together, Exc.2, Exc.3, Exc.4, and Exc.5 together,
# Exc.6, Exc.8, and Exc.9 together, Inh.1 and Inh.6 together,
# Inh.2 and Inh.4 together, Inh.5 and Inh.7 together
# There are no VLMCs or Pericytes in Leng

# Mathys has more original subclasses than the mapped ones but they match fairly well:
if (reference_dataset == "mathys") {
  mapping <- c(
    "Ex0" = "Exc.1",
    "Ex1" = "Exc.4",
    "Ex11" = "Exc.10",
    "Ex12" = "Exc.9",
    "Ex14" = "Exc.6",
    "Ex2" = "Exc.1", # smaller amount to Exc.2 and Exc.4
    "Ex3" = "Exc.7",
    "Ex4" = "Exc.1",
    "Ex5" = "Exc.2", # smaller amount to Exc.5
    "Ex6" = "Exc.1",
    "Ex7" = "Exc.1", # Equal split between Exc.1 and Exc.3, is the only type that maps to Exc.3
    "Ex8" = "Exc.1",
    "Ex9" = "Exc.8",
    "In0" = "Inh.1",
    "In1" = "Inh.2",
    "In10" = "Inh.4",
    "In11" = "Inh.3",
    "In2" = "Inh.3",
    "In3" = "Inh.4",
    "In4" = "Inh.3",
    "In5" = "Inh.5",
    "In6" = "Inh.2", # smaller amount to Inh.4
    "In7" = "Inh.1",
    "In8" = "Inh.7",
    "In9" = "Inh.6"
  )
}
# Merge Exc.1 and Exc.3 together, Exc.2 and Exc.5 together for mathys
# There are no VLMCs in mathys

# seaRef has mostly good matches to the mapped subclasses:
if (reference_dataset == "seaRef") {
  mapping <- c(
    "Chandelier" = "Inh.6",
    "L2/3 IT" = "Exc.1",
    "L4 IT" = "Exc.2", # some also map to Exc.1, Exc.3, and Ex.4
    "L5 ET" = "Exc.1", # some also map to Exc.4
    "L5 IT" = "Exc.4", # some also map to Exc.5
    "L5/6 NP" = "Exc.6",
    "L6 CT" = "Exc.9",
    "L6 IT" = "Exc.7",
    "L6 IT Car3" = "Exc.10",
    "L6b" = "Exc.8",
    "Lamp5" = "Inh.5", # Some also match to Exc.4
    "Lamp5 Lhx6" = "Inh.7",
    "Pax6" = "Inh.4",
    "Pvalb" = "Inh.1", # small amount also match to Inh.3
    "Sncg" = "Inh.4",
    "Sst" = "Inh.3", # small amount also match to Inh.1
    "Sst Chodl" = "Inh.3",
    "Vip" = "Inh.2" # some also match to Inh.4
  )
}
# Merge Exc.4 and Exc.5 together, Exc.3 and Exc.2 together
# seaRef doesn't have any Pericyte cells

sce$original_sub_class <- as.character(sce$original_sub_class)
sce$original_sub_class <- str_replace(sce$original_sub_class, "Ast.*", "Astrocyte")
sce$original_sub_class <- str_replace(sce$original_sub_class, "Oli.*", "Oligodendrocyte")
sce$original_sub_class <- str_replace(sce$original_sub_class, "Mic.*", "Microglia")
sce$original_sub_class <- str_replace(sce$original_sub_class, "Opc.*", "OPC")
sce$original_sub_class <- str_replace(sce$original_sub_class, "End.*", "Endothelial")
sce$original_sub_class <- str_replace(sce$original_sub_class, "Per.*", "Pericyte")

mapping <- c(mapping, "Astrocyte" = "Astrocyte", "Oligodendrocyte" = "Oligodendrocyte",
             "Microglia" = "Microglia", "OPC" = "OPC", "Endothelial" = "Endothelial",
             "Pericyte" = "Pericyte", "VLMC" = "VLMC")

sce$original_sub_class <- mapping[sce$original_sub_class]

### Leng only
if (reference_dataset == "leng") {
  sce$sub_class[sce$sub_class == "Exc.10"] <- "Exc.1"
  sce$sub_class[sce$sub_class %in% c("Exc.2", "Exc.3", "Exc.5")] <- "Exc.4"
  sce$sub_class[sce$sub_class %in% c("Exc.6", "Exc.9")] <- "Exc.8"
  sce$sub_class[sce$sub_class == "Inh.6"] <- "Inh.1"
  sce$sub_class[sce$sub_class == "Inh.4"] <- "Inh.2"
  sce$sub_class[sce$sub_class == "Inh.7"] <- "Inh.5"
  sce$sub_class <- factor(sce$sub_class)
}
###

### Mathys only
if (reference_dataset == "mathys") {
  sce$sub_class[sce$sub_class == "Exc.3"] <- "Exc.1"
  sce$sub_class[sce$sub_class == "Exc.5"] <- "Exc.2"
  sce$sub_class <- factor(sce$sub_class)
}
###

### SeaRef only
if (reference_dataset == "seaRef") {
  sce$sub_class[sce$sub_class == "Exc.5"] <- "Exc.4"
  sce$sub_class[sce$sub_class == "Exc.3"] <- "Exc.2"
  sce$sub_class <- factor(sce$sub_class)
}

sim <- table(sce$sub_class, sce$original_sub_class)

### Leng only
if (reference_dataset == "leng") {
  sim <- cbind(cbind(sim, rep(0, nrow(sim))), rep(0, nrow(sim)))
  colnames(sim)[(ncol(sim)-1):ncol(sim)] <- c("Pericyte", "VLMC")
}
###

### Mathys only
if (reference_dataset == "mathys") {
  sim <- cbind(sim, rep(0, nrow(sim)))
  colnames(sim)[ncol(sim)] <- "VLMC"
}
###

### seaRef only
if (reference_dataset == "seaRef") {
  sim <- cbind(sim, rep(0, nrow(sim)))
  colnames(sim)[ncol(sim)] <- "Pericyte"
}
###

print(sim)

# Things that were originally assigned cell type X and are not mapped to cell type X
false_neg <- sapply(rownames(sim), function(N) {
  (sum(sim[, N]) - sim[N, N]) / sum(sim[, N])
})
print(false_neg)

# Things that were mapped to cell type X that were not originally assigned cell type X
false_pos <- sapply(rownames(sim), function(N) {
  (sum(sim[N, ]) - sim[N, N]) / sum(sim[N, ])
})
print(false_pos)

