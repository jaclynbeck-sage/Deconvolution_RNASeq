library(Matrix)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)

source("Filenames.R")

dataset <- "morabito"
datatype <- "training"
algorithm = "deconRNASeq"

est_fields = list("dtangle" = "estimates",
                  "music" = "Est.prop.weighted", # TODO analyze Est.prop.allgene too
                  "music2" = "Est.prop",
                  "hspe" = "estimates",
                  "deconRNASeq" = "Est.prop")
est_field = est_fields[[algorithm]]

errs <- readRDS(file = file.path(dir_output, paste0("errors_", dataset, "_",
                                                    datatype, "_broad.rds")))

se <- readRDS(file.path(dir_pseudobulk,
                               paste0("pseudobulk_", dataset, "_", datatype,
                                      "_broadcelltypes.rds")))
props <- colData(se)

# Not implemented yet
#errs_test <- readRDS(file = file.path(dir_output,
#                                      paste0("errors_", dataset, "_broad_ROSMAP.rds")))

err <- errs[[algorithm]][["means"]]
ests <- readRDS(file = file.path(dir_params_lists,
                                 paste0(algorithm, "_list_", dataset, "_",
                                        datatype, "_broad.rds")))

bests <- lapply(colnames(err), FUN = function(X) {
  if (length(grep("cor", X)) > 0) {
    ord <- order(err[, X], decreasing = TRUE)
  }
  else {
    ord <- order(err[, X], decreasing = FALSE)
  }
  ord <- ord[1:min(3, length(ord))]
  data.frame(name = rownames(err[ord,]), group = X)
})

bests <- do.call(rbind, bests)
#bests <- distinct(bests)

gof_rows <- grepl("gof", bests$group)

par(mfcol = c(3, 6))

for (G in unique(bests$group[!gof_rows])) {
  params <- subset(bests, group == G)$name
  best_ests <- lapply(ests[params], FUN = function(X) {
    X[[est_field]][rownames(props), colnames(props)]
  })

  best_errs <- err[names(best_ests),]
  for (N in names(best_ests)) {
    matplot(props, best_ests[[N]], xlim = c(0,1), ylim=c(0,1),
            xlab="Truth", ylab="Estimates", main = paste(G, round(best_errs[N,G], digits = 3)))
    abline(a = 0, b = 1)
  }
}

# Goodness of fit
#par(mfcol = c(3, 4))
#for (G in unique(bests$group)[6:9]) {
#  params <- subset(bests, group == G)$name
#  best_ests <- lapply(ests[params], FUN = function(X) {
#    X$estimates[rownames(props), colnames(props)]
#  })
#
#  best_errs <- err[names(best_ests),]
#  for (N in names(best_ests)) {
#    matplot(props, best_ests[[N]], xlim = c(0,1), ylim=c(0,1),
#            xlab="Truth", ylab="Estimates", main = paste(G, round(best_errs[N,G], digits = 3)))
#    abline(a = 0, b = 1)
#  }
#}

par(mfrow = c(1,1))
errs_means <- lapply(errs, '[[', "means")
errs_means <- lapply(names(errs_means), FUN = function(X) {
  errs_means[[X]]$method <- X
  errs_means[[X]]
})

errs_df <- melt(do.call(rbind, errs_means))
colnames(errs_df) <- c("method", "error.type", "value")

errs1 <- subset(errs_df, !grepl("gof", error.type))
#errs2 <- subset(errs_df, grepl("gof", error.type))

ggplot(errs1, aes(x = method, y = value, fill = method)) +
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~error.type, scales = "free")

#ggplot(errs2, aes(x = method, y = value, fill = method)) +
#  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  facet_wrap(~error.type, scales = "free")

# How well is each cell type characterized?

errs_by_celltype <- errs[[algorithm]][["by_celltype"]] #lapply(errs, "[[", "by_celltype")
errs_by_celltype <- lapply(names(errs_by_celltype), FUN = function(X) {
  errs_by_celltype[[X]]$name <- X
  errs_by_celltype[[X]]$celltype <- rownames(errs_by_celltype[[X]])
  errs_by_celltype[[X]]
})
errs_df <- melt(do.call(rbind, errs_by_celltype))
colnames(errs_df) <- c("name", "celltype", "error.type", "value")

errs1 <- subset(errs_df, !grepl("gof", error.type) & name %in% bests$name)
#errs2 <- subset(errs_df, grepl("gof", error.type))

ggplot(errs1, aes(x = celltype, y = value, fill = celltype)) +
  geom_violin() + geom_jitter(size = 0.5, width = 0.1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~error.type, scales = "free")


##### Donors only
# How well was each subject characterized? -- not as useful
errs_by_subject <- errs[[algorithm]][["by_subject"]] #lapply(errs, "[[", "by_subject")
errs_by_subject <- lapply(names(errs_by_subject), FUN = function(X) {
  errs_by_subject[[X]]$name <- X
  errs_by_subject[[X]]$subject <- rownames(errs_by_subject[[X]])
  errs_by_subject[[X]]
})
errs_df <- melt(do.call(rbind, errs_by_subject))
colnames(errs_df) <- c("name", "subject", "error.type", "value")

errs1 <- subset(errs_df, !grepl("gof", error.type) & name %in% bests$name)
ggplot(errs1, aes(x = subject, y = value, fill = subject)) +
  geom_violin() + geom_jitter(size = 0.5, width = 0.1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~error.type, scales = "free")

##### Training data only
# How well is each cell type characterized by proportion?
props_melt <- as.data.frame(props)
props_melt$subject <- rownames(props_melt)
props_melt <- melt(props_melt)
colnames(props_melt) <- c("subject", "celltype", "prop_truth")

ests_melt <- lapply(ests, "[[", est_field)
ests_melt <- lapply(names(ests_melt), FUN = function(X) {
  ests_melt[[X]] <- as.data.frame(ests_melt[[X]])
  ests_melt[[X]]$name <- X
  ests_melt[[X]]$subject <- rownames(ests_melt[[X]])
  ests_melt[[X]]
})
ests_melt <- do.call(rbind, ests_melt)
ests_melt <- melt(ests_melt)
colnames(ests_melt) <- c("name", "subject", "celltype", "prop_est")

ests_df <- merge(props_melt, ests_melt, by = c("subject", "celltype"))

ests1 <- subset(ests_df, name %in% bests$name)

ggplot(ests1, aes(x = prop_truth, y = prop_est, color = name)) +
  geom_jitter(size = 0.5) + geom_abline(slope = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", aspect.ratio = 1) +
  facet_wrap(~celltype, scales = "free", ncol = 4)

# Zoom in on the crowded area < 30%
ests2 <- subset(ests1, prop_est <= 0.3 | prop_truth <= 0.3)
ggplot(ests2, aes(x = prop_truth, y = prop_est, color = name)) +
  geom_jitter(size = 0.5) + geom_abline(slope = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", aspect.ratio = 1) +
  facet_wrap(~celltype, scales = "free", ncol = 4)

# Lines instead of points
ests3 <- ests1 %>% mutate(prop_truth = round(prop_truth, digits = 2)) %>%
          group_by(name, prop_truth, celltype) %>%
          mutate(avg_est = mean(prop_est)) %>%
          select(-prop_est) %>% distinct()
ggplot(ests3, aes(x = prop_truth, y = avg_est, color = name)) +
  geom_line() + geom_abline(slope = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", aspect.ratio = 1) +
  facet_wrap(~celltype, scales = "free", ncol = 4)

ests4 <- subset(ests3, prop_truth <= 0.3 | avg_est <= 0.3)
ggplot(ests4, aes(x = prop_truth, y = avg_est, color = name)) +
  geom_line() + geom_abline(slope = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", aspect.ratio = 1) +
  facet_wrap(~celltype, scales = "free", ncol = 4)
