library(Matrix)
library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)
library(reshape2)
library(patchwork)

source(file.path("functions", "Plotting_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step11_Error_HelperFunctions.R"))
source(file.path("functions", "FileIO_HelperFunctions.R"))

options(scipen = 999)

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

est_fields = list("CibersortX" = "estimates",
                  "DeconRNASeq" = "out.all",
                  "Dtangle" = "estimates",
                  "DWLS" = "estimates",
                  "HSPE" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "Scaden" = "estimates",
                  "Baseline" = "estimates")

algorithms <- names(est_fields)

granularity <- c("broad_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors <- Get_AllBestErrorsAsDf(bulk_datasets, granularity)

bulk_dataset <- "Mayo"
dataset <- "leng"
ests_alg <- Get_AllEstimatesAsDf(dataset, bulk_dataset, algorithms,
                                 granularity, best_errors, est_fields)

ct <- "Excitatory"
metric <- "^cor"

set.seed(12345)
samps <- sample(unique(ests_alg$sample), size = 20, replace = FALSE)
ests_sub <- subset(ests_alg, celltype == ct & sample %in% samps)
ests_sub <- ests_sub[grepl(metric, ests_sub$metrics),]
ests_sub$sample_numeric <- as.numeric(factor(ests_sub$sample))

ggplot(ests_sub, aes(x = sample, y = pct_est, fill = algorithm)) +
  geom_boxplot() +
  theme(legend.position = "bottom")

tmp <- ests_sub %>% group_by(algorithm, sample) %>% mutate(rel_pct_est = pct_est / mean(pct_est))
ggplot(tmp, aes(x = sample, y = rel_pct_est, fill = algorithm)) +
  geom_boxplot() +
  theme(legend.position = "bottom")

tmp2 <- ests_sub %>% subset(sample == samps[1]) %>% group_by(algorithm) %>%
  summarize(mean_pct = mean(pct_est)) %>% merge(ests_sub, by = "algorithm")
tmp2$rel_pct_est <- tmp2$pct_est / tmp2$mean_pct
ggplot(tmp2, aes(x = sample, y = rel_pct_est, fill = algorithm)) +
  geom_boxplot() +
  theme(legend.position = "bottom")
