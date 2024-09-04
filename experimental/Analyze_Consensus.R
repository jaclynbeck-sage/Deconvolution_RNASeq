library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)
#library(foreach)
#library(doParallel)

source(file.path("functions", "Step11_Error_HelperFunctions.R"))

#cores <- 6
#cl <- makeCluster(cores, type = "FORK", outfile = "")
#registerDoParallel(cl)

reference_datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE", "Music",
                "Scaden", "Baseline")

granularity <- c("sub_class")

bulk_datasets <- c("Mayo")#, "MSBB", "ROSMAP")

combined_metadata <- lapply(bulk_datasets, function(B) {
  bulk <- Load_BulkData(B)
  meta <- as.data.frame(colData(bulk)) %>% select(sample, diagnosis, tissue)
})
combined_metadata <- do.call(rbind, combined_metadata)

best_errors_list <- Get_AllBestErrorsAsDf(bulk_datasets, granularity)
best_errors <- best_errors_list$errors
#best_estimates <- best_errors_list$estimates

# We are analyzing each tissue separately, so we will use the best parameter
# sets for each tissue instead of the best overall
best_errors <- subset(best_errors, tissue != "All")

best_signature <- best_errors %>%
  subset(algorithm != "Baseline") %>%
  group_by(tissue, param_id) %>%
  summarize(best_cor = signature[which.max(cor)],
            best_rMSE = signature[which.min(rMSE)],
            best_mAPE = signature[which.min(mAPE)],
            .groups = "drop")

get_best_signature <- function(best_cor, best_rMSE, best_mAPE) {
  tmp <- table(c(best_cor, best_rMSE, best_mAPE))
  names(tmp)[which.max(tmp)]
}

best_signature <- best_signature %>%
  group_by(tissue) %>%
  summarize(best_sig = get_best_signature(best_cor, best_rMSE, best_mAPE),
            .groups = "drop")

best_errors <- merge(best_errors, best_signature,
                     by = c("tissue"))
best_errors <- subset(best_errors, signature == best_sig) %>%
  select(-best_sig)

# The baseline "zeros" entries have NA correlation so they need to return
# NA for best_cor since which.max would return 0 entries
best_params <- best_errors %>%
  group_by(tissue, reference_data_name, test_data_name, reference_input_type,
           normalization, regression_method, algorithm) %>%
  summarize(best_cor = ifelse(all(is.na(cor)), NA, param_id[which.max(cor)]),
            best_rMSE = param_id[which.min(rMSE)],
            best_mAPE = param_id[which.min(mAPE)],
            .groups = "drop") %>%
  melt(measure.vars = c("best_cor", "best_rMSE", "best_mAPE"),
       value.name = "param_id") %>%
  select(tissue, param_id) %>%
  subset(!is.na(param_id)) %>%
  distinct()

best_params_all <- unique(best_params$param_id)

best_errors <- merge(best_params, best_errors, by = c("tissue", "param_id"),
                     all.y = FALSE)

best_estimates <- Get_AllBestEstimatesAsDf(bulk_datasets, granularity,
                                           combined_metadata, best_params)

gc()

#best_errors <- subset(best_errors, normalization %in% c("counts", "cpm", "log_cpm") &
#                        regression_method == "edger")

mean_props_all <- lapply(unique(best_errors$tissue), function(T) {
  errs_tmp <- subset(best_errors, tissue == T)
  bulk_dataset <- unique(errs_tmp$test_data_name) # Tissues are unique to a single bulk dataset

  bulk_se <- Load_BulkData(bulk_dataset)
  metadata <- as.data.frame(colData(bulk_se)) %>%
    subset(tissue == T) %>%
    select(sample, diagnosis, tissue)

  props_b <- lapply(unique(errs_tmp$reference_data_name), function(ref_dataset) {
    print(paste(bulk_dataset, T, ref_dataset))

    param_ids <- errs_tmp %>%
      subset(reference_data_name == ref_dataset) %>%
      select(param_id) %>%
      unlist()

    ests_ad <- subset(best_estimates, param_id %in% param_ids &
                        sample %in% metadata$sample)

    # TODO make functions to do this and move these files to a better place
    ests_ad <- merge(ests_ad, metadata, by = "sample") #%>%
                  #subset(diagnosis %in% c("CT", "AD")) %>%
                  #subset(celltype %in% levels(ests_ad$celltype)) # Gets rid of added cell types from ROSMAP IHC

    significant <- lapply(unique(ests_ad$param_id), function(param_set) { #foreach(param_set = unique(ests_ad$param_id)) %dopar% {
      ests_param <- subset(ests_ad, param_id == param_set)

      anov <- aov(percent ~ diagnosis*celltype, data = ests_param)
      summ <- summary(anov)[[1]]
      tuk <- TukeyHSD(anov, "diagnosis:celltype")

      comparisons <- paste0("CT:", levels(ests_ad$celltype), "-AD:", levels(ests_ad$celltype))
      tuk <- as.data.frame(tuk[[1]][comparisons,])
      tuk$p_adj <- tuk$`p adj`
      tuk$celltype <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,3]
      tuk$tissue <- T
      tuk$significant <- tuk$p_adj <= 0.05
      tuk$anova_significant <- summ["diagnosis:celltype", "Pr(>F)"] <= 0.05
      tuk$param_id <- param_set
      return(tuk)
    })

    significant <- do.call(rbind, significant) %>%
      select(celltype, tissue, p_adj, significant, param_id) #%>%
                      #merge(n_subjects, by = c("celltype", "tissue"))

    mean_props <- ests_ad %>% subset(diagnosis %in% c("CT", "AD")) %>%
                    group_by(param_id, diagnosis, tissue, celltype) %>%
                    summarise(mean_pct = mean(percent),
                              sd_pct = sd(percent),
                              count = n(),
                              .groups = "drop")


    # Make one column for AD and one for CT for each of mean_pct, sd_pct, and count,
    # calculate fold-change between AD and CT means
    # TODO how to deal with cases where mean_pct_CT = 0?
    mean_props <- pivot_wider(mean_props, names_from = "diagnosis",
                              values_from = c("mean_pct", "sd_pct", "count")) %>%
                    mutate(fc = mean_pct_AD / mean_pct_CT,
                           log2_fc = log2(fc)) # equivalent to log2(AD)-log2(CT)
                           #sd_fc = calc_sd_fc(mean_pct_AD,))

    mean_props <- merge(mean_props, significant, by = c("param_id", "celltype", "tissue"))
    return(mean_props)
  })

  return(do.call(rbind, props_b))
})

mean_props_all <- do.call(rbind, mean_props_all) %>%
                    subset(mean_pct_AD != 0 & mean_pct_CT != 0) %>%
                    mutate(log_p = log(p_adj + 1e-8)) # pseudocount necessary to avoid log(0)

# TODO have to remove baseline from combined stats & fisher

# TODO a heatmap of "direction" would be a nice way to display agreement
# across runs -- closer to 1 or -1 means stronger agreement
combined_stats <- mean_props_all %>% group_by(celltype, tissue) %>%
                    summarise(mean_fc = mean(fc),
                              sd_fc = sd(fc),
                              mean_log2fc = mean(log2_fc),
                              sd_log2fc = sd(log2_fc),
                              direction = mean(sign(log2_fc)),
                              mean_p = mean(p_adj),
                              geo_mean_p = exp(mean(log_p)),
                              .groups = "drop")

print(subset(combined_stats, geo_mean_p <= 0.05 & abs(mean_log2fc) >= 0.5))


fisher <- mean_props_all %>% group_by(celltype, tissue) %>%
            summarise(fisher_score = sum(-log_p),
                      min_p = min(p_adj),
                      max_p = max(p_adj))

pval <- sapply(fisher$fisher_score, function(x) {
  ## calculate Fisher p-value
  pchisq(2 * x, lower.tail=FALSE, 2 * length(fisher$fisher_score) )
})
qval <- p.adjust(pval, method="BH")

stat <- sapply(fisher$fisher_score, function(x) {2*x})

fisher$pooled_p <- pval
fisher$pooled_FDR <- qval


# Tried seeing what an ANOVA would look like, treating the multiple runs like
# repeated measures, but it requires over a terabyte of memory for broad class
# cell types, probably way more for subclass. So this averages the top parameter
# sets from each algorithm and then does an ANOVA on the averages

metadata <- lapply(c("Mayo", "MSBB", "ROSMAP"), function(bulk_dataset) {
  bulk <- Load_BulkData(bulk_dataset)
  meta <- as.data.frame(colData(bulk)) %>%
    select(sample, tissue, diagnosis)
  rm(bulk)
  return(meta)
})
metadata <- do.call(rbind, metadata)
gc()

best_estimates$algorithm <- factor(str_replace(best_estimates$param_id, "_.*", ""))

mean_ests <- best_estimates %>%
  group_by(algorithm, celltype, sample) %>%
  summarize(mean_pct = mean(percent),
            sd_pct = sd(percent),
            rel_sd_pct = sd_pct / mean_pct,
            .groups = "drop")

mean_ests <- merge(mean_ests, metadata, by = "sample")

res <- lapply(levels(mean_ests$tissue), function(tiss) {
  ests_df <- subset(mean_ests, tissue == tiss)

  #anov <- aov(percent ~ diagnosis*celltype + Error(sample/param_id), data = ests_df)
  fit <- lmer(mean_pct ~ diagnosis*celltype*algorithm + (1|sample), data = ests_df)
  anov <- car::Anova(fit)
  emm <- emmeans::emmeans(fit, ~ diagnosis * celltype)

  comparisons <- paste0("AD ", levels(ests_df$celltype), " - CT ", levels(ests_df$celltype))
  pw <- pairs(emm)
  pw <- as.data.frame(subset(pw, contrast %in% comparisons))

  pw$celltype <- str_split(pw$contrast, pattern = " ", simplify = TRUE)[,2]
  pw$tissue <- tiss
  pw$significant <- pw$p.value <= 0.05
  pw$anova_significant <- anov["diagnosis:celltype", "Pr(>Chisq)"] < 0.05

  return(pw)
})

res <- do.call(rbind, res)

stopCluster(cl)

#library(meta)

mean_astro <- subset(mean_props_all, celltype == "Astrocyte")
mean_astro_tcx <- subset(mean_astro, tissue == "TCX")
mean_astro_cbe <- subset(mean_astro, tissue == "CBE")

res_astro_tcx <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = param_id,
                          data = mean_astro_tcx)

res_astro_cbe <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = param_id,
                          data = mean_astro_cbe,
                          sm = "SMD")


mean_props2 <- mean_props_all
tmp = str_split(mean_props2$param_id, "_", simplify = TRUE)[,1:2]
mean_props2$group <- paste(tmp[,1], tmp[,2], sep = "_")

mean_ratios <- mean_props2 %>% group_by(group, tissue, celltype) %>%
                summarise(mean_ratio = mean(mean_pct_AD / mean_pct_CT),
                          sd_ratio = sd(mean_pct_AD / mean_pct_CT),
                          n = mean(count_AD) + mean(count_CT),
                          .groups = "drop")

mean_astro <- subset(mean_ratios, celltype == "Astrocyte")
mean_astro_tcx <- subset(mean_astro, tissue == "TCX")
mean_astro_cbe <- subset(mean_astro, tissue == "CBE")

res_astro_tcx <- metamean(n = n, mean = mean_ratio, sd = sd_ratio,
                          studlab = group, data = mean_astro_tcx,
                          null.effect = 1)

res_astro_cbe <- metamean(n = n, mean = mean_ratio, sd = sd_ratio,
                          studlab = group, data = mean_astro_cbe,
                          null.effect = 1)

##### This section seems to make the most sense #####
# But we either need to map all datasets into the same space or keep the
# results on the diff datasets separate -- large heterogeneity between the
# datasets because the cell type definitions are not the same

# Removing ROSMAP IHC bests for now
ok <- sapply(str_split(best_params_all$metrics, ", "), function(M) {
  any(M %in% c("cor", "rMSE", "mAPE"))
})

#best_params_sub <- best_params_all[ok,]
#to_keep <- unique(subset(errs_tissue, param_id %in% best_params_sub$param_id)$param_id)

#mean_props_sub <- subset(mean_props_all, param_id %in% best_errors)
mean_props_sub <- mean_props_all

mean_astro <- subset(mean_props_sub, celltype == "Astrocyte")
mean_astro_tcx <- subset(mean_astro, tissue == "TCX")
mean_astro_cbe <- subset(mean_astro, tissue == "CBE")

res_astro_tcx <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = param_id,
                          data = mean_astro_tcx,
                          sm = "ROM") # Ratio of means

res_astro_cbe <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = param_id,
                          data = mean_astro_cbe,
                          sm = "ROM")


##### End #####

reference <- str_split(mean_astro_cbe$param_id, "_", simplify = TRUE)[,2]

res_astro_cbe_meta <- metareg(res_astro_cbe, ~reference)

mean_props_sub$reference <- str_split(mean_props_sub$param_id, "_", simplify = TRUE)[,2]

mean_astro_cain_cbe <- subset(mean_props_sub, celltype == "Astro" & tissue == "CBE" & reference == "cain")
res_astro_cain_cbe <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = param_id,
                          data = mean_astro_cain_cbe,
                          sm = "ROM")

mean_astro_cain_tcx <- subset(mean_props_sub, celltype == "Astro" & tissue == "TCX" & reference == "cain")
res_astro_cain_tcx <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                               n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                               studlab = param_id,
                               data = mean_astro_cain_tcx,
                               sm = "ROM")


metas <- lapply(levels(mean_props_sub$celltype), function(ct) {
  res_tissue <- lapply(levels(mean_props_sub$tissue), function(tiss) {
    mean_props_tmp <- subset(mean_props_sub, celltype == ct & tissue == tiss)
    tryCatch({
      res <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                      n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                      studlab = param_id,
                      data = mean_props_tmp,
                      sm = "ROM")
      return(data.frame(celltype = ct, tissue = tiss,
                        effect.common = res$TE.common, pval.common = res$pval.common,
                        effect.random = res$TE.random, pval.random = res$pval.random,
                        I2 = res$I2, Q = res$Q, pval.Q = res$pval.Q))
    },
    error = function(err) {
      print(paste(ct, tiss))
      print(err)
      return(NULL)
    })
  })
  res_tissue <- do.call(rbind, res_tissue)
  return(res_tissue)
})

metas <- do.call(rbind, metas)
# pval.adjust ??
