library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)

source(file.path("functions", "Error_HelperFunctions.R"))

reference_datasets <- c("cain", "lau", "leng", "mathys", "seaRef")
est_fields = list("Dtangle" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "HSPE" = "estimates",
                  "DeconRNASeq" = "out.all",
                  "DWLS" = "estimates")

algorithms <- c(names(est_fields), "Random")

granularity <- c("broad_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors <- Get_AllBestErrorsAsDf(bulk_datasets, granularity)

#best_params_all <- Get_AllBestParamsAsDf(reference_datasets, granularity)

#errs_all <- Get_AllErrorsAsDf(bulk_datasets, reference_datasets, algorithms,
#                              granularity, best_params_all$params)

#errs_all <- melt(errs_all) %>% dplyr::rename(error_type = "variable")

errs_melt <- melt(best_errors) %>% dplyr::rename(error_type = "variable")

errs_total <- subset(errs_melt, tissue == "All")

errs_tissue <- subset(errs_melt, tissue != "All")
errs_tissue$tissue <- paste(errs_tissue$test_data_name, errs_tissue$tissue)

file_params <- best_errors %>%
  select(reference_data_name, test_data_name, algorithm,
         reference_input_type, normalization,
         regression_method) %>%
  distinct()

mean_props_all <- lapply(reference_datasets, function(ref_dataset) {
  props_b <- lapply(bulk_datasets, function(bulk_dataset) {
    bulk_se <- Load_BulkData(bulk_dataset)
    metadata <- as.data.frame(colData(bulk_se))

    ests_alg <- Get_AllEstimatesAsDf(ref_dataset, bulk_dataset, algorithms,
                                     granularity, file_params)

    # TODO make functions to do this and move these files to a better place
    ests_ad <- merge(ests_alg, metadata, by = "sample") %>%
                  subset(diagnosis %in% c("CT", "AD")) %>%
                  subset(celltype %in% levels(ests_alg$celltype)) # Gets rid of added cell types from ROSMAP IHC
    ests_ad$celltype <- factor(ests_ad$celltype)

    significant <- lapply(unique(ests_ad$name), function(param_set) {
      ests_param <- subset(ests_ad, name == param_set)
      anov <- aov(pct_est ~ diagnosis*celltype*tissue, data = ests_param)
      tuk <- TukeyHSD(anov, "diagnosis:celltype:tissue")

      ct_v_tissue <- expand.grid(levels(ests_ad$celltype), levels(ests_ad$tissue))
      ct_v_tissue <- paste(ct_v_tissue$Var1, ct_v_tissue$Var2, sep = ":")

      comparisons <- paste0("CT:", ct_v_tissue, "-AD:", ct_v_tissue)
      tuk <- as.data.frame(tuk[[1]][comparisons,])
      tuk$p_adj <- tuk$`p adj`
      tuk$celltype <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,2]
      tuk$tissue <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,5]
      tuk$significant <- tuk$p_adj <= 0.05
      tuk$name <- param_set
      return(tuk)
    })

    significant <- do.call(rbind, significant) %>%
                      select(p_adj, celltype, tissue, significant, name) #%>%
                      #merge(n_subjects, by = c("celltype", "tissue"))

    mean_props <- ests_ad %>%
                    group_by(name, algorithm, diagnosis, tissue, celltype) %>%
                    summarise(mean_pct = mean(pct_est),
                              sd_pct = sd(pct_est),
                              count = n(),
                              .groups = "drop")


    # Make one column for AD and one for CT for each of mean_pct, sd_pct, and count,
    # calculate fold-change between AD and CT means
    mean_props <- pivot_wider(mean_props, names_from = "diagnosis",
                              values_from = c("mean_pct", "sd_pct", "count")) %>%
                    mutate(fc = mean_pct_AD / mean_pct_CT,
                           log2_fc = log2(fc))
                           #sd_fc = calc_sd_fc(mean_pct_AD,))

    mean_props <- merge(mean_props, significant, by = c("name", "celltype", "tissue"))
    return(mean_props)
  })

  return(do.call(rbind, props_b))
})

mean_props_all <- do.call(rbind, mean_props_all) %>%
                    subset(mean_pct_AD != 0 & mean_pct_CT != 0) %>%
                    mutate(log_p = log(p_adj + 1e-8)) # pseudocount necessary to avoid log(0)

combined_stats <- mean_props_all %>% group_by(celltype, tissue) %>%
                    summarise(mean_fc = mean(fc),
                              sd_fc = sd(fc),
                              mean_log2fc = mean(log2_fc),
                              sd_log2fc = sd(log2_fc),
                              direction = mean(sign(log2_fc)),
                              mean_p = mean(p_adj),
                              geo_mean_p = exp(mean(log_p)),
                              .groups = "drop")

print(subset(combined_stats, geo_mean_p <= 0.05))


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

library(meta)

mean_astro <- subset(mean_props_all, celltype == "Astro")
mean_astro_tcx <- subset(mean_astro, tissue == "TCX")
mean_astro_cbe <- subset(mean_astro, tissue == "CBE")

res_astro_tcx <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = name,
                          data = mean_astro_tcx)

res_astro_cbe <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = name,
                          data = mean_astro_cbe,
                          sm = "SMD")


mean_props2 <- mean_props_all
tmp = str_split(mean_props2$name, "_", simplify = TRUE)[,1:2]
mean_props2$group <- paste(tmp[,1], tmp[,2], sep = "_")

mean_ratios <- mean_props2 %>% group_by(group, tissue, celltype) %>%
                summarise(mean_ratio = mean(mean_pct_AD / mean_pct_CT),
                          sd_ratio = sd(mean_pct_AD / mean_pct_CT),
                          n = mean(count_AD) + mean(count_CT),
                          .groups = "drop")

mean_astro <- subset(mean_ratios, celltype == "Astro")
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
#to_keep <- unique(subset(errs_tissue, param_id %in% best_params_sub$param_id)$name)

#mean_props_sub <- subset(mean_props_all, name %in% best_errors)
mean_props_sub <- mean_props_all

mean_astro <- subset(mean_props_sub, celltype == "Astrocyte")
mean_astro_tcx <- subset(mean_astro, tissue == "TCX")
mean_astro_cbe <- subset(mean_astro, tissue == "CBE")

res_astro_tcx <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = name,
                          data = mean_astro_tcx,
                          sm = "ROM") # Ratio of means

res_astro_cbe <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = name,
                          data = mean_astro_cbe,
                          sm = "ROM")


##### End #####

reference <- str_split(mean_astro_cbe$name, "_", simplify = TRUE)[,2]

res_astro_cbe_meta <- metareg(res_astro_cbe, ~reference)

mean_props_sub$reference <- str_split(mean_props_sub$name, "_", simplify = TRUE)[,2]

mean_astro_cain_cbe <- subset(mean_props_sub, celltype == "Astro" & tissue == "CBE" & reference == "cain")
res_astro_cain_cbe <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                          n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                          studlab = name,
                          data = mean_astro_cain_cbe,
                          sm = "ROM")

mean_astro_cain_tcx <- subset(mean_props_sub, celltype == "Astro" & tissue == "TCX" & reference == "cain")
res_astro_cain_tcx <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                               n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                               studlab = name,
                               data = mean_astro_cain_tcx,
                               sm = "ROM")


metas <- lapply(levels(mean_props_sub$celltype), function(ct) {
  res_tissue <- lapply(levels(mean_props_sub$tissue), function(tiss) {
    mean_props_tmp <- subset(mean_props_sub, celltype == ct & tissue == tiss)
    tryCatch({
      res <- metacont(n.e = count_AD, mean.e = mean_pct_AD, sd.e = sd_pct_AD,
                      n.c = count_CT, mean.c = mean_pct_CT, sd.c = sd_pct_CT,
                      studlab = name,
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
