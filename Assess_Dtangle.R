library(dtangle)
library(Matrix)
library(SummarizedExperiment)
library(dplyr)
library(Metrics)
library(stringr)
library(ggplot2)

source("Filenames.R")

####for each bulk data set, pick the best parameter set with those genes, based on correlation with pseudobulk "ground truth"###

datasets = list("mathys")#,"cain","lau","morabito","lengSFG","lengEC")

for (dataset in datasets) { #(datasource in c("msbb","mayotcx","mayocbe","mayotcxcbe","rosmap")) {
  dtangle_list <- readRDS(file.path(dir_output, paste0("dtangle_list_", dataset, "_training_logcpm_broad.rds"))) #paste0("dtangle_list_", dataset, "_donors_broad.rds")))
  se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", dataset, "_broadcelltypes.rds"))) #paste0("pseudobulk_", dataset, "_bydonor_broadcelltypes.rds"))) ###pseudobulk ground truth and expression table

  propval <- as.matrix(colData(se))

  # Correlation broken down by cell type
  cor_celltype <- list()

  # Correlation broken down by subject/donor
  cor_subject <- list()

  rMSE <- list()
  meanAE <- list()
  meanAPE <- list()

  ###step 1: find correlations between predictions and pseudobulk proportions
  for (ii in names(dtangle_list)) {
    if (!(any(is.na(dtangle_list[[ii]]$estimates)))) {  ###exclude all predictions that return any NA values on the pseudobulk
      tmp <- dtangle_list[[ii]]
      tmp$estimates <- tmp$estimates[rownames(propval), colnames(propval)]
      cor_celltype[[ii]]<- diag(cor(tmp$estimates,propval))
      cor_subject[[ii]] <- diag(cor(t(tmp$estimates), t(propval)))
      rMSE[[ii]] <- rmse(propval, tmp$estimates)
      meanAE[[ii]] <- mae(propval, tmp$estimates)
      meanAPE[[ii]] <- smape(propval, tmp$estimates)
    }
  }

  cor_celltype <- do.call(rbind, cor_celltype)
  cor_subject <- do.call(rbind, cor_subject)
  errs <- cbind(do.call(rbind, rMSE), do.call(rbind, meanAE), do.call(rbind, meanAPE))
  colnames(errs) <- c("rMSE", "meanAE", "meanAPE")

  # Mean correlation
  celltype_mean <- rowSums(cor_celltype) / ncol(cor_celltype)
  subject_mean <- rowSums(cor_subject) / ncol(cor_subject)

  best_celltype <- top_n(as.data.frame(celltype_mean),
                         wt = celltype_mean, n = 10) %>%
                   arrange(desc(celltype_mean))

  best_subject <- top_n(as.data.frame(subject_mean),
                        wt = subject_mean, n = 10) %>%
                  arrange(desc(subject_mean))

  best_rmse <- top_n(as.data.frame(errs), wt = -rMSE, n = 10) %>% arrange(rMSE)
  best_meanae <- top_n(as.data.frame(errs), wt = -meanAE, n = 10) %>% arrange(meanAE)
  best_meanape <- top_n(as.data.frame(errs), wt = -meanAPE, n = 10) %>% arrange(meanAPE)

  for (K in rownames(best_celltype[1])) {
    matplot(propval, dtangle_list[[K]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
            xlab="Truth", ylab="Estimates", main = K)
    abline(a = 0, b = 1)
  }

  for (K in rownames(best_subject)[1]) {
    matplot(propval, dtangle_list[[K]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
            xlab="Truth", ylab="Estimates", main = K)
    abline(a = 0, b = 1)
  }

  for (K in rownames(best_rmse)[1]) {
    matplot(propval, dtangle_list[[K]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
            xlab="Truth", ylab="Estimates", main = K)
    abline(a = 0, b = 1)
  }

  for (K in rownames(best_meanae)[1]) {
    matplot(propval, dtangle_list[[K]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
            xlab="Truth", ylab="Estimates", main = K)
    abline(a = 0, b = 1)
  }

  for (K in rownames(best_meanape)[1]) {
    matplot(propval, dtangle_list[[K]]$estimates[rownames(propval),], xlim = c(0,1), ylim=c(0,1),
            xlab="Truth", ylab="Estimates", main = K)
    abline(a = 0, b = 1)
  }


  # Error by cell proportion
  meta <- sapply(names(dtangle_list), FUN = str_split, pattern = "_")
  meta <- as.data.frame(do.call(rbind, meta))
  meta <- meta[, c(1, 2, 4, 6, 8, 10, 12)]
  colnames(meta) <- c("Dataset", "Fineness", "Method", "Gamma",
                      "Summary.Function", "N.Markers", "Norm.Method")


  tmp <- dtangle_list[[14]]
  tmp$estimates <- tmp$estimates[rownames(propval), colnames(propval)]

  ct = "Astro"

  samps <- grep(ct, rownames(tmp$estimates), value = TRUE)
  props <- unique(unlist(lapply(str_split(samps, "_"), '[[', 2)))
  props_list <- lapply(props, FUN = function(X) {
    grep(paste0("_", X, "_"), samps, value = TRUE)
  })
  names(props_list) <- props

  stats <- list()
  ses <- list()
  sapes <- list()
  subs <- list()

  for (P in names(props_list)) {
    ests <- tmp$estimates[props_list[[P]], ]
    truth <- propval[props_list[[P]], ]

    # TODO should we only consider the cell type column, or use all columns?
    stats[[P]] <- c(
      rmse(truth, ests),
      mae(truth, ests),
      smape(truth, ests),
      mean(cor_subject[14, props_list[[P]]])
    )

    ses[[P]] <- sapply(1:nrow(truth), FUN = function(X) {mse(truth[X,], ests[X,])})
    sapes[[P]] <- sapply(1:nrow(truth), FUN = function(X) {smape(truth[X,], ests[X,])})
    #2 * abs(truth[X,] - ests[X,]) / (abs(truth[X,]) + abs(ests[X,]))
    subs[[P]] <- cor_subject[14, props_list[[P]]]
  }

  stats <- as.data.frame(do.call(rbind, stats))
  colnames(stats) <- c("RMSE", "mAE", "SMAPE", "Cor")
  stats$Proportion <- rownames(stats)

  p <- ggplot(stats, aes(x = Proportion, y = Cor)) + geom_point()
  p

  sapes <- melt(sapes)
  colnames(sapes) <- c("RPE", "Proportion")
  sapes$Type <- rep(c(rep(1,10), rep(2,10)), nrow(sapes) / 20)

  p <- ggplot(sapes, aes(x = Proportion, y = RPE, color = Type)) + geom_point()
  p

  subs <- melt(subs)
  colnames(subs) <- c("Cor", "Proportion")
  subs$Type <- rep(c(rep(1,10), rep(2,10)), nrow(sapes) / 20)

  p <- ggplot(subs, aes(x = Proportion, y = Cor, color = Type)) + geom_point()
  p



  keeppars_realbulk=gsub("pseudo","bulk",keeppars) ###best parameter sets, but now applied to bulk

  # I think this is for the datasets with unknown proportions, using best parameter sets?
  ###plot predicted distributions for cell types over all real bulk samples, using the top 10 best parameter sets
  pdf(paste0("proportion_boxplots_broad_",datasource,"_mathys.pdf"),width=16)
  for (ii in 1:10) {
    tempmat=dtangle_list[[keeppars_realbulk[ii]]]$estimates
    if (!(any(is.na(tempmat)))) {
      boxplot(c(tempmat)~rep(colnames(dtangle_list[[keeppars_realbulk[ii]]]$estimates),each=nrow(dtangle_list[[keeppars_realbulk[ii]]]$estimates)),
              xlab="",ylab="Predicted Proportion",main=keeppars_realbulk[ii])
    }
  }
  dev.off()
}



###fine predictions###
for (datasource in c("msbb","mayotcx","mayocbe","mayotcxcbe","rosmap")) {
  load(paste0("dtangle_list_fine_",datasource,"_mathys.rda"))
  load("pseudobulk_mathys_finecelltypes_30percentlimit.rda")
  cormat=matrix(0,nrow=length(dtangle_list),ncol=ncol(dtangle_list[[1]]$estimates))
  colnames(cormat)=colnames(dtangle_list[[1]]$estimates)
  rownames(cormat)=names(dtangle_list)
  pseudobulk_elements=grep("pseudo",names(dtangle_list),val=T)
  for (ii in pseudobulk_elements) {
    if (!(any(is.na(dtangle_list[[ii]])))) {
      cormat[ii,]=diag(cor(dtangle_list[[ii]]$estimates,propval))
    }
  }

  corrows=apply(cormat,1,function(x){!(any(is.na(x)))})
  corrows=intersect(which(corrows),grep("pseudo",rownames(cormat)))
  minval=apply(cormat[corrows,],1,min)  ###find minimum correlation over all cell types within each parameter set's predictions
  minord=order(-minval)
  keeppars=rownames(cormat)[corrows[minord]] ###parameter set names ordered by minimum correlation value with ground truth


  keeppars_realbulk=gsub("pseudo","bulk",keeppars) ###best parameter sets, but now applied to bulk
  ###plot predicted distributions for cell types over all real bulk samples, using the top 10 best parameter sets
  pdf(paste0("proportion_boxplots_fine_",datasource,"_mathys.pdf"),width=16)
  for (ii in 1:10) {
    tempmat=dtangle_list[[keeppars_realbulk[ii]]]$estimates
    if (!(any(is.na(tempmat)))) {
      boxplot(c(tempmat)~rep(colnames(dtangle_list[[keeppars_realbulk[ii]]]$estimates),each=nrow(dtangle_list[[keeppars_realbulk[ii]]]$estimates)),
              xlab="",ylab="Predicted Proportion",main=keeppars_realbulk[ii])
    }
  }
  dev.off()
}


###Output tables of predictions with best prediction outputs###
load("dtangle_list_broad_msbb_mathys.rda")
outmat=dtangle_list[["broadbulk_0.1_auto_p.value_counts_counts_log"]]$estimates
load("dtangle_list_fine_msbb_mathys.rda")
outmat2=dtangle_list[["broadbulk_0.2_auto_ratio_counts_counts_log"]]$estimates
fullmat=cbind(outmat,outmat2)
write.csv(fullmat,file="dtangle_estimates_msbb_mathys.csv")

load("dtangle_list_broad_rosmap_mathys.rda")
outmat=dtangle_list[["broadbulk_0.15_auto_p.value_counts_counts_log"]]$estimates
load("dtangle_list_fine_rosmap_mathys.rda")
outmat2=dtangle_list[["broadbulk_0.2_auto_ratio_counts_counts_log"]]$estimates
fullmat=cbind(outmat,outmat2)
write.csv(fullmat,file="dtangle_estimates_rosmap_mathys.csv")

load("dtangle_list_broad_mayotcx_mathys.rda")
outmat=dtangle_list[["broadbulk_0.03_auto_ratio_cpm_cpm_log"]]$estimates
load("dtangle_list_fine_mayotcx_mathys.rda")
outmat2=dtangle_list[["broadbulk_0.03_auto_ratio_cpm_cpm_log"]]$estimates
fullmat=cbind(outmat,outmat2)
write.csv(fullmat,file="dtangle_estimates_mayotcx_mathys.csv")

load("dtangle_list_broad_mayocbe_mathys.rda")
outmat=dtangle_list[["broadbulk_0.02_auto_diff_counts_counts_log"]]$estimates
load("dtangle_list_fine_mayocbe_mathys.rda")
outmat2=dtangle_list[["broadbulk_0.15_auto_ratio_cpm_cpm_log"]]$estimates
fullmat=cbind(outmat,outmat2)
write.csv(fullmat,file="dtangle_estimates_mayocbe_mathys.csv")


###output dtangle models for best paramater set - can extract gene sets from these models for other purposes
load("dtangle_list_broad_msbb_mathys.rda")
dtangle_model=dtangle_list[["broadbulk_0.1_auto_p.value_counts_counts_log"]]
save(dtangle_model,file="dtangle_model_broadtypes_msbb_mathys.rda")
load("dtangle_list_fine_msbb_mathys.rda")
dtangle_model=dtangle_list[["broadbulk_0.2_auto_ratio_counts_counts_log"]]
save(dtangle_model,file="dtangle_model_finetypes_msbb_mathys.rda")

load("dtangle_list_broad_rosmap_mathys.rda")
dtangle_model=dtangle_list[["broadbulk_0.15_auto_p.value_counts_counts_log"]]
save(dtangle_model,file="dtangle_model_broadtypes_rosmap_mathys.rda")
load("dtangle_list_fine_rosmap_mathys.rda")
dtangle_model=dtangle_list[["broadbulk_0.2_auto_ratio_counts_counts_log"]]
save(dtangle_model,file="dtangle_model_finetypes_rosmap_mathys.rda")

load("dtangle_list_broad_mayotcx_mathys.rda")
dtangle_model=dtangle_list[["broadbulk_0.03_auto_ratio_cpm_cpm_log"]]
save(dtangle_model,file="dtangle_model_broadtypes_mayoctx_mathys.rda")
load("dtangle_list_fine_mayotcx_mathys.rda")
dtangle_model=dtangle_list[["broadbulk_0.03_auto_ratio_cpm_cpm_log"]]
save(dtangle_model,file="dtangle_model_finetypes_mayoctx_mathys.rda")

load("dtangle_list_broad_mayocbe.rda")
dtangle_model=dtangle_list[["broadbulk_0.02_auto_diff_counts_counts_log"]]
save(dtangle_model,file="dtangle_model_broadtypes_mayocbe.rda")
load("dtangle_list_fine_mayocbe.rda")
dtangle_model=dtangle_list[["broadbulk_0.15_auto_ratio_cpm_cpm_log"]]
save(dtangle_model,file="dtangle_model_finetypes_mayocbe.rda")
