library(Metrics)
library(nnls)
library(lme4)
library(Seurat)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "morabito")#, "seaRef") #, "seaAD")

reference_data_name <- "cain"
test_data_name <- "donors" #"ROSMAP"
granularity <- "broad"

DESeq2::DESeq()
DESeq2::estimateSizeFactorsForMatrix()

data <- Load_AlgorithmInputData(reference_data_name, test_data_name,
                                granularity,
                                reference_input_type = "singlecell",
                                output_type = "counts")

sig_full <- Load_SignatureMatrix(reference_data_name, granularity)
#markers <- Load_Markers(reference_data_name, granularity, "dtangle", "diff", "pseudobulk")
expr_obs <- assay(data$test, "counts")
sig_full <- sig_full[intersect(rownames(sig_full), rownames(expr_obs)),]

sig_filt <- FilterSignature(sig_full, 3, reference_data_name, granularity,
                            5, "dtangle", "diff", "pseudobulk")

expr_obs <- expr_obs[intersect(rownames(sig_filt), rownames(expr_obs)),]

# ROSMAP stuff
if (test_data_name == "ROSMAP") {
  ihc_props <- as.matrix(read.csv(file_rosmap_ihc_proportions, row.names = 1))
  #rownames(ihc_props) <- make.names(rownames(ihc_props))
  A <- Load_AvgLibSize(reference_data_name, granularity)
  A2 <- A
  # TODO this isn't quite right and the A matrix needs to be re-processed
  A2["Neuro"] <- mean(c(A2["Exc"], A2["Inh"]))
  A2["Oligo"] <- mean(c(A2["Oligo"], A2["OPC"]))
  A2 <- A2[colnames(ihc_props)]
  A2 <- A2 / sum(A2)
  ihc_pct <- ConvertPropCellsToPctRNA(ihc_props, A2)

  specimen2projid <- subset(colData(data$test), donor %in% rownames(ihc_props) &
                              specimenID %in% colnames(expr_obs) &
                              tissue == "DLPFC")

  expr_obs <- expr_obs[,as.character(specimen2projid$specimenID)]

  meta <- read.table("data/input/rosmap_raw/ROSMAP_Covariates_ages_censored.tsv",
                     header = TRUE)
  meta$specimenID <- make.names(meta$specimenID)
  rownames(meta) <- meta$specimenID

  meta <- meta[colnames(expr_obs),]
  meta$final_batch <- factor(meta$final_batch)
}

p_init <- matrix(runif(ncol(expr_obs) * ncol(sig_full)), nrow = ncol(expr_obs))
rownames(p_init) <- colnames(expr_obs)
colnames(p_init) <- colnames(sig_full)

p_init <- sweep(p_init, 1, rowSums(p_init), "/")

celltypes <- colnames(p_init)
p <- p_init

#S <- matrix(nrow = nrow(expr_obs), ncol = ncol(p_init))
#rownames(S) <- rownames(expr_obs)
#colnames(S) <- colnames(p_init)
S <- sig_filt

tolerance <- 1e-8
last_error <- 10000
continue <- TRUE
iter <- 1

# This loop doesn't out-perform one run of nnls(S,expr_obs) on either broad or
# fine cell types for rmse(p, truth). It does, however, result in lower error
# for rmse(expr_obs, S %*% t(p)). There also seems to be fewer 0 values in the
# final result.
while(continue == TRUE & iter <= 1000) {

  # Primal

  # Diff between nnls and lm seems trivial for p
  for (d in colnames(expr_obs)) {
    res_p <- nnls(S, expr_obs[,d]) # Could also include other technical covariates
    p[d,] <- res_p$x / sum(res_p$x) # Norm to 1 may not matter much?
  }
  #res_p <- lm(expr_obs ~ 0+., data = as.data.frame(S))
  #p <- t(coef(res_p))
  #p[p<0] <- 0
  #p <- sweep(p, 1, rowSums(p), "/")

  # Dual

  for (g in rownames(expr_obs)) {
    res_sig <- nnls(p, expr_obs[g,])
    S[g,] <- res_sig$x
  }

  # Letting S go negative seems to perform better than nnls -- but not for fine
  # cell types. Fine cell types produces some NAs.
  #res_sig <- lm(t(expr_obs) ~ 0+., data = as.data.frame(p))
  #S <- t(coef(res_sig))

  error <- rmse(expr_obs, S %*% t(p))

  print(paste("**** ITER", iter, ", ERROR: ", error, "****"))
  print(p[1:3,])
  print(error - last_error)
  if (abs(error - last_error) < tolerance || error > last_error) {
    continue <- FALSE
  }
  last_error <- error
  iter <- iter + 1
}

#truth <- data$test@metadata[["pctRNA"]]

# RMSE for ihc_props is usually lower than ihc_pct... hmm.
truth <- ihc_props[as.character(specimen2projid$donor),]

p2 <- as.data.frame(p) %>% mutate(Neuro = Exc + Inh, Endo = Endo + Peri, Oligo = Oligo + OPC) %>%
  select(-Exc,-Inh,-OPC,-Peri) %>% as.matrix()

rmse(p2, truth)


# With other covariates

# Initial p estimate
props <- table(data$reference$donor, data$reference$broadcelltype)
props <- sweep(props, 1, rowSums(props), "/")

noise <- sapply(1:ncol(sig_filt), function(X) {rnorm(ncol(expr_obs), sd = sd(props[,X])*2)})

p_init <- abs(sweep(noise, 2, colMeans(props), "+"))
p_init <- sweep(p_init, 1, rowSums(p_init), "/")
colnames(p_init) <- colnames(sig_filt)
rownames(p_init) <- colnames(expr_obs)

p <- p_init

tolerance <- 1e-8
last_error <- 10000
continue <- TRUE
iter <- 1

libSize <- colSums(assay(data$test, "counts"))
libSize <- libSize[colnames(expr_obs)]
#libSize <- libSize / max(libSize)

while(continue == TRUE & iter <= 1000) {

  # Primal

  # Including final batch seems to mess up Peri coefficient (NA)
  #covar_cols <- c("final_batch", "pmi", "RIN", "RnaSeqMetrics_PCT_CODING_BASES")
  covar_cols <- c("rnaBatch", "pmi", "medianUMIs")
  #m2 <- meta %>% select(all_of(covar_cols))
  m2 <- sc_meta %>% select(all_of(covar_cols))
  m2 <- cbind(m2, p)
  m2$pmi <- m2$pmi / 24 # values between 0 and 1
  #m2$RIN <- m2$RIN / 10 # Same
  m2$medianUMIs <- m2$medianUMIs / max(m2$medianUMIs)

  #form <- paste0("y ~ 1 + (1|final_batch) + ",
  #               paste(setdiff(colnames(m2), "final_batch"), collapse = " + "))
  form <- paste0("y ~ 0 + (1|rnaBatch) + ",
                 paste(setdiff(colnames(m2), "rnaBatch"), collapse = " + "))

  cc <- apply(expr_obs, 1, function(row) {
    y <- as.numeric(row)
    #res_glm <- glm(row ~ 0 + ., family = poisson(), data = m2, offset = libSize)
    res_glm <- glmer(as.formula(form), family = poisson, data = m2, offset = log(libSize))
    #res_glm <- lmer(as.formula(form), data = m2)
    #coef(res_glm)
    #rand_effects <- as.numeric(t(ranef(res_glm)$final_batch))
    rand_effects <- as.numeric(t(ranef(res_glm)$rnaBatch))
    #names(rand_effects) <- paste0("final_batch", rownames(ranef(res_glm)$final_batch))
    names(rand_effects) <- paste0("rnaBatch", rownames(ranef(res_glm)$rnaBatch))

    c(rand_effects, fixef(res_glm))
  })

  #res_cov <- lm(t(expr_obs) ~ 0 + ., data = m2)
  #coefs <- t(coef(res_cov))
  coefs <- t(cc)

  mod <- model.matrix(~ 0 + ., data = m2 %>% select(-all_of(celltypes)))
  adjust <- t(mod %*% t(coefs[,colnames(mod)])) + log(libSize)

  expr_new <- as.matrix(expr_obs - exp(adjust))
  #expr_new <- as.matrix(log1p(expr_obs) - t(adjust))
  #expr_new <- as.matrix(expr_obs - t(adjust))
  S <- coefs[,colnames(p)]

  # Dual

  # Diff between nnls and lm seems trivial for p
  for (d in colnames(expr_obs)) {
    res_p <- nnls(S, expr_new[,d])
    p[d,] <- res_p$x / sum(res_p$x) # Norm to 1 may not matter much?
  }

  error <- rmse(expr_new, S %*% t(p))

  print(paste("**** ITER", iter, ", ERROR: ", error, "****"))
  print(p[1:3,])
  print(error - last_error)
  if (abs(error - last_error) < tolerance || error > last_error) {
    continue <- FALSE
  }
  last_error <- error
  iter <- iter + 1
}

truth <- ihc_props[as.character(specimen2projid$donor),]
truth2 <- ihc_pct[as.character(specimen2projid$donor),]

p2 <- as.data.frame(p) %>% mutate(Neuro = Exc + Inh, Endo = Endo + Peri, Oligo = Oligo + OPC) %>%
  select(-Exc,-Inh,-OPC,-Peri) %>% as.matrix()
p2 <- p2[,colnames(ihc_props)]

rmse(p2, truth)
rmse(p2, truth2)


# Try covariates outside the loop

covar_cols <- c("rnaBatch", "pmi", "medianUMIs")
#m2 <- meta %>% select(all_of(covar_cols))
m2 <- sc_meta %>% select(all_of(covar_cols))
m2$pmi <- m2$pmi / 24 # values between 0 and 1
#m2$RIN <- m2$RIN / 10 # Same
m2$medianUMIs <- m2$medianUMIs / max(m2$medianUMIs)

#form <- paste0("y ~ 1 + (1|final_batch) + ",
#               paste(setdiff(colnames(m2), "final_batch"), collapse = " + "))
form <- paste0("y ~ 1 + (1|rnaBatch) + ",
               paste(setdiff(colnames(m2), "rnaBatch"), collapse = " + "))

cc <- apply(expr_obs, 1, function(row) {
  y <- as.numeric(row)
  #res_glm <- glm(row ~ 0 + ., family = poisson(), data = m2)
  res_glm <- glmer(as.formula(form), family = poisson, data = m2, offset = log(libSize))
  #res_glm <- lmer(as.formula(form), data = m2)
  #coef(res_glm)
  #rand_effects <- as.numeric(t(ranef(res_glm)$final_batch))
  rand_effects <- as.numeric(t(ranef(res_glm)$rnaBatch))
  #names(rand_effects) <- paste0("final_batch", rownames(ranef(res_glm)$final_batch))
  names(rand_effects) <- paste0("rnaBatch", rownames(ranef(res_glm)$rnaBatch))

  c(rand_effects, fixef(res_glm))
})

#res_cov <- lm(t(expr_obs) ~ 0 + ., data = m2)
#coefs <- t(coef(res_cov))
coefs <- t(cc)

mod <- model.matrix(~ 0 + ., data = m2)
mod <- model.matrix(~ 1 + ., data = m2)
#adjust <- t(mod %*% t(coefs[,colnames(mod)]))
adjust <- t(mod %*% t(coefs[,colnames(mod)])) + log(libSize)

expr_new <- as.matrix(expr_obs - exp(adjust))
expr_new <- sweep(expr_new, 2, libSize, "/") * 1e6
#expr_new <- expm1(expr_new)
#expr_new <- as.matrix(expr_obs - t(adjust))

p <- p_init

tolerance <- 1e-8
last_error <- 10000
continue <- TRUE
iter <- 1

libSize <- colSums(assay(data$test, "counts"))
libSize <- libSize[colnames(expr_obs)]
#libSize <- libSize / max(libSize)

while(continue == TRUE & iter <= 1000) {

  # Primal

  #for (g in rownames(expr_new)) {
  #  res_sig <- nnls(p, expr_new[g,])
  #  S[g,] <- res_sig$x
  #}

  # Letting S go negative seems to perform better than nnls -- but not for fine
  # cell types. Fine cell types produces some NAs.
  res_sig <- lm(t(expr_new) ~ 0+., data = as.data.frame(p))
  S <- t(coef(res_sig))


  # Dual

  # Diff between nnls and lm seems trivial for p
  for (d in colnames(expr_obs)) {
    res_p <- nnls(S, expr_new[,d])
    p[d,] <- res_p$x / sum(res_p$x) # Norm to 1 may not matter much?
  }

  error <- rmse(expr_new, S %*% t(p))

  print(paste("**** ITER", iter, ", ERROR: ", error, "****"))
  print(p[1:3,])
  print(error - last_error)
  if (abs(error - last_error) < tolerance || error > last_error) {
    continue <- FALSE
  }
  last_error <- error
  iter <- iter + 1
}

truth <- data$test@metadata$pctRNA
rmse(truth, p)

truth <- ihc_props[as.character(specimen2projid$donor),]
truth2 <- ihc_pct[as.character(specimen2projid$donor),]

p2 <- as.data.frame(p) %>% mutate(Neuro = Exc + Inh, Endo = Endo + Peri, Oligo = Oligo + OPC) %>%
  select(-Exc,-Inh,-OPC,-Peri) %>% as.matrix()
p2 <- p2[,colnames(ihc_props)]

rmse(p2, truth)
rmse(p2, truth2)


library(edgeR)

expr_obs <- assay(data$test, "counts")
expr_obs <- expr_obs[,as.character(specimen2projid$specimenID)]
y <- DGEList(expr_obs)
y <- calcNormFactors(y)

mod <- model.matrix(~ 1 + ., data = m2)

y <- estimateDisp(y, mod)
fit <- glmQLFit(y, mod)
coefs <- fit$coefficients
fitted <- fit$fitted.values

# This equals 'fitted.values'
offset <- log(y$samples$lib.size * y$samples$norm.factors)
adjust <- t(mod %*% t(coefs))
adjust <- sweep(adjust, 2, offset, "+")
adjust <- exp(adjust)

expr_new <- expr_obs - adjust
expr_new <- as.matrix(expr_new[rownames(sig_filt),])

S_adjusted <- sweep(sig_filt, 1, rowMeans(sig_filt), "-")

p <- p_init
S <- S_adjusted

# Diff between nnls and lm seems trivial for p
for (d in colnames(expr_obs)) {
  res_p <- nnls(S, expr_new[,d])
  p[d,] <- res_p$x / sum(res_p$x) # Norm to 1 may not matter much?
}


res_sig <- lm(t(expr_new) ~ 0+., data = as.data.frame(p))
S <- t(coef(res_sig))

tolerance <- 1e-8
last_error <- 10000
continue <- TRUE
iter <- 1

while(continue == TRUE & iter <= 1000) {

  # Primal

  # Letting S go negative seems to perform better than nnls -- but not for fine
  # cell types. Fine cell types produces some NAs.
  res_sig <- lm(t(expr_new) ~ 0+., data = as.data.frame(p))
  S <- t(coef(res_sig))


  # Dual

  # Diff between nnls and lm seems trivial for p
  for (d in colnames(expr_obs)) {
    res_p <- nnls(S, expr_new[,d])
    p[d,] <- res_p$x / sum(res_p$x) # Norm to 1 may not matter much?
  }

  error <- rmse(expr_new, S %*% t(p))

  print(paste("**** ITER", iter, ", ERROR: ", error, "****"))
  print(p[1:3,])
  print(error - last_error)
  if (abs(error - last_error) < tolerance || error > last_error) {
    continue <- FALSE
  }
  last_error <- error
  iter <- iter + 1
}

truth <- ihc_props[as.character(specimen2projid$donor),]
truth2 <- ihc_pct[as.character(specimen2projid$donor),]

p2 <- as.data.frame(p) %>% mutate(Neuro = Exc + Inh, Endo = Endo + Peri, Oligo = Oligo + OPC) %>%
  select(-Exc,-Inh,-OPC,-Peri) %>% as.matrix()
p2 <- p2[,colnames(ihc_props)]

rmse(p2, truth)
rmse(p2, truth2)

