library(Matrix)
library(stringr)
library(dplyr)
library(scuttle)
library(Metrics)

set.seed(12345)
sig <- matrix(c(runif(3, min = 7000, max = 10000), # celltype 1, genes 1-3 are markers
                runif(6, min = 0, max = 1000), # low expression
                runif(3, min = 2000, max = 5000), # not markers for anything

                runif(3, min = 0, max = 2000), # celltype 2, low expression 
                runif(3, min = 7000, max = 10000), # genes 4-6 are markers
                runif(3, min = 0, max = 1000), # low expression 
                runif(3, min = 2000, max = 5000), # not markers
              
                runif(6, min = 0, max = 1000), # celltype 3, low expression 
                runif(3, min = 7000, max = 10000), # genes 7-9 are markers
                runif(3, min = 2000, max = 5000)), # not markers
              nrow = 12) 
              
sig <- calculateCPM(sig) / 1e3 # More manageable scale
colnames(sig) <- paste("CT", 1:3, sep = "_")
rownames(sig) <- paste("G", 1:12, sep = "_")

set.seed(6789)
pcts <- matrix(c(runif(5, min = 0.1, max = 0.4), # celltype 1, 5 samples
                 runif(5, min = 0, max = 0.1), # celltype 2
                 runif(5, min = 0.5, max = 1)), # celltype 3
               nrow = 5)

pcts <- sweep(pcts, 1, rowSums(pcts), "/")
colnames(pcts) <- colnames(sig)
rownames(pcts) <- paste("S", 1:5, sep = "_")

# ground "truth"
meas_expr <- sig %*% t(pcts)

# pretend estimated signature
set.seed(123)
noise <- matrix(runif(36, min = -100, max = 100), nrow = 12)
sig_mod <- sig + noise
sig_mod[sig_mod < 0] <- 0
sig_mod <- calculateCPM(sig_mod) / 1e3

# pretend estimated percents
set.seed(456)
pcts_est <- matrix(c(runif(5, min = 0, max = 0.6), # wider ranges
                     runif(5, min = 0, max = 0.3), 
                     runif(5, min = 0, max = 1)),
                   nrow = 5)

pcts_est <- sweep(pcts_est, 1, rowSums(pcts_est), "/")
colnames(pcts_est) <- colnames(sig)
rownames(pcts_est) <- paste("S", 1:5, sep = "_")

# bad estimate
set.seed(789)
pcts_bad <- matrix(c(runif(5, min = 0, max = 0.1), # bad ranges
                     runif(5, min = 0.8, max = 1), 
                     runif(5, min = 0, max = 0.3)),
                   nrow = 5)

pcts_bad <- sweep(pcts_bad, 1, rowSums(pcts_bad), "/")
colnames(pcts_bad) <- colnames(sig)
rownames(pcts_bad) <- paste("S", 1:5, sep = "_")

print_errors <- function(truth, est) {
  df <- data.frame("cor" = mean(diag(cor(log2(truth+1),
                                         log2(est+1), 
                                         use = "na.or.complete"))),
                   "rmse" = rmse(truth, est),
                   "mape" = smape(truth, est) / 2)
  print(df)
}

# Error calculations -- as-is
print_errors(meas_expr, sig %*% t(pcts)) # 1 0 0
print_errors(meas_expr, sig_mod %*% t(pcts)) # 0.82 22.2 0.21
print_errors(meas_expr, sig %*% t(pcts_est)) # 0.9 24.2 0.13
print_errors(meas_expr, sig_mod %*% t(pcts_est)) # 0.77 25.6 0.19
print_errors(meas_expr, sig %*% t(pcts_bad)) # -0.54 104.5 0.4
print_errors(meas_expr, sig_mod %*% t(pcts_bad)) # -0.17 82.8 0.43

# lm - linear scale
fit1 <- lm(t(meas_expr) ~ 0 + pcts)
print(max(abs(sig - t(coef(fit1))))) # close to 0
print_errors(meas_expr, t(fitted(fit1))) # 1 ~0 ~0 

fit2 <- lm(t(meas_expr) ~ 0 + pcts_est)
print(max(abs(sig - t(coef(fit2))))) # 528
print_errors(meas_expr, t(fitted(fit2))) # 0.996 3.9 0.03

fit3 <- lm(t(meas_expr) ~ 0 + pcts_bad)
print(max(abs(sig - t(coef(fit3))))) # 576
print_errors(meas_expr, t(fitted(fit3))) # 0.98 7.0 0.05

# lm - log scale
fit1 <- lm(t(log2(meas_expr+1)) ~ 0 + pcts)
print_errors(meas_expr, t(2^fitted(fit1)-1)) # 0.999 0.68 0.005

fit2 <- lm(t(log2(meas_expr+1)) ~ 0 + pcts_est)
print_errors(meas_expr, t(2^fitted(fit2)-1)) # 0.996 4.0 0.03

fit3 <- lm(t(log2(meas_expr+1)) ~ 0 + pcts_bad) 
print_errors(meas_expr, t(2^fitted(fit3)-1)) # 0.98 7.0 0.05


# Add noise to measurement
set.seed(12)
noise <- matrix(rnorm(60, sd = 10), nrow = 12)
meas_noise <- meas_expr + noise
meas_noise[meas_noise < 0] <- 0
meas_noise <- calculateCPM(meas_noise) / 1e3

print_errors(meas_noise, sig %*% t(pcts)) # 0.96 8.6 0.09
print_errors(meas_noise, sig_mod %*% t(pcts)) # 0.78 23.1 0.22
print_errors(meas_noise, sig %*% t(pcts_est)) # 0.83 26.8 0.18
print_errors(meas_noise, sig_mod %*% t(pcts_est)) # 0.72 27.4 0.22
print_errors(meas_noise, sig %*% t(pcts_bad)) # -0.51 105.9 0.51
print_errors(meas_noise, sig_mod %*% t(pcts_bad)) # -0.17 84.8 0.43

# lm - linear scale
fit1 <- lm(t(meas_noise) ~ 0 + pcts)
print(max(abs(sig - t(coef(fit1))))) # 232
print_errors(meas_noise, t(fitted(fit1))) # 0.983 6.4 0.0768

fit2 <- lm(t(meas_noise) ~ 0 + pcts_est)
print(max(abs(sig - t(coef(fit2))))) # 493
print_errors(meas_noise, t(fitted(fit2))) # 0.985 6.5 0.0772

fit3 <- lm(t(meas_noise) ~ 0 + pcts_bad)
print(max(abs(sig - t(coef(fit3))))) # 642
print_errors(meas_noise, t(fitted(fit3))) # 0.95 10.3 0.10

# lm - log scale
fit1 <- lm(t(log2(meas_noise+1)) ~ 0 + pcts)
print_errors(meas_noise, t(2^fitted(fit1)-1)) # 0.978 7.1 0.084

fit2 <- lm(t(log2(meas_noise+1)) ~ 0 + pcts_est)
print_errors(meas_noise, t(2^fitted(fit2)-1)) # 0.984 6.8 0.085

fit3 <- lm(t(log2(meas_noise+1)) ~ 0 + pcts_bad) 
print_errors(meas_noise, t(2^fitted(fit3)-1)) # 0.941 11.5 0.116


# Add noise to measurement in a batch-dependent manner
set.seed(34)
noise <- cbind(matrix(rnorm(36, sd = 5), nrow = 12),
               matrix(rnorm(24, sd = 20), nrow = 12))
meas_batch <- meas_expr + noise
meas_batch[meas_batch < 0] <- 0
meas_batch <- calculateCPM(meas_batch) / 1e3

batch <- factor(paste("B", c(1, 1, 1, 2, 2), sep = "_"))

print_errors(meas_batch, sig %*% t(pcts)) # 0.965 13.1 0.0994
print_errors(meas_batch, sig_mod %*% t(pcts)) # 0.786 26.9 0.228
print_errors(meas_batch, sig %*% t(pcts_est)) # 0.901 24.1 0.175
print_errors(meas_batch, sig_mod %*% t(pcts_est)) # 0.750 27.3 0.227
print_errors(meas_batch, sig %*% t(pcts_bad)) # -0.571 104.8 0.524
print_errors(meas_batch, sig_mod %*% t(pcts_bad)) # -0.208 82.3 0.449

# lm - linear scale w/o batch
fit1 <- lm(t(meas_batch) ~ 0 + pcts)
print(max(abs(sig - t(coef(fit1))))) # 250
print_errors(meas_batch, t(fitted(fit1))) # 0.981 9.2 0.0933

fit2 <- lm(t(meas_batch) ~ 0 + pcts_est)
print(max(abs(sig - t(coef(fit2))))) # 765
print_errors(meas_batch, t(fitted(fit2))) # 0.987 7.97 0.0791

fit3 <- lm(t(meas_batch) ~ 0 + pcts_bad)
print(max(abs(sig - t(coef(fit3))))) # 560
print_errors(meas_batch, t(fitted(fit3))) # 0.969 11.7 0.107

# lm - linear scale with batch
fit1 <- lm(t(meas_batch) ~ 0 + pcts + batch)
print_errors(meas_batch, t(fitted(fit1))) # 0.983 4.8 0.0658

fit2 <- lm(t(meas_batch) ~ 0 + pcts_est + batch)
print_errors(meas_batch, t(fitted(fit2))) # 0.995 3.6 0.0422

fit3 <- lm(t(meas_batch) ~ 0 + pcts_bad + batch)
print_errors(meas_batch, t(fitted(fit3))) # 0.977 9.5 0.0783

# lm - log scale w/o batch
fit1 <- lm(t(log2(meas_batch+1)) ~ 0 + pcts)
print_errors(meas_batch, t(2^fitted(fit1)-1)) # 0.982 9.4 0.107

fit2 <- lm(t(log2(meas_batch+1)) ~ 0 + pcts_est)
print_errors(meas_batch, t(2^fitted(fit2)-1)) # 0.985 8.3 0.0935

fit3 <- lm(t(log2(meas_batch+1)) ~ 0 + pcts_bad)
print_errors(meas_batch, t(2^fitted(fit3)-1)) # 0.968 11.92 0.123

# lm - log scale with batch
fit1 <- lm(t(log2(meas_batch+1)) ~ 0 + pcts + batch)
print_errors(meas_batch, t(2^fitted(fit1)-1)) # 0.983 5.8 0.0723

fit2 <- lm(t(log2(meas_batch+1)) ~ 0 + pcts_est + batch)
print_errors(meas_batch, t(2^fitted(fit2)-1)) # 0.995 3.7 0.0422

fit3 <- lm(t(log2(meas_batch+1)) ~ 0 + pcts_bad + batch)
print_errors(meas_batch, t(2^fitted(fit3)-1)) # 0.976 9.7 0.0905

# lm - Using the signature to estimate
do_lm <- function(sig1, pcts1, meas) {
  est <- sig1 %*% t(pcts1)
  ff <- sapply(rownames(meas), function(G) {
    fit <- lm(meas[G,] ~ est[G,])
    return(fitted(fit))
  })
  ff <- t(ff)
  print_errors(meas, ff)
}

do_lm(sig, pcts, meas_expr) # 1 ~0 ~0
do_lm(sig_mod, pcts, meas_expr) # 0.994 2.3 0.0266
do_lm(sig, pcts_est, meas_expr) # 0.978 9.3 0.0758
do_lm(sig_mod, pcts_est, meas_expr) # 0.978 9.1 0.0698
do_lm(sig, pcts_bad, meas_expr) # 0.982 8.6 0.0713  # too good
do_lm(sig_mod, pcts_bad, meas_expr) # 0.982 8.5 0.0660 # too good

do_lm(sig, pcts, meas_noise) # 0.9687 7.5 0.0890
do_lm(sig_mod, pcts, meas_noise) # 0.9691 8.0 0.0975
do_lm(sig, pcts_est, meas_noise) # 0.945 11.9 0.124
do_lm(sig_mod, pcts_est, meas_noise) # 0.9397 11.5 0.118
do_lm(sig, pcts_bad, meas_noise) # 0.931 12.3 0.124
do_lm(sig_mod, pcts_bad, meas_noise) # 0.936 11.9 0.122

do_lm_log <- function(sig1, pcts1, meas) {
  est <- log2(sig1 %*% t(pcts1)+1)
  ff <- sapply(rownames(meas), function(G) {
    fit <- lm(log2(meas[G,]+1) ~ est[G,])
    return(fitted(fit))
  })
  ff <- t(2^ff-1)
  print_errors(meas, ff)
}


do_lm_log(sig, pcts, meas_expr) # 1 ~0 ~0 
do_lm_log(sig_mod, pcts, meas_expr) # 0.993 2.3 0.0256
do_lm_log(sig, pcts_est, meas_expr) # 0.9785 9.3 0.0751
do_lm_log(sig_mod, pcts_est, meas_expr) # 0.9793 9.2 0.0693
do_lm_log(sig, pcts_bad, meas_expr) # 0.983 8.6 0.0692  # Too good
do_lm_log(sig_mod, pcts_bad, meas_expr) # 0.9829 8.4 0.0642  # too good

do_lm_log(sig, pcts, meas_noise) # 0.971 7.5 0.0963
do_lm_log(sig_mod, pcts, meas_noise) # 0.968 8.2 0.0994
do_lm_log(sig, pcts_est, meas_noise) # 0.948 12.5 0.139
do_lm_log(sig_mod, pcts_est, meas_noise) # 0.931 12.5 0.133
do_lm_log(sig, pcts_bad, meas_noise) # 0.930 12.8 0.132  # too good
do_lm_log(sig_mod, pcts_bad, meas_noise) # 0.934 12.4 0.131  # too good

# batch
do_lm_batch <- function(sig1, pcts1, meas, batch1) {
  est <- sig1 %*% t(pcts1)
  ff <- sapply(rownames(meas), function(G) {
    fit <- lm(meas[G,] ~ est[G,] + batch1)
    return(fitted(fit))
  })
  ff <- t(ff)
  print_errors(meas, ff)
}

do_lm_log_batch <- function(sig1, pcts1, meas, batch1) {
  est <- log2(sig1 %*% t(pcts1)+1)
  ff <- sapply(rownames(meas), function(G) {
    fit <- lm(log2(meas[G,]+1) ~ est[G,] + batch1)
    return(fitted(fit))
  })
  ff <- t(2^ff-1)
  print_errors(meas, ff)
}

do_lm(sig, pcts, meas_batch) # 0.979 10.2 0.09416
do_lm(sig_mod, pcts, meas_batch) # 0.972 10.3 0.094196
do_lm(sig, pcts_est, meas_batch) # 0.969 12.9 0.10936
do_lm(sig_mod, pcts_est, meas_batch) # 0.971 12.8 0.108
do_lm(sig, pcts_bad, meas_batch) # 0.968 12.9 0.1094  # too good
do_lm(sig_mod, pcts_bad, meas_batch) # 0.968 12.6 0.106  # too good

do_lm_batch(sig, pcts, meas_batch, batch) # 0.980 6.5 0.0758
do_lm_batch(sig_mod, pcts, meas_batch, batch) # 0.977 6.7 0.0839
do_lm_batch(sig, pcts_est, meas_batch, batch) # 0.9726 10.8 0.102
do_lm_batch(sig_mod, pcts_est, meas_batch, batch) # 0.9716 10.8 0.0921
do_lm_batch(sig, pcts_bad, meas_batch, batch) # 0.9736 11.4 0.106  # too good
do_lm_batch(sig_mod, pcts_bad, meas_batch, batch) # 0.9739 11.1 0.09898  # too good


do_lm_log(sig, pcts, meas_batch) # 0.979 10.2 0.108
do_lm_log(sig_mod, pcts, meas_batch) # 0.971 10.5 0.105
do_lm_log(sig, pcts_est, meas_batch) # 0.969 13.1 0.128
do_lm_log(sig_mod, pcts_est, meas_batch) # 0.968 13.1 0.126
do_lm_log(sig, pcts_bad, meas_batch) # 0.966 13.2 0.1297  # too good
do_lm_log(sig_mod, pcts_bad, meas_batch) # 0.965 12.9 0.125  # too good

do_lm_log_batch(sig, pcts, meas_batch, batch) # 0.976 7.5 0.0858
do_lm_log_batch(sig_mod, pcts, meas_batch, batch) # 0.967 7.7 0.0972
do_lm_log_batch(sig, pcts_est, meas_batch, batch) # 0.971 11.1 0.115
do_lm_log_batch(sig_mod, pcts_est, meas_batch, batch) # 0.968 11.1 0.108
do_lm_log_batch(sig, pcts_bad, meas_batch, batch) # 0.9727 11.6 0.112  # too good
do_lm_log_batch(sig_mod, pcts_bad, meas_batch, batch) # 0.9737 11.3 0.107  # too good
