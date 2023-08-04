library(Metrics)
library(nnls)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "morabito")#, "seaRef") #, "seaAD")


reference_data_name <- "cain"
test_data_name <- "donors"
granularity <- "broad"

data <- Load_AlgorithmInputData(reference_data_name, test_data_name,
                                granularity,
                                reference_input_type = "singlecell",
                                output_type = "cpm")

sig_full <- Load_SignatureMatrix(reference_data_name, granularity)

markers <- Load_Markers(reference_data_name, granularity, "dtangle", "diff", "pseudobulk")
sig_filt <- sig_full[unique(unlist(markers)),]

p_init <- runif(ncol(sig_filt))
p_init <- p_init / sum(p_init)
names(p_init) <- colnames(sig_filt)

celltypes <- colnames(sig_filt)
p <- p_init

# Gibbs sampling

tolerance <- 1e-8
last_error <- 0
continue <- TRUE
iter <- 1

while(continue == TRUE & iter <= 1000) {
  for (ind_vary in 1:(length(celltypes)-1)) {
    for (ind_adjust in setdiff(1:length(celltypes), ind_vary)) {
      #print(paste("Varying", ind_vary, ", adjusting", ind_adjust))
      ct_vary <- celltypes[ind_vary]
      ct_adjust <- celltypes[ind_adjust]
      ct_constants <- setdiff(celltypes, c(ct_vary, ct_adjust))

      p_const <- p[ct_constants]

      p_hat <- c(1-sum(p_const), p_const)
      names(p_hat)[1] <- ct_adjust

      if (p[ct_vary] == 0 & p[ct_adjust] == 0) {
        next
      }

      S_hat <- sig_filt[, c(ct_adjust, ct_constants)] %*% p_hat
      S_ct <- as.matrix(sig_filt[,ct_vary] - sig_filt[,ct_adjust])

      #M_obs <- assay(data$test, "counts")[rownames(sig_filt),1]

      res <- nnls(S_ct, as.matrix(M_obs - S_hat))

      # constraints
      new_p <- res$x
      new_p <- min(res$x, 1-sum(p_const))
      new_p <- max(new_p, 0)

      p[ct_vary] <- new_p
      p[ct_adjust] <- 1 - sum(p_const) - new_p
      #print(p)
    }
  }
  error <- rmse(M_obs, (sig_filt %*% p)[,1])

  print(paste("**** ITER", iter, ", ERROR: ", error, "****"))
  print(p)
  print(error - last_error)
  if (abs(error - last_error) < tolerance) {
    continue <- FALSE
  }
  last_error <- error
  iter <- iter + 1
}


truth <- data$test@metadata[["pctRNA"]][1,]

rmse(truth, p_init)
rmse(truth, p)

library(WGCNA)
plotDendroAndColors(cl2, as.numeric(factor(cl2$labels)), cex.dendroLabels = 0.3)
