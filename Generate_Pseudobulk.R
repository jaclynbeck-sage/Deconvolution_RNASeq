#### Systematic generation of pseudobulks from clusters ###
# TODO turn the two for loops that generate random broad/fine samples into a
# single function. Can make a new metadata DF that contains "cell" and "type"
# that is populated by meta$broadcelltype or meta$subcluster as appropriate,
# and pass it into this function.
# Also pass in: fullmat, ints, numreps, numcells
# Can return pseudobulk and propval as a list, or just save to file
# Make a separate pseudobulk matrix for full data and donor data, since
# this don't change based on broad/fine cell types. Then it can cbind to the
# pb matrix returned by the function.
# Sacrificed readability for speed. Doing matrix multiplication is >10x faster
# than using rowSums().

library(Matrix)
library(stringr)

input_dir <- file.path("DeconvolutionData", "input")
output_dir <- file.path("DeconvolutionData", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- list("mathys") #, "cain", "lau", "leng_SFG", "leng_EC", "lau", "morabito")

for (dataset in datasets) {
  load(file.path(input_dir, paste0(dataset,"_counts.rda")))
  meta <- read.csv(file.path(input_dir, paste0(dataset,"_metadata.csv")), as.is = T)
  meta <- meta[,-1]
  rownames(meta) <- meta$cellid

  # Ensure that the cells in metadata and count matrix are the same and are
  # in the same order
  keep <- intersect(meta$cellid, colnames(fullmat))
  meta <- meta[keep,]
  fullmat <- fullmat[, keep]

  ### top level - broad cell types ###
  broadtypes <- unique(meta$broadcelltype)
  donors <- unique(meta$donor)

  pure_samples_broad <- list()
  for (ii in broadtypes) {
    pure_samples_broad[[ii]] <- which(meta$broadcelltype==ii)
  }

  numreps <- 10
  numcells <- 20000 # TODO -- or just use the number of cells available per cell type? Or number of cells in data set?

  ints <- seq(from = 0, to = 1, by = 0.1)

  pseudobulk <- matrix(0, nrow = nrow(fullmat),
                       ncol = length(donors) + length(donors)*length(broadtypes) + numreps * 2 * length(ints) * length(pure_samples_broad))
  rownames(pseudobulk) <- rownames(fullmat)
  colnames(pseudobulk) <- 1:ncol(pseudobulk)

  propval <- matrix(0, nrow = ncol(pseudobulk), ncol = length(pure_samples_broad))
  colnames(propval) <- names(pure_samples_broad)

  meta$donor = factor(meta$donor)

  # This is SIGNIFICANTLY faster than calling rowSums on individual donor sets

  y = model.matrix(~0 + donor, data = meta)
  inds = 1:(ncol(y) - 1)
  pseudobulk[, inds] <- as.array(fullmat %*% y)
  colnames(pseudobulk)[inds] = colnames(y)

  tab1 <- table(meta$donor, meta$broadcelltype)
  propval[inds, colnames(tab1)] <- tab1 / rowSums(tab1)

  index <- max(inds) + 1


  # Pure samples by donor
  meta$celltypedonor <- factor(paste(meta$broadcelltype, meta$donor, sep = "_"))
  y = model.matrix(~0 + celltypedonor, data = meta)
  inds = index:(index + ncol(y) - 1)
  pseudobulk[, inds] <- as.array(fullmat %*% y)
  colnames(pseudobulk)[inds] = str_replace(colnames(y), "celltypedonor", "puresample_")

  tab1 <- table(meta$celltypedonor, meta$broadcelltype)
  propval[inds, colnames(tab1)] <- tab1 / rowSums(tab1)

  index <- max(inds) + 1

  # Now add randomly-sampled data sets to fill out the rest of pseudobulk

  for (ii in 1:length(pure_samples_broad)) {
    cellkeep <- pure_samples_broad[[ii]]
    othercells <- setdiff(1:nrow(meta), cellkeep)

    for (jj in ints) {
      groups <- matrix(0, nrow = nrow(meta), ncol = numreps * 2)

      colnames(groups) <- paste(names(pure_samples_broad)[ii], jj, 1:(numreps*2), sep = "_")
      rownames(groups) <- rownames(meta)

      inds = index:(index + ncol(groups) - 1)

      for (kk in 1:numreps) {
        set.seed(ii*10000 + jj*100 + kk)
        numtokeep <- round(jj*numcells)
        keepinds <- sample(cellkeep, numtokeep, replace = T)

        ###remaining cells###
        numtofill <- numcells - numtokeep
        keepinds2 <- sample(othercells, numtofill, replace = T)

        keepinds <- c(keepinds, keepinds2)
        groups[sort(unique(keepinds)), kk] = table(keepinds)

        tab1 <- table(meta$broadcelltype[keepinds])
        propval[index, names(tab1)] <- tab1 / sum(tab1)

        index = index + 1
      }

      for (kk in 1:numreps) {
        set.seed(ii*10000 + jj*100 + kk)
        numtokeep <- round(jj*numcells)
        keepinds <- sample(cellkeep, numtokeep, replace = T)

        ###remaining cells###
        numtofill <- numcells - numtokeep
        props <- sample(1:10, ncol(propval)-1, replace = TRUE)
        props <- (props / sum(props)) * (1-jj) # Make proportions out of 100%, accounting for the celltype already used
        names(props) <- names(pure_samples_broad)[setdiff(1:length(pure_samples_broad), ii)]

        for (ct in names(props)) {
          numtofill2 <- round(numcells * props[ct])

          # The last one gets slightly different math to account for rounding error
          if (ct == names(props)[length(props)]) {
            numtofill2 <- numcells - length(keepinds)
          }

          keepinds2 <- sample(pure_samples_broad[[ct]], numtofill2, replace = T)
          keepinds <- c(keepinds, keepinds2)
        }

        groups[sort(unique(keepinds)), kk + numreps] = table(keepinds)

        tab1 <- table(meta$broadcelltype[keepinds])
        propval[index, names(tab1)] <- tab1 / sum(tab1)

        index = index + 1
      }

      pseudobulk[, inds] <- as.array(fullmat %*% groups)
      colnames(pseudobulk)[inds] = colnames(groups)

      print(c(dataset, index-1))
    }
  }

  rownames(propval) <- colnames(pseudobulk)

  # Some donors didn't have specific cell types, so this filters out those cases
  # where pure sample by donor resulted in 0 counts
  pseudobulk = pseudobulk[, which(colSums(pseudobulk) > 0)]
  propval = propval[colnames(pseudobulk), ]

  save(pseudobulk, propval, file = file.path(output_dir, paste0("pseudobulk_", dataset, "_broadcelltypes.rda")))

} # Artifically closing for loop to avoid running code below



  ###finer subtypes##
  finetypes <- unique(meta$subcluster)
  pure_samples_fine <- list()
  for (ii in finetypes) {
    pure_samples_fine[[ii]] <- which(meta$subcluster == ii)
  }

  numreps <- 10
  numcells <- 1000 # TODO adjust this number since the max fine subtype might be smaller than 20000

  ints <- seq(from = 0, to = 0.3,by = 0.02)

  pseudobulk <- matrix(0, nrow = nrow(fullmat),
                       ncol = numreps * length(ints) * length(pure_samples_fine))
  index <- 1
  propval <- matrix(0, nrow = ncol(pseudobulk), ncol = length(pure_samples_fine))
  colnames(propval) <- names(pure_samples_fine)

  for (ii in 1:length(pure_samples_fine)) {
    cellkeep <- pure_samples_fine[[ii]]
    othercells <- setdiff(1:nrow(meta), cellkeep)

    for (jj in ints) {
      for (kk in 1:numreps) {
        set.seed(ii*10000+jj*100+kk)
        numtokeep <- round(jj*numcells)
        keepinds=sample(cellkeep,numtokeep, replace = T)

        ###remaining cells###
        numtofill <- numcells - numtokeep
        keepinds2 <- sample(othercells, numtofill, replace = T)

        keepinds <- c(keepinds, keepinds2)
        pseudobulk[, index] <- rowSums(fullmat[, keepinds])
        colnames(pseudobulk)[index] <- paste(names(pure_samples_broad)[ii], jj, kk, sep = "_")

        tab1 <- table(meta$subcluster[keepinds])
        propval[index, names(tab1)] = tab1 / sum(tab1)
        index <- index + 1
        print(c(dataset, "fine", index - 1))
      }
    }
  }

  rownames(propval) <- colnames(pseudobulk)

  save(pseudobulk, propval, file = file.path(output_dir, paste0("pseudobulk_",dataset,"_finecelltypes_30percentlimit.rda")))
}
