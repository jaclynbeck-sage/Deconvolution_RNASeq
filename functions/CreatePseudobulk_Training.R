CreatePseudobulk_Training <- function(singlecell_counts, metadata, main_celltype, proportion, numreps) {
  celltypes <- levels(metadata$broad.cell.type)
  celltype_ind <- which(celltypes == main_celltype)

  set.seed(celltype_ind*10000 + proportion*100)

  cellkeep <- which(metadata$broad.cell.type == main_celltype)
  othercells <- setdiff(1:nrow(metadata), cellkeep)

  # Create the model matrix by randomly assigning cells to groups
  groups <- matrix(0, nrow = nrow(metadata), ncol = numreps)

  colnames(groups) <- paste(main_celltype, proportion, 1:numreps, sep = "_")
  rownames(groups) <- rownames(metadata)

  # Dummy first column so we can cbind to it below
  probs <- Matrix(0, nrow = nrow(metadata), sparse = FALSE)

  # Cells that are not the desired cell type will have proportion of (1-proportion)%.
  # Dividing this evenly between those cells results in each cell getting a
  # (1-proportion)/length(othercells) chance of being drawn. Cells that are of the
  # desired celltype will have a proportion of (proportion)%. Dividing this evenly
  # between these cells results in each desired cell getting a proportion / length(cellkeep)
  # chance of being drawn.

  # The first 10 samples all have the same probabilities, so this gets copied
  # across first 10 rows of "probs". For some reason doing cbind is way faster
  # than pre-allocating the matrix and assigning at certain indices.
  #pp <- rep((1-proportion)/length(othercells), nrow(metadata))
  #pp[cellkeep] <- proportion / length(cellkeep)
  #for (kk in 1:numreps) {
  #  probs <- cbind(probs, pp)
  #}

  # 10 samples: main cell type is <proportion>% of population, the remainder
  # of the cell types are randomly assigned a proportion of the population and
  # sampled proportionally. This allows the chance for rarer cell types to make
  # up a larger amount of the sample than would happen if sampling from a pool
  # of all remaining cells.
  for (kk in 1:numreps) { #(numreps+1):(numreps*2)) {
    tmp <- sample.int(100, length(celltypes)-1, replace = TRUE)
    tmp <- c(proportion, (tmp / sum(tmp)) * (1-proportion))
    names(tmp) <- c(main_celltype, setdiff(celltypes, main_celltype))

    pp <- rep(0, nrow(metadata))

    # Same concept as above except each cell type gets its own proportion
    for (N in names(tmp)) {
      keep = metadata$broad.cell.type == N
      pp[keep] <- tmp[N] / sum(keep)
    }

    probs <- cbind(probs, pp)
  }

  # Get rid of dummy first column
  probs <- probs[,-1]

  propval <- matrix(0, nrow = ncol(groups), ncol = length(celltypes))
  colnames(propval) <- celltypes
  rownames(propval) <- colnames(groups)

  for (kk in 1:ncol(probs)) {
    # Will sample the cells with the probabilities above. We will get *approximately*
    # the assigned proportions, but not exactly due to randomness.
    keepinds <- sample(metadata$cell.id, size = numcells, replace = TRUE, prob = probs[,kk])
    totals <- table(keepinds)
    groups[names(totals), kk] <- totals

    # Actual cell type proportions
    tab1 <- table(metadata[keepinds, "broad.cell.type"])
    propval[kk, names(tab1)] <- tab1 / sum(tab1)
  }

  counts = singlecell_counts %*% groups

  return(list("counts" = counts, "propval" = propval))
}
