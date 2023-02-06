library(reshape2)

CalculatePercentRNA_old <- function(sce, donors, celltypes) {
  pctRNA <- sapply(levels(donors), function(dn) {
    cells1 <- donors == dn
    cellsize <- sum(counts(sce)[, cells1])

    # Percent of RNA contributed to the sample's total RNA by each cell type
    sapply(levels(celltypes), function(ct) {
      cells2 <- cells1 & celltypes == ct
      count <- counts(sce)[,cells2]
      return(sum(count) / cellsize)
    })
  })
  pctRNA <- t(pctRNA)
}

# Using model.matrix and matrix multiplication to sum counts is way faster
# than iterating over cell type/donor combos and using sum().
# singlecell_counts = matrix of counts (rows = genes, cols = cells)
# samples = vector of sample IDs that is the same length as ncol(singlecell_counts).
#           Must either be a factor (contains multiple sample IDs) or a
#           vector with the same sample ID repeated (for all cells belonging to
#           a single sample).
# celltypes = vector of cell type assignments that is the same length as ncol(sce).
#             Must be a factor.
CalculatePercentRNA <- function(singlecell_counts, samples, celltypes) {
  # Sums all counts for each donor/celltype combination
  y <- model.matrix(~0 + samples:celltypes)

  count_sums <- singlecell_counts %*% y

  # Sum over all genes for each donor/celltype combination
  count_sums <- colSums(count_sums)

  col_info <- str_split(names(count_sums), ":", simplify = TRUE)
  sample_list <- unique(col_info[,1])

  pctRNA <- sapply(sample_list, function (samp) {
    # All entries for this donor = 1 of each cell type
    cols <- which(col_info[,1] == samp)
    pct <- count_sums[cols] / sum(count_sums[cols])

    # Remove extra labels added by model.matrix
    names(pct) <- str_replace(names(pct), "samples.*celltypes", "")
    return(pct)
  })

  colnames(pctRNA) <- str_replace(colnames(pctRNA), "samples", "")
  return(t(pctRNA))
}

CalculateA <- function(sce, samples, celltypes) {
  # Sums all counts for each donor/celltype combination
  y <- model.matrix(~0 + samples:celltypes)

  count_sums <- counts(sce) %*% y

  # Sum over all genes for each donor/celltype combination
  count_sums <- colSums(count_sums)

  A_d <- sapply(levels(samples), function(d) {
    sapply(levels(celltypes), function(ct) {
      name <- paste0("samples", d, ":celltypes", ct)
      cells <- celltypes == ct & samples == d

      # Average library size = sum over all genes for this donor/celltype
      # divided by number of cells for this donor/celltype
      return(as.numeric(count_sums[name] / sum(cells)))
    })
  })

  # Normalize each donor row to sum to 1
  A_d <- t(sweep(A_d, 2, colSums(A_d, na.rm = TRUE), "/"))

  # Take the means but exclude entries for donors who don't have a certain cell type
  A_d[A_d == 0] = NA
  A = colMeans(A_d, na.rm = TRUE)
  A = A / sum(A) # Enforce sum to 1
}


CalculateSignature <- function(sc.cpm, samples, celltypes) {
  # Sums all counts for each donor/celltype combination
  y <- model.matrix(~0 + samples:celltypes)

  count_sums <- sc.cpm %*% y

  n_cells <- table(samples, celltypes)
  n_cells <- melt(n_cells) %>% mutate(name = paste0("samples", samples, ":celltypes", celltypes))
  rownames(n_cells) <- n_cells$name
  n_cells <- n_cells[colnames(count_sums),]

  # Divide total counts for each donor/celltype by the number of cells in
  # that combo
  count_sums <- sweep(count_sums, 2, n_cells$value, "/")

  # Get the mean over all donors, for each gene and cell type
  sig <- sapply(levels(celltypes), function(ct) {
    cols <- which(n_cells$celltypes == ct)
    rowMeans(count_sums[,cols], na.rm = TRUE)
  })

  return(sig)
}


ConvertPropCellsToPctRNA <- function(propCells, A) {
  pct <- sweep(propCells, 2, A, "*")
  pct <- sweep(pct, 1, rowSums(pct), "/")
  return(pct)
}
