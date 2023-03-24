# This function creates some pseudobulk samples by randomly sampling cells
# from the single cell dataset and adding the sampled cells' counts together.
# Cells are sampled with replacement, so cells can be included multiple times.
# To mimic how real bulk RNA sequencing works, the counts are not normalized
# prior to being added.
#
# Unlike the other two CreatePseudobulk functions, this function does not
# write to a file because it is intended to be used to generate one piece
# of a larger data set. Instead, it returns the pieces it creates.
#
# We sample cells in such a way that we force each cell type to be a certain
# proportion of the sample:
#   The "main_celltype" will have a known proportion (argument "proportion"),
#   while all the other cell types will be assigned a random proportion to make
#   up the rest of the sample. Each individual cell in the dataset is then
#   assigned a probability of being drawn by sample(), based on its cell type's
#   assigned proportion, such that the proportions of the sample are
#   approximately what we wanted +/- some randomness. Doing it this way ensures
#   that rare cell types have an opportunity to be sampled at higher proportions
#   than they exist at in the sample when they are randomly assigned a
#   proportion. It also ensures that we are NOT creating training samples that
#   all have the same relative proportions between non-main cell types.
#
# This script sacrifices readability for speed, as matrix multiplication is
# >10x faster than using rowSums() at this scale. Below is how things work:
#
# 1. We create the model matrix "groups" by randomly assigning cells to each
#    sample as described above. For each sample, the number of times each cell
#    came up is put in the row/col for that cell/sample. For example:
#            Samp1 Samp2 Samp3 ...
#     Cell1  0     1     2
#     Cell2  4     1     2
#     Cell3  1     3     1
#
#    indicates that sampling 5 cells from the population of {Cell1, Cell2, Cell3}
#    returned Cell1 0 times, Cell2 4 times, and Cell3 1 time for Samp1, etc.
#
# 2. Multiplying this group matrix with the counts matrix adds the counts of the
#    the sampled cells to create each sample:
#        Cell1 Cell2 Cell3            Samp1 Samp2 Samp3          Samp1 Samp2 Samp3
# Gene1  0     1     1         Cell1  0     1     2       Gene1  5     4     3
# Gene2  1     10    5     x   Cell2  4     1     2     = Gene2  45    26    27
# Gene3  10    5     0         Cell3  1     3     1       Gene3  20    15    30
# Gene4  2     0     2                                    Gene4  2     8     6
#
# Arguments:
#   singlecell_counts - a gene x cell matrix of counts
#   cell_assigns - a factored vector the same length as the number of cells.
#                  Each entry should be the cell type assignment of the
#                  corresponding cell in singlecell_counts
#   main_celltype - the cell type that has a pre-assigned proportion
#   proportion - the assigned proportion for main_celltype
#   num_cells - the number of cells to sample
#   num_samples - the number of pseudobulk samples to generate using this
#                 main_celltype / proportion combination
#
# Returns: a list with entries for:
#           "counts" - gene x num_samples pseudobulk counts matrix
#           "propCells" - num_samples x n_celltypes matrix of cell proportions
#                         in each sample
#           "pctRNA" - num_samples x n_celltypes matrix of percent RNA in each
#                      sample

CreatePseudobulk_Training <- function(singlecell_counts, cell_assigns,
                                      main_celltype, proportion, num_cells,
                                      num_samples) {

  # Ensure we always get the same result from the same data
  celltype_ind <- which(levels(cell_assigns) == main_celltype)
  set.seed(celltype_ind*10000 + proportion*100)

  n_celltypes <- length(levels(cell_assigns))
  names(cell_assigns) <- colnames(singlecell_counts)

  # The groups matrix will be filled in as we sample
  groups <- matrix(0, nrow = ncol(singlecell_counts), ncol = num_samples)

  colnames(groups) <- paste(main_celltype, proportion, 1:num_samples, sep = "_")
  rownames(groups) <- colnames(singlecell_counts)

  # Probability of each cell being drawn
  probs <- list()

  # Main cell type is <proportion>% of population, the remainder of the cell
  # types are randomly assigned a proportion of the remaining population.
  for (samp in 1:num_samples) {
    # Randomly pick a percent for each non-main cell type, normalize it to be
    # out of (1-proportion) instead of 100.
    tmp <- sample.int(100, n_celltypes-1, replace = TRUE)
    tmp <- c(proportion, (tmp / sum(tmp)) * (1-proportion))
    names(tmp) <- c(main_celltype, setdiff(levels(cell_assigns), main_celltype))

    probs_samp <- rep(0, ncol(singlecell_counts))

    # Each cell type gets its own proportion. That proportion is divided evenly
    # across all cells of that cell type so that each of those cells has a
    # probability of being drawn equal to <proportion>/n_cells. Assigning
    # probabilities this way for each cell type in the data set will result in
    # cells being drawn in approximately the correctly-assigned proportions.
    for (N in names(tmp)) {
      cells <- (cell_assigns == N)
      probs_samp[cells] <- tmp[N] / sum(cells)
    }

    probs[[samp]] <- probs_samp
  }

  # Turn into matrix
  probs <- do.call(cbind, probs)

  propCells <- matrix(0, nrow = ncol(groups), ncol = n_celltypes)
  colnames(propCells) <- levels(cell_assigns)
  rownames(propCells) <- colnames(groups)

  pctRNA <- matrix(0, nrow = ncol(groups), ncol = n_celltypes)
  colnames(pctRNA) <- levels(cell_assigns)
  rownames(pctRNA) <- colnames(groups)

  for (samp in 1:ncol(probs)) {
    # Will sample the cells with the probabilities above. We will get *approximately*
    # the assigned proportions, but not exactly due to randomness.
    cells <- sample(colnames(singlecell_counts), size = num_cells,
                    replace = TRUE, prob = probs[,samp])
    totals <- table(cells)
    groups[names(totals), samp] <- totals

    # Actual cell type proportions
    tab1 <- table(cell_assigns[cells])
    propCells[samp, names(tab1)] <- tab1 / sum(tab1)

    # Percent RNA
    pct <- CalculatePercentRNA(singlecell_counts[,cells],
                               rep(samp, length(cells)),
                               cell_assigns[cells])
    pctRNA[samp, colnames(pct)] <- pct
  }

  # This matrix multiplication adds the counts of the randomly-selected cells
  # for each created sample
  counts <- singlecell_counts %*% groups

  return(list("counts" = counts, "propCells" = propCells, "pctRNA" = pctRNA))
}
