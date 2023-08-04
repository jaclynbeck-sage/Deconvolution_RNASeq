# mkrs is fine cell type markers
# normalized_counts are DESeq2::vst() of bulk data

mm <- lapply(names(mkrs), function(N) {
  data.frame(gene = mkrs[[N]], celltype = N)
})
mm_fine <- do.call(rbind, mm)
mm_fine$major_celltype <- str_replace("\\..*", "")

cc_fine = cor(normalized_counts[,mm_fine$gene])

a3 <- subset(mm_fine, celltype == "Astr.3")

edges <- cc_fine[a3$gene, a3$gene]

edges[edges < 0] <- 0
tt <- melt(edges)
tt <- subset(tt, Var1 != Var2)
tt <- subset(tt, value > 0)
tmp <- paste(tt$Var1, tt$Var2, sep = ",")
tmp <- paste(tmp, collapse = ",")
gg <- igraph::make_undirected_graph(str_split(tmp, ",")[[1]])
cliq <- igraph::largest_cliques(gg)

# There may be multi-way ties
cliq_means <- sapply(cliq, function(cl) {
  mean(cc_fine[names(cl), (names(cl))])
})

best <- which.max(cliq_means)
if (length(best) > 1) {
  best <- best[1]
}

new_genes <- names(cliq[[best]])
ord <- order(rowMeans(cc_fine[new_genes, new_genes]), decreasing = TRUE)

new_genes[ord]

corrplot::corrplot(cc[ast$gene, ast$gene], tl.cex = 0.5, addgrid.col = NA, order = "FPC")
corrplot::corrplot(cc[new_genes[ord], new_genes[ord]], tl.cex = 0.5, addgrid.col = NA)
