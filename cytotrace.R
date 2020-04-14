library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("CytoTRACE")

combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/combined.Rds")
epi      <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/epi.rds")
jpool    <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/jpool.rds")
tj       <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")

combined$dataset <- "mesenchyme"
epi$dataset      <- "epithelium"
jpool$dataset    <- "jaw"
tj$dataset       <- "tooth"
all_obj <- list(combined, epi, jpool, tj)

# Top 50 Genes (Least Differentiated)
# # all_top_genes <- lapply(0:length(all_obj), function(x) c())
# for (i in 1:length(all_obj)) {
#   obj <- all_obj[[i]]
#   obj_str <- unique(obj$dataset)
#   results <- CytoTRACE(as.matrix(obj@assays$RNA@counts))
#   genes <- rownames(obj@assays$RNA@counts)
#   
#   this_result <- data.frame()
#   for (gene in genes) {
#     expr <- FetchData(object = obj, vars = gene, slot = "counts")
#     gene_cells <- colnames(obj@assays$RNA@counts[, which(x = expr > 0)])
#     if (length(gene_cells) > 30) {
#       gene_p <- "NA"
#       try(gene_p <- t.test(results$CytoTRACE[gene_cells], results$CytoTRACE)$p.value, silent = TRUE)
#       this_result <- rbind(this_result, t(c( obj_str, gene, length(gene_cells), mean(results$CytoTRACE[gene_cells]), gene_p  )))
#     }
#   } # end gene for
#   colnames(this_result) <- c("dataset", "gene", "num_cells", "mean", "p")
#   this_result <- this_result[order(this_result$mean),]
#   this_result <- this_result[1:50,]
#   # write.table(this_result, paste("/nv/hp10/ggruenhagen3/scratch/d_tooth/", obj_str, "_top_50.tsv", sep=""), ssep = "\t", quote = FALSE)
# } # end object for

# Celsr1 by Cluster
all_celsr1 <- c()
all_results <- c()
all_results_cluster <- data.frame()
for (obj in all_obj) {
  obj_str <- unique(obj$dataset)
  results <- CytoTRACE(as.matrix(obj@assays$RNA@counts))
  all_results <- c(all_results, results$CytoTRACE)

  if (obj_str == "mesenchyme" || obj_str == "epithelium") {
    expr <- FetchData(object = obj, vars = "Celsr1", slot = "counts")
  } else {
    expr <- FetchData(object = obj, vars = "celsr1a", slot = "counts")
  }
  celsr1_cells <- colnames(obj@assays$RNA@counts[, which(x = expr > 0)])

  conds <- unique(obj$cond)
  num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
  for (cond in conds) {
    Idents(obj) <- "cond"
    this_cond_cells <- WhichCells(obj, idents = c(cond))
    Idents(obj) <- "seurat_clusters"
    for (cluster in 0:num_clusters) {
      this_cluster_cells <- WhichCells(obj, idents = c(cluster))
      this_cells <- this_cond_cells[which(this_cond_cells %in% this_cluster_cells)]
      celsr1_cluster <- celsr1_cells[which(celsr1_cells %in% this_cells)]
      celsr1_cluster_p <- "NA"
      try(celsr1_cluster_p <- t.test(results$CytoTRACE[celsr1_cluster], results$CytoTRACE[this_cells])$p.value, silent = TRUE)
      all_results_cluster <- rbind(all_results_cluster, t(c( obj_str, cond, cluster, mean(results$CytoTRACE[celsr1_cluster]), mean(results$CytoTRACE[this_cells]), celsr1_cluster_p )))
    }
  }
  # # plotCytoTRACE(results, emb = obj@reductions$umap@cell.embeddings[,1:2], outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")
}
colnames(all_results_cluster) <- c("dataset", "cond", "cluster", "mean_of_celsr1_in_cluster", "mean_of_cluster", "p")
all_results_cluster$p_sig_greater <- as.vector(all_results_cluster$mean_of_celsr1_in_cluster) > as.vector(all_results_cluster$mean_of_cluster) & as.vector(all_results_cluster$p) < 0.05
write.table(all_results_cluster, "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/celsr1_cluster.tsv", sep="\t", quote = FALSE, row.names = FALSE)