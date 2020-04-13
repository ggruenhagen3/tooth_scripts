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

all_obj <- c(combined, epi, jpool, tj)
all_celsr1 <- c()
all_results <- c()
all_results_cluster <- data.frame()
for (obj in all_obj) {
  obj_str <- obj$dataset
  results <- CytoTRACE(as.matrix(obj@assays$RNA@counts))
  all_results <- c(all_results, results$CytoTRACE)
  
  expr <- FetchData(object = jpool, vars = "celsr1a", slot = "counts")
  if (obj_str == "combined" || obj_str == "epi") {
    expr <- FetchData(object = jpool, vars = "Celsr1", slot = "counts")
  }
  celsr1_cells <- colnames(jpool@assays$RNA@counts[, which(x = expr > 0)])
  
  Idents(obj) <- "seurat_clusters"
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  for (cluster in 0:num_clusters) {
    this_cells <- WhichCells(obj, idents = c(cluster))
    celsr1_cluster <- celsr1_cells[which(duplicated(c(celsr1_cells), this_cells))]
    celsr1_cluster_p <- t.test(results$CytoTRACE[celsr1_cluster], results$CytoTRACE[this_cells])$p.value
    all_results_cluster <- rbind(all_results_cluster, t(c( obj_str, cluster, mean(results$CytoTRACE[celsr1_cluster]), mean(results$CytoTRACE[this_cells]), celsr1_cluster_p )))
  }
  # plotCytoTRACE(results, emb = obj@reductions$umap@cell.embeddings[,1:2], outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")
}
colnames(all_results_cluster) <- c("dataset", "cluster", "mean_of_celsr1_in_cluster", "mean_of_cluster", "p")
write.table(all_results_cluster, "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/celsr1_cluster.tsv", sep="\t")