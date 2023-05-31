## RunCCA ##
test <- RunCCA(tj, jaw, renormalize = TRUE, rescale = TRUE)
tj_jaw <- FindVariableFeatures(test)
tj_jaw <- RunPCA(tj_jaw, npcs = 30, verbose = FALSE)
tj_jaw <- RunUMAP(tj_jaw, reduction = "pca", dims = 1:12) 
tj_jaw <- FindNeighbors(tj_jaw, reduction = "umap", dims = 1:2)
tj_jaw <- FindClusters(tj_jaw, resolution = 0.25)
DimPlot(tj_jaw, reduction = "umap", split.by = "cond", label = TRUE)

test <- RunCCA(mz, combined, renormalize = TRUE, rescale = TRUE)
tj_jaw <- FindVariableFeatures(test)
tj_jaw <- RunPCA(tj_jaw, npcs = 50, verbose = FALSE)
tj_jaw <- FindNeighbors(tj_jaw, reduction = "pca", dims = 1:50)
tj_jaw <- FindClusters(tj_jaw, resolution = 1.86)
tj_jaw <- RunUMAP(tj_jaw, dims = 1:50)
DimPlot(tj_jaw, reduction = "umap", split.by = "cond", label = TRUE)

old_cluster <- colnames(tj_jaw)
for (i in 0:40) {
  old_cluster[which(old_cluster %in% substr(WhichCells(b1b2c1mz, idents = i), 6, 21))] <- i
  old_cluster[which(old_cluster %in% WhichCells(b1b2c1mz, idents = i))] <- i
}
old_cluster <- as.numeric(old_cluster)
old_cluster[is.na(old_cluster)] <- "Removed"
old_cluster <- factor(old_cluster, levels = c(0:40, "Removed"))
tj_jaw$old_cluster <- old_cluster
Idents(tj_jaw) <- tj_jaw$old_cluster
DimPlot(tj_jaw, reduction = "umap", split.by = "cond", label = TRUE, pt.size = 1.2) + ggtitle("Zack - Master Clusters")
removed_cells <-WhichCells(tj_jaw, idents = "Removed") 
Idents(tj_jaw) <- "seurat_clusters"
DimPlot(tj_jaw[,removed_cells], reduction = "umap", split.by = "cond", label = TRUE, pt.size = 1.2) + ggtitle("Removed Cells")
DimPlot(tj_jaw, reduction = "umap", split.by = "cond", label = TRUE, pt.size = 1.2) + ggtitle("New Clusters")
## ClustMap ##
# Pre-analysis
library(ClusterMap)
jaw_deg <- read.table("C:/Users/miles/Downloads/d_tooth/results/jpool_deg.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
tj_deg <- read.table("C:/Users/miles/Downloads/d_tooth/results/tj_deg.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
tj_num <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
jaw_num <- as.numeric(tail(levels(jaw@meta.data$seurat_clusters), n=1))
cell_num_list <- list()
tj_out <- data.frame()
n_keep <- 50
for (i in 0:tj_num) {
  this_rows <- tj_deg[which(tj_deg$cluster == i),c("cluster", "gene")]
  # tj_out <- rbind(tj_out, this_rows[1:n_keep,])
  tj_out <- rbind(tj_out, this_rows)
  cell_num_list[["tj"]] <- c(cell_num_list[["tj"]], length(WhichCells(tj, idents = i)))
}
jaw_out <- data.frame()
for (i in 0:jaw_num) {
  this_rows <- jaw_deg[which(jaw_deg$cluster == i),c("cluster", "gene")]
  # jaw_out <- rbind(jaw_out, this_rows[1:n_keep,])
  jaw_out <- rbind(jaw_out, this_rows)
  cell_num_list[["jaw"]] <- c(cell_num_list[["jaw"]], length(WhichCells(jaw, idents = i)))
}
write.csv(jaw_out, "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv/jaw.csv")
write.csv(tj_out,  "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv/tj.csv")

file_list <- list()
file_list[["tj"]]  <- "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv/tj.csv"
file_list[["jaw"]] <- "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv/jaw.csv"
single_obj_list <- list()
single_obj_list[["tj"]]  <- tj
single_obj_list[["jaw"]] <- jaw
cluster_map(file_list,
            output = "C:/Users/miles/Downloads/d_tooth/data/clustermap_output/",
            cell_num_list = cell_num_list,
            single_obj_list = single_obj_list)

# All ClusterMap
data_list <- list()
data_list[["b1"]] <- "C:/Users/miles/Downloads/brain/data/B1-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/"
data_list[["b2"]] <- "C:/Users/miles/Downloads/brain/data/B2-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/"
data_list[["c1"]] <- "C:/Users/miles/Downloads/brain/data/C1-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/"
data_list[["mz"]] <- "C:/Users/miles/Downloads/brain/data/MZ-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/"
new_obj_list <- list()
new_obj_list[["b1"]] <- make_single_obj(data_list[["b1"]], is.10X = TRUE, output = "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/b1")
new_obj_list[["b2"]] <- make_single_obj(data_list[["b2"]], is.10X = TRUE, output = "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/b2")
new_obj_list[["c1"]] <- make_single_obj(data_list[["c1"]], is.10X = TRUE, output = "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/c1")
new_obj_list[["mz"]] <- make_single_obj(data_list[["mz"]], is.10X = TRUE, output = "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/mz")
new_cell_list <- list()
for (str in c("b1", "b2", "c1", "mz")) {
  print(str)
  obj <- new_obj_list[[str]]
  num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
  new_v <- c()
  for (i in 0:num_clusters) {
    new_v <- c(new_v, length(WhichCells(obj, idents = i)))
  }
  new_cell_list[[str]] <- new_v
}
new_file_list <- list()
new_file_list[["b1"]] <- "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/b1.markers.csv"
new_file_list[["b2"]] <- "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/b2.markers.csv"
new_file_list[["c1"]] <- "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/c1.markers.csv"
new_file_list[["mz"]] <- "C:/Users/miles/Downloads/d_tooth/data/clustermap_csv_2/mz.markers.csv"
res <- cluster_map(new_file_list,
            output = "C:/Users/miles/Downloads/d_tooth/data/clustermap_output_2/",
            cell_num_list = new_cell_list,
            single_obj_list = new_obj_list)
cm_comb <- make_comb_obj(data_list, is.10X = TRUE)
# cm_comb$b1 <- colnames(cm_comb) %in% colnames(new_obj_list[["b1"]])
sample_id <- colnames(cm_comb)
sample_id[which(sample_id %in% paste0("b1-", colnames(new_obj_list[["b1"]])))] <- "b1"
sample_id[which(sample_id %in% paste0("b2-", colnames(new_obj_list[["b2"]])))] <- "b2"
sample_id[which(sample_id %in% paste0("c1-", colnames(new_obj_list[["c1"]])))] <- "c1"
sample_id[which(sample_id %in% paste0("mz-", colnames(new_obj_list[["mz"]])))] <- "mz"
cm_comb$sample <- sample_id
unique(cm_comb$sample)
DimPlot(cm_comb, label = TRUE, reduction = "umap", split.by = "sample")
