library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("CytoTRACE")

nrow(combined@reductions$umap@cell.embeddings[,1:2])
ncol(as.matrix(combined@assays$RNA@counts))

non_empty <- combined@assays$RNA@counts[which(rowSums(combined@assays$RNA@counts != 0)),]
nrow(non_empty)

combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/combined.Rds")
results <- CytoTRACE(as.matrix(combined@assays$RNA@counts))
plotCytoTRACE(results, phenotype = as.matrix(combined@assays$RNA@counts), emb = combined@reductions$umap@cell.embeddings[,1:2], outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")