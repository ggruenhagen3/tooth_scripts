library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("CytoTRACE")

combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/combined.Rds")
results <- CytoTRACE(as.matrix(combined@assays$RNA@counts))
plotCytoTRACE(results, phenotype = as.matrix(combined@assays$RNA@counts), emb = combined@reductions$umap[,1:2], outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")