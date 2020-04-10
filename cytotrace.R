library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("CytoTRACE")

combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/combined.Rds")
epi      <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/epi.rds")
jpool    <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/jpool.rds")
tj       <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")

obj <- jpool
results <- CytoTRACE(as.matrix(obj@assays$RNA@counts))
plotCytoTRACE(results, emb = obj@reductions$umap@cell.embeddings[,1:2], outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")