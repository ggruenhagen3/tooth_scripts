# # *HUMAN MOLAR* #
# library("Seurat")
# library("Matrix")
# library("reticulate")
# library("stringr")
# library("cowplot")
# library("ggplot2")
# 
# rna_path <- "C:/Users/miles/Downloads/d_tooth/"
# 
# data = read.table(paste0(rna_path, "data/GSM4365601_counts_human_germinectomy/matrix.mtx"))
# germ = CreateSeuratObject(counts = data, project = "GERM")
# germ$sample = "germ"
# 
# data = read.table(paste0(rna_path, "data/GSM4365602_counts_human_molar_pulp/matrix.mtx"))
# plpy = CreateSeuratObject(counts = data, project = "PLPY")
# plpy$sample = "plpy"
# 
# data = read.table(paste0(rna_path, "data/GSM4365603_counts_human_molar_pulp_old/matrix.mtx"))
# plpo = CreateSeuratObject(counts = data, project = "PLPO")
# plpo$sample = "plpo"
# 
# data = read.table(paste0(rna_path, "data/GSM4365608_counts_human_molar_healthy_1/matrix.mtx"))
# h1 = CreateSeuratObject(counts = data, project = "H1")
# h1$sample = "h1"
# 
# data = read.table(paste0(rna_path, "data/GSM4365607_counts_human_molar_healthy_2/matrix.mtx"))
# h2 = CreateSeuratObject(counts = data, project = "H2")
# h2$sample = "h2"
# 
# data = read.table(paste0(rna_path, "data/GSM4365609_counts_human_germ_molar_apical_papilla_female_15yo/matrix.mtx"))
# y15 = CreateSeuratObject(counts = data, project = "Y15")
# y15$sample = "y15"
# 
# data = read.table(paste0(rna_path, "data/GSM4365610_counts_No5_24_yo_healthy_retained/matrix.mtx"))
# y24 = CreateSeuratObject(counts = data, project = "Y24")
# y24$sample = "y24"
# 
# hm = merge(x=y24, y = list(germ, plpy, plpo, h1, h2, y15), add.cell.ids=c("y24", "germ", "plpy", "plpo", "h1", "h2", "y15"))
# 
# their_clusters = read.table(paste0(rna_path, "data/annotation_human.txt"))
# hm$their_clusters = their_clusters[,2]
# 
# hm <- NormalizeData(hm, normalization.method = "LogNormalize", scale.factor = 100000)
# hm = subset(hm, subset = nFeature_RNA > 500 & nFeature_RNA < 7500)
# hm = FindVariableFeatures(object = hm, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
# hm <- ScaleData(object = hm, vars.to.regress = NULL)
# hm <- RunPCA(hm, npcs = 30, verbose = FALSE)
# # t-SNE and Clustering
# hm <- RunUMAP(hm, reduction = "pca", dims = 1:30)
# hm <- FindNeighbors(hm, reduction = "umap", dims = 1:2)
# hm <- FindClusters(hm, resolution = 0.12)
# hm$seurat_clusters = factor(as.numeric(as.vector(hm$seurat_clusters)) + 1)
# 
# 
# Idents(hm) = hm$seurat_clusters
# png()
# DimPlot(hm, reduction = 'umap', label=T)
# dev.off()
# system(paste("rclone moveto ", "~/scratch/d_tooth/Rplot001.png", " dropbox:BioSci-Streelman/George/tmp/cluster.png", sep=""))
# 
# Idents(hm) = hm$sample
# png()
# DimPlot(hm, reduction = 'umap', label=T)
# dev.off()
# system(paste("rclone moveto ", "~/scratch/d_tooth/Rplot001.png", " dropbox:BioSci-Streelman/George/tmp/cluster_by_sample.png", sep=""))
# 
# 
# Idents(hm) = hm$their_clusters
# png()
# DimPlot(hm, reduction = 'umap', label=T)
# dev.off()
# system(paste("rclone moveto ", "~/scratch/d_tooth/Rplot001.png", " dropbox:BioSci-Streelman/George/tmp/cluster_theirs.png", sep=""))
# 
# Idents(hm) = hm$seurat_clusters
# png()
# FeaturePlot(hm, features="CELSR1", label = T, pt.size = 1)
# dev.off()
# system(paste("rclone moveto ", "~/scratch/d_tooth/Rplot001.png", " dropbox:BioSci-Streelman/George/tmp/celsr1.png", sep=""))
# 
# hm = readRDS("~/scratch/d_tooth/data/hm.rds")
# Idents(hm) = hm$their_clusters
# hm.deg = FindAllMarkers(hm, only.pos = T)
# hm.deg = hm.deg[which(hm.deg$p_val_adj < 0.05),]


########################################
# Huge Matrix: Human v Mouse v Cichlid #
########################################
# Load Libraries and Scripts
library("Seurat")
library("Matrix")
library("stringr")
library("cowplot")
library("ggplot2")
library("limma")
source("~/scratch/brain/brain_scripts/all_f.R")

# Load Data
incsr = readRDS("~/scratch/d_tooth/data/igor_incsr.rds")
im    = readRDS("~/scratch/d_tooth/data/igor_incsr_molar.rds")
hm    = readRDS("~/scratch/d_tooth/data/hm.rds")
tj    = readRDS("~/scratch/d_tooth/data/tj.rds")
jaw   = readRDS("~/scratch/d_tooth/data/jpool.rds")

# Find DEGs and Put them in HGNC format for cross-species comparison
# Idents(incsr) = incsr$annot
# incsr.all = FindAllMarkers(incsr)
# saveRDS(incsr.all, "~/scratch/d_tooth/data/incsr_deg_all_unfiltered.rds")
# incsr.deg = incsr.all[which(incsr.all$avg_logFC > 0),]
# incsr.deg = incsr.deg[which(incsr.deg$p_val_adj < 0.05),]
# incsr.deg = convertMouseDataFrameToHgnc(incsr.deg, 7)

Idents(im) = im$annot
im.all = FindAllMarkers(im)
saveRDS(im.all, "~/scratch/d_tooth/data/im_deg_all_unfiltered.rds")
im.deg = im.all[which(im.all$avg_logFC > 0),]
im.deg = im.deg[which(im.deg$p_val_adj < 0.05),]

Idents(hm$their_clusters)
hm.all = FindAllMarkers(hm)
saveRDS(hm.all, "~/scratch/d_tooth/data/hm_deg_all_unfiltered.rds")
hm.deg = hm.all[which(hm.all$avg_logFC > 0),]
hm.deg = hm.deg[which(hm.deg$p_val_adj < 0.05),]
# hm_deg = readRDS("~/scratch/d_tooth/data/hm_deg.rds")

tj.all = FindAllMarkers(tj)
saveRDS(tj.all, "~/scratch/d_tooth/data/tj_deg_all_unfiltered.rds")
tj.deg = tj.all[which(tj.all$avg_logFC > 0),]
tj.deg = tj.deg[which(tj.deg$p_val_adj < 0.05),]
tj.deg = hgncMzebraInPlace(tj.deg, 7, rownames(tj))

jaw.all = FindAllMarkers(jaw)
saveRDS(jaw.all, "~/scratch/d_tooth/data/jaw_deg_all_unfiltered.rds")
jaw.deg = jaw.all[which(jaw.all$avg_logFC > 0),]
jaw.deg = jaw.deg[which(jaw.deg$p_val_adj < 0.05),]
jaw.deg = hgncMzebraInPlace(jaw.deg, 7, rownames(jaw))

dfs = list(incsr.deg, im.deg, hm_deg, tj.deg, jaw.deg)
samples = samples=c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth", "Cichlid Jaw")

heatmapComparisonMulti(dfs = dfs, samples=samples, filename="mvhvc", "~/scratch/d_tooth/results/mvhvc/")

## Random from PACE ##
tj = readRDS("~/scratch/d_tooth/data/tj.rds")
tj.deg = FindAllMarkers(tj, only.pos = T)
tj.deg = tj.deg[which(tj.deg$p_val_adj < 0.05),]
tj.deg = hgncMzebraInPlace(tj.deg, 7, rownames(tj))

  source("~/scratch/brain/brain_scripts/all_f.R")
  hgncMzebra <- function(genes, gene_names) {
    pat <- read.table("~/scratch/m_zebra_ref/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
    valid_genes <- validGenes(genes, gene_names)
    all_hgnc <- convertToHgnc(gene_names)
    ind <- match(genes,pat$V2)
    ind <- ind[! is.na(ind)]
    found_names <- as.vector(pat$V7[ind])
    found_names <- found_names[!is.na(found_names)]
    found_names_hgnc <- as.vector(pat$V8[ind])
    found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
    
    df1 <- setNames(as.data.frame(found_names_hgnc), c("hgnc"))
    found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")
    
    pseudo_hgnc <- toupper(genes)
    df1 <- setNames(as.data.frame(pseudo_hgnc), c("hgnc"))
    found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")
    
    found_mzebra <- found_mzebra[,2:1]
    found_names_hgnc <- found_names_hgnc[,2:1]
    good_df <- rbind(all_hgnc, setNames(found_names, names(all_hgnc)), setNames(found_mzebra, names(all_hgnc)), setNames(found_names_hgnc, names(all_hgnc)))
    good_df <- unique(good_df)
    return(good_df)
  }
