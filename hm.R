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

rownames(test) = gsub("-", ".", rownames(test))
rownames(test)[which(startsWith(rownames(test), "human_germinectomy"))] = paste0("germ_", rownames(test)[which(startsWith(rownames(test), "human_germinectomy"))])
rownames(test)[which(startsWith(rownames(test), "No5"))] = paste0("y24_", rownames(test)[which(startsWith(rownames(test), "No5"))])
rownames(test)[which(startsWith(rownames(test), "human_molar_healthy_2_one_"))] = paste0("h2_", rownames(test)[which(startsWith(rownames(test), "human_molar_healthy_2_one_"))])
rownames(test)[which(startsWith(rownames(test), "human_molar_healthy_1_one_"))] = paste0("h1_", rownames(test)[which(startsWith(rownames(test), "human_molar_healthy_1_one_"))])
rownames(test)[which(startsWith(rownames(test), "human_molar_pulp_old_10_one_"))] = paste0("plpo_", rownames(test)[which(startsWith(rownames(test), "human_molar_pulp_old_10_one_"))])
rownames(test)[which(startsWith(rownames(test), "human_molar_pulp_4_one_"))] = paste0("plpy_", rownames(test)[which(startsWith(rownames(test), "human_molar_pulp_4_one_"))])
rownames(test)[which(startsWith(rownames(test), "human_germ_molar_apical_papilla_female_15yo_one_"))] = paste0("y15_", rownames(test)[which(startsWith(rownames(test), "human_germ_molar_apical_papilla_female_15yo_one_"))])

length(rownames(test) %in% colnames(hm))
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
incsr.all = readRDS("~/scratch/d_tooth/data/incsr_deg_all_unfiltered.rds")
# incsr.deg = incsr.all[which(incsr.all$avg_logFC > 0),]
incsr.deg = incsr.all
incsr.deg = incsr.deg[which(incsr.deg$p_val_adj < 0.05),]
incsr.deg$orig = 0
for (cluster in unique(incsr.deg$cluster)) {
  this_cluster_rows = nrow(incsr.deg[which(incsr.deg$cluster == cluster),])
  incsr.deg$orig[which(incsr.deg$cluster == cluster)] = this_cluster_rows
}
incsr.deg = convertMouseDataFrameToHgnc(incsr.deg, 7)
incsr.deg$new = 0
for (cluster in unique(incsr.deg$cluster)) {
  this_cluster_rows = nrow(incsr.deg[which(incsr.deg$cluster == cluster),])
  incsr.deg$new[which(incsr.deg$cluster == cluster)] = this_cluster_rows
}
incsr.deg$correction_factor = incsr.deg$orig/incsr.deg$new

# Idents(im) = im$annot
# im.all = FindAllMarkers(im)
# saveRDS(im.all, "~/scratch/d_tooth/data/im_deg_all_unfiltered.rds")
im.all = readRDS("~/scratch/d_tooth/data/im_deg_all_unfiltered.rds")
# im.deg = im.all[which(im.all$avg_logFC > 0),]
im.deg = im.all
im.deg = im.deg[which(im.deg$p_val_adj < 0.05),]
im.deg$correction_factor = 1

# Idents(hm) = hm$their_clusters
# hm.all = FindAllMarkers(hm)
# saveRDS(hm.all, "~/scratch/d_tooth/data/hm_deg_all_unfiltered.rds")
hm.all = readRDS("~/scratch/d_tooth/data/hm_deg_all_unfiltered.rds")
# hm.deg = hm.all[which(hm.all$avg_logFC > 0),]
hm.deg = hm.all
hm.deg = hm.deg[which(hm.deg$p_val_adj < 0.05),]
hm.deg$correction_factor = 1
# hm_deg = readRDS("~/scratch/d_tooth/data/hm_deg.rds")

# tj.all = FindAllMarkers(tj)
# saveRDS(tj.all, "~/scratch/d_tooth/data/tj_deg_all_unfiltered.rds")
tj.all = readRDS("~/scratch/d_tooth/data/tj_deg_all_unfiltered.rds")
# tj.deg = tj.all[which(tj.all$avg_logFC > 0),]
tj.deg = tj.all
tj.deg = tj.deg[which(tj.deg$p_val_adj < 0.05),]
tj.deg$orig = 0
for (cluster in unique(tj.deg$cluster)) {
  this_cluster_rows = nrow(tj.deg[which(tj.deg$cluster == cluster),])
  tj.deg$orig[which(tj.deg$cluster == cluster)] = this_cluster_rows
}
tj.deg = hgncMzebraInPlace(tj.deg, 7, rownames(tj), onPACE = T)
tj.deg$new = 0
for (cluster in unique(tj.deg$cluster)) {
  this_cluster_rows = nrow(tj.deg[which(tj.deg$cluster == cluster),])
  tj.deg$new[which(tj.deg$cluster == cluster)] = this_cluster_rows
}
tj.deg$correction_factor = tj.deg$orig/tj.deg$new

# jaw.all = FindAllMarkers(jaw)
# saveRDS(jaw.all, "~/scratch/d_tooth/data/jaw_deg_all_unfiltered.rds")
jaw.all = readRDS("~/scratch/d_tooth/data/jaw_deg_all_unfiltered.rds")
# jaw.deg = jaw.all[which(jaw.all$avg_logFC > 0),]
jaw.deg = jaw.all
jaw.deg = jaw.deg[which(jaw.deg$p_val_adj < 0.05),]
jaw.deg$orig = 0
for (cluster in unique(jaw.deg$cluster)) {
  this_cluster_rows = nrow(jaw.deg[which(jaw.deg$cluster == cluster),])
  jaw.deg$orig[which(jaw.deg$cluster == cluster)] = this_cluster_rows
}
jaw.deg = hgncMzebraInPlace(jaw.deg, 7, rownames(jaw), onPACE = T)
jaw.deg$new = 0
for (cluster in unique(jaw.deg$cluster)) {
  this_cluster_rows = nrow(jaw.deg[which(jaw.deg$cluster == cluster),])
  jaw.deg$new[which(jaw.deg$cluster == cluster)] = this_cluster_rows
}
jaw.deg$correction_factor = jaw.deg$orig/jaw.deg$new

incsr.deg = read.table("~/scratch/d_tooth/results/igor/incsr_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)
im.deg = read.table("~/scratch/d_tooth/results/igor/im_sig_degs.tsv", sep="\t", header = T, stringsAsFactors = F)
hm.deg = read.table("~/scratch/d_tooth/results/igor/hm_sig_degs.tsv", sep="\t", header = T, stringsAsFactors = F)
tj.deg = read.table("~/scratch/d_tooth/results/igor/tj_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)
jaw.deg = read.table("~/scratch/d_tooth/results/igor/jaw_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)

dfs = list(incsr.deg, im.deg, hm.deg, tj.deg, jaw.deg)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth", "Cichlid Jaw")

rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

system(paste0("rclone copy ~/scratch/brain/results/bb_j_sig_cor.csv dropbox:BioSci-Streelman/George/Brain/bb/results/coexp/"))

incsr.deg2 = na.omit(plyr::ldply(unique(incsr.deg$cluster), function(x) my_var = incsr.deg[which(incsr.deg$cluster == x)[1:100],]))
im.deg2    = na.omit(plyr::ldply(unique(im.deg$cluster),    function(x) my_var = im.deg[which(im.deg$cluster       == x)[1:100],]))
hm.deg2    = na.omit(plyr::ldply(unique(hm.deg$cluster),    function(x) my_var = hm.deg[which(hm.deg$cluster       == x)[1:100],]))
tj.deg2    = na.omit(plyr::ldply(unique(tj.deg$cluster),    function(x) my_var = tj.deg[which(tj.deg$cluster       == x)[1:100],]))
jaw.deg2   = na.omit(plyr::ldply(unique(jaw.deg$cluster),   function(x) my_var = jaw.deg[which(jaw.deg$cluster     == x)[1:100],]))
dfs2 = list(incsr.deg2, im.deg2, hm.deg2, tj.deg2, jaw.deg2)
source("~/scratch/brain/brain_scripts/all_f.R")
rs = heatmapComparisonMulti(dfs = dfs2, samples=samples,  filename="mvhvc100", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc100_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# my_var = incsr.deg[which(incsr.deg$cluster == x),]; my_var[1:50,]
# Single Comparisons
# TJ vs INCSR
rs = heatmapComparison(incsr.deg, tj.deg, "Mouse Incisor", "Cichlid Tooth",         "mi_v_ct", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/mi_v_ct_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_ct_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_ct_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ vs IM
rs = heatmapComparison(im.deg, tj.deg, "Mouse Incisor+Molar", "Cichlid Tooth",         "im_v_ct", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/im_v_ct_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_ct_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_ct_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ vs HM
rs = heatmapComparison(hm.deg, tj.deg, "Human Molar", "Cichlid Tooth",         "hm_v_ct", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/hm_v_ct_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_ct_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_ct_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Jaw vs INCSR
rs = heatmapComparison(incsr.deg, jaw.deg, "Mouse Incisor", "Cichlid Jaw",         "mi_v_cj", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/mi_v_cj_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_cj_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_cj_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_cj_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Jaw vs IM
rs = heatmapComparison(im.deg, jaw.deg, "Mouse Incisor+Molar", "Cichlid Jaw",         "im_v_cj", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/im_v_cj_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_cj_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_cj_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_cj_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Jaw vs HM
rs = heatmapComparison(hm.deg, jaw.deg, "Human Molar", "Cichlid Jaw",         "hm_v_cj", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/hm_v_cj_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_cj_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_cj_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_cj_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

#
rs = heatmapComparison(incsr.deg, hm.deg, "Mouse Incisor", "Human Molar",           "mi_v_hm", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth/mvhvc/mi_v_hm_genes.txt", sep="\t", quote = F)
rs = heatmapComparison(incsr.deg, im.deg, "Mouse Incisor", "Mouse Incisor & Molar", "mi_v_im", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth/mvhvc/mi_v_im_genes.txt", sep="\t", quote = F)

# rs = heatmapComparison(incsr.deg, tj.deg,  "Mouse Incisor", "Cichlid Tooth",       "mi_v_ct", "C:/Users/miles/Downloads/d_tooth/results/mvhvc/")
# rs = heatmapComparison(incsr.deg, jaw.deg, "Mouse Incisor", "Cichlid Jaw",         "mi_v_cj", "C:/Users/miles/Downloads/d_tooth/results/mvhvc/")
# write.table(rs[[2]], file="C:/Users/miles/Downloads/d_tooth/results/mvhvc/mi_v_cj_genes.txt", sep="\t", quote = F)
# write.table(rs[[2]], file="C:/Users/miles/Downloads/d_tooth/results/mvhvc/mi_v_ct_genes.txt", sep="\t", quote = F)
#######################
# Expression Heatmaps #
#######################
incsr_hgnc = convertToHgncObj(incsr, "mouse")
tj_hgnc = convertToHgncObj(tj, "mzebra")
jaw_hgnc = convertToHgncObj(jaw, "mzebra")

saveRDS(incsr_hgnc, "~/scratch/d_tooth/data/incsr_hgnc.rds")
saveRDS(tj_hgnc,    "~/scratch/d_tooth/data/tj_hgnc.rds")
saveRDS(jaw_hgnc,   "~/scratch/d_tooth/data/jaw_hgnc.rds")

incsr_hgnc = readRDS("~/scratch/d_tooth/data/incsr_hgnc.rds")
tj_hgnc    = readRDS("~/scratch/d_tooth/data/tj_hgnc.rds")
jaw_hgnc   = readRDS("~/scratch/d_tooth/data/jaw_hgnc.rds")

incsr_hgnc$project = "Mouse Incisor"
im$project = "Mouse Incisor+Molar"
hm$project = "Human Molar"
tj_hgnc$project = "Cichlid Tooth"
jaw_hgnc$project = "Cichlid Jaw"

incsr_hgnc$my_seurat_clusters = incsr_hgnc$seurat_clusters
im$my_seurat_clusters = im$seurat_clusters
hm$my_seurat_clusters = hm$seurat_clusters

incsr_hgnc$seurat_clusters = incsr_hgnc$annot
im$seurat_clusters = im$annot
hm$seurat_clusetrs = hm$their_clusters

source("~/scratch/brain/brain_scripts/all_f.R")
objs = list(incsr_hgnc, im, hm, tj_hgnc, jaw_hgnc)
expressionDend(objs = objs)

#############
# CytoTRACE #
#############
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("CytoTRACE")
library("ggplot2")
library("RColorBrewer")
library("tidyverse")

incsr_results = readRDS("C:/Users/miles/Downloads/d_tooth/data/incsr_cyto.rds")
im_results    = readRDS("C:/Users/miles/Downloads/d_tooth/data/im_cyto.rds")
hm_results    = readRDS("C:/Users/miles/Downloads/d_tooth/data/hm_cyto.rds")

incsr_data = readRDS("C:/Users/miles/Downloads/d_tooth/data/incisor_data.rds")
incsr_umap = incsr_data$emb # this is actually tsne
incsr_data = "hi"

im_data = readRDS("C:/Users/miles/Downloads/d_tooth/data/incisor_molar_data.rds")
im_umap = im_data$UMAP
im_data = "hi"

hm_data = readRDS("C:/Users/miles/Downloads/d_tooth/data/human_data.rds")
hm_umap = hm_data$umap
hm_data = "hi"

rownames(im_umap) = colnames(im)

hm_umap = readRDS("~/scratch/d_tooth/data/hm_umap.rds")
goodnames = substr(colnames(hm), 5, 1000L)
new_hm_umap = data.frame()
for (i in 1:nrow(hm_umap)) {
  cell = rownames(hm_umap)[i]
  new_cell = str_replace(cell, "-", ".")
  if (grepl("24_yo", new_cell)) {
    rownames(hm_umap)[i] = paste0("y24", "_", new_cell)
  } else if (grepl("15yo", new_cell)) {
    rownames(hm_umap)[i] = paste0("y15", "_", new_cell)
  } else if (grepl("germinectomy", new_cell)) {
    rownames(hm_umap)[i] = paste0("germ", "_", new_cell)
  } else if (grepl("human_molar_pulp_4", new_cell)) {
    rownames(hm_umap)[i] = paste0("plpy", "_", new_cell)
  } else if (grepl("human_molar_pulp_old", new_cell)) {
    rownames(hm_umap)[i] = paste0("plpo", "_", new_cell)
  } else if (grepl("human_molar_healthy_1", new_cell)) {
    rownames(hm_umap)[i] = paste0("h1", "_", new_cell)
  } else if (grepl("human_molar_healthy_2", new_cell)) {
    rownames(hm_umap)[i] = paste0("h2", "_", new_cell)
  }
  # if (new_cell %in% goodnames) {
  #   newRow = data.frame(UMAP_1 = hm_umap[i,1], UMAP_2 = hm_umap[i,2])
  #   new_hm_umap = rbind(new_hm_umap, newRow)
  # }
}
# test = hm_umap[match(rownames(hm_umap), substr(head(colnames(hm)), 5, 1000L)),]
hm_umap = hm_umap[which(rownames(hm_umap) %in% colnames(hm)),]

plotCytoTRACE(incsr_results, emb = incsr_umap, outputDir = "C:/Users/miles/Downloads/d_tooth/results/")
plotCytoTRACE(im_results, emb = im_umap, outputDir = "C:/Users/miles/Downloads/d_tooth/results/")
plotCytoTRACE(hm_results, emb = new_hm_umap, outputDir = "C:/Users/miles/Downloads/d_tooth/results/")
plotCytoTRACE(results, emb = new_hm_umap, outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")

incsr_epi = subset(incsr, cells = names(epi_data$clusters))
incsr_epi$their_clusters = epi_data$clusters
Idents(incsr_epi) = incsr_epi$their_clusters
incsr_epi@reductions$umap@cell.embeddings = as.matrix(epi_data_emb[,2:3])
colnames(incsr_epi@reductions$umap@cell.embeddings) = c("UMAP_1", "UMAP_2")
rownames(incsr_epi@reductions$umap@cell.embeddings) = colnames(incsr_epi)
plotCyto(incsr_epi)
cytoScoreByIdent(incsr_epi) + ggtitle("Igor Incisor Epithelium")
res = cytoBINdeg(incsr_epi)
res[[1]] + ggtitle("Number of Cells in CytoBINs per Cluster - Igor Incisor Epi")
res[[2]] + ggtitle("Percent of Cells in CytoBINs per Cluster - Igor Incisor Epi")
write.table(res[[3]], "C:/Users/miles/Downloads/d_tooth/results/igor_epi_cytoBIN_deg.tsv", sep="\t", row.names = F, quote = F)


#########
# Conos #
#########
library("SeuratWrappers")
library("Seurat")
library("conos")
quickPlot = function(p, width = 800, height = 800, res = 100, suffix = "test") {
  png(paste0("/mnt/c/Users/miles/Downloads/", suffix, ".png"), width = width, height = height, res = res)
  print(p)
  dev.off()
}
# incsr = readRDS("/mnt/c/Users/miles/Downloads/d_tooth/data/igor_incsr.rds")
# incsr = readRDS("/mnt/c/Users/miles/Downloads/incsr_hgnc.rds")
im = readRDS("/mnt/c/Users/miles/Downloads/d_tooth/data/igor_incsr_molar.rds")
hm = readRDS("/mnt/c/Users/miles/Downloads/d_tooth/data/hm.rds")
tj = readRDS("/mnt/c/Users/miles/Downloads/tj_hgnc.rds")
jaw = readRDS("/mnt/c/Users/miles/Downloads/jaw_hgnc.rds")

# incsr$sample = "INCSR"
im$sample = "Mouse Incisor + Molar"
hm$sample = "Human Molar"
tj$sample = "Cichlid Tooth"
jaw$sample = "Cichlid Jaw"
tj$annot = tj$seurat_clusters
jaw$annot = jaw$seurat_clusters
# incsr$sample_annot = paste(incsr$sample, incsr$annot)
im$sample_annot = paste(im$sample, im$annot)
hm$sample_annot = paste(hm$sample, hm$annot)

all_list = list()
# all_list[["INCSR"]] =  NormalizeData(incsr) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
all_list[["IM"]] =  NormalizeData(im) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
all_list[["HM"]] =  NormalizeData(hm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
all_list[["TJ"]] =  NormalizeData(tj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
all_list[["JAW"]] =  NormalizeData(jaw) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
all <- Conos$new(all_list)
all$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 25, n.odgenes = 2000, matching.method = "mNN", metric = "angular", score.component.variance = TRUE, verbose = TRUE)
all$findCommunities()
all$embedGraph()
all_s = as.Seurat(all)

p = DimPlot(all_s, reduction = "largeVis", split.by = "sample")
quickPlot(p, width = 1600)
Idents(all_s) = all_s$sample_annot
p = DimPlot(all_s, reduction = "largeVis", split.by = "sample", label = T)
quickPlot(p, width = 1600, suffix = "test_sample_annot")

Idents(all_s) = all_s$annot
p = DimPlot(all_s, reduction = "largeVis", split.by = "sample", label = T)
quickPlot(p, width = 1600, suffix = "test_annot")


# Incsr vs Molar
im$anatomy.annot = paste(im$anatomy, im$annot)
Idents(im) = im$anatomy.annot
all_deg = data.frame()
for (this_annot in unique(im$annot)) {
  this_deg = FindMarkers(im, ident.1 = paste("Incisor", this_annot), ident.2 = paste("Molar", this_annot))
  this_deg$gene = rownames(this_deg)
  this_deg$annot = this_annot
  all_deg = rbind(all_deg, this_deg)
}
all_deg_sig = all_deg[which(all_deg$p_val_adj < 0.05),]

all_deg_sig$rank = 0
for (this_annot in unique(im$annot)) {
  print(this_annot)
  all_deg_sig$rank[which(all_deg_sig$annot == this_annot)] = 1:length(which(all_deg_sig$annot == this_annot))
  # all_deg_sig$rank[which(all_deg_sig$annot == this_annot)] = 1
}
write.csv(all_deg_sig, "C:/Users/miles/Downloads/d_tooth/results/im_incsr_vs_molar_sig_cluster.csv")

# Incisor Epi Markers in Our Datasets
incsr_epi_markers = data.frame(Cell_Type = "SI_progenitors", Genes = toupper(c("Cdh6", "Lrp11", "Cpne5")), Color = "red")
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "SR", Genes = toupper(c("Vat1l", "Fam19a4", "Hey2")), Color = "red"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Preameloblasts", Genes = toupper(c("Col22a1", "Vwde", "Kif5c")), Color = "red"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Am_secr", Genes = toupper(c("Enam", "Amelx", "Ctnna2")), Color = "yellow"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Am_ryr2", Genes = toupper(c("Ryr2", "Mylk", "Sox5")), Color = "blue"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Am_mat", Genes = toupper(c("Klk4", "Gpr155", "Slc34a2")), Color = "blue"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Am_post_mat", Genes = toupper(c("Gm17660", "Slc5a8", "Ptpn22")), Color = "blue"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "SI", Genes = toupper(c("Psmb10", "C1qb", "Ibsp")), Color = "green"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "SI", Genes = toupper(c("Rab3il1", "Pmch", "Cyp2s1")), Color = "green"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Cuboidal", Genes = toupper(c("Thbd", "Gnrh1", "Jph4")), Color = "yellow"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "OEE", Genes = toupper(c("Slc04a1", "Th", "Amer1")), Color = "green"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Lgr5_SC", Genes = toupper(c("Lrig1", "Disc1", "Fez1")), Color = "orange"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "Sfrp5_prog", Genes = toupper(c("Sfrp5", "Grp", "Rhoc")), Color = "orange"))
incsr_epi_markers = rbind(incsr_epi_markers, 
                          data.frame(Cell_Type = "OEE_prog", Genes = toupper(c("Fos", "Egr1", "Vrtn")), Color = "orange"))

tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
jaw <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/pool_shiny/data/jpool.rds")
tj_cyto  <- readRDS("C:/Users/miles/Downloads/d_tooth/data/tooth_cyto.rds")
jaw_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/jaw_cyto.rds")
tj$cyto = tj_cyto$CytoTRACE
jaw$cyto = jaw_cyto$CytoTRACE

obj = jaw
temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
pal = colorRampPalette(temp)
cell_type_df = data.frame()
for (i in 1:nrow(incsr_epi_markers)) {
  hgnc = as.vector(incsr_epi_markers[i, "Genes"])
  cell_type = as.vector(incsr_epi_markers[i, "Cell_Type"])
  color = as.vector(incsr_epi_markers[i, "Color"])
  gene = gene_info$ens[match(hgnc, gene_info$human)]
  gene_pos = c()
  if (hgnc %in% rownames(obj))
    gene = hgnc
  if (tolower(hgnc) %in% rownames(obj))
    gene = tolower(hgnc)
  if (! is.na(gene) & gene %in% rownames(obj) ) 
    gene_pos = colnames(obj)[which(obj@assays$RNA@counts[gene,] > 0)]
  if ( length(gene_pos) > 0 ) {
    newRow = data.frame(CytoTRACE = as.numeric(as.vector(obj$cyto[gene_pos])), Genes = gene, Cell_Type = cell_type, Color = color)
    cell_type_df = rbind(cell_type_df, newRow) 
  }
}
rownames(cell_type_df) <- NULL

pdf("C:/Users/miles/Downloads/d_tooth/results/incsr_epi_markers_in_jaw.pdf", width = 5.50, height = 7.50)
print(ggplot(cell_type_df, aes(x=Genes, y=CytoTRACE)) + geom_boxplot(alpha=0.8, aes(fill=Cell_Type)) + ylim(0,1) + geom_jitter(position=position_dodge2(width=0.5), alpha=0.5, aes(color = CytoTRACE)) + coord_flip() + theme_classic() + scale_color_gradientn(colors = pal(50)) + guides(color = F) + ggtitle("Igor Incsr Epi Markers in Cichlid Jaw") + scale_fill_manual(values = as.vector(sapply(levels(cell_type_df$Cell_Type), function(x) cell_type_df$Color[which(cell_type_df$Cell_Type == x)[1]]))))
dev.off()

# Paint Genes in Igor Incsr Mes
igor_incsr_mes = readRDS("C:/Users/miles/Downloads/d_tooth/data/igor_incsr_mes.rds")
FeaturePlot(igor_incsr_mes, "Celsr1", order = T, label = T, pt.size = 1.5)
FeaturePlot(igor_incsr_mes, "Gli1", order = T, label = T, pt.size = 1.5)
FeaturePlot(igor_incsr_mes, "Thy1", order = T, label = T, pt.size = 1.5)
FeaturePlot(igor_incsr_mes, "Ng2", order = T, label = T, pt.size = 1.1)

incsr_results = readRDS("C:/Users/miles/Downloads/d_tooth/data/incsr_cyto.rds")

# I'm suppose to re-run it, not take the original
# igor_incsr_mes$cyto = incsr_results$CytoTRACE[match(colnames(igor_incsr_mes), names(incsr_results$CytoTRACE))]
# 
# igor_incsr_mes$bin <- igor_incsr_mes$cyto
# igor_incsr_mes$bin[which(igor_incsr_mes$cyto < quantile(igor_incsr_mes$cyto, 0.33))] <- "relative_low"
# igor_incsr_mes$bin[which(igor_incsr_mes$cyto > quantile(igor_incsr_mes$cyto, 0.33) & igor_incsr_mes$cyto < quantile(igor_incsr_mes$cyto, 0.66))] <- "relative_medium"
# igor_incsr_mes$bin[which(igor_incsr_mes$cyto > quantile(igor_incsr_mes$cyto, 0.66))] <- "relative_high"
# 
# Idents(igor_incsr_mes) <- igor_incsr_mes$bin
# igor_incsr_mes_cyto_deg <- FindAllMarkers(igor_incsr_mes)
# write.csv(igor_incsr_mes_cyto_deg, "C:/Users/miles/Downloads/d_tooth/results/igor_incsr_mes_cyto_deg_raw.csv")
# igor_incsr_mes_cyto_deg_sig = igor_incsr_mes_cyto_deg[which(igor_incsr_mes_cyto_deg$p_val_adj < 0.05),]
# write.csv(igor_incsr_mes_cyto_deg_sig, "C:/Users/miles/Downloads/d_tooth/results/igor_incsr_mes_cyto_deg_sig.csv")
# FeaturePlot(igor_incsr_mes, "cyto", order = T, label = T) + scale_color_gradientn(colors = pal(50))

igor_incsr_mes_cyto2 = CytoTRACE::CytoTRACE(mat = as.matrix(igor_incsr_mes@assays$RNA@counts))
igor_incsr_mes$cyto2 = igor_incsr_mes_cyto2$CytoTRACE
FeaturePlot(igor_incsr_mes, "cyto2", order = T, label = T) + scale_color_gradientn(colors = pal(50))

igor_incsr_mes$bin2 <- igor_incsr_mes$cyto2
igor_incsr_mes$bin2[which(igor_incsr_mes$cyto2 <= quantile(igor_incsr_mes$cyto2, 0.33))] <- "relative_low"
igor_incsr_mes$bin2[which(igor_incsr_mes$cyto2 > quantile(igor_incsr_mes$cyto2, 0.33) & igor_incsr_mes$cyto2 < quantile(igor_incsr_mes$cyto2, 0.66))] <- "relative_medium"
igor_incsr_mes$bin2[which(igor_incsr_mes$cyto2 >= quantile(igor_incsr_mes$cyto2, 0.66))] <- "relative_high"

Idents(igor_incsr_mes) <- igor_incsr_mes$bin2
igor_incsr_mes_cyto_deg <- FindAllMarkers(igor_incsr_mes)
write.csv(igor_incsr_mes_cyto_deg, "C:/Users/miles/Downloads/d_tooth/results/igor_incsr_mes_cyto_deg_raw.csv")
igor_incsr_mes_cyto_deg_sig = igor_incsr_mes_cyto_deg[which(igor_incsr_mes_cyto_deg$p_val_adj < 0.05),]
write.csv(igor_incsr_mes_cyto_deg_sig, "C:/Users/miles/Downloads/d_tooth/results/igor_incsr_mes_cyto_deg_sig.csv")

#===================================================================================================
# Quiescence 03/25/2021 ============================================================================
#===================================================================================================
# Using quiescent CytoTRACE range (0.6-0.7) from CytoTRACE paper to find quiescent cells
obj_list = list()
obj_list[[1]] = incsr
obj_list[[2]] = tj
obj_list[[3]] = jaw
sig_qui_deg_list = list()
for (i in 1:length(obj_list)) {
  obj = obj_list[[i]]
  obj$qui = "none"
  obj$qui[which(obj$cyto >= 0.6)] = "qui"
  obj$qui[which(obj$cyto > 0.7)] = "high"
  
  Idents(obj) = obj$qui
  qui_deg = FindMarkers(obj, ident.1 = "qui", ident.2 = "high")
  qui_deg$gene = rownames(qui_deg)
  sig_qui_deg_list[[i]] = qui_deg[which(qui_deg$p_val_adj < 0.05),] 
}

write.csv(sig_qui_deg_list[[1]], "~/scratch/d_tooth/results/incsr_qui_60_70_v_high.csv")
write.csv(sig_qui_deg_list[[2]], "~/scratch/d_tooth/results/tj_qui_60_70_v_high.csv")

# 03/30/2021
# Take Hoxb5+ and Mki67- population (See CytoTRACE paper Fig 3). Find the CytoTRACE Range.
qui_pop = colnames(hm)[which(hm@assays$RNA@counts["HOXB5",] > 0 & hm@assays$RNA@counts["MKI67",] == 0)]
hm$cyto = hm_results$CytoTRACE
hm$qui_pop = factor(colnames(hm) %in% qui_pop)
hm$qui_pop = plyr::revalue(hm$qui_pop, replace = c("TRUE" = "Qui", "FALSE" = "Other"))
Idents(hm) = hm$qui_pop
cytoScoreByIdent(hm)
range(hm$cyto[which(hm$qui_pop == "Qui")])
# df = data.frame(cell = colnames(hm), qui_pop = colnames(hm) %in% qui_pop, cyto = hm$cyto)

# Split Igor's Human Molar into Epithelium and Mesenchyme. Then do that, see how it matches in Igor's Mouse Incisor, then split our data into epi and mes, and finally identify a quiescent population.
hm_epi = hm$

#===================================================================================================
# Growing vs Adult =================================================================================
#===================================================================================================
hm$age = hm$sample
hm$age[which(hm$age %in% c("germ", "y15"))] = "grow"
hm$age[which(hm$age != "grow")] = "adult"

Idents(hm) = hm$age
age_deg = FindMarkers(hm, ident.1 = "grow", ident.2 = "adult")
age_deg = age_deg[which(age_deg$p_val_adj < 0.05),]
age_deg$gene = rownames(age_deg)
write.csv(age_deg, "C:/Users/miles/Downloads/d_tooth/results/hm_grow_v_adult_deg_sig.csv")

Idents(hm) = hm$annot
for (i in 1:100) {
  gene = age_deg$gene[i]
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/hm_grow_v_adult/", gene, ".png"), width = 1000, height = 500, res = 80)
  print(myFeaturePlot(hm, gene, my.split.by = "age", my.title = paste0(gene, " - #", i)))
  dev.off()
}

# As a control test other samples vs each other
hm$ctrl_age = hm$sample
hm$ctrl_age[which(hm$ctrl_age %in% c("h1", "plpy"))] = "ctrl_grow"
hm$ctrl_age[which(hm$ctrl_age %in% c("h2", "plpo"))] = "ctrl_adult"
Idents(hm) = hm$ctrl_age
ctrl_deg = FindMarkers(hm, ident.1 = "ctrl_grow", ident.2 = "ctrl_adult")
ctrl_deg = ctrl_deg[which(ctrl_deg$p_val_adj < 0.05),]

hm$ctrl_age2 = hm$sample
hm$ctrl_age2[which(hm$ctrl_age2 %in% c("germ"))] = "ctrl_grow"
hm$ctrl_age2[which(hm$ctrl_age2 %in% c("y15"))] = "ctrl_adult"
Idents(hm) = hm$ctrl_age2
ctrl_deg2 = FindMarkers(hm, ident.1 = "ctrl_grow", ident.2 = "ctrl_adult")
ctrl_deg2 = ctrl_deg2[which(ctrl_deg2$p_val_adj < 0.05),]

#===================================================================================================
# Neural Recruitment Genes in Paul's Mes and Cluster Proportion ====================================
#===================================================================================================
# BoxPlots of Genes
mes <- readRDS("C:/Users/miles/Downloads/rna/data/combined.rds")
stem = c("Celsr1", "Gli1", "Sox2", "Sox9")
neuro_genes = c("Tnfrsf1b", "Ntf5", "Ntf3", "Ntrk1", "Ngf", "Bdnf", "Gdnf", "Sema3f", "Cxcr4", "Cxcl12", "Epha4", "Ret", "Gnao1", "Gnai2", "Stx7", "Rtn1", "Strap", "Nefh")
for (gene in stem) {
  values = mes@assays$RNA@data[gene,]
  df = data.frame(values = values, gene = rep(gene, length(values)), cond = mes$cond, cluster = mes$seurat_clusters, sample = mes$cond)
  
  col_pal = c("#F2444A","#0077b6")
  p = ggplot(df, aes(x=1, y = values, fill=sample, color=sample)) + geom_split_violin(alpha=0.6) + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) + scale_color_manual(values=col_pal, name = "Condition") + scale_fill_manual(values=col_pal, name = "Condition") + ylab("Normalized Expression") + xlab("") + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(), axis.ticks=element_blank())
  p1 = ggplot(df, aes(x=cluster, y = values, fill=sample, color=sample)) + geom_split_violin(alpha=0.6) + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) + scale_color_manual(values=col_pal, name = "Condition") + scale_fill_manual(values=col_pal, name = "Condition") + ylab("Normalized Expression") + xlab("Cluster") + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5))
  p2 = myFeaturePlot(mes, gene, my.split.by = "cond", my.col.pal = pal, na.blank = T, my.title = gene)
  
  pdf(paste0("C:/Users/miles/Downloads/d_tooth/results/painting/stem/", gene, "_cond.pdf"), width = 6, height = 7.5)
  print(p)
  dev.off()
  
  pdf(paste0("C:/Users/miles/Downloads/d_tooth/results/painting/stem/", gene, "_cond_cluster.pdf"), width = 15, height = 5)
  print(p1)
  dev.off()
  
  pdf(paste0("C:/Users/miles/Downloads/d_tooth/results/painting/stem/", gene, ".pdf"), width = 10, height = 5)
  print(p2)
  dev.off()
}

#Volcano plot
res = FindMarkers(mes, ident.1 = "CTRL", ident.2 = "CLIPP", min.pct = 0.01, logfc.threshold = 0.1)
library("scales")
sig_deg = read.csv("C:/Users/miles/Downloads/d_tooth/results/paul_mes_ctrl_v_clipp_deg_sig.csv")
res$plot_p = -log(res$p_val_adj, base = 10)
res$label = rownames(res)
res$sig = as.character(res$label %in% sig_deg$gene)
res[stem, "sig"] = "Stem"
res$sig = plyr::revalue(res$sig, replace = c("TRUE" = "sig", "FALSE" = "non-sig", "Stem" = "Stem"))
res = res[order(res$sig, decreasing = F),]
ggplot(res, aes(x = avg_logFC, y = plot_p, color = sig)) + geom_point(alpha = 1, pt.size = 1.5) + geom_text_repel(data = res[stem,],aes(label = label), color = "black") + ylab("-Log Adjusted P") + scale_color_manual(values = rev(hue_pal()(3)))

# Find the Pct Diff and Average Log Change for Every Gene between CTRL and CLIPP
clipp_cells = colnames(mes)[which(mes$cond == "CLIPP")]
ctrl_cells = colnames(mes)[which(mes$cond == "CTRL")]
df_dif_logFC = pct_dif_avg_logFC(mes, cells.1 = ctrl_cells, cells.2 = clipp_cells)

df_dif_logFC$list = "non_sig"
df_dif_logFC$list[which(df_dif_logFC$genes %in% sig_deg$gene)] = "sig"
df_dif_logFC[neuro2, "list"] = "list_axon"
df_dif_logFC = df_dif_logFC[order(df_dif_logFC$list, decreasing = T),]
ggplot(df_dif_logFC, aes(x = pct_dif, y = avg_logFC, color=list)) + geom_point(alpha = 0.7) + geom_text_repel(data = df_dif_logFC[neuro2,], aes(label = genes), color = "black") + ggtitle("Axon") + labs(color = "Color")

growth = c("Tmod2", "Pacs1", "Rtn1", "Stx7", "Gnai2", "Gnao1", "Fabp7", "Cotl1", "Cap1", "Capzb", "Sept2", "Strap", "Clptm1", "Cyfip1", "Crmp1", "Farp2", "Cxcl12")
neuro2 = c("Bdnf", "Gdnf", "Ngf", "Sema3f", "Ntf3", "Ntf5", "Ntrk1", "Tnfrsf1b")
colnames(df_dif_logFC)[c(3,4)] = c("CTRL", "CLIPP")
plot_df = melt(df_dif_logFC[neuro2,c("genes", "CTRL", "CLIPP")], id = "genes")
colnames(plot_df)[2] = "Condition"
ggplot(plot_df, aes(x=genes, y = value, color = Condition, fill = Condition)) + geom_bar(stat = "identity", alpha = 0.8, position = "dodge") + ylab("% Cells")

# Find Differences in Cluster Proportion in Paul's Mes
df_prop = data.frame(cluster = levels(mes$seurat_clusters), clipp = 0, ctrl = 0, clipp_of_cluster = 0, ctrl_of_cluster = 0, cluster_of_clipp = 0, cluster_of_ctrl = 0)
for (cluster in levels(mes$seurat_clusters)) {
  this_cells = colnames(mes)[which(mes$seurat_clusters == cluster)]
  this_clipp = this_cells[which(this_cells %in% clipp_cells)]
  this_ctrl = this_cells[which(this_cells %in% ctrl_cells)]
  
  clipp_of_cluster = length(this_clipp) / length(this_cells)
  ctrl_of_cluster = length(this_ctrl) / length(this_cells)
  cluster_of_clipp = length(this_clipp) / length(clipp_cells)
  cluster_of_ctrl = length(this_ctrl) / length(ctrl_cells)
  
  df_prop[which(df_prop$cluster == cluster),] = c(cluster, length(this_clipp), length(this_ctrl), 
                                                  clipp_of_cluster, ctrl_of_cluster, 
                                                  cluster_of_clipp, cluster_of_ctrl)
}

plot_prop = melt(df_prop[,c("cluster", "clipp_of_cluster", "ctrl_of_cluster")], id = "cluster")
plot_prop$value = as.numeric(as.vector(plot_prop$value))
plot_prop$cluster = factor(plot_prop$cluster, levels = levels(mes$seurat_clusters))
ggplot(plot_prop, aes(x = cluster, y = value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") + ylab("% of Cluster")
plot_prop = melt(df_prop[,c("cluster", "cluster_of_clipp", "cluster_of_ctrl")], id = "cluster")
plot_prop$value = as.numeric(as.vector(plot_prop$value))
plot_prop$cluster = factor(plot_prop$cluster, levels = levels(mes$seurat_clusters))
ggplot(plot_prop, aes(x = cluster, y = value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") + ylab("% of Condition")

#==================================================================================
# Attempt to find Epi Stem Markers in Our Epi 03/31/2021 ==========================
#==================================================================================
epi_stem_Lgr5 = c("Pknox2", "Spock1", "Sfrp5", "Grp", "Arhgef33", "Ccdc80", "Pcp4", "Pla2g4a", "Cd27", "Vsig2", "Trpm5", "Lgr5", "Gm1110", "Clcnkb", "Fez1", "Icam2", "Prob1", "Disc1", "Zscan10", "Slc35f3")
epi_stem_Sox2 = c("Gjb3", "Sox2", "Ccdc112", "Trpm4", "Gstm2", "AA465934", "Mex3a", "Myh10", "Cks1b", "Moxd1", "Bex1", "Mdk", "Ccdc34", "Dlgap5", "Uhrf1", "Rfc4", "Rrm2", "Gins2", "Rad54b", "Dtl")
epi_stem = c(epi_stem_Lgr5, epi_stem_Sox2, "Lrig1")
incsr_epi = readRDS("C:/Users/miles/Downloads/d_tooth/data/igor_incsr_epi.rds")
incsr_epi$stem_score = colSums(incsr_epi@assays$RNA@data[epi_stem,])
FeaturePlot(incsr_epi, "stem_score", order = T, label = T, pt.size = 3) + scale_color_gradientn(colors = pal(50))
myFeaturePlot(incsr_epi, "stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epithelial Stem Markers in Igor Incisor Epithelium")
blind_epi_stem_cells = colnames(incsr_epi)[which(incsr_epi$stem_score >= quantile(incsr_epi$stem_score, 0.95))]
DimPlot(incsr_epi, cells.highlight = blind_epi_stem_cells, label = T, pt.size = 3, sizes.highlight = 3)

mz_epi_stem = convertHgncDataFrameToMzebra(data.frame(toupper(epi_stem)), 1, rownames(tj), na.rm = T, return_vect = T)
tj$epi_stem_score = colSums(tj@assays$RNA@data[mz_epi_stem,])
jaw$epi_stem_score = colSums(jaw@assays$RNA@data[mz_epi_stem,])
FeaturePlot(tj, "epi_stem_score", order = T, label = T, pt.size = 2) + scale_color_gradientn(colors = pal(50))
FeaturePlot(jaw, "epi_stem_score", order = T, label = T, pt.size = 2) + scale_color_gradientn(colors = pal(50))
myFeaturePlot(tj, "epi_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epithelial Stem Markers in Cichlid Jaw")
myFeaturePlot(jaw, "epi_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epithelial Stem Markers in Cichlid Tooth")
tj_blind_epi_stem_cells = colnames(tj)[which(tj$epi_stem_score >= quantile(tj$epi_stem_score, 0.95))]
jaw_blind_epi_stem_cells = colnames(jaw)[which(jaw$epi_stem_score >= quantile(jaw$epi_stem_score, 0.95))]
DimPlot(tj,  cells.highlight = tj_blind_epi_stem_cells,  label = T, pt.size = 2, sizes.highlight = 2)
DimPlot(jaw, cells.highlight = jaw_blind_epi_stem_cells, label = T, pt.size = 2, sizes.highlight = 2)

tj_epi = subset(tj, idents = c(2,3,4))
tj_epi$epi_stem_score = colSums(tj@assays$RNA@data[mz_epi_stem,])
FeaturePlot(tj_epi, "epi_stem_score", order = T, label = T, pt.size = 2) + scale_color_gradientn(colors = pal(50))
myFeaturePlot(tj_epi, "epi_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epithelial Stem Markers in Cichlid Tooth Epithelium")
tj_epi$epi_stem_score2 = colSums(tj@assays$RNA@counts[mz_epi_stem,])
myFeaturePlot(tj_epi, "epi_stem_score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Number of Epithelial Stem Markers Expressed in Cichlid Tooth Epithelium")

tj_jaw$epi_stem_score = colSums(tj_jaw@assays$RNA@data[mz_epi_stem,])
mat = tj_jaw@assays$RNA@counts[mz_epi_stem,]
mat[which(mat > 0)] = 1
tj_jaw$epi_stem_score2 = colSums(mat)
myFeaturePlot(tj_jaw, "epi_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epithelial Stem Markers in Cichlid Tooth and Jaw")
myFeaturePlot(tj_jaw, "epi_stem_score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Number of Epithelial Stem Markers Expressed in Cichlid Tooth and Jaw")
tj_jaw_epi_stem_cells = colnames(tj_jaw)[which(tj_jaw$epi_stem_score >= quantile(tj_jaw$epi_stem_score, 0.95))]
DimPlot(tj_jaw, cells.highlight = setNames(list(jaw_blind_epi_stem_cells), "95th"), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle("Cells in 95th quantile of expression of Epi Stem Markers")

# Epi + Mes Stem
mz_mes_stem = c("celsr1a", "gli1", "ENSMZEG00005018875")
mz_epi_mes_stem = c(mz_epi_stem, mz_mes_stem)
mat = tj_jaw@assays$RNA@counts[mz_epi_mes_stem,]
mat[which(mat > 0)] = 1
tj_jaw$epi_mes_stem_score = colSums(tj_jaw@assays$RNA@data[mz_epi_mes_stem,])
tj_jaw$epi_mes_stem_score2 = colSums(mat)
myFeaturePlot(tj_jaw, "epi_mes_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epi+Mes Stem Markers in Cichlid Tooth and Jaw")
myFeaturePlot(tj_jaw, "epi_mes_stem_score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Number of Epi+Mes Stem Markers Expressed in Cichlid Tooth and Jaw")
tj_jaw_epi_mes_stem_cells = colnames(tj_jaw)[which(tj_jaw$epi_mes_stem_score >= quantile(tj_jaw$epi_mes_stem_score, 0.95))]
DimPlot(tj_jaw, cells.highlight = setNames(list(tj_jaw_epi_mes_stem_cells), "95th"), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle("Cells in 95th quantile of expression of Epi+Mes Stem Markers")

# Find Mes Markers that are similar in expression to Gli1 and Thy1
igor_incsr_mes = readRDS("C:/Users/miles/Downloads/d_tooth/data/igor_incsr_mes.rds")
mes_stem = c("Celsr1", "Gli1", "Thy1")
igor_incsr_mes$mes_stem_score = colSums(igor_incsr_mes@assays$RNA@data[mes_stem,])
top_cells = colnames(igor_incsr_mes)[which(igor_incsr_mes$mes_stem_score > quantile(igor_incsr_mes$mes_stem_score, 0.95))]
myFeaturePlot(igor_incsr_mes, "mes_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Celsr1, Gli1 and Thy1 in Incisor Mes")
DimPlot(igor_incsr_mes, cells.highlight = setNames(list(top_cells), "95th"), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle("95th quantile of expression of Celsr1, Gli1 and Thy1")

igor_incsr_mes$mes_stem = "Other"
igor_incsr_mes$mes_stem[top_cells] = "Stem"
Idents(igor_incsr_mes) = igor_incsr_mes$mes_stem
mes_stem_deg = FindMarkers(igor_incsr_mes, ident.1 = "Stem", ident.2 = "Other")
mes_stem_deg_sig = mes_stem_deg[which(mes_stem_deg$p_val_adj < 0.05),]

# Figure of Coexpression
mes_stem_deg_sig = read.csv("C:/Users/miles/Downloads/d_tooth/results/igor_incsr_mes_gli1_thy1_cells_deg_sig.csv", stringsAsFactors = F)
mes_stem_deg_sig_pos = as.vector(mes_stem_deg_sig$gene[which(mes_stem_deg_sig$avg_logFC > 0)])
mat = as.matrix(igor_incsr_mes@assays$RNA@data[mes_stem_deg_sig_pos,])
mat = mat[, order(igor_incsr_mes$mes_stem_score, decreasing = T)]

igor_incsr_mes = ScaleData(igor_incsr_mes, features = mes_stem_deg_sig_pos)
mat2 = as.matrix(igor_incsr_mes@assays$RNA@scale.data[mes_stem_deg_sig_pos,])
mat2 = mat2[, order(igor_incsr_mes$mes_stem_score, decreasing = T)]

test_genes = as.vector(mes_stem_deg_sig$gene[which(mes_stem_deg_sig$pct.2 < 0.55 & mes_stem_deg_sig$avg_logFC > 0)])
pheatmap::pheatmap(mat[1:10,], cluster_rows = F, cluster_cols = F, show_colnames = F)
pheatmap::pheatmap(mat[mes_stem2,], cluster_rows = F, cluster_cols = F, show_colnames = F, main = "Markers Used")
pheatmap::pheatmap(mat[test_genes[1:10],], cluster_rows = F, cluster_cols = F, show_colnames = F)
pheatmap::pheatmap(mat[cor_genes[1:10],], cluster_rows = F, cluster_cols = F, show_colnames = F, main = "Markers Most Correlated with Gli1+Thy1")

pheatmap::pheatmap(mat2[1:10,], cluster_rows = F, cluster_cols = F, show_colnames = F)
pheatmap::pheatmap(mat2[mes_stem2,], cluster_rows = F, cluster_cols = F, show_colnames = F)
pheatmap::pheatmap(mat2[test_genes[1:10],], cluster_rows = F, cluster_cols = F, show_colnames = F)

test_cor = c()
for (gene in test_genes) {
  test_cor = c(test_cor, cor(igor_incsr_mes@assays$RNA@data[gene,], igor_incsr_mes$mes_stem_score))
}
names(test_cor) = test_genes
cor_genes = names(test_cor)[order(test_cor, decreasing = T)]

# common_genes = mes_stem_deg_sig$gene[which(mes_stem_deg_sig$gene %in% mes_stem_deg_sig_backup$gene)]
# test = data.frame(common = common_genes)
# mes_stem_deg_sig_backup$rank = 1:nrow(mes_stem_deg_sig_backup)
# mes_stem_deg_sig$rank = 1:nrow(mes_stem_deg_sig)
# test$no_celsr1 = mes_stem_deg_sig[common_genes, "rank"]
# test$w_celsr1 = mes_stem_deg_sig_backup[common_genes, "rank"]
# test$avg_rank = test$no_celsr1 + test$w_celsr1
# test$avg_rank = test$avg_rank / 2
# mes_stem2 = as.vector(test$common[which(test$no_celsr1 < 20 & test$w_celsr1 < 20)])

mes_stem2 = mes_stem_deg_sig[which(mes_stem_deg_sig$avg_logFC > 0)[1:10],"gene"]
mes_stem2 = c("Thy1", "Gli1", "Foxd1", "Bcl2", "Fgf10", "Rimbp2", "Gem", "Hs3st6", "Etv5", "Ptn")
mz_mes_stem2 = convertHgncDataFrameToMzebra(data.frame(toupper(mes_stem2)), 1, gene_names = rownames(tj), return_vect = T, na.rm = T)

for (gene in mes_stem2) {
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_mes_stem/", gene, ".png"), width = 600, height = 485)
  print(FeaturePlot(igor_incsr_mes, gene, order = T, label = T, pt.size = 3))
  dev.off()
}

mz_epi_mes_stem = c(mz_epi_stem, mz_mes_stem2)
mat = tj_jaw@assays$RNA@counts[mz_epi_mes_stem,]
mat[which(mat > 0)] = 1
tj_jaw$epi_mes_stem_score = colSums(tj_jaw@assays$RNA@data[mz_epi_mes_stem,])
tj_jaw$epi_mes_stem_score2 = colSums(mat)
myFeaturePlot(tj_jaw, "epi_mes_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epi+Mes Stem Markers in Cichlid Tooth and Jaw")
myFeaturePlot(tj_jaw, "epi_mes_stem_score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Number of Epi+Mes Stem Markers Expressed in Cichlid Tooth and Jaw")
tj_jaw_epi_mes_stem_cells = colnames(tj_jaw)[which(tj_jaw$epi_mes_stem_score >= quantile(tj_jaw$epi_mes_stem_score, 0.95))]
DimPlot(tj_jaw, cells.highlight = setNames(list(tj_jaw_epi_mes_stem_cells), "95th"), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle("Cells in 95th quantile of expression of Epi+Mes Stem Markers")


library("circlize")
test = read.csv("C:/Users/miles/Downloads/test.csv")
chord_df = test[, c("df1_cluster", "df2_cluster", "ovlp")]
colnames(chord_df) = c("from", "to", "value")
# chord_df = chord_df[which( (startsWith(chord_df$from, "(CT)") | startsWith(chord_df$from, "(HM)")) & (startsWith(chord_df$to, "(CT)") | startsWith(chord_df$to, "(HM)")) ),]
# chord_df = chord_df[which( (startsWith(chord_df$from, "(CT)")) & (startsWith(chord_df$to, "(HM)")) ),]
# chord_df = chord_df[which( startsWith(chord_df$from, "(CT)") & ! startsWith(chord_df$to, "(CT)") ),]

# TJ vs MIM
chord_df = test[, c("df1_cluster", "df2_cluster", "ovlp")]
colnames(chord_df) = c("from", "to", "value")
chord_df = chord_df[which( (startsWith(chord_df$from, "(CT)")) & (startsWith(chord_df$to, "(MIM)")) ),]
paste(chord_df$to, collapse = "'='gray', '")
grid.col = c("(CT) Mesenchymal" = "red", "(CT) Endothelial" = "blue", "(CT) Epithelial" = "green", "(CT) Glia" = "orange", "(CT) Immune" = "purple", '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells' = 'gray' )
lty_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MIM) Pulp cells", "(MIM) Endothelial cells", "(MIM) Epithelial cells", "(MIM) Pulp cells", "(MIM) Macrophages"), c(2,2,2, 2, 2))
lwd_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MIM) Pulp cells", "(MIM) Endothelial cells", "(MIM) Epithelial cells", "(MIM) Pulp cells", "(MIM) Macrophages"), c(2, 2, 2, 2, 2))
border_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MIM) Pulp cells", "(MIM) Endothelial cells", "(MIM) Epithelial cells", "(MIM) Pulp cells", "(MIM) Macrophages"), c(1, 1, 1, 1, 1))
chordDiagram(chord_df, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df, link.border = border_df)

# TJ vs HM
chord_df = test[, c("df1_cluster", "df2_cluster", "ovlp")]
colnames(chord_df) = c("from", "to", "value")
chord_df = chord_df[which( (startsWith(chord_df$from, "(CT)")) & (startsWith(chord_df$to, "(HM)")) ),]
paste(chord_df$to, collapse = "'='gray', '")
grid.col = c("(CT) Mesenchymal" = "red", "(CT) Endothelial" = "blue", "(CT) Epithelial" = "green", "(CT) Glia" = "orange", "(CT) Immune" = "purple", '(HM) Pulp cells'='gray', '(HM) Perivascular cells'='gray', '(HM) Endothelial cells'='gray', '(HM) Glial cells'='gray', '(HM) Odontoblasts'='gray', '(HM) Immune cells'='gray', '(HM) PDL'='gray', '(HM) Pulp cells'='gray', '(HM) Perivascular cells'='gray', '(HM) Endothelial cells'='gray', '(HM) Glial cells'='gray', '(HM) Odontoblasts'='gray', '(HM) Immune cells'='gray', '(HM) PDL'='gray', '(HM) Pulp cells'='gray', '(HM) Perivascular cells'='gray', '(HM) Endothelial cells'='gray', '(HM) Glial cells'='gray', '(HM) Odontoblasts'='gray', '(HM) Immune cells'='gray', '(HM) PDL'='gray', '(HM) Pulp cells'='gray', '(HM) Perivascular cells'='gray', '(HM) Endothelial cells'='gray', '(HM) Glial cells'='gray', '(HM) Odontoblasts'='gray', '(HM) Immune cells'='gray', '(HM) PDL'='gray', '(HM) Pulp cells'='gray', '(HM) Perivascular cells'='gray', '(HM) Endothelial cells'='gray', '(HM) Glial cells'='gray', '(HM) Odontoblasts'='gray', '(HM) Immune cells'='gray', '(HM) PDL' = 'gray' )
lty_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(HM) Pulp cells", "(HM) Endothelial cells", "(HM) Endothelial cells", "(HM) Pulp cells", "(HM) Immune cells"), c(2,2,2, 2, 2))
lwd_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(HM) Pulp cells", "(HM) Endothelial cells", "(HM) Endothelial cells", "(HM) Pulp cells", "(HM) Immune cells"), c(2, 2, 2, 2, 2))
border_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(HM) Pulp cells", "(HM) Endothelial cells", "(HM) Endothelial cells", "(HM) Pulp cells", "(HM) Immune cells"), c(1, 1, 1, 1, 1))
png("C:/Users/miles/Downloads/ct_v_hm_chord.png", width = 1500, height = 1500, res = 150)
chordDiagram(chord_df, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df, link.border = border_df)
dev.off()
circos.clear()

# TJ vs MIM
chord_df = test[, c("df1_cluster", "df2_cluster", "ovlp")]
colnames(chord_df) = c("from", "to", "value")
chord_df = chord_df[which( (startsWith(chord_df$from, "(CT)")) & (startsWith(chord_df$to, "(MIM)")) ),]
paste(chord_df$to, collapse = "'='gray', '")
grid.col = c("(CT) Mesenchymal" = "red", "(CT) Endothelial" = "blue", "(CT) Epithelial" = "green", "(CT) Glia" = "orange", "(CT) Immune" = "purple", '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells'='gray', '(MIM) Pulp cells'='gray', '(MIM) PDL'='gray', '(MIM) Macrophages'='gray', '(MIM) Pericytes'='gray', '(MIM) Endothelial cells'='gray', '(MIM) Glia'='gray', '(MIM) Epithelial cells' = 'gray' )
lty_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MIM) Pulp cells", "(MIM) Endothelial cells", "(MIM) Epithelial cells", "(MIM) Pulp cells", "(MIM) Macrophages"), c(2,2,2, 2, 2))
lwd_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MIM) Pulp cells", "(MIM) Endothelial cells", "(MIM) Epithelial cells", "(MIM) Pulp cells", "(MIM) Macrophages"), c(2, 2, 2, 2, 2))
border_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MIM) Pulp cells", "(MIM) Endothelial cells", "(MIM) Epithelial cells", "(MIM) Pulp cells", "(MIM) Macrophages"), c(1, 1, 1, 1, 1))
png("C:/Users/miles/Downloads/ct_v_mim_chord.png", width = 1500, height = 1500, res = 150)
chordDiagram(chord_df, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df, link.border = border_df)
dev.off()
circos.clear()

# TJ vs MI
chord_df = test[, c("df1_cluster", "df2_cluster", "ovlp")]
colnames(chord_df) = c("from", "to", "value")
chord_df = chord_df[which( (startsWith(chord_df$from, "(CT)")) & (startsWith(chord_df$to, "(MI)")) ),]
paste(chord_df$to, collapse = "'='gray', '")
grid.col = c("(CT) Mesenchymal" = "red", "(CT) Endothelial" = "blue", "(CT) Epithelial" = "green", "(CT) Glia" = "orange", "(CT) Immune" = "purple", '(MI) Endothelial'='gray', '(MI) Maturing pulp'='gray', '(MI) Glia'='gray', '(MI) Lymphocytes'='gray', '(MI) Distal pulp'='gray', '(MI) Perivascular'='gray', '(MI) Apical pulp'='gray', '(MI) Lyve1 Macrophages'='gray', '(MI) Macrophages'='gray', '(MI) Pre-odontoblasts'='gray', '(MI) SI + SR'='gray', '(MI) Dental follicle 1'='gray', '(MI) Dental follicle 2'='gray', '(MI) OEE'='gray', '(MI) Alveolar osteo.'='gray', '(MI) Ameloblasts'='gray', '(MI) Innate leukocytes'='gray', '(MI) Endothelial'='gray', '(MI) Maturing pulp'='gray', '(MI) Glia'='gray', '(MI) Lymphocytes'='gray', '(MI) Distal pulp'='gray', '(MI) Perivascular'='gray', '(MI) Apical pulp'='gray', '(MI) Lyve1 Macrophages'='gray', '(MI) Macrophages'='gray', '(MI) Pre-odontoblasts'='gray', '(MI) SI + SR'='gray', '(MI) Dental follicle 1'='gray', '(MI) Dental follicle 2'='gray', '(MI) OEE'='gray', '(MI) Alveolar osteo.'='gray', '(MI) Ameloblasts'='gray', '(MI) Innate leukocytes'='gray', '(MI) Endothelial'='gray', '(MI) Maturing pulp'='gray', '(MI) Glia'='gray', '(MI) Lymphocytes'='gray', '(MI) Distal pulp'='gray', '(MI) Perivascular'='gray', '(MI) Apical pulp'='gray', '(MI) Lyve1 Macrophages'='gray', '(MI) Macrophages'='gray', '(MI) Pre-odontoblasts'='gray', '(MI) SI + SR'='gray', '(MI) Dental follicle 1'='gray', '(MI) Dental follicle 2'='gray', '(MI) OEE'='gray', '(MI) Alveolar osteo.'='gray', '(MI) Ameloblasts'='gray', '(MI) Innate leukocytes'='gray', '(MI) Endothelial'='gray', '(MI) Maturing pulp'='gray', '(MI) Glia'='gray', '(MI) Lymphocytes'='gray', '(MI) Distal pulp'='gray', '(MI) Perivascular'='gray', '(MI) Apical pulp'='gray', '(MI) Lyve1 Macrophages'='gray', '(MI) Macrophages'='gray', '(MI) Pre-odontoblasts'='gray', '(MI) SI + SR'='gray', '(MI) Dental follicle 1'='gray', '(MI) Dental follicle 2'='gray', '(MI) OEE'='gray', '(MI) Alveolar osteo.'='gray', '(MI) Ameloblasts'='gray', '(MI) Innate leukocytes'='gray', '(MI) Endothelial'='gray', '(MI) Maturing pulp'='gray', '(MI) Glia'='gray', '(MI) Lymphocytes'='gray', '(MI) Distal pulp'='gray', '(MI) Perivascular'='gray', '(MI) Apical pulp'='gray', '(MI) Lyve1 Macrophages'='gray', '(MI) Macrophages'='gray', '(MI) Pre-odontoblasts'='gray', '(MI) SI + SR'='gray', '(MI) Dental follicle 1'='gray', '(MI) Dental follicle 2'='gray', '(MI) OEE'='gray', '(MI) Alveolar osteo.'='gray', '(MI) Ameloblasts'='gray', '(MI) Innate leukocytes'='gray' )
lty_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MI) Maturing pulp", "(MI) Endothelial", "(MI) OEE", "(MI) Maturing pulp", "(MI) Macrophages"), c(2,2,2, 2, 2))
lwd_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MI) Maturing pulp", "(MI) Endothelial", "(MI) OEE", "(MI) Maturing pulp", "(MI) Macrophages"), c(2, 2, 2, 2, 2))
border_df = data.frame(c("(CT) Mesenchymal", "(CT) Endothelial", "(CT) Epithelial", "(CT) Glia", "(CT) Immune"), c("(MI) Maturing pulp", "(MI) Endothelial", "(MI) OEE", "(MI) Maturing pulp", "(MI) Macrophages"), c(1, 1, 1, 1, 1))
png("C:/Users/miles/Downloads/ct_v_mi_chord.png", width = 2000, height = 2000, res = 150)
chordDiagram(chord_df, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df, link.border = border_df)
dev.off()
circos.clear()

# test_mat = acast(test, df1_cluster ~ df2_cluster, value.var = "pct_same_dir")
# test_mat2 = acast(chord_df, to ~ from)
# test_mat2 = scale(t(test_mat2)) ** 2
# chordDiagram(test_mat2)

png("C:/Users/miles/Downloads/chord_diagram.png", width = 1000, height = 1000)
chordDiagram(chord_df)
dev.off()
circos.clear()

scaleValues = function(values) {
  values_norm = (values - min(values)) / (max(values) - min(values))
  # col_fun <- colorRamp(viridis(100))
  col_fun = colorRamp(rev(brewer.pal(11, "RdYlBu")))
  cols <- col_fun(values_norm)
  cols = rgb(cols[,1], cols[,2], cols[,3], maxColorValue = 255)
  return(cols)
}

bb_df = data.frame(cluster = 0:14, col = convert15$col[match(0:14, convert15$old)], new = convert15$new.full[match(0:14, convert15$old)], num_cells = aggregate(nCount_RNA ~ seuratclusters15, bb@meta.data, length)[,2], ieg = aggregate(ieg_score ~ seuratclusters15, bb@meta.data, mean)[,2])
bb_df = bb_df[which(bb_df$num_cells > 100),]
bb_df$ieg_col = scaleValues(bb_df$ieg)
circos.par("gap.degree" = 0, cell.padding = c(0, 0, 0, 0))
circos.initialize(bb_df$cluster, xlim = c(0, 1), sector.width = bb_df$num_cells)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  # print(xlim[0])
  print(CELL_META$sector.index)
  circos.rect(0, 0, 1, 1, col = bb_df$col[which(bb_df$cluster == as.numeric(CELL_META$sector.index))])
}, bg.border = 1, track.height = 0.15)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  # print(xlim[0])
  print(CELL_META$sector.index)
  circos.rect(0, 0, 1, 1, col = bb_df$ieg_col[which(bb_df$cluster == as.numeric(CELL_META$sector.index))])
}, bg.border = 1, track.height = 0.15)
circos.clear()

circos.initialize(bb_df$cluster, xlim = c(0, 1), sector.width = bb_df$num_cells)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  # print(xlim[0])
  print(CELL_META$sector.index)
  circos.rect(0, 0, 1, 1, col = bb_df$col[which(bb_df$cluster == as.numeric(CELL_META$sector.index))])
}, bg.border = 1, track.height = 0.15)
circos.clear()

bb = ScaleData(bb, features = c("LOC106675461", "egr1"))
bb$aroma = bb@assays$RNA@scale.data["LOC106675461",]
bb_df = data.frame(sample = unique(bb$sample), col = c("#9d0208", "#d00000", "#dc2f02", "#e85d04", "#f4b906", "#03045e", "#023e8a", "#0077b6", "#0096c7", "#00b4d8"), cond = c(rep("BHVE", 5), rep("CTRL", 5)), num_cells = aggregate(nCount_RNA ~ sample, bb@meta.data, length)[,2], depth = aggregate(depth ~ sample, bb@meta.data, mean)[,2], aroma = aggregate(aroma ~ sample, bb@meta.data, mean)[,2] )
bb_df$prop13 = unlist(sapply(bb_df$sample, function(x) length(which(bb$seuratclusters53 == 13 & bb$sample == x)) ))
bb_df$prop13 = unlist(sapply(bb_df$sample, function(x) bb_df$prop13[which(bb_df$sample == x)] /length(which(bb$seuratclusters15 == 0 & bb$sample == x)) ))
# bb_df$prop13 = bb_df$prop13 / bb_df$num_cells
bb_df$sample = factor(bb_df$sample, levels = c("c4", "c5", "b5", "b4", "b3", "b2", "b1", "c1", "c2", "c3"))
bb_df$aroma_col = scaleValues(bb_df$aroma)
# bb_df$depth[which(bb_df$cond == "CTRL")] = F

pdf("C:/Users/miles/Downloads/test2.pdf", width =  10, height = 10)
circos.par("gap.degree" = 0, cell.padding = c(0, 0, 0, 0), "start.degree" = -18.25)
circos.initialize(bb_df$sample, xlim = c(0, 1))

# Depth
track1_breaks = pretty(bb_df$depth, n = 3)
circos.track(ylim = c(0, max(track1_breaks)), bg.col = NA, bg.border = NA, track.height = 0.15, panel.fun = function(x, y) {
  for (tb in track1_breaks) {
    circos.segments(0, tb, 1, tb, col = "gray90")
  }
  
  value = bb_df$depth[which(bb_df$sample == CELL_META$sector.index)]
  circos.barplot(value, CELL_META$xcenter, col = bb_df$col[which(bb_df$sample == CELL_META$sector.index)], bar_width = 0.2, border = NA)
})
circos.yaxis(at = track1_breaks, sector.index = "b1", track.index = 1, side = "right")

# Proportion
track2_breaks = pretty(bb_df$prop13, n = 1)
circos.track(ylim = c(0, max(track2_breaks)), bg.col = NA, bg.border = NA, track.height = 0.15, panel.fun = function(x, y) {
  for (tb in track2_breaks) {
    circos.segments(0, tb, 1, tb, col = "gray90")
  }
  
  value = bb_df$prop13[which(bb_df$sample == CELL_META$sector.index)]
  circos.barplot(value, CELL_META$xcenter, col = bb_df$col[which(bb_df$sample == CELL_META$sector.index)], bar_width = 0.2, border = NA)

})
circos.yaxis(at=track2_breaks, sector.index = "b1", track.index = 2, side = "right")

# Aromatase
# track3_breaks = pretty(bb_df$aroma, n = 3)
# circos.track(ylim = c(min(track3_breaks), max(track3_breaks)), bg.col = NA, bg.border = NA, track.height = 0.15, panel.fun = function(x, y) {
#   # Negative zone
#   circos.rect(0, min(track3_breaks), 1, 0, col = "gray70", border = NULL)
#   
#   for (tb in track3_breaks) {
#     circos.segments(0, tb, 1, tb, col = "gray90")
#   }
#   
#   value = bb_df$aroma[which(bb_df$sample == CELL_META$sector.index)]
#   circos.barplot(value, CELL_META$xcenter, col = bb_df$aroma_col[which(bb_df$sample == CELL_META$sector.index)], bar_width = 0.2, border = NA)
# })
# circos.yaxis(at=track3_breaks, sector.index = "b1", track.index = 3, side = "right")

# Sample
circos.track(ylim = c(0, 1), bg.col = NA, bg.border = NA, track.height = 0.15, panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  circos.rect(0, 0, 1, 1, col = bb_df$col[which(bb_df$sample == CELL_META$sector.index)], border = NA)
})
dev.off()
