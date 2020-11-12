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

dfs = list(incsr.deg, im.deg, hm.deg, tj.deg, jaw.deg)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth", "Cichlid Jaw")
source("~/scratch/brain/brain_scripts/all_f.R")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc", "~/scratch/d_tooth/mvhvc/")

# Single Comparisons
rs = heatmapComparison(incsr.deg, tj.deg, "Mouse Incisor", "Cichlid Tooth",         "mi_v_ct", "~/scratch/d_tooth/mvhvc/")
write.table(rs[[2]], file="~/scratch/d_tooth/mvhvc/mi_v_ct_genes.txt", sep="\t", quote = F)
rs = heatmapComparison(incsr.deg, hm.deg, "Mouse Incisor", "Human Molar",           "mi_v_hm", "~/scratch/d_tooth/mvhvc/")
write.table(rs[[2]], file="~/scratch/d_tooth/mvhvc/mi_v_hm_genes.txt", sep="\t", quote = F)
rs = heatmapComparison(incsr.deg, im.deg, "Mouse Incisor", "Mouse Incisor & Molar", "mi_v_im", "~/scratch/d_tooth/mvhvc/")
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
