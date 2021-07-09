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

#tj.all = read.csv("~/scratch/d_tooth/results/igor/tj_annot_cluster_deg_sig.csv")
#colnames(tj.all)[3] = "avg_logFC"
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
# tj.deg = tj.deg[which(!is.na(tj.deg$hgnc)),]
# tj.deg$gene = tj.deg$hgnc
tj.deg$new = 0
for (cluster in unique(tj.deg$cluster)) {
  this_cluster_rows = nrow(tj.deg[which(tj.deg$cluster == cluster),])
  tj.deg$new[which(tj.deg$cluster == cluster)] = this_cluster_rows
}
tj.deg$correction_factor = tj.deg$orig/tj.deg$new

# jaw.all = read.csv("~/scratch/d_tooth/results/igor/jaw_annot_cluster_deg_sig.csv")
# colnames(jaw.all)[3] = "avg_logFC"
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
# jaw.deg = jaw.deg[which(!is.na(jaw.deg$hgnc)),]
# jaw.deg$gene = jaw.deg$hgnc
jaw.deg$new = 0
for (cluster in unique(jaw.deg$cluster)) {
  this_cluster_rows = nrow(jaw.deg[which(jaw.deg$cluster == cluster),])
  jaw.deg$new[which(jaw.deg$cluster == cluster)] = this_cluster_rows
}
jaw.deg$correction_factor = jaw.deg$orig/jaw.deg$new

# tj_jaw.deg = read.table("~/scratch/d_tooth/results/igor/tj_jaw_sig_description_hgnc.tsv", sep = "\t", header = T)
tj_jaw.deg = jaw.deg[which(tj_jaw.deg$p_val_adj < 0.05),]
tj_jaw.deg$orig = 0
for (cluster in unique(tj_jaw.deg$cluster)) {
  this_cluster_rows = nrow(tj_jaw.deg[which(tj_jaw.deg$cluster == cluster),])
  tj_jaw.deg$orig[which(tj_jaw.deg$cluster == cluster)] = this_cluster_rows
}
tj_jaw.deg = tj_jaw.deg[which(! is.na(tj_jaw.deg$hgnc)),]
tj_jaw.deg$new = 0
for (cluster in unique(tj_jaw.deg$cluster)) {
  this_cluster_rows = nrow(tj_jaw.deg[which(tj_jaw.deg$cluster == cluster),])
  tj_jaw.deg$new[which(tj_jaw.deg$cluster == cluster)] = this_cluster_rows
}
tj_jaw.deg$correction_factor = tj_jaw.deg$orig/tj_jaw.deg$new

incsr.deg = read.table("~/scratch/d_tooth/results/igor/incsr_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)
im.deg = read.table("~/scratch/d_tooth/results/igor/im_sig_degs.tsv", sep="\t", header = T, stringsAsFactors = F)
hm.deg = read.table("~/scratch/d_tooth/results/igor/hm_sig_degs.tsv", sep="\t", header = T, stringsAsFactors = F)
tj.deg = read.table("~/scratch/d_tooth/results/igor/tj_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)
jaw.deg = read.table("~/scratch/d_tooth/results/igor/jaw_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)
tj_jaw.deg = read.table("~/scratch/d_tooth/results/igor/tj_jaw_sig_degs_hgnc.tsv", sep="\t", header = T, stringsAsFactors = F)

# Big Heatmap - Tooth and Jaw Separate
dfs = list(incsr.deg, im.deg, hm.deg, tj.deg, jaw.deg)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth", "Cichlid Jaw")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Big Heatmap - Tooth and Jaw Separate - Remove Mouse Incsior
dfs = list(im.deg, hm.deg, tj.deg, jaw.deg)
samples = c("Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth", "Cichlid Jaw")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc5", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc5_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

heat_df = rs[[1]]
heat_df = heat_df[which( startsWith(as.vector(heat_df$df1_cluster), "Cichlid") & (startsWith(as.vector(heat_df$df2_cluster), "Human") | startsWith(as.vector(heat_df$df2_cluster), "Mouse")) ),]
heat_mat = reshape2::acast(heat_df, df2_cluster ~ df1_cluster, value.var = "pct_same_dir")
rownames(heat_mat) = str_replace(rownames(heat_mat), "Mouse Incisor\\+Molar", "(MIM)")
rownames(heat_mat) = str_replace(rownames(heat_mat), "Mouse Incisor", "(MI)")
rownames(heat_mat) = str_replace(rownames(heat_mat), "Human Molar", "(HM)")
pdf("~/scratch/d_tooth/results/igor/test_col.pdf", width = 10, height = 10)
pheatmap::pheatmap(heat_mat, scale = "column", cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 20, angle_col = "45", color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100))
dev.off()
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/test_col.pdf dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
pdf("~/scratch/d_tooth/results/igor/test2_col.pdf", width = 10, height = 10)
pheatmap::pheatmap(heat_mat, scale = "column", cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 20, angle_col = "45", color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100))
dev.off()
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/test2_col.pdf dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Big Heatmap - Tooth and Jaw Combined
dfs = list(incsr.deg, im.deg, hm.deg, tj_jaw.deg)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth+Jaw")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc2", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc2_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc2_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Big Heatmap Unique - Tooth and Jaw Combined
incsr_unique = read.table("~/scratch/d_tooth/results/igor/incsr_unique_100.tsv", sep = "\t", header = T, stringsAsFactors = F)
incsr.deg.unique = data.frame(gene = toupper(unlist(incsr_unique)))
# incsr.deg.unique$cluster = substr(rownames(incsr.deg.unique),1,nchar(rownames(incsr.deg.unique))-1)
incsr.deg.unique$cluster = gsub('[[:digit:]]+', '', rownames(incsr.deg.unique))
incsr.deg.unique$avg_logFC = 1
im_unique = read.table("~/scratch/d_tooth/results/igor/im_unique_100.tsv", sep = "\t", header = T, stringsAsFactors = F)
im.deg.unique = data.frame(gene = unlist(im_unique))
# im.deg.unique$cluster = substr(rownames(im.deg.unique),1,nchar(rownames(im.deg.unique))-1)
im.deg.unique$cluster = gsub('[[:digit:]]+', '', rownames(im.deg.unique))
im.deg.unique$avg_logFC = 1
hm_unique = read.table("~/scratch/d_tooth/results/igor/hm_unique_100.tsv", sep = "\t", header = T, stringsAsFactors = F)
hm.deg.unique = data.frame(gene = unlist(hm_unique))
# hm.deg.unique$cluster = substr(rownames(hm.deg.unique),1,nchar(rownames(hm.deg.unique))-1)
hm.deg.unique$cluster = gsub('[[:digit:]]+', '', rownames(hm.deg.unique))
hm.deg.unique$avg_logFC = 1
tj_jaw.deg.pos = tj_jaw.deg[which(tj_jaw.deg$avg_logFC > 0),]
dfs = list(incsr.deg.unique, im.deg.unique, hm.deg.unique, tj_jaw.deg.pos)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth+Jaw")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc2_unique", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc2_unique_ovlp.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc2_unique_pct.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Big Heatmap - TJ and Jaw Annot Separate and Unique
tj.annot.deg.pos = read.delim("~/scratch/d_tooth/results/igor/tj_annot_unique_long_hgnc_100.tsv", sep = "\t", header = T)
tj.annot.deg.pos$avg_logFC = 1
tj.annot.deg.pos$description = NULL
jaw.annot.deg.pos = read.delim("~/scratch/d_tooth/results/igor/jaw_annot_unique_long_hgnc_100.tsv", sep = "\t", header = T)
jaw.annot.deg.pos$avg_logFC = 1
jaw.annot.deg.pos$description = NULL
dfs = list(incsr.deg.unique, im.deg.unique, hm.deg.unique, jaw.annot.deg.pos, tj.annot.deg.pos)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Jaw", "Cichlid Tooth")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc3_unique", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc3_unique_ovlp.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc3_unique_pct.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# Big Heatmap - TJ and Jaw Annot Separate
tj.annot.deg  = read.csv("~/scratch/d_tooth/results/igor/tj_annot_cluster_deg_sig_clean_2.csv")
jaw.annot.deg = read.csv("~/scratch/d_tooth/results/igor/jaw_annot_cluster_deg_sig_clean_no_tb_no_pigmented_2.csv")
incsr.deg = incsr.deg[order(incsr.deg$cluster),]
im.deg = im.deg[order(im.deg$cluster),]
hm.deg = hm.deg[order(hm.deg$cluster),]
jaw.annot.deg = jaw.annot.deg[order(jaw.annot.deg$cluster),]
tj.annot.deg = tj.annot.deg[order(tj.annot.deg$cluster),]
dfs = list(incsr.deg, im.deg, hm.deg, jaw.annot.deg, tj.annot.deg)
samples = c("(MI)", "(MIM)", "(HM)", "(CJ)", "(CT)")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc_figure", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc_figure_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/annot/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc_figure_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/annot/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc_figure_pct_same_dir.pdf dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/annot/"))

# Find common DEGs across Cichlid Tooth, Mouse Incisor + Molar, Human Molar
tj.annot.deg.pos$gene[which(tj.annot.deg.pos$cluster == "Mesenchymal" & tj.annot.deg.pos$gene %in% im.deg.pos$gene[which(im.deg.pos$cluster %in% c("Pulp cell","PDL") )] & tj.annot.deg.pos$gene %in% hm.deg.pos$gene[which(hm.deg.pos$cluster %in% c("Pulp cell","PDL","Odontoblasts") )] )]
tj.annot.deg.pos$gene[which(tj.annot.deg.pos$cluster == "Endothelial" & tj.annot.deg.pos$gene %in% im.deg.pos$gene[which(im.deg.pos$cluster %in% c("Endothelial cells") )] & tj.annot.deg.pos$gene %in% hm.deg.pos$gene[which(hm.deg.pos$cluster %in% c("Endothelial cells") )] )]
tj.annot.deg.pos$gene[which(tj.annot.deg.pos$cluster == "Glia" & tj.annot.deg.pos$gene %in% im.deg.pos$gene[which(im.deg.pos$cluster %in% c("Glia") )] & tj.annot.deg.pos$gene %in% hm.deg.pos$gene[which(hm.deg.pos$cluster %in% c("Glial cells") )] )]

heat_df = rs[[1]]
heat_mat = reshape2::acast(heat_df, df2_cluster ~ df1_cluster, value.var = "pct_same_dir")
heat_mat = heat_mat[,which( startsWith(colnames(heat_mat), "(CT)") |  startsWith(colnames(heat_mat), "(CJ)") )]
heat_mat_mi = heat_mat[which( startsWith(rownames(heat_mat), "(MI)") ),]
heat_mat_mi = heat_mat_mi[order(rownames(heat_mat_mi), decreasing = T),]
p = pheatmap::pheatmap(heat_mat_mi, scale = "column", border_color = NA, cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 20, angle_col = "45", color = colorRampPalette(viridis(n = 11))(100))
p_df = melt(p$gtable$grobs[[1]]$children[[1]]$gp$fill)
p_df$Var1 = factor(p_df$Var1 , levels = unique(p_df$Var1)[order(unique(p_df$Var1), decreasing = T)])
p2 = ggplot(p_df, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pdf("~/scratch/d_tooth/results/igor/heat_figure_mi.pdf", width = 10, height = 10)
print(p2)
dev.off()
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/heat_figure_mi.pdf dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/annot/"))
heat_mat_hm = heat_mat[which( startsWith(rownames(heat_mat), "(HM)") ),]
heat_mat_hm = heat_mat_hm[order(rownames(heat_mat_hm), decreasing = T),]
p = pheatmap::pheatmap(heat_mat_hm, scale = "column", border_color = NA, cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 20, angle_col = "45", color = colorRampPalette(viridis(n = 11))(100))
p_df = melt(p$gtable$grobs[[1]]$children[[1]]$gp$fill)
p_df$Var1 = factor(p_df$Var1 , levels = unique(p_df$Var1)[order(unique(p_df$Var1), decreasing = T)])
p2 = ggplot(p_df, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pdf("~/scratch/d_tooth/results/igor/heat_figure_hm.pdf", width = 10, height = 5)
print(p2)
dev.off()
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/heat_figure_hm.pdf dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/annot/"))
heat_mat_mim = heat_mat[which( startsWith(rownames(heat_mat), "(MIM)") ),]
heat_mat_mim = heat_mat_mim[order(rownames(heat_mat_mim), decreasing = T),]
p = pheatmap::pheatmap(heat_mat_mim, scale = "column", border_color = NA, cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 20, angle_col = "45", color = colorRampPalette(viridis(n = 11))(100))
p_df = melt(p$gtable$grobs[[1]]$children[[1]]$gp$fill)
p_df$Var1 = factor(p_df$Var1 , levels = unique(p_df$Var1)[order(unique(p_df$Var1), decreasing = T)])
p2 = ggplot(p_df, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pdf("~/scratch/d_tooth/results/igor/heat_figure_mim.pdf", width = 10, height = 5)
print(p2)
dev.off()
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/heat_figure_mim.pdf dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/annot/"))

# Big Heatmap - TJ and Jaw Annot Separate
dfs = list(incsr.deg, im.deg, hm.deg, tj.deg, jaw.deg)
samples = c("Mouse Incisor", "Mouse Incisor+Molar", "Human Molar", "Cichlid Tooth", "Cichlid Jaw")
rs = heatmapComparisonMulti(dfs = dfs, samples=samples,  filename="mvhvc4", "~/scratch/d_tooth/results/igor/")
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc4_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mvhvc4_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))


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
rs = heatmapComparison(incsr.deg, tj.deg, "Mouse Incisor", "Cichlid Tooth",         "mi_v_ct2", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/mi_v_ct2_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_ct2_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_ct2_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

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

# TJ+Jaw vs HM
rs = heatmapComparison(hm.deg, tj_jaw.deg, "Human Molar", "Cichlid Tooth+Jaw", "hm_v_tj_jaw", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/hm_v_tj_jaw_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_tj_jaw_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_tj_jaw_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_tj_jaw_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ+Jaw vs IM
rs = heatmapComparison(im.deg, tj_jaw.deg, "Mouse Incisor+Molar", "Cichlid Tooth+Jaw", "im_v_tj_jaw", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/im_v_tj_jaw_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_tj_jaw_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_tj_jaw_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_tj_jaw_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ+Jaw vs INCSR
rs = heatmapComparison(incsr.deg, tj_jaw.deg, "Mouse Incisor", "Cichlid Tooth+Jaw", "mi_v_tj_jaw", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/mi_v_tj_jaw_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_tj_jaw_ovlp_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_tj_jaw_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/mi_v_tj_jaw_pct_same_dir.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ+Jaw vs HM - Unique
rs = heatmapComparison(hm.deg.unique, tj_jaw.deg.pos, "Human Molar", "Cichlid Tooth+Jaw", "hm_v_tj_jaw_unique", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/hm_v_tj_jaw_unique_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_tj_jaw_unique_ovlp.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_tj_jaw_unique_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/hm_v_tj_jaw_unique_pct.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ+Jaw vs IM - Unique
rs = heatmapComparison(im.deg.unique, tj_jaw.deg.pos, "Mouse Incisor+Molar", "Cichlid Tooth+Jaw", "im_v_tj_jaw_unique", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/im_v_tj_jaw_unique_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_tj_jaw_unique_ovlp.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_tj_jaw_unique_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/im_v_tj_jaw_unique_pct.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))

# TJ+Jaw vs Incisor - Unique
rs = heatmapComparison(incsr.deg.unique, tj_jaw.deg.pos, "Mouse Incisor", "Cichlid Tooth+Jaw", "incsr_v_tj_jaw_unique", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth//results/igor/incsr_v_tj_jaw_unique_genes.txt", sep="\t", quote = F)
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/incsr_v_tj_jaw_unique_ovlp.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/incsr_v_tj_jaw_unique_genes.txt dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))
system(paste0("rclone copy ~/scratch/d_tooth/results/igor/incsr_v_tj_jaw_unique_pct.png dropbox:BioSci-Streelman/George/Tooth/igor/results/heatmaps/"))


#
rs = heatmapComparison(incsr.deg, hm.deg, "Mouse Incisor", "Human Molar",           "mi_v_hm", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth/mvhvc/mi_v_hm_genes.txt", sep="\t", quote = F)
rs = heatmapComparison(incsr.deg, im.deg, "Mouse Incisor", "Mouse Incisor & Molar", "mi_v_im", "~/scratch/d_tooth/results/igor/")
write.table(rs[[2]], file="~/scratch/d_tooth/mvhvc/mi_v_im_genes.txt", sep="\t", quote = F)

# rs = heatmapComparison(incsr.deg, tj.deg,  "Mouse Incisor", "Cichlid Tooth",       "mi_v_ct", "C:/Users/miles/Downloads/d_tooth/results/mvhvc/")
# rs = heatmapComparison(incsr.deg, jaw.deg, "Mouse Incisor", "Cichlid Jaw",         "mi_v_cj", "C:/Users/miles/Downloads/d_tooth/results/mvhvc/")
# write.table(rs[[2]], file="C:/Users/miles/Downloads/d_tooth/results/mvhvc/mi_v_cj_genes.txt", sep="\t", quote = F)
# write.table(rs[[2]], file="C:/Users/miles/Downloads/d_tooth/results/mvhvc/mi_v_ct_genes.txt", sep="\t", quote = F)
#=======================================================================================================
# Expression Heatmaps ==================================================================================
#=======================================================================================================
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
neuro2 = c("Bdnf", "Gdnf", "Ngf", "Sema3f", "Ntf3", "Ntf5", "Ntrk1", "Tnfrsf1b", "Artn", "Il3")
colnames(df_dif_logFC)[c(3,4)] = c("CTRL", "CLIPP")
plot_df = melt(df_dif_logFC[neuro2,c("genes", "CTRL", "CLIPP")], id = "genes")
colnames(plot_df)[2] = "Condition"
ggplot(plot_df, aes(x=genes, y = value, color = Condition, fill = Condition)) + geom_bar(stat = "identity", alpha = 0.8, position = "dodge") + ylab("% Cells")

# Markers From Zack
raw_angio = "TIAM1/SEMA5B/EPHA10/CNTN4/CNTN6/ADCY1/OLFM1/NRXN1/NPTX1/PLXNA1/GAB1/RELN/ROBO3/CHN1/CCK/NOTCH2/LRP1/PRKCA/LAMA5/CDKL5/NCAM1/TBR1/PTPRS/LIMK1/SIPA1L1/EFNA1/SYNGAP1/RAP1GAP/GFRA2/MAP2K1/SLITRK5/LINGO1/BDNF/BRSK2/NPTN/PAK1/UST/SEMA3E/SEMA6A/RNF165/EPHB3/FZD3/TRIO/ADARB1/GFRA4/GPC1/LRP4/CDH11/MARK2/ARHGAP4/EFNB1/ABL1/DSCAML1/FGF13/KLF7/SLITRK3/DCLK1/NUMBL/SPTBN4/L1CAM/GAB2/SH3KBP1/VLDLR/SEMA6B/APBB2/TIAM2/FYN/SPTBN1/TAOK2/ULK1/BAIAP2/PITPNA/EXT1/EFNA2"
raw_axon = "SEMA5B/OLFM1/LRP1/CDKL5/PTPRS/LIMK1/BDNF/PAK1/SEMA3E/SEMA6A/ARHGAP4/ABL1/DCLK1/L1CAM/SEMA6B/ULK1"
angio = str_to_title( str_split(raw_angio, "/")[[1]] )
axon = str_to_title( str_split(raw_axon, "/")[[1]] )
angio_hgnc = toupper(angio)
axon_hgnc = toupper(axon)
axon_hgnc[which(axon_hgnc == "NTF5")] = "NTF4"

# Find Correlations
clpp_r = data.frame()
ctrl_r = data.frame()
grow_r = data.frame()
adult_r = data.frame()
for (i in 1:length(angio)) {
  angio_m_i = angio[i]
  angio_h_i = angio_hgnc[i]
  for (j in 1:length(axon)) {
    axon_m_j = axon[j]
    axon_h_j = axon_hgnc[j]
    clpp_r = rbind(clpp_r, t(c(angio_m_i, axon_m_j, paul_mes_clipp_r[angio_m_i, axon_m_j])))
    ctrl_r = rbind(ctrl_r, t(c(angio_m_i, axon_m_j, paul_mes_ctrl_r[angio_m_i, axon_m_j])))
    grow_r = rbind(grow_r, t(c(angio_h_i, axon_h_j, hm_grow_r[angio_h_i, axon_h_j])))
    adult_r = rbind(adult_r, t(c(angio_h_i, axon_h_j, hm_adult_r[angio_h_i, axon_h_j])))
  }
}
colnames(clpp_r) = colnames(ctrl_r) = colnames(grow_r) = colnames(adult_r) = c("angio", "axon", "cor")
clpp_r$cor     = as.numeric(as.vector(clpp_r$cor))
ctrl_r$cor     = as.numeric(as.vector(ctrl_r$cor))
grow_r$cor     = as.numeric(as.vector(grow_r$cor))
adult_r$cor     = as.numeric(as.vector(adult_r$cor))

clpp_r$cor[which(clpp_r$angio == clpp_r$axon)] = NA
ctrl_r$cor[which(ctrl_r$angio == ctrl_r$axon)] = NA
grow_r$cor[which(grow_r$angio == grow_r$axon)] = NA
adult_r$cor[which(adult_r$angio == adult_r$axon)] = NA

png("~/scratch/d_tooth/results/angio_axon_cor_clpp.png", width = 1000, height = 2000, res = 120)
print(ggplot(clpp_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Axon Extension Markers") + ylab("Angiogenesis Markers") + ggtitle("Correlation b/w Angiogenesis and Axon Extension Markers in Paul Mes CLPP") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
png("~/scratch/d_tooth/results/angio_axon_cor_ctrl.png", width = 1000, height = 2000, res = 120)
print(ggplot(ctrl_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Axon Extension Markers") + ylab("Angiogenesis Markers") + ggtitle("Correlation b/w Angiogenesis and Axon Extension Markers in Paul Mes CTRL") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
png("~/scratch/d_tooth/results/angio_axon_cor_grow.png", width = 1000, height = 2000, res = 120)
print(ggplot(grow_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Axon Extension Markers") + ylab("Angiogenesis Markers") + ggtitle("Correlation in Growing Human Molar") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
png("~/scratch/d_tooth/results/angio_axon_cor_adult.png", width = 1000, height = 2000, res = 120)
print(ggplot(adult_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Axon Extension Markers") + ylab("Angiogenesis Markers") + ggtitle("Correlation in Adult Human Molar") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()

clpp_r$cor_dif = clpp_r$cor - ctrl_r$cor
png("~/scratch/d_tooth/results/angio_axon_cor_clpp_ctrl_dif.png", width = 1000, height = 2000, res = 120)
print(ggplot(clpp_r, aes(x=axon, y = angio, fill = cor_dif)) + geom_tile() + scale_fill_viridis_c() + xlab("Axon Extension Markers") + ylab("Angiogenesis Markers") + ggtitle("R Difference b/w Angiogenesis and Axon Extension Markers in Paul Mes") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()

grow_r$cor_dif = grow_r$cor - adult_r$cor
png("~/scratch/d_tooth/results/angio_axon_cor_grow_adult_dif.png", width = 1000, height = 2000, res = 120)
print(ggplot(grow_r, aes(x=axon, y = angio, fill = cor_dif)) + geom_tile() + scale_fill_viridis_c() + xlab("Axon Extension Markers") + ylab("Angiogenesis Markers") + ggtitle("R Difference b/w Growing and Adult Human Molar") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()

# Growth and Neuro recruit
png("~/scratch/d_tooth/results/growth_neuro_cor_clpp.png", width = 1000, height = 2000, res = 120)
print(ggplot(clpp_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Neuro Recruitment Markers") + ylab("Neuro Growth Markers") + ggtitle("Correlations in Paul Mes CLPP") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
png("~/scratch/d_tooth/results/growth_neuro_cor_ctrl.png", width = 1000, height = 2000, res = 120)
print(ggplot(ctrl_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Neuro Recruitment Markers") + ylab("Neuro Growth Markers") + ggtitle("Correlations in Paul Mes CTRL") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
png("~/scratch/d_tooth/results/growth_neuro_cor_grow.png", width = 1000, height = 2000, res = 120)
print(ggplot(grow_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Neuro Recruitment Markers") + ylab("Neuro Growth Markers") + ggtitle("Correlations in Growing Human Molar") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
png("~/scratch/d_tooth/results/growth_neuro_cor_adult.png", width = 1000, height = 2000, res = 120)
print(ggplot(adult_r, aes(x=axon, y = angio, fill = cor)) + geom_tile() + scale_fill_viridis_c() + xlab("Neuro Recruitment Markers") + ylab("Neuro Growth Markers") + ggtitle("Correlations in Adult Human Molar") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()

clpp_r$cor_dif = clpp_r$cor - ctrl_r$cor
png("~/scratch/d_tooth/results/growth_neuro_cor_clpp_ctrl_dif.png", width = 1000, height = 2000, res = 120)
print(ggplot(clpp_r, aes(x=axon, y = angio, fill = cor_dif)) + geom_tile() + scale_fill_viridis_c() + xlab("Neuro Recruitment Markers") + ylab("Angiogenesis Markers") + ggtitle("R Difference b/w CLPP and CTRL") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()

grow_r$cor_dif = grow_r$cor - adult_r$cor
png("~/scratch/d_tooth/results/growth_neuro_cor_grow_adult_dif.png", width = 1000, height = 2000, res = 120)
print(ggplot(grow_r, aes(x=axon, y = angio, fill = cor_dif)) + geom_tile() + scale_fill_viridis_c() + xlab("Neuro Recruitment Markers") + ylab("Angiogenesis Markers") + ggtitle("R Difference b/w Growing and Adult Human Molar") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()


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

# Epi + Mes Stem (OLD)
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

# Epi + Stem (NEW)
mz_epi_mes_stem = c(mz_epi_stem, mz_mes_stem2)
mat = tj_jaw@assays$RNA@counts[mz_epi_mes_stem,]
mat[which(mat > 0)] = 1
tj_jaw$epi_mes_stem_score = colSums(tj_jaw@assays$RNA@data[mz_epi_mes_stem,])
tj_jaw$epi_mes_stem_score2 = colSums(mat)
myFeaturePlot(tj_jaw, "epi_mes_stem_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Epi+Mes Stem Markers in Cichlid Tooth and Jaw")
myFeaturePlot(tj_jaw, "epi_mes_stem_score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Number of Epi+Mes Stem Markers Expressed in Cichlid Tooth and Jaw")
tj_jaw_epi_mes_stem_cells = colnames(tj_jaw)[which(tj_jaw$epi_mes_stem_score >= quantile(tj_jaw$epi_mes_stem_score, 0.95))]
DimPlot(tj_jaw, cells.highlight = setNames(list(tj_jaw_epi_mes_stem_cells), "95th"), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle("Cells in 95th quantile of expression of Epi+Mes Stem Markers")

# Find Celsr1 Quiescent Markers 04/06/2021
igor_incsr_mes = readRDS("~/research/tooth/data/igor_incsr_mes.rds")
igor_incsr_epi = readRDS("~/research/tooth/data/igor_incsr_epi.rds")

findCoexpressedGenes = function(obj, gene, num_top_genes = 10, my.scale = "none") {
  obj$gene = "Other"
  obj$gene[which(obj@assays$RNA@counts[gene,] > 0)] = gene
  Idents(obj) = obj$gene
  deg = FindMarkers(obj, ident.1 = gene, ident.2 = "Other")
  deg$gene = rownames(deg)
  
  gene_deg_sig_pos = as.vector(deg$gene[which(deg$avg_log2FC > 0 & deg$p_val_adj < 0.05)])
  mat = as.matrix(obj@assays$RNA@data[gene_deg_sig_pos,])
  mat = mat[, order(obj@assays$RNA@data[gene,], decreasing = T)]
  p = pheatmap::pheatmap(mat[gene_deg_sig_pos[1:num_top_genes],], scale = my.scale, cluster_rows = F, cluster_cols = F, show_colnames = F, main = paste("Top", gene, "Markers"))
  
  return(list(deg, p))
}

celsr1_res = findCoexpressedGenes("Celsr1")
celsr1_deg = celsr1_res[[1]]
celsr1_res[[2]]

celsr1_res = findCoexpressedGenes(igor_incsr_epi, "Celsr1")
celsr1_deg = celsr1_res[[1]]
celsr1_res[[2]]

gli1_res = findCoexpressedGenes("Gli1", my.scale = "row")
gli1_deg = gli1_res[[1]]
gli1_res[[2]]

thy1_res = findCoexpressedGenes("Thy1")
thy1_deg = thy1_res[[1]]
thy1_res[[2]]

thy1_sig_pos = c("Thy1", "Fgf3", "Ptn", "Wnt10a", "Cdca7", "Foxd1", "Fgf10", "5930403L14Rik", "Etv5", "Ecel1")
mz_thy1_sig_pos = convertHgncDataFrameToMzebra(data.frame(toupper(thy1_sig_pos)), gene_column = 1, gene_names = rownames(tj_jaw), na.rm = T, return_vect = T)
mz_thy1_res = expressionOfList(tj_jaw, mz_thy1_sig_pos, "Igor Mes Thy1 Markers")
mz_thy1_res[[1]]
mz_thy1_res[[2]]
mz_thy1_res[[3]]
write.csv(res_sig, "~/research/tooth/results/")

# Celsr1 markers in TJ_JAW
tj_jaw = readRDS("~/research/tooth/data/tj_jaw.RDS")
tj     = readRDS("~/research/tooth/data/tj.rds")
jaw    = readRDS("~/research/tooth/data/jpool.rds")
celsr1_makers = c("Celsr1", "Notum", "Sall1", "Bmp8a", "D430041D05Rik", "Wnt6", "Nmnat2", "Gsc", "Sms", "Otof")
mz_celsr1_makers = convertHgncDataFrameToMzebra(data.frame(toupper(celsr1_makers)), 1, gene_names = rownames(tj_jaw), return_vect = T, na.rm = T)

mat = tj_jaw@assays$RNA@counts[mz_celsr1_makers,]
mat[which(mat > 0)] = 1
tj_jaw$celsr1_score = colSums(tj_jaw@assays$RNA@data[mz_celsr1_makers,])
tj_jaw$celsr1_score2 = colSums(mat)
myFeaturePlot(tj_jaw, "celsr1_score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Expression of Celsr1 Markers in Cichlid Tooth and Jaw")
myFeaturePlot(tj_jaw, "celsr1_score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle("Number of Celsr1 Markers Expressed in Cichlid Tooth and Jaw")
tj_jaw_celsr1_markers_cells = colnames(tj_jaw)[which(tj_jaw$celsr1_score >= quantile(tj_jaw$celsr1_score, 0.95))]
DimPlot(tj_jaw, cells.highlight = setNames(list(tj_jaw_celsr1_markers_cells), "95th"), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle("Cells in 95th quantile of expression of Celsr1 Markers")

# TAC
tj_jaw_cc = cellCycle(tj_jaw)
tj_jaw$cc = tj_jaw_cc
Idents(tj_jaw) = tj_jaw$cc
DimPlot(tj_jaw, label = T)
DimPlot(tj_jaw, label = T, split.by = "cc")
FeaturePlot(tj_jaw, "ENSMZEG00005009812", order = T) + ggtitle("ENSMZEG00005009812 (Mki67)")
FeaturePlot(tj_jaw, "ENSMZEG00005001142", order = T)

igor_tac_markers = c("Mki67", "Vat1l", "Fam19a4", "Hey2", "Cdh6", "Lrp11", "Cpne5")
igor_epi_tac_markers = c("Mki67", "Cenpf", "1190002F15Rik", "Fam64a", "Ube2c", "Nuf2", "2810417H13Rik", "Birc5", "Anln", "Cdca3")
mz_tac_markers = convertHgncDataFrameToMzebra(data.frame(toupper(igor_tac_markers)), 1, gene_names = rownames(tj_jaw), return_vect = T, na.rm = T)
mz_epi_tac_markers = convertHgncDataFrameToMzebra(data.frame(toupper(igor_epi_tac_markers)), 1, gene_names = rownames(tj_jaw), return_vect = T, na.rm = T)

tac_tj_jaw_res = expressionOfList(tj_jaw, mz_tac_markers, "Igor Epi TAC Markers")
tac_tj_jaw_res[[1]]
tac_tj_jaw_res[[2]]
tac_tj_jaw_res[[3]]

tac_tj_jaw_res = expressionOfList(tj_jaw, mz_epi_tac_markers, "Igor Epi True TAC Markers")
tac_tj_jaw_res[[1]]
tac_tj_jaw_res[[2]]
tac_tj_jaw_res[[3]]

incsr_mes_tac_res = expressionOfList(igor_incsr_mes, c("Mki67"), "Mki67")
igor_incsr_mes_mki67_res = findCoexpressedGenes(igor_incsr_mes, "Mki67")
igor_incsr_mes_mki67_markers = c("Mki67", "Birc5", "Ska1", "Cenpf", "Ccnb1", "Kif23", "Ncapg", "Casc5", "Kif11", "Cdca5")
igor_incsr_mes_mki67_deg_sig_pos = as.vector(igor_incsr_mes_mki67_res[[1]]$gene[which(igor_incsr_mes_mki67_res[[1]]$avg_log2FC > 0 & igor_incsr_mes_mki67_res[[1]]$p_val_adj < 0.05)])[1:10]
mz_mes_tac_markers = convertHgncDataFrameToMzebra(data.frame(toupper(igor_incsr_mes_mki67_deg_sig_pos)), 1, gene_names = rownames(tj_jaw), return_vect = T, na.rm = T)
mz_mes_tac_res = expressionOfList(tj_jaw, mz_mes_tac_markers, "Igor Mes TAC Markers")
mz_mes_tac_res[[2]]
mz_mes_tac_res[[1]]
mz_mes_tac_res[[3]]

mz_all_tac = unique(c(mz_epi_tac_markers, mz_mes_tac_markers))
mz_all_tac_res = expressionOfList(tj_jaw, mz_all_tac, "Igor All TAC Markers")
mz_all_tac_res[[2]]
mz_all_tac_res[[1]]
mz_all_tac_res[[3]]

expressionOfList = function(obj, gene_list, list_name = "Gene List", my_q = 95) {
  mat = obj@assays$RNA@counts[gene_list,]
  mat[which(mat > 0)] = 1
  if (length(gene_list) == 1 ) {
    obj$score = obj@assays$RNA@data[gene_list,]
    obj$score2 = mat
  } else {
    obj$score = colSums(obj@assays$RNA@data[gene_list,])
    obj$score2 = colSums(mat)
  }
  p1 = myFeaturePlot(obj, "score", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle(paste0("Expression of ", list_name))
  p2 = myFeaturePlot(obj, "score2", my.col.pal = pal, na.blank = T, my.pt.size = 3) + ggtitle(paste0("Number of ", list_name, " Expressed"))
  print(quantile(obj$score, 0.95))
  cells = colnames(obj)[which(obj$score >= quantile(obj$score, my_q/100))]
  p3 = DimPlot(obj, cells.highlight = setNames(list(cells), as.character(my_q)), label = T, pt.size = 2, sizes.highlight = 2) + ggtitle(paste0("Cells in ", my_q,"th quantile of expression"))
  return(list(p1, p2, p3, cells))
}

# Tj + JAW State
tj_jaw_qui_cells = tj_jaw_celsr1_markers_cells[which(tj_jaw$seurat_clusters[tj_jaw_celsr1_markers_cells] == 0)]
mz_mes_tac_cells = colnames(tj_jaw)[which(colSums(tj_jaw@assays$RNA@counts[mz_mes_tac_markers,]) > 0)]

tj_jaw$state = "None"
tj_jaw$state[tj_jaw_qui_cells] = "Quiescent"
tj_jaw$state[tj_jaw_epi_mes_stem_cells] = "Stem"
tj_jaw$state[mz_mes_tac_cells] = "TAC"

library("scales")
Idents(tj_jaw) = tj_jaw$state
cytoScoreByIdent(tj_jaw, my_idents = c("Quiescent", "Stem", "TAC"), pt.alpha = 0.5) + ggtitle("CytoTRACE by Cell State in Cichlid Tooth + Jaw")
DimPlot(tj_jaw, order = T, pt.size = 2) + scale_color_manual(values = c("lightgray", gc.ramp <- hue_pal()(3)))

# TJ State
tj$state = "None"
tj$state[tj_jaw_qui_cells[which(tj_jaw_qui_cells %in% colnames(tj))]] = "Quiescent"
tj$state[tj_jaw_epi_mes_stem_cells[which(tj_jaw_epi_mes_stem_cells %in% colnames(tj))]] = "Stem"
tj$state[mz_mes_tac_cells[which(mz_mes_tac_cells %in% colnames(tj))]] = "TAC"

Idents(tj) = tj$state
cytoScoreByIdent(tj, my_idents = c("Stem", "TAC"), pt.alpha = 0.8) + ggtitle("CytoTRACE by Cell State in Cichlid Tooth")

# JAW State
jaw$state = "None"
jaw$state[tj_jaw_qui_cells[which(tj_jaw_qui_cells %in% colnames(jaw))]] = "Quiescent"
jaw$state[tj_jaw_epi_mes_stem_cells[which(tj_jaw_epi_mes_stem_cells %in% colnames(jaw))]] = "Stem"
jaw$state[mz_mes_tac_cells[which(mz_mes_tac_cells %in% colnames(jaw))]] = "TAC"

Idents(jaw) = jaw$state
cytoScoreByIdent(jaw, my_idents = c("Quiescent", "Stem", "TAC"), pt.alpha = 0.8) + ggtitle("CytoTRACE by Cell State in Cichlid Jaw")

# Investigate Igor Mes
my_igor_mes = readRDS("~/research/tooth/data/igor_incsr_mes_my.RDS")
my_igor_mes$celsr = "None"
my_igor_mes$celsr[celsr_cells] = "Positive"
my_igor_mes$celsr[mat_celsr_cells] = "Maturing Celsr1"
my_igor_mes$celsr[pre_celsr_cells] = "Arm Celsr1"
Idents(my_igor_mes) = my_igor_mes$celsr
DimPlot(my_igor_mes)
celsr_res = FindMarkers(my_igor_mes, ident.1 = "Maturing Celsr1", ident.2 = "Arm Celsr1")
celsr_res_sig = celsr_res[which(celsr_res$p_val_adj < 0.05),]

my_igor_mes$celsr2 = "None"
my_igor_mes$celsr2[celsr_cells] = "Positive"
my_igor_mes$celsr2[celsr_cells[which(my_igor_mes@reductions$umap@cell.embeddings[celsr_cells,1] < 2 & my_igor_mes@reductions$umap@cell.embeddings[celsr_cells,1] > -2.5 & my_igor_mes@reductions$umap@cell.embeddings[celsr_cells,2] < -0.5)]] = "Bot Celsr1"
my_igor_mes$celsr2[celsr_cells[which(my_igor_mes$arm[celsr_cells] == "Celsr1")]] = "Arm Celsr1"
Idents(my_igor_mes) = my_igor_mes$celsr2
DimPlot(my_igor_mes)
celsr_res2 = FindMarkers(my_igor_mes, ident.1 = "Bot Celsr1", ident.2 = "Arm Celsr1")
celsr_res_sig2 = celsr_res2[which(celsr_res2$p_val_adj < 0.05),]

thy1_pos = colnames(my_igor_mes)[which(my_igor_mes@assays$RNA@counts["Thy1",] > 0)]
thy1_celsr1 = thy1_pos[which(thy1_pos %in% celsr_cells)]
DimPlot(my_igor_mes, cells.highlight = thy1_celsr1)
my_igor_mes$thy1_celsr1 = "None"
my_igor_mes$thy1_celsr1[thy1_celsr1] = "Both"
Idents(my_igor_mes) = my_igor_mes$thy1_celsr1
thy1_celsr1_res = FindMarkers(my_igor_mes, ident.1 = "Both", ident.2 = "None")
thy1_celsr1_res_sig = thy1_celsr1_res[which(thy1_celsr1_res$p_val_adj < 0.05),]

my_igor_mes$thy1_celsr12 = "None"
my_igor_mes$thy1_celsr12[thy1_celsr1[which(my_igor_mes@reductions$umap@cell.embeddings[thy1_celsr1,2] < 2.5)]] = "Lower"
my_igor_mes$thy1_celsr12[thy1_celsr1[which(my_igor_mes@reductions$umap@cell.embeddings[thy1_celsr1,2] > 2.5)]] = "Upper"
Idents(my_igor_mes) = my_igor_mes$thy1_celsr12
DimPlot(my_igor_mes, order = T)
thy1_celsr1_res2 = FindMarkers(my_igor_mes, ident.1 = "Lower", ident.2 = "None")
thy1_celsr1_res_sig2 = thy1_celsr1_res2[which(thy1_celsr1_res2$p_val_adj < 0.05),]
thy1_celsr1_res3 = FindMarkers(my_igor_mes, ident.1 = "Upper", ident.2 = "None")
thy1_celsr1_res_sig3 = thy1_celsr1_res3[which(thy1_celsr1_res3$p_val_adj < 0.05),]
thy1_celsr1_res4 = FindMarkers(my_igor_mes, ident.1 = "Upper", ident.2 = "Lower")
thy1_celsr1_res_sig4 = thy1_celsr1_res4[which(thy1_celsr1_res4$p_val_adj < 0.05),]

# Markers + Cyto in Igor Epi
igor_epi_tac_markers = c("Mki67", "Cenpf", "1190002F15Rik", "Fam64a", "Ube2c", "Nuf2", "2810417H13Rik", "Birc5", "Anln", "Cdca3")
igor_mes_tac_markers = c("Mki67", "Birc5", "Ska1", "Cenpf", "Ccnb1", "Kif23", "Ncapg", "Casc5", "Kif11", "Cdca5")
igor_tac_markers = unique(c(igor_epi_tac_markers, igor_mes_tac_markers))
epi_stem_Lgr5 = c("Pknox2", "Spock1", "Sfrp5", "Grp", "Arhgef33", "Ccdc80", "Pcp4", "Pla2g4a", "Cd27", "Vsig2", "Trpm5", "Lgr5", "Gm1110", "Clcnkb", "Fez1", "Icam2", "Prob1", "Disc1", "Zscan10", "Slc35f3")
epi_stem_Sox2 = c("Gjb3", "Sox2", "Ccdc112", "Trpm4", "Gstm2", "AA465934", "Mex3a", "Myh10", "Cks1b", "Moxd1", "Bex1", "Mdk", "Ccdc34", "Dlgap5", "Uhrf1", "Rfc4", "Rrm2", "Gins2", "Rad54b", "Dtl")
epi_stem = c(epi_stem_Lgr5, epi_stem_Sox2, "Lrig1")
mes_stem2 = c("Thy1", "Gli1", "Foxd1", "Bcl2", "Fgf10", "Rimbp2", "Gem", "Hs3st6", "Etv5", "Ptn")
igor_stem = unique(c(epi_stem, mes_stem2))
celsr1_makers = c("Celsr1", "Notum", "Sall1", "Bmp8a", "D430041D05Rik", "Wnt6", "Nmnat2", "Gsc", "Sms", "Otof")
cols = gc.ramp <- hue_pal()(3)

igor_epi_tac_res = expressionOfList(igor_incsr_epi, gene_list = igor_tac_markers, list_name = "All TAC Markers")
igor_epi_tac_cells = igor_epi_tac_res[[4]]
igor_epi_stem_res = expressionOfList(igor_incsr_epi, gene_list = igor_stem, list_name = "All Stem Markers")
igor_epi_stem_cells = igor_epi_stem_res[[4]]
igor_epi_qui_res = expressionOfList(igor_incsr_epi, gene_list = celsr1_makers, list_name = "Celsr1 Markers")
igor_epi_qui_cells = igor_epi_qui_res[[4]]
all_cells = c(igor_epi_tac_cells, igor_epi_stem_cells, igor_epi_qui_cells)
table(table(all_cells))

igor_incsr_epi$state = "None"
igor_incsr_epi$state[igor_epi_tac_cells] = "TAC"
igor_incsr_epi$state[igor_epi_stem_cells] = "Stem"
igor_incsr_epi$state[igor_epi_qui_cells] = "Quiescent"
igor_incsr_epi$state = factor(igor_incsr_epi$state, levels = c("None", "TAC", "Stem", "Quiescent"))
Idents(igor_incsr_epi) = igor_incsr_epi$state
DimPlot(igor_incsr_epi, order = T, pt.size = 1.5) + scale_color_manual(values = c("lightgray", cols[2], cols[1], cols[3]))
cytoScoreByIdent(igor_incsr_epi, my_idents = c("TAC", "Stem", "Quiescent"), pt.alpha = 0.8) + ggtitle("Igor Incisor Epithelium - State")
t.test(x=igor_incsr_epi$cyto[which(igor_incsr_epi$state == "Quiescent")], y=igor_incsr_epi$cyto[which(igor_incsr_epi$state == "TAC")])
t.test(x=igor_incsr_epi$cyto[which(igor_incsr_epi$state == "Quiescent")], y=igor_incsr_epi$cyto[which(igor_incsr_epi$state == "Stem")])
t.test(x=igor_incsr_epi$cyto[which(igor_incsr_epi$state == "Stem")], y=igor_incsr_epi$cyto[which(igor_incsr_epi$state == "TAC")])

igor_mes_tac_res = expressionOfList(igor_incsr_mes, gene_list = igor_tac_markers, list_name = "All TAC Markers")
igor_mes_tac_cells = igor_mes_tac_res[[4]]
igor_mes_stem_res = expressionOfList(igor_incsr_mes, gene_list = igor_stem, list_name = "All Stem Markers")
igor_mes_stem_cells = igor_mes_stem_res[[4]]
igor_mes_qui_res = expressionOfList(igor_incsr_mes, gene_list = celsr1_makers, list_name = "Celsr1 Markers")
igor_mes_qui_cells = igor_mes_qui_res[[4]]
all_cells = c(igor_mes_tac_cells, igor_mes_stem_cells, igor_mes_qui_cells)
table(table(all_cells))

igor_incsr_mes$state = "None"
igor_incsr_mes$state[igor_mes_tac_cells] = "TAC"
igor_incsr_mes$state[igor_mes_stem_cells] = "Stem"
igor_incsr_mes$state[igor_mes_qui_cells] = "Quiescent"
igor_incsr_mes$state = factor(igor_incsr_mes$state, levels = c("None", "TAC", "Stem", "Quiescent"))
Idents(igor_incsr_mes) = igor_incsr_mes$state
DimPlot(igor_incsr_mes, order = T, pt.size = 1.5) + scale_color_manual(values = c("lightgray", cols[2], cols[1], cols[3]))
cytoScoreByIdent(igor_incsr_mes, my_idents = c("TAC", "Stem", "Quiescent"), pt.alpha = 0.8) + ggtitle("Igor Incisor Mesenchyme - State")
t.test(x=igor_incsr_mes$cyto[which(igor_incsr_mes$state == "Quiescent")], y=igor_incsr_mes$cyto[which(igor_incsr_mes$state == "TAC")])
t.test(x=igor_incsr_mes$cyto[which(igor_incsr_mes$state == "Quiescent")], y=igor_incsr_mes$cyto[which(igor_incsr_mes$state == "Stem")])
t.test(x=igor_incsr_mes$cyto[which(igor_incsr_mes$state == "Stem")], y=igor_incsr_mes$cyto[which(igor_incsr_mes$state == "TAC")])

igor_tac_res = expressionOfList(igor_incsr, gene_list = igor_tac_markers, list_name = "All TAC Markers")
igor_tac_cells = igor_tac_res[[4]]
igor_stem_res = expressionOfList(igor_incsr, gene_list = igor_stem, list_name = "All Stem Markers")
igor_stem_cells = igor_stem_res[[4]]
igor_qui_res = expressionOfList(igor_incsr, gene_list = celsr1_makers, list_name = "Celsr1 Markers")
igor_qui_cells = igor_qui_res[[4]]
all_cells = c(igor_tac_cells, igor_stem_cells, igor_qui_cells)
table(table(all_cells))

igor_incsr$state = "None"
igor_incsr$state[igor_tac_cells] = "TAC"
igor_incsr$state[igor_stem_cells] = "Stem"
igor_incsr$state[igor_qui_cells] = "Quiescent"
igor_incsr$state = factor(igor_incsr$state, levels = c("None", "TAC", "Stem", "Quiescent"))
Idents(igor_incsr) = igor_incsr$annot
DimPlot(igor_incsr, cells.highlight = setNames(list(igor_tac_cells, igor_stem_cells, igor_qui_cells), c("TAC", "Stem", "Quiescent")), cols.highlight = cols, label = T, pt.size = 1, sizes.highlight = 1)
Idents(igor_incsr) = igor_incsr$state
cytoScoreByIdent(igor_incsr, my_idents = c("TAC", "Stem", "Quiescent"), pt.alpha = 0.8) + ggtitle("Igor Incisor - State")

# Combined
tac_res = findCoexpressedGenes(igor_incsr, "Mki67")
tac_res_deg = tac_res[[1]]
tac_res_deg$hgnc = toupper(tac_res_deg$gene)
tac_res_deg[3, "hgnc"] = "PCLAF"
combined_tac = c("MKI67", "BIRC5", "PCLAF", "NUF2", "CCNB1", "CASC5", "CDCA3", "SKA1", "CDCA8", "CCNB2")
combined_tac_mouse =  c("Mki67", "Birc5", "2810417H13Rik", "Nuf2", "Ccnb1", "Casc5", "Cdca3", "Ska1", "Cdca8", "Ccnb2")

combined_stem = igor_stem

qui_res = findCoexpressedGenes(igor_incsr, "Celsr1")
qui_res_deg = qui_res[[1]]
qui_res_deg$hgnc = toupper(qui_res_deg$gene)
# qui_res_deg[3, "hgnc"] = "PCLAF"
combined_qui = c("CELSR1", "PITX2", "ISL1", "TRP63", "PROM2", "FXYD3", "TBX1", "CDH1", "FERMT1", "SFN")


igor_tac_res = expressionOfList(igor_incsr, gene_list = combined_tac_mouse, list_name = "Mki67 Markers")
igor_tac_cells = igor_tac_res[[4]]
igor_stem_res = expressionOfList(igor_incsr, gene_list = igor_stem, list_name = "All Stem Markers")
igor_stem_cells = igor_stem_res[[4]]
igor_qui_res = expressionOfList(igor_incsr, gene_list = str_to_title(combined_qui), list_name = "Celsr1 Markers")
igor_qui_cells = igor_qui_res[[4]]
igor_incsr$state = "None"
igor_incsr$state[igor_tac_cells] = "TAC"
igor_incsr$state[igor_stem_cells] = "Stem"
igor_incsr$state[igor_qui_cells] = "Quiescent"
igor_incsr$state = factor(igor_incsr$state, levels = c("None", "TAC", "Stem", "Quiescent"))
Idents(igor_incsr) = igor_incsr$annot
DimPlot(igor_incsr, cells.highlight = setNames(list(igor_tac_cells, igor_stem_cells, igor_qui_cells), c("TAC", "Stem", "Quiescent")), cols.highlight = cols, label = T, pt.size = 1, sizes.highlight = 1)
Idents(igor_incsr) = igor_incsr$state
cytoScoreByIdent(igor_incsr, my_idents = c("TAC", "Stem", "Quiescent"), pt.alpha = 0.8) + ggtitle("Igor Incisor - State")

# Quiescence
# ccAF
tj_jaw$cc = read.csv("~/research/tooth/results/tj_jaw_ccAf.txt", header = F)[,1]
tj_jaw$cc = factor(tj_jaw$cc, levels = c(unique(tj_jaw$cc)[-c(5)], "Neural G0"))
Idents(tj_jaw) = tj_jaw$cc
DimPlot(tj_jaw, order = T, label = T, pt.size = 2)
DimPlot(tj_jaw, order = T, label = T, pt.size = 2, cols = c("lightgray", "lightgray", "lightgray", "lightgray", "lightgray", "lightgray", "lightgray", "red")) + ggtitle("ccAF")
DimPlot(tj_jaw, order = T, label = T, pt.size = 2) + ggtitle("ccAF")

jaw$cc = read.csv("~/research/tooth/results/jaw_ccAf.txt", header = F)[,1]
jaw$cc = factor(jaw$cc, levels = c(unique(jaw$cc)[-c(4)], "Neural G0"))
Idents(jaw) = jaw$cc
DimPlot(jaw, order = T, label = T, pt.size = 2)
DimPlot(jaw, order = T, label = T, pt.size = 2, cols = c("lightgray", "lightgray", "lightgray", "lightgray", "lightgray", "lightgray", "lightgray", "red")) + ggtitle("ccAF")
DimPlot(jaw, order = T, label = T, pt.size = 2) + ggtitle("ccAF")

# RNA Velocity
jaw_celsr = read.csv("~/research/tooth/results/jaw_celsr1_velo.csv")
library("parallel")
numCores = detectCores()
jaw_celsr$idx = unlist(mclapply(jaw_celsr[,1], function(x) grep(x, colnames(jaw))[1], mc.cores = numCores))
jaw_celsr$cell = colnames(jaw)[jaw_celsr$idx]
jaw_celsr$exp = jaw@assays$RNA@data["celsr1a",match(jaw_celsr$cell, colnames(jaw))]
jaw_celsr$raw = jaw@assays$RNA@counts["celsr1a",match(jaw_celsr$cell, colnames(jaw))]
jaw_celsr$celsr1a_at_pos = NA
jaw_celsr$celsr1a_at_pos[which(jaw_celsr$raw >0)] = jaw_celsr$celsr1a[which(jaw_celsr$raw > 0)]

jaw$celsr1a_velo = jaw_celsr$celsr1a[match(colnames(jaw), jaw_celsr$cell)]
Idents(jaw) = jaw$seurat_clusters
FeaturePlot(jaw, "celsr1a_velo", order = T, label = T) + scale_color_gradientn(colors = pal(50))
myFeaturePlot(jaw, "celsr1a_velo", my.col.pal = pal, my.pt.size = 2) + ggtitle("Celsr1a Velocity")
myFeaturePlot(jaw, "celsr1a_velo", my.col.pal = pal, my.pt.size = 2) + ggtitle("Quiescent on Bottom")

jaw$celsr1a_velo2= jaw_celsr$celsr1a_at_pos[match(colnames(jaw), jaw_celsr$cell)]
myFeaturePlot(jaw, "celsr1a_velo2", my.col.pal = pal, my.pt.size = 2) + ggtitle("Velocity at Celsr1+ Cells")

# All RNA velo
jaw_velo = read.csv("~/research/tooth/results/jaw_velo.csv")
jaw_velo$idx = jaw_celsr$idx
jaw_velo$cell = jaw_celsr$cell
jaw$all_velo = jaw_velo$velocity_self_transition[match(colnames(jaw), jaw_velo$cell)]
myFeaturePlot(jaw, "all_velo", my.col.pal = pal, my.pt.size = 2) + ggtitle("All Velocity")

# Talha Qui 
myFeaturePlot(jaw, "cdk6", my.col.pal = pal, my.pt.size = 2, na.blank = T)
myFeaturePlot(jaw, "cdk4", my.col.pal = pal, my.pt.size = 2, na.blank = T)
jaw$cdks = colSums(jaw@assays$RNA@data[c("cdk4", "cdk6"),])
myFeaturePlot(jaw, "cdks", my.col.pal = pal, my.pt.size = 2, na.blank = T) + ggtitle("cdk4 + cdk6")
myFeaturePlot(jaw, "celsr1a", my.col.pal = pal, my.pt.size = 2, na.blank = T)

markerCellPerCluster(jaw, c("cdk4")) + ggtitle("Normalized Number of Cells Expressing Cdk4")
markerCellPerCluster(jaw, c("cdk6")) + ggtitle("Normalized Number of Cells Expressing Cdk6")
markerCellPerCluster(jaw, c("cdk4", "cdk6")) + ggtitle("Normalized Number of Cells Expressing Cdk4 and Cdk6")

markerLogFC(jaw, c("cdk4"))[[1]] + ggtitle("LogFC for Cdk4")
markerLogFC(jaw, c("cdk6"))[[1]] + ggtitle("LogFC for Cdk6")
markerLogFC(jaw, c("cdk4", "cdk6"))[[1]] + ggtitle("LogFC for Cdk4 and Cdk6")

new_qui = c("cdkn1bb", "rb1", "tp53", "cdkn1cb")
expressionOfList(jaw, gene_list = new_qui, list_name = "New Quiescent Markers")

expressionOfList(jaw, gene_list = s.genes, list_name = "S Markers")

all_cc = unique(s.genes, g2m.genes)
expressionOfList(jaw, gene_list = all_cc, list_name = "All Cell Cycle Markers")
top_long$hgnc = jaw_annot_deg_sig$hgnc[match(top_long$value, jaw_annot_deg_sig$gene)]



# High Quality Plot for figures
tj_jaw$seurat_clusters = factor(tj_jaw$seurat_clusters, levels = 0:10)
Idents(tj_jaw) = tj_jaw$seurat_clusters
tj_jaw = RenameIdents(tj_jaw, '0' = "Epithelial", '1' = "Epithelial", '2' = "Mesenchymal", '3' = "Immune", '4' = "Immune", "5" = "Epithelial", "6" = "Epithelial", "7" = "Shared", "8"= "Pigmented", "9" = "Mature TB", "10" = "Immature TB")
svg("~/research/tooth/results/tj_jaw_annot.svg")
DimPlot(tj_jaw, pt.size = 1.5, label = F) + coord_fixed()
dev.off()

mz_thy1  = "ENSMZEG00005018875"
mz_eng   = "ENSMZEG00005019543"
mz_nt5e  = "nt5e"
mz_mki67 = "ENSMZEG00005009812"
mz_ccl2 = "ENSMZEG00005000472"

#==================================================================================================
# Celsr1 by CytoBIN ===============================================================================
#==================================================================================================
# Load Single Cell Datasets
incsr = readRDS("~/research/tooth/data/igor_incsr.rds")
im = readRDS("~/research/tooth/data/igor_incsr_molar.rds")
hm = readRDS("~/research/tooth/data/hm.rds")
igor_incsr_mes = readRDS("~/research/tooth/data/igor_incsr_mes.rds")
igor_incsr_epi = readRDS("~/research/tooth/data/igor_incsr_epi.rds")
tj_jaw = readRDS("~/research/tooth/data/tj_jaw.RDS")
tj     = readRDS("~/research/tooth/data/tj.rds")
jaw    = readRDS("~/research/tooth/data/jpool.rds")

# Load CytoTRACE data
incsr_cyto = readRDS("~/research/tooth/data/incsr_cyto.rds")
im_cyto = readRDS("~/research/tooth/data/im_cyto.rds")
hm_cyto = readRDS("~/research/tooth/data/hm_cyto.rds")
igor_incsr_epi_cyto = readRDS("~/research/tooth/data/igor_incsr_epi_cyto.rds")
igor_incsr_mes_cyto = readRDS("~/research/tooth/data/igor_incsr_mes_cyto.rds")
tj_jaw_cyto = CytoTRACE(as.matrix(tj_jaw@assays$RNA@counts))
tj_cyto = CytoTRACE(as.matrix(tj@assays$RNA@counts))
jaw_cyto = CytoTRACE(as.matrix(jaw@assays$RNA@counts))

# Add CytoTRACE score as metadata
incsr$cyto = incsr_cyto$CytoTRACE
im$cyto = im_cyto$CytoTRACE
hm$cyto = hm_cyto$CytoTRACE
igor_incsr_epi$cyto = igor_incsr_epi_cyto$CytoTRACE
igor_incsr_mes$cyto = igor_incsr_mes_cyto$CytoTRACE
tj_jaw$cyto = tj_jaw_cyto$CytoTRACE
tj$cyto = tj_cyto$CytoTRACE
jaw$cyto = jaw_cyto$CytoTRACE

incsr_celsr1_bin_res = splitGenebyCytoBIN(incsr, "Celsr1")
im_celsr1_bin_res = splitGenebyCytoBIN(im, "CELSR1")
hm_celsr1_bin_res = splitGenebyCytoBIN(hm, "CELSR1")
igor_incsr_epi_celsr1_bin_res = splitGenebyCytoBIN(igor_incsr_epi, "Celsr1")
jaw_celsr1_bin_res = splitGenebyCytoBIN(jaw, "celsr1a")
tj_jaw_celsr1_bin_res = splitGenebyCytoBIN(tj_jaw, "celsr1a")

incsr_celsr1_bin_deg = incsr_celsr1_bin_res[[3]]
im_celsr1_bin_deg = im_celsr1_bin_res[[3]]
hm_celsr1_bin_deg = hm_celsr1_bin_res[[3]]
igor_incsr_epi_celsr1_bin_deg = igor_incsr_epi_celsr1_bin_res[[3]]
igor_incsr_mes_celsr1_bin_deg = igor_incsr_mes_celsr1_bin_res[[3]]
jaw_celsr1_bin_deg = jaw_celsr1_bin_res[[3]]
tj_jaw_celsr1_bin_deg = tj_jaw_celsr1_bin_res[[3]]

jaw_celsr1_bin_deg = degWithHgncAndDescription(jaw_celsr1_bin_deg, rownames(tj))

write.csv(incsr_celsr1_bin_deg, "~/research/tooth/results/incsr_celsr1_cytoBIN_deg.csv")
write.csv(im_celsr1_bin_deg, "~/research/tooth/results/im_celsr1_cytoBIN_deg.csv")
write.csv(hm_celsr1_bin_deg, "~/research/tooth/results/hm_celsr1_cytoBIN_deg.csv")
write.csv(igor_incsr_epi_celsr1_bin_deg, "~/research/tooth/results/igor_incsr_epi_celsr1_cytoBIN_deg.csv")
write.csv(igor_incsr_mes_celsr1_bin_deg, "~/research/tooth/results/igor_incsr_mes_celsr1_cytoBIN_deg.csv")
write.csv(jaw_celsr1_bin_deg, "~/research/tooth/results/jaw_celsr1_cytoBIN_deg.csv")

incsr_celsr1_bin_deg_pos = incsr_celsr1_bin_deg[which(incsr_celsr1_bin_deg$avg_logFC > 0),]
im_celsr1_bin_deg_pos = im_celsr1_bin_deg[which(im_celsr1_bin_deg$avg_logFC > 0),]
hm_celsr1_bin_deg_pos = hm_celsr1_bin_deg[which(hm_celsr1_bin_deg$avg_logFC > 0),]
igor_incsr_epi_celsr1_bin_deg_pos = igor_incsr_epi_celsr1_bin_deg[which(igor_incsr_epi_celsr1_bin_deg$avg_logFC > 0),]
igor_incsr_mes_celsr1_bin_deg_pos = igor_incsr_mes_celsr1_bin_deg[which(igor_incsr_mes_celsr1_bin_deg$avg_logFC > 0),]
jaw_celsr1_bin_deg_pos = jaw_celsr1_bin_deg[which(jaw_celsr1_bin_deg$avg_logFC > 0),]
tj_jaw_celsr1_bin_deg_pos = tj_jaw_celsr1_bin_deg[which(tj_jaw_celsr1_bin_deg$avg_logFC > 0),]

write.csv(incsr_celsr1_bin_deg_pos, "~/research/tooth/results/incsr_celsr1_cytoBIN_deg_pos.csv")
write.csv(im_celsr1_bin_deg_pos, "~/research/tooth/results/im_celsr1_cytoBIN_deg_pos.csv")
write.csv(hm_celsr1_bin_deg_pos, "~/research/tooth/results/hm_celsr1_cytoBIN_deg_pos.csv")
write.csv(jaw_celsr1_bin_deg_pos, "~/research/tooth/results/jaw_celsr1_cytoBIN_deg_pos.csv")
write.csv(tj_jaw_celsr1_bin_deg_pos, "~/research/tooth/results/tj_jaw_celsr1_cytoBIN_deg_pos.csv")
write.csv(igor_incsr_epi_celsr1_bin_deg_pos, "~/research/tooth/results/igor_incsr_epi_celsr1_cytoBIN_deg_pos.csv")
write.csv(igor_incsr_mes_celsr1_bin_deg_pos, "~/research/tooth/results/igor_incsr_mes_celsr1_cytoBIN_deg_pos.csv")

jaw_celsr1_bin_deg_pos$im_cluster = im_celsr1_bin_deg_pos$cluster[match(jaw_celsr1_bin_deg_pos$hgnc, im_celsr1_bin_deg_pos$gene)]
jaw_celsr1_bin_deg_pos$incsr_cluster = incsr_celsr1_bin_deg_pos$cluster[match(jaw_celsr1_bin_deg_pos$hgnc, toupper(incsr_celsr1_bin_deg_pos$gene))]
jaw_celsr1_bin_deg_pos$igor_incsr_epi_cluster = igor_incsr_epi_celsr1_bin_deg_pos$cluster[match(jaw_celsr1_bin_deg_pos$hgnc, toupper(igor_incsr_epi_celsr1_bin_deg_pos$gene))]
jaw_celsr1_bin_deg_pos$igor_incsr_mes_cluster = igor_incsr_mes_celsr1_bin_deg_pos$cluster[match(jaw_celsr1_bin_deg_pos$hgnc, toupper(igor_incsr_mes_celsr1_bin_deg_pos$gene))]

tj_jaw_celsr1_bin_deg_pos$im_cluster = im_celsr1_bin_deg_pos$cluster[match(tj_jaw_celsr1_bin_deg_pos$hgnc, im_celsr1_bin_deg_pos$gene)]
tj_jaw_celsr1_bin_deg_pos$incsr_cluster = incsr_celsr1_bin_deg_pos$cluster[match(tj_jaw_celsr1_bin_deg_pos$hgnc, toupper(incsr_celsr1_bin_deg_pos$gene))]
tj_jaw_celsr1_bin_deg_pos$igor_incsr_epi_cluster = igor_incsr_epi_celsr1_bin_deg_pos$cluster[match(tj_jaw_celsr1_bin_deg_pos$hgnc, toupper(igor_incsr_epi_celsr1_bin_deg_pos$gene))]
tj_jaw_celsr1_bin_deg_pos$igor_incsr_mes_cluster = igor_incsr_mes_celsr1_bin_deg_pos$cluster[match(tj_jaw_celsr1_bin_deg_pos$hgnc, toupper(igor_incsr_mes_celsr1_bin_deg_pos$gene))]

incsr_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/incsr_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
im_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/im_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
jaw_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/jaw_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
tj_jaw_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/tj_jaw_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
igor_incsr_epi_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/igor_incsr_epi_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")

jaw_high_toppgene$inIncsr = jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & jaw_high_toppgene$Name %in% incsr_high_toppgene$Name
jaw_high_toppgene$inIM = jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & jaw_high_toppgene$Name %in% im_high_toppgene$Name

tj_jaw_high_toppgene$inIncsr = tj_jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & tj_jaw_high_toppgene$Name %in% incsr_high_toppgene$Name
tj_jaw_high_toppgene$inIM = tj_jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & tj_jaw_high_toppgene$Name %in% im_high_toppgene$Name
tj_jaw_high_toppgene$inEpi = tj_jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & tj_jaw_high_toppgene$Name %in% igor_incsr_epi_high_toppgene$Name

length(which(jaw_high_toppgene$inIncsr))
length(which(jaw_high_toppgene$inIM))
length(which(jaw_high_toppgene$inIncsr & jaw_high_toppgene$inIM))

length(which(tj_jaw_high_toppgene$inIncsr))
length(which(tj_jaw_high_toppgene$inIM))
length(which(tj_jaw_high_toppgene$inIncsr & tj_jaw_high_toppgene$inIM))


# Do the same, except remove the low CytoBIN
incsr_celsr1_bin_no_low_res = splitGenebyCytoBIN(incsr, "Celsr1", rm.low = T)
im_celsr1_bin_no_low_res = splitGenebyCytoBIN(im, "CELSR1", rm.low = T)
hm_celsr1_bin_no_low_res = splitGenebyCytoBIN(hm, "CELSR1", rm.low = T)
igor_incsr_epi_celsr1_bin_no_low_res = splitGenebyCytoBIN(igor_incsr_epi, "Celsr1", rm.low = T)
igor_incsr_mes_celsr1_bin_no_low_res = splitGenebyCytoBIN(igor_incsr_mes, "Celsr1", rm.low = T)
jaw_celsr1_bin_no_low_res = splitGenebyCytoBIN(jaw, "celsr1a", rm.low = T)

incsr_celsr1_bin_no_low_deg = incsr_celsr1_bin_no_low_res[[3]]
im_celsr1_bin_no_low_deg = im_celsr1_bin_no_low_res[[3]]
hm_celsr1_bin_no_low_deg = hm_celsr1_bin_no_low_res[[3]]
igor_incsr_epi_celsr1_bin_no_low_deg = igor_incsr_epi_celsr1_bin_no_low_res[[3]]
igor_incsr_mes_celsr1_bin_no_low_deg = igor_incsr_mes_celsr1_bin_no_low_res[[3]]
jaw_celsr1_bin_no_low_deg = jaw_celsr1_bin_no_low_res[[3]]

incsr_celsr1_bin_no_low_deg_pos = incsr_celsr1_bin_no_low_deg[which(incsr_celsr1_bin_no_low_deg$avg_logFC > 0),]
im_celsr1_bin_no_low_deg_pos = im_celsr1_bin_no_low_deg[which(im_celsr1_bin_no_low_deg$avg_logFC > 0),]
hm_celsr1_bin_no_low_deg_pos = hm_celsr1_bin_no_low_deg[which(hm_celsr1_bin_no_low_deg$avg_logFC > 0),]
igor_incsr_epi_celsr1_bin_no_low_deg_pos = igor_incsr_epi_celsr1_bin_no_low_deg[which(igor_incsr_epi_celsr1_bin_no_low_deg$avg_logFC > 0),]
igor_incsr_mes_celsr1_bin_no_low_deg_pos = igor_incsr_mes_celsr1_bin_no_low_deg[which(igor_incsr_mes_celsr1_bin_no_low_deg$avg_logFC > 0),]
jaw_celsr1_bin_no_low_deg_pos = jaw_celsr1_bin_no_low_deg[which(jaw_celsr1_bin_no_low_deg$avg_logFC > 0),]

write.csv(incsr_celsr1_bin_no_low_deg_pos, "~/research/tooth/results/incsr_celsr1_cytoBIN_no_low_deg_pos.csv")
write.csv(im_celsr1_bin_no_low_deg_pos, "~/research/tooth/results/im_celsr1_cytoBIN_no_low_deg_pos.csv")
write.csv(hm_celsr1_bin_no_low_deg_pos, "~/research/tooth/results/hm_celsr1_cytoBIN_no_low_deg_pos.csv")
write.csv(igor_incsr_epi_celsr1_bin_no_low_deg_pos, "~/research/tooth/results/igor_incsr_epi_celsr1_cytoBIN_no_low_deg_pos.csv")
write.csv(igor_incsr_mes_celsr1_bin_no_low_deg_pos, "~/research/tooth/results/igor_incsr_mes_celsr1_cytoBIN_no_low_deg_pos.csv")
write.csv(jaw_celsr1_bin_no_low_deg_pos, "~/research/tooth/results/jaw_celsr1_cytoBIN_no_low_deg_pos.csv")


splitGenebyCytoBIN = function(obj, gene, rm.low = FALSE) {
  #' Split Gene+ cells by cells with relative high, medium, and low CytoTRACE scores.
  #' Find DEGs between the BINs.
  #' @param obj seurat object
  #' @param gene gene to subset by
  #' @param rm.low remove low CytoBIN? ie only compare medium to high CytoBINs?
  #' @return 
  
  gene_pos_cells = colnames(obj)[which(obj@assays$RNA@counts[gene,] > 0)]
  gene_obj = subset(obj, cells = gene_pos_cells)

  # Define the CytoBINs by the quantile of CytoTRACE score
  gene_obj$bin <- gene_obj$cyto
  gene_obj$bin[which(gene_obj$cyto <= quantile(gene_obj$cyto, 0.33))] <- "relative_low"
  gene_obj$bin[which(gene_obj$cyto > quantile(gene_obj$cyto, 0.33) & gene_obj$cyto <= quantile(gene_obj$cyto, 0.66))] <- "relative_medium"
  gene_obj$bin[which(gene_obj$cyto > quantile(gene_obj$cyto, 0.66))] <- "relative_high"
  
  # Remove low CytoBIN if necessary
  Idents(gene_obj) = gene_obj$bin
  if (rm.low) {
    gene_obj = subset(gene_obj, idents = "relative_low", invert = T)
  }
  
  # Find DEGs between CytoBINS
  bin_deg = FindAllMarkers(gene_obj)
  bin_deg_sig = bin_deg[which(bin_deg$p_val_adj < 0.05),]
  colnames(bin_deg_sig)[which(colnames(bin_deg_sig) == "avg_log2FC")] = "avg_logFC"
  
  # Visualize Gene and Cyto
  obj$gene_cyto = NA
  obj$gene_cyto[gene_pos_cells] = obj$cyto[gene_pos_cells]
  p = myFeaturePlot(obj, "gene_cyto", my.col.pal = pal, na.blank = T) + ggtitle("CytoTRACE in Celsr1+ Cells")
  p1 = cytoScoreByIdent(gene_obj, pt.alpha = 0.2)
  
  return(list(p, p1, bin_deg_sig))
}