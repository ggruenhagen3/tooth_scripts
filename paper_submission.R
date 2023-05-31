# Load Packages
library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("biomaRt")
library("stringr")
library("dplyr")
library("CytoTRACE")
library("ggplot2")
library('RColorBrewer')

# Load Seurat Objects
combined <- readRDS("C:/Users/miles/Downloads/rna/data/combined.rds")
tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
jpool <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/pool_shiny/data/jpool.rds")
mes <- readRDS("C:/Users/miles/Downloads/rna/data/combined.rds")
epi <- readRDS("C:/Users/miles/Downloads/d_tooth/data/epi_full.rds")
epi$cond[is.na(epi$cond)] <- "INJR"
jaw <- jpool
mes <- combined

# Load CytoTRACE Data
mes_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/mesenchyme_cyto.rds")
epi_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/epithelium_cyto.rds")
tj_cyto  <- readRDS("C:/Users/miles/Downloads/d_tooth/data/tooth_cyto.rds")
jaw_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/jaw_cyto.rds")

# Save CytoTRACE Data As Metadata in the Seurat Object
mes$cyto <- mes_cyto$CytoTRACE
epi$cyto <- epi_cyto$CytoTRACE
tj$cyto  <- tj_cyto$CytoTRACE
jaw$cyto <- jaw_cyto$CytoTRACE

# Figure Scratch Folder
fig_path = "C:/Users/miles/Downloads/rna/paper_submission/figures/"
scratch  = paste0(fig_path, "scratch/")

############
# Figure 1 #
############
# 4 columns Mouse mes (combo treatments), mouse epi (combo treatments), cichlid jaw, cichlid tooth.
# 4 rows Clusters, YAP paint, piezo1 paint, Cytotrace
png_w = 1000
png_h = 600
png_res = 100
lbl_size = 9
text_size = 20

# Row 1 - Clusters
all_p <- list()
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_clusters_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = DimPlot(ob.list[[x]], label = TRUE, pt.size = 1.5) + labs(title = paste("Mouse Mesenchyme -", cond.list[[x]])) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'))
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_clusters.png"), width = png_w, height = png_h, res = png_res)
p <- DimPlot(mes, label = TRUE, pt.size = 1.5, label.size = lbl_size) + labs(title = paste("Mouse Mesenchyme")) + theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_clusters_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p <- DimPlot(ob.list[[x]], label = TRUE, pt.size = 1.5) + labs(title = paste("Mouse Epithelium -", cond.list[[x]])) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'))
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_clusters.png"), width = png_w, height = png_h, res = png_res)
p <- DimPlot(epi, label = TRUE, pt.size = 1.5, label.size = lbl_size) + labs(title = paste("Mouse Epithelium")) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_clusters.png"), width = png_w, height = png_h, res = png_res)
tj_p <- DimPlot(tj, label = TRUE, pt.size = 2.5, label.size = lbl_size) + labs(title = paste("Cichlid Tooth")) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(tj_p)
dev.off()

png(paste0(scratch, "jaw_clusters.png"), width = png_w, height = png_h, res = png_res)
jaw_p <- DimPlot(jaw, label = TRUE, pt.size = 2.5, label.size = lbl_size) + labs(title = paste("Cichlid Jaw")) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(jaw_p)
dev.off()

all_p[[length(all_p)+1]] <- tj_p
all_p[[length(all_p)+1]] <- jaw_p

# Row 2 - Yap
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_yap_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], feature = "Yap1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_yap.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(mes, feature = "Yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_yap_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p <- FeaturePlot(ob.list[[x]], features = "Yap1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_yap.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(epi, feature = "Yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_yap.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(tj, features = "yap1", label = TRUE, pt.size = 2.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "jaw_yap.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(jaw, features = "yap1", label = TRUE, pt.size = 2.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# Row 3 - Piezo
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_piezo_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], feature = "Piezo1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_piezo.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(mes, feature = "Piezo1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_piezo_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p <- FeaturePlot(ob.list[[x]], features = "Piezo1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_piezo.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(epi, feature = "Piezo1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_piezo.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(tj, features = "piezo1", label = TRUE, pt.size = 2.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "jaw_piezo.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(jaw, features = "piezo1", label = TRUE, pt.size = 2.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# Row 4 - CytoTRACE
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_cyto_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank())
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(mes, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_cyto_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank())
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(epi, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(tj, features = "cyto", reduction = "umap", pt.size = 2.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "jaw_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(jaw, features = "cyto", reduction = "umap", pt.size = 2.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(file = paste0(fig_path, "fig_1_2.png"), width = png_w*4, height = png_h*4, res = png_res)
p = plot_grid(plotlist=all_p, ncol = 4)
print(p)
dev.off()

##################
# Figure 1 - Alt #
##################
# Using tj_jaw combo instead of separate
test <- RunCCA(tj, jaw, renormalize = TRUE, rescale = TRUE)
tj_jaw <- FindVariableFeatures(test)
tj_jaw <- RunPCA(tj_jaw, npcs = 30, verbose = FALSE)
tj_jaw <- RunUMAP(tj_jaw, reduction = "pca", dims = 1:12) 
tj_jaw <- FindNeighbors(tj_jaw, reduction = "umap", dims = 1:2)
tj_jaw <- FindClusters(tj_jaw, resolution = 0.25)
DimPlot(tj_jaw, reduction = "umap", split.by = "cond", label = TRUE)

tj_jaw_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/tj_jaw_cyto.rds")
tj_jaw$cyto <- tj_jaw_cyto$CytoTRACE

text_size = 25

# Row 1 - Clusters
all_p <- list()
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_clusters_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = DimPlot(ob.list[[x]], label = TRUE, pt.size = 1.5) + labs(title = paste("Mouse Mesenchyme -", cond.list[[x]])) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'))
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_clusters.png"), width = png_w, height = png_h, res = png_res)
p <- DimPlot(mes, label = TRUE, pt.size = 1.5, label.size = lbl_size) + labs(title = paste("Mouse Mesenchyme")) + theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_clusters_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p <- DimPlot(ob.list[[x]], label = TRUE, pt.size = 1.5) + labs(title = paste("Mouse Epithelium -", cond.list[[x]])) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'))
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_clusters.png"), width = png_w, height = png_h, res = png_res)
p <- DimPlot(epi, label = TRUE, pt.size = 1.5, label.size = lbl_size) + labs(title = paste("Mouse Epithelium")) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_jaw_clusters.png"), width = png_w, height = png_h, res = png_res)
p <- DimPlot(tj_jaw, label = TRUE, pt.size = 2.5, label.size = lbl_size) + labs(title = paste("Cichlid Tooth")) +  theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# Row 2 - Yap
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_yap_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], feature = "Yap1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_yap.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(mes, feature = "Yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_yap_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p <- FeaturePlot(ob.list[[x]], features = "Yap1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_yap.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(epi, feature = "Yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_jaw_yap.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(tj_jaw, features = "yap1", label = TRUE, pt.size = 2.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# Row 3 - Piezo
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_piezo_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], feature = "Piezo1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_piezo.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(mes, feature = "Piezo1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_piezo_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p <- FeaturePlot(ob.list[[x]], features = "Piezo1", label = TRUE, pt.size = 1.5, order = T) + theme(plot.title = element_blank())
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_piezo.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(epi, feature = "Piezo1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_jaw_piezo.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(tj_jaw, features = "piezo1", label = TRUE, pt.size = 2.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# Row 4 - CytoTRACE
# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Clipped"
# ob.list <- SplitObject(mes, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "mes_cyto_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank())
#   all_p[[length(all_p)+1]] <- p
#   print(p)
#   dev.off()
# }
png(paste0(scratch, "mes_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(mes, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

# cond.list <- list()
# cond.list[[1]] <- "Control"
# cond.list[[2]] <- "Injury"
# ob.list <- SplitObject(epi, split.by = "orig.ident")
# for (x in 1:length(ob.list)) {
#   png(paste0(scratch, "epi_cyto_", cond.list[[x]], ".png"), width = png_w, height = png_h, res = png_res)
#   p = FeaturePlot(ob.list[[x]], features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank())
#   print(p)
#   all_p[[length(all_p)+1]] <- p
#   dev.off()
# }
png(paste0(scratch, "epi_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(epi, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_jaw_cyto.png"), width = png_w, height = png_h, res = png_res)
p = FeaturePlot(tj_jaw, features = "cyto", reduction = "umap", pt.size = 2.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(file = paste0(fig_path, "fig_1_3.png"), width = png_w*4, height = png_h*4, res = png_res)
p = plot_grid(plotlist=all_p, ncol = 3)
print(p)
dev.off()

############
# Figure 2 #
############
# maybe try a figure for YAP like you made for celsr 1, 
# using only YAP1 cells and looking at expression and 
# cytotrace side by side across the data sets

# Find yap positive cells
mes_yap_pos_cells <- names(which(mes@assays$RNA@data["Yap1",] > 0))
epi_yap_pos_cells <- names(which(epi@assays$RNA@data["Yap1",] > 0))
tj_yap_pos_cells  <- names(which(tj@assays$RNA@data["yap1",] > 0))
jaw_yap_pos_cells <- names(which(jaw@assays$RNA@data["yap1",] > 0))

# Subset Seurat objects by yap positive cells
mes_yap <- subset(mes, cells = mes_yap_pos_cells)
epi_yap <- subset(epi, cells = epi_yap_pos_cells)
tj_yap  <- subset(tj,  cells = tj_yap_pos_cells)
jaw_yap <- subset(jaw, cells = jaw_yap_pos_cells)

# Recluster the objects
mes_yap$orig.clust <- mes$seurat_clusters[mes_yap_pos_cells]
mes_yap <- FindVariableFeatures(object = mes_yap, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
mes_yap <- ScaleData(object = mes_yap, vars.to.regress = NULL)
mes_yap <- RunPCA(mes_yap, npcs = 30, verbose = FALSE)
mes_yap <- RunUMAP(mes_yap, reduction = "pca", dims = 1:12)
mes_yap <- FindNeighbors(mes_yap, reduction = "umap", dims = 1:2)
mes_yap <- FindClusters(mes_yap, resolution = 0.30)

epi_yap$orig.clust <- epi$seurat_clusters[epi_yap_pos_cells]
epi_yap <- FindVariableFeatures(object = epi_yap, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
epi_yap <- ScaleData(object = epi_yap, vars.to.regress = NULL)
epi_yap <- RunPCA(epi_yap, npcs = 30, verbose = FALSE)
epi_yap <- RunUMAP(epi_yap, reduction = "pca", dims = 1:12)
epi_yap <- FindNeighbors(epi_yap, reduction = "umap", dims = 1:2)
epi_yap <- FindClusters(epi_yap, resolution = 0.30)

tj_yap$orig.clust <- tj$seurat_clusters[tj_yap_pos_cells]
tj_yap <- FindVariableFeatures(object = tj_yap, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
tj_yap <- ScaleData(object = tj_yap, vars.to.regress = NULL)
tj_yap <- RunPCA(tj_yap, npcs = 30, verbose = FALSE)
tj_yap <- RunUMAP(tj_yap, reduction = "pca", dims = 1:12)
tj_yap <- FindNeighbors(tj_yap, reduction = "umap", dims = 1:2)
tj_yap <- FindClusters(tj_yap, resolution = 0.30)

jaw_yap$orig.clust <- jaw$seurat_clusters[jaw_yap_pos_cells]
jaw_yap <- FindVariableFeatures(object = jaw_yap, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
jaw_yap <- ScaleData(object = jaw_yap, vars.to.regress = NULL)
jaw_yap <- RunPCA(jaw_yap, npcs = 30, verbose = FALSE)
jaw_yap <- RunUMAP(jaw_yap, reduction = "pca", dims = 1:12)
jaw_yap <- FindNeighbors(jaw_yap, reduction = "umap", dims = 1:2)
jaw_yap <- FindClusters(jaw_yap, resolution = 0.30)

# Make the plots
all_p <- list()
text_size = 18

# Row 1 - Yap Expressions
p <- FeaturePlot(mes_yap, features = "Yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

p <- FeaturePlot(epi_yap, features = "Yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

p <- FeaturePlot(tj_yap, features = "yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

p <- FeaturePlot(jaw_yap, features = "yap1", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

# Row 2 - CytoTRACE
p = FeaturePlot(mes_yap, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

p = FeaturePlot(epi_yap, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

p = FeaturePlot(tj_yap, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

p = FeaturePlot(jaw_yap, features = "cyto", reduction = "umap", pt.size = 1.5, label=T, order = T, label.size = lbl_size) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + theme(plot.title = element_blank(), text=element_text(size=text_size))
all_p[[length(all_p)+1]] <- p

png(filename = paste0(fig_path, "fig_2_1.png"), width=600*4, height=png_h*2, res = png_res)
p = plot_grid(plotlist = all_p, ncol = 4)
print(p)
dev.off()

############
# Figure 3 #
############
# cartoon of cichlid jaw vs tooth samples, clusters for each dissection, 
# then paint keratin 5+15 in both samples
fig_3_w = 700
fig_3_h = 400
fig_3_res = 100
lbl_size = 6

tj_krt_5   <- names(which(tj@assays$RNA@data["krt5",] > 0))
jaw_krt_5  <- names(which(jaw@assays$RNA@data["krt5",] > 0))
tj_krt_15  <- names(which(tj@assays$RNA@data["krt15",] > 0))
jaw_krt_15 <- names(which(jaw@assays$RNA@data["krt15",] > 0))
tj_krt5_15  <- tj_krt_5[which(tj_krt_5 %in% tj_krt_15)]
jaw_krt5_15 <- jaw_krt_5[which(jaw_krt_5 %in% jaw_krt_15)]

tj$krt5_15  <- tj@assays$RNA@data["krt5",]  + tj@assays$RNA@data["krt15",]
jaw$krt5_15 <- jaw@assays$RNA@data["krt5",] + jaw@assays$RNA@data["krt15",]
tj$krt5_15[which(! names(tj$krt5_15) %in% tj_krt5_15)] = 0
jaw$krt5_15[which(! names(jaw$krt5_15) %in% jaw_krt5_15)] = 0

all_p <- list()
png(paste0(scratch, "tj_clusters_fig_3.png"), width = fig_3_w, height = fig_3_h, res = fig_3_res)
p <- DimPlot(tj, label = TRUE, pt.size = 1.5, label.size = lbl_size) + theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "jaw_clusters_fig_3.png"), width = fig_3_w, height = fig_3_h, res = fig_3_res)
p <- DimPlot(jaw, label = TRUE, pt.size = 1.5, label.size = lbl_size) + theme(plot.title = element_text(hjust = 0.5, face = 'plain'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "tj_krt5_krt15.png"), width = fig_3_w, height = fig_3_h, res = fig_3_res)
p <- FeaturePlot(tj, features = "krt5_15", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank())
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(paste0(scratch, "jaw_krt5_krt15.png"), width = png_w, height = png_h, res = png_res)
p <- FeaturePlot(jaw, features = "krt5_15", label = TRUE, pt.size = 1.5, order = T, label.size = lbl_size) + theme(plot.title = element_blank())
print(p)
all_p[[length(all_p)+1]] <- p
dev.off()

png(file = paste0(fig_path, "fig_3_1.png"), width = fig_3_w*2, height = fig_3_h*2, res = 120)
p = plot_grid(plotlist=all_p)
print(p)
dev.off()
