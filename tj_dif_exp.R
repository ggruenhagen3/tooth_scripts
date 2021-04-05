library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("ggplot2")

rna_path <- "C:/Users/miles/Downloads/d_tooth/"

tj.data <- b1.data <- Read10X(data.dir = paste(rna_path, "data/TJ/outs/filtered_feature_bc_matrix/", sep=""))

tj <- CreateSeuratObject(counts = tj.data, project = "TJ")
tj$cond <- "TJ"
tj$sample <- rep("TJ", ncol(tj))
tj <- NormalizeData(tj, normalization.method = "LogNormalize", scale.factor = 100000)

FeatureScatter(tj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
tj <- subset(tj, subset = nFeature_RNA > 125) # changed from default in order to not remove too many cells
tj <- subset(tj, subset = nFeature_RNA < 2500)
tj <- ScaleData(object = tj, vars.to.regress = NULL)
tj <- FindVariableFeatures(object = tj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
tj <- RunPCA(tj, npcs = 30, verbose = FALSE, features = tj@assays$RNA@var.features)
# JackStraw Values: 5.969996e-102  3.050209e-41  3.754158e-18  2.218550e-18  1.471733e-03  1.482937e-11 1.000000e+00  5.555104e-20  7.345818e-02  2.480356e-01  1.000000e+00  4.107480e-02 1.000000e+00  1.334201e-01  1.334201e-01  2.322131e-02  1.804553e-07  1.323612e-02 4.793902e-01  4.793902e-01
# tj <- JackStraw(tj, num.replicate = 100)
# tj <- ScoreJackStraw(tj, dims = 1:20)
# JackStrawPlot(tj, dims = 1:20)

tj <- RunUMAP(tj, reduction = "pca", dims = 1:12) 
tj <- FindNeighbors(tj, reduction = "umap", dims = 1:2)
tj <- FindClusters(tj, resolution = 0.25)
DimPlot(tj, reduction = "umap", split.by = "cond", label = TRUE)
# saveRDS(tj, "C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")

# Label Cell Types
# Idents(tj) <- "seurat_clusters"
# new.cluster.ids <- c("Stem Progenitor", "Stem Progenitor", "Stem Progenitor", "MSC", "MSC", "TAC", "TAC","MSC", "TAC", "Immune", "Immune", "TAC", "Immune", "Immune", "Immune", "Immune", "Neural Crest", "Immune", "Immune", "Immune", "Epithelial")
# names(new.cluster.ids) <- levels(tj)
# tj <- RenameIdents(tj, new.cluster.ids)
# png(filename = paste(rna_path, "results/", "tj.png", sep=""), width = 2400, height = 1800, unit="px", res=300)
# DimPlot(tj, reduction = "umap", cols = c("pink", "cadetblue1", "chartreuse", "gold", "purple", "orange"))
# dev.off()


# DEG Across Clusters
rna_path <- "C:/Users/miles/Downloads/d_tooth/"
tj.nk.markers <- FindAllMarkers(tj)
tj.nk.markers <- tj.nk.markers[,c(ncol(tj.nk.markers)-1, ncol(tj.nk.markers), 1:(ncol(tj.nk.markers)-2)),]
tj.nk.markers <- tj.nk.markers[which(tj.nk.markers$p_val_adj < 0.05),]
write.table(tj.nk.markers, file = paste(rna_path, "/results/tj_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

## DEG Against Jaw
# First merge the objects
tj$project  <- "TJ"
jaw$project <- "JAW"
Idents(tj)  <- tj$project
Idents(jaw) <- jaw$project
tj_jaw <- merge(tj, jaw, merge.data = FALSE)
tj_jaw <- NormalizeData(tj_jaw, normalization.method = "LogNormalize", scale.factor = 100000)
tj_jaw$project <- Idents(tj_jaw)
tj_jaw <- ScaleData(object = tj_jaw, vars.to.regress = NULL)
tj_jaw <- FindVariableFeatures(object = tj_jaw, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
tj_jaw <- RunPCA(tj_jaw, npcs = 30, verbose = FALSE, features = tj_jaw@assays$RNA@var.features)
tj_jaw <- RunUMAP(tj_jaw, reduction = "pca", dims = 1:12) 
tj_jaw <- FindNeighbors(tj_jaw, reduction = "umap", dims = 1:2)
tj_jaw <- FindClusters(tj_jaw, resolution = 0.25)
DimPlot(tj_jaw, reduction = "umap", split.by = "cond", label = TRUE)
Idents(tj_jaw) <- "cond"
VlnPlot(tj_jaw, features = "fgfr2", split.by = "cond", slot = "data")

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
Idents(tj_jaw) <- tj_jaw$project
tj.markers  <- FindMarkers(tj_jaw, ident.1 = "TJ", ident.2 = "JAW", only.pos = TRUE)
jaw.markers <- FindMarkers(tj_jaw, ident.1 = "JAW", ident.2 = "TJ", only.pos = TRUE)
tj.markers$gene  <- rownames(tj.markers)
jaw.markers$gene <- rownames(jaw.markers)
tj.markers  <- tj.markers[,c(ncol(tj.markers), 1:(ncol(tj.markers)-1))]
jaw.markers <- jaw.markers[,c(ncol(jaw.markers), 1:(ncol(jaw.markers)-1))]
write.table(tj.markers, file = paste(rna_path, "/results/tj_vs_jaw_tj_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)
write.table(jaw.markers, file = paste(rna_path, "/results/tj_vs_jaw_jaw_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

num_jaw_clusters <- as.numeric(tail(levels(jaw@meta.data$seurat_clusters), n=1))
num_tj_clusters  <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
Idents(tj)  <- tj$seurat_clusters
Idents(jaw) <- jaw$seurat_clusters
for (jaw_clust in 0:num_jaw_clusters) {
  tj_jaw <- SetIdent(tj_jaw, cells=WhichCells(jaw, idents = jaw_clust), value=paste("jaw", jaw_clust, sep="_"))
}
for (tj_clust in 0:num_tj_clusters) {
  tj_jaw <- SetIdent(tj_jaw, cells=WhichCells(tj, idents = tj_clust), value=paste("tj", tj_clust, sep="_"))
}
DimPlot(tj_jaw, reduction = "umap", split.by = "project", label = TRUE)

#####################
# Cell Trajectories #
#####################
# .libPaths(c("C:/Users/miles/Documents/my_R_packages", "C:/Users/miles/Documents/R/win-library/3.6", "C:/Program Files/R/R-3.6.1/library"))
tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
jpool <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/pool_shiny/data/jpool.rds")
jaw <- jpool
library("slingshot")
library("SingleCellExperiment")
library("RColorBrewer")
library("ggplot2")
tj_jaw_sce <- as.SingleCellExperiment(tj_jaw, assay = "RNA")
sds <- slingshot(tj_jaw_sce, clusterLabels = tj_jaw$seurat_clusters, reducedDim = "UMAP", start.clus = 7)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
plotcol2 <- colors[cut(1:ncol(tj_jaw), breaks=100)]

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
png(paste(rna_path, "results/trajectory/tj_jaw_slingshot_pseduotime_start_7.png", sep=""), width = 3600, height=2400, unit = "px", res = 300)
sds_ds <- SlingshotDataSet(sds)
plot(reducedDims(sds)$UMAP, col = brewer.pal(11,'Paired')[sds$seurat_clusters], pch=16, asp = 1)
legend("bottomleft", legend = 0:10, fill = brewer.pal(11,'Paired'), bty = "n", cex = 2)
for (curve in sds_ds@curves) {
  x <- curve$s[curve$ord,1]
  y <- curve$s[curve$ord,2]
  segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=plotcol2, lwd=3)
}
dev.off()

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
png(paste(rna_path, "results/trajectory/tj_jaw_slingshot_dot_start_7.png", sep=""), width = 3600, height=2400, unit = "px", res = 300)
plot(reducedDims(sds)$UMAP, col = brewer.pal(11,'Paired')[sds$seurat_clusters], pch=16, asp = 1)
legend("bottomleft", legend = 0:12, fill = brewer.pal(11,'Paired'), bty = "n", cex = 1.8)
lines(SlingshotDataSet(sds), lwd=3, type = 'lineages', col = 'black')
forest <- slingAdjacency(sds)
clusters <- rownames(forest)
nclus <- nrow(forest)
x <- reducedDim(sds_ds)
centers <- t(vapply(clusters,function(clID){
  w <- slingClusterLabels(sds_ds)[,clID]
  return(apply(x, 2, weighted.mean, w = w))
}, rep(0,ncol(reducedDim(sds_ds)))))
rownames(centers) <- clusters
for (lineage in sds_ds@lineages) {
  for (i in 1:(length(lineage)-1)) {
    cluster_start <- as.character(lineage[i])
    cluster_end <- as.character(lineage[i+1])
    arrows(centers[cluster_start,1],centers[cluster_start,2], centers[cluster_end,1], centers[cluster_end,2], lwd = 3)
  }
}
dev.off()

# Tooth Cell Trajectories
tj_sce <- as.SingleCellExperiment(tj, assay = "RNA")
sds <- slingshot(tj_sce, clusterLabels = tj$seurat_clusters, reducedDim = "UMAP", start.clus=6)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
plotcol[is.na(plotcol)] <- "grey"
plotcol2 <- colors[cut(1:739, breaks=100)]

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
png(paste(rna_path, "results/trajectory/tj_slingshot_pseduotime_start_6.png", sep=""), width = 3600, height=2400, unit = "px", res = 300)
sds_ds <- SlingshotDataSet(sds)
plot(reducedDims(sds)$UMAP, col = brewer.pal(8,'Paired')[sds$seurat_clusters], pch=16, asp = 1)
legend("bottomleft", legend = 0:7, fill = brewer.pal(8,'Paired'), bty = "n", cex = 2)
for (curve in sds_ds@curves) {
  x <- curve$s[curve$ord,1]
  y <- curve$s[curve$ord,2]
  segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=plotcol2, lwd=3)
}
dev.off()

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
png(paste(rna_path, "results/trajectory/tj_slingshot_dot_start_6.png", sep=""), width = 3600, height=2400, unit = "px", res = 300)
plot(reducedDims(sds)$UMAP, col = brewer.pal(8,'Paired')[sds$seurat_clusters], pch=16, asp = 1)
legend("bottomleft", legend = 0:7, fill = brewer.pal(8,'Paired'), bty = "n", cex = 1.8)
lines(SlingshotDataSet(sds), lwd=3, type = 'lineages', col = 'black')
forest <- slingAdjacency(sds)
clusters <- rownames(forest)
nclus <- nrow(forest)
x <- reducedDim(sds_ds)
centers <- t(vapply(clusters,function(clID){
  w <- slingClusterLabels(sds_ds)[,clID]
  return(apply(x, 2, weighted.mean, w = w))
}, rep(0,ncol(reducedDim(sds_ds)))))
rownames(centers) <- clusters
for (lineage in sds_ds@lineages) {
  for (i in 1:(length(lineage)-1)) {
    cluster_start <- as.character(lineage[i])
    cluster_end <- as.character(lineage[i+1])
    arrows(centers[cluster_start,1],centers[cluster_start,2], centers[cluster_end,1], centers[cluster_end,2], lwd = 3)
  }
}
dev.off()

# Jaw Cell Trajectories
jaw_sce <- as.SingleCellExperiment(jaw, assay = "RNA")
sds <- slingshot(jaw_sce, clusterLabels = jaw$seurat_clusters, reducedDim = "UMAP", start.clus = 6)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
plotcol[is.na(plotcol)] <- "grey"
plotcol2 <- colors[cut(1:ncol(jaw), breaks=100)]

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
png(paste(rna_path, "results/trajectory/jaw_slingshot_pseduotime_start_6.png", sep=""), width = 3600, height=2400, unit = "px", res = 300)
sds_ds <- SlingshotDataSet(sds)
plot(reducedDims(sds)$UMAP, col = brewer.pal(12,'Paired')[sds$seurat_clusters], pch=16, asp = 1)
legend("bottomleft", legend = 0:11, fill = brewer.pal(12,'Paired'), bty = "n", cex = 2)
for (curve in sds_ds@curves) {
  x <- curve$s[curve$ord,1]
  y <- curve$s[curve$ord,2]
  segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=plotcol2, lwd=3)
}
dev.off()

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
png(paste(rna_path, "results/trajectory/jaw_slingshot_dot_start_6.png", sep=""), width = 3600, height=2400, unit = "px", res = 300)
plot(reducedDims(sds)$UMAP, col = brewer.pal(11,'Paired')[sds$seurat_clusters], pch=16, asp = 1)
legend("bottomleft", legend = 0:11, fill = brewer.pal(12,'Paired'), bty = "n", cex = 1.8)
lines(SlingshotDataSet(sds), lwd=3, type = 'lineages', col = 'black')
x <- reducedDim(sds_ds)
forest <- slingAdjacency(sds)
clusters <- rownames(forest)
nclus <- nrow(forest)
centers <- t(vapply(clusters,function(clID){
  w <- slingClusterLabels(sds_ds)[,clID]
  return(apply(x, 2, weighted.mean, w = w))
}, rep(0,ncol(reducedDim(sds_ds)))))
rownames(centers) <- clusters
for (lineage in sds_ds@lineages) {
  for (i in 1:(length(lineage)-1)) {
    cluster_start <- as.character(lineage[i])
    cluster_end <- as.character(lineage[i+1])
    arrows(centers[cluster_start,1],centers[cluster_start,2], centers[cluster_end,1], centers[cluster_end,2], lwd = 3)
  }
}
dev.off()

# Painting
rna_path <- "C:/Users/miles/Downloads/d_tooth/"
gene_names <- rownames(tj)
average <- TRUE
Idents(tj) <- tj$seurat_clusters
num_clusters <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
poster_genes <- c(gene_names[which(startsWith(gene_names, "bmp"))], gene_names[which(startsWith(gene_names, "fgf"))], gene_names[which(startsWith(gene_names, "wnt"))], gene_names[which(startsWith(gene_names, "tnf"))])
col <- c(gene_names[which(startsWith(gene_names, "col"))])
mesenchyme <- c("serpinf1","pcolcea","pcolceb","msx1a","col1a1a","col1a1b","MFAP4","ctsk","NBL1","cxcl12a","cxcl12b")
clust_3 <- c("shha", "bmp4", "edar", "pitx2", "ENSMZEG00005009006")
immune <- c("ccr7", "ENSMZEG00005018861", "ptpra", "ptprc", "PTPRB (1 of many)", "PTPRB (1 of many).1", "ptprfa", "ptprfb", "fcer1g", "c1qb")
paul <- c("edar", "shha", "celsr1a")
cur_list <- adrenergic_mz
for (gene in cur_list) {
  print(gene)
  
  total_tj <- c()
  expr1 <- FetchData(object = tj, vars = gene)
  pos_cells <- tryCatch({
    colnames(tj[, which(x = expr1 > 1)])
  }, error = function(e) {
    c()
  })
  for (i in 0:num_clusters) {
    tj_cells_in_cluster <- 0
    this_cells <- WhichCells(tj, idents = i)
    try(tj_cells_in_cluster <- length(pos_cells[pos_cells %in% this_cells]), silent=TRUE)
    all_tj_cells_in_cluster <- length(this_cells)
    
    if (average == TRUE) {
      total_tj <- c(total_tj, round( (tj_cells_in_cluster/all_tj_cells_in_cluster) * 100, 2))
    } else {
      total_tj <- c(total_tj, tj_cells_in_cluster)
    }
  }
  df <- data.frame(condition <- c(rep("TJ", length(total_tj))),cluster_num <- c(0:num_clusters),
                   value <- c(total_tj))
  colnames(df) <- c("condition", "cluster_num", "value")
  my_title <- paste("Number of Cells Expressing", paste(gene, collapse = ' and '), "per Cluster")
  my_ylab <- "Number of Cells"
  if (average == TRUE) {
    my_title <- paste("% Cells Expressing", paste(gene, collapse = ' and '), "per Cluster")
    my_ylab <- "% Cells"
  }
  p2 <- ggplot(df, aes(fill=condition, x=cluster_num, y=value)) +
    geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    ggtitle(my_title) +
    xlab("Cluster") +
    ylab(my_ylab) +
    scale_x_continuous(breaks = 0:40) +
    geom_text(aes(label=value), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), plot.title=element_text(size=18,face="bold"))
  theme_minimal()
  
  png(paste(rna_path, "/results/painting/tj/tj_", gene, "_umap.png", sep=""), width = 800, height = 750, unit = "px")
  p1 <- FeaturePlot(tj, features = gene, reduction = "umap", split.by = "cond", pt.size = 2, label=TRUE, order = TRUE)
  print(plot_grid(p1,p2,ncol=1))
  dev.off()
}

# Regress out Cell Cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- hgncGood(s.genes, rownames(tj@assays$RNA@counts))
g2m.genes <- hgncGood(g2m.genes, rownames(tj@assays$RNA@counts))
tj <- CellCycleScoring(tj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
tj <- RunPCA(tj, features = c(s.genes, g2m.genes))
tj <- ScaleData(tj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(tj))
tj <- RunPCA(tj, features = VariableFeatures(tj), nfeatures.print = 10)
tj <- RunPCA(tj, features = c(s.genes, g2m.genes))
DimPlot(tj)
tj <- RunUMAP(tj, reduction = "pca", dims = 1:30)
tj <- FindNeighbors(tj, reduction = "umap", dims = 1:2)
tj <- FindClusters(tj, resolution = 0.1)
DimPlot(tj, reduction = "umap", split.by = "cond", label = TRUE)

adrenergic_mz <- c("adra2b", "adra1d", "adra2c", "adrb3a", "adra1aa", "adra1ab", "adra1bb", "adrb2a", "adra2a")
val_adr_mz <- c("adrb3a", "src", "pnpla2", "ABHD5",  "srebf1", "adcy2a", "adcy2b", "adcy7",  "snap23.1", "snap29", "stx3a", "vamp1", "vamp2",  "VAMP3",  "vamp8")
for (gene in val_adr) {
  print(gene)
  print(grep(gene, rownames(tj), value = TRUE, ignore.case = TRUE))
}

# Exploring Transcriptional Diversity and Expression Level in CytoBINs
mat = tj@assays$RNA@counts
mat[which(mat > 0)] = 1
tj$div = colSums(mat)
df <- data.frame(names(tj$div), tj$div, tj$nCount_RNA)
ggplot(df, aes(tj.div, tj.nCount_RNA)) + geom_point(alpha = 0.1)
df <- data.frame(names(tj$div), tj$div, tj$nCount_RNA, tj$cyto)
df <- df[order(df$tj.cyto, decreasing = T),]
ggplot(df, aes(tj.div, tj.nCount_RNA, color = tj.cyto)) + geom_point(alpha = 0.5) + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1)) + xlab("Number of Genes Expressed (Transcription Diversity)") + ylab("Number of Transcripts (Expression Level)") + labs(title="Transcriptional Diversity and Expression Level in Cichlid Tooth by Cell with CytoTRACE Score")

tj$bin <- tj$cyto
tj$bin[which(tj$cyto <= 0.33)] = "low"
tj$bin[which(tj$cyto > 0.33 & tj$cyto <= 0.66)] = "medium"
tj$bin[which(tj$cyto > 0.66)] = "high"
df <- data.frame(names(tj$div), tj$div, tj$nCount_RNA, tj$cyto, tj$bin)
ggplot(df, aes(tj.bin, tj.div, fill = tj.bin, color = tj.bin)) + geom_boxplot(alpha=0.3) + geom_point(position=position_jitterdodge(jitter.width = 1.8), alpha = 0.25) + NoLegend() + xlab("CytoBIN") + ylab("Transcriptional Diversity per Cell") + labs(title = "Transcriptional Diversity in CytoBINs")
ggplot(df, aes(tj.bin, tj.nCount_RNA, fill = tj.bin, color = tj.bin)) + geom_boxplot(alpha=0.3) + geom_point(position=position_jitterdodge(jitter.width = 1.8), alpha = 0.25) + NoLegend() + xlab("CytoBIN") + ylab("Expression Level per Cell") + labs(title = "Expression Level in CytoBINs")

jaw$bin <- jaw$cyto
jaw$bin[which(jaw$cyto <= 0.33)] = "low"
jaw$bin[which(jaw$cyto > 0.33 & jaw$cyto <= 0.66)] = "medium"
jaw$bin[which(jaw$cyto > 0.66)] = "high"

customPaint = function(obj, gene, obj_name, average=T) {
  rna_path = "C:/Users/miles/Downloads/d_tooth/"
  total_tj <- c()
  num_clusters = max(unique(as.numeric(as.vector(Idents(obj)))))
  expr1 <- FetchData(object = obj, vars = gene, slot="counts")
  pos_cells <- tryCatch({
    colnames(obj[, which(x = expr1 > 0)])
  }, error = function(e) {
    c()
  })
  for (i in 0:num_clusters) {
    tj_cells_in_cluster <- 0
    this_cells <- WhichCells(obj, idents = i)
    try(tj_cells_in_cluster <- length(pos_cells[pos_cells %in% this_cells]), silent=TRUE)
    all_tj_cells_in_cluster <- length(this_cells)
    
    if (average == TRUE) {
      total_tj <- c(total_tj, round( (tj_cells_in_cluster/all_tj_cells_in_cluster) * 100, 2))
    } else {
      total_tj <- c(total_tj, tj_cells_in_cluster)
    }
  }
  df <- data.frame(condition <- c(rep("TJ", length(total_tj))),cluster_num <- c(0:num_clusters),
                   value <- c(total_tj))
  colnames(df) <- c("condition", "cluster_num", "value")
  my_title <- paste("Number of Cells Expressing", paste(gene, collapse = ' and '), "per Cluster")
  my_ylab <- "Number of Cells"
  if (average == TRUE) {
    my_title <- paste("% Cells Expressing", paste(gene, collapse = ' and '), "per Cluster")
    my_ylab <- "% Cells"
  }
  p2 <- ggplot(df, aes(fill=condition, x=cluster_num, y=value)) +
    geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    ggtitle(my_title) +
    xlab("Cluster") +
    ylab(my_ylab) +
    scale_x_continuous(breaks = 0:40) +
    geom_text(aes(label=value), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), plot.title=element_text(size=18,face="bold")) + guides(fill = FALSE)
  theme_minimal()
  
  png(paste(rna_path, "/results/painting/", obj_name, "/", obj_name, "_", gene, "_umap.png", sep=""), width = 800, height = 750, unit = "px")
  p1 <- FeaturePlot(obj, features = gene, reduction = "umap", split.by = "cond", pt.size = 2, label=TRUE, order = TRUE)
  print(plot_grid(p1,p2,ncol=1))
  dev.off()
}

EnhancedVolcano(tj_deg, lab = "gene", x="avg_logFC", y="p_val_adj")

# tal = readxl::read_excel("C:/Users/miles/Downloads/Specific and enriched Dental Genes.xlsx")
# tal = tal[,which(! startsWith(colnames(tal), "..."))]
# tal = readxl::read_excel("C:/Users/miles/Downloads/d_tooth/results/incsr_unique_100.xlsx")
# tal = readxl::read_excel("C:/Users/miles/Downloads/d_tooth/results/im_unique_100.xlsx", skip = 1)
tal = readxl::read_excel("C:/Users/miles/Downloads/d_tooth/results/hm_unique_100.xlsx", skip = 1)

for (i in 1:ncol(tal)) {
  name = str_replace(colnames(tal)[i], "/", "_")
  name_ns = str_replace(name, " ", "_")
  print(name)
  
  cur_list = hgncGood(as.data.frame(tal[which(! is.na(tal[,i])), i])[,1], rownames(tj), as_df = T)
  write.table(cur_list, paste0("C:/Users/miles/Downloads/d_tooth/data/markers/igor/hm/", name_ns, ".tsv"), sep="\t", quote = F, row.names = F)
  cur_list = read.table(paste0("C:/Users/miles/Downloads/d_tooth/data/markers/igor/hm/", name_ns, ".tsv"), sep="\t", header=T)

  res = markersInDEGUpDown(tj_deg_all, cur_list$mzebra)
  res_pct = markersInDEGUpDown(tj_deg_all, cur_list$mzebra, pct=T)[[2]]
  dot_res = myDotPlot(tj, cur_list$mzebra)
  dot_plot = dot_res[[1]]
  dot_bin = dot_res[[2]]
  dot_bin_weight = dot_res[[3]]

  # UMAP plot of Expression of Markers per Cell
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_exp_per_cell.png"), width = 500, height = 400)
  print(markerExpPerCell(tj, cur_list$mzebra) + ggtitle(paste0("Expression of ", name, " Markers per Cell in Cichlid Tooth")))
  dev.off()

  # Boxplot of Genes per Cell per Cluster - Corrected
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_avg_genes_per_cell_per_cluster_correct.png"), width = 600, height = 400)
  print(markerExpPerCellPerCluster(tj, cur_list$mzebra, correct = T)[[1]] + ggtitle(paste0("Average Number of ", name, " Markers per Cell per Cluster in Cichlid Tooth - Corrected")))
  dev.off()

  # Boxplot of Expression of Genes per Cell per Cluster - Corrected
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_exp_per_cell_per_cluster_correct.png"), width = 600, height = 400)
  print(markerExpPerCellPerCluster(tj, cur_list$mzebra, correct = T, n_markers = F)[[1]] + ggtitle(paste0("Average Expression of ", name, " Markers per Cell per Cluster in Cichlid Tooth - Corrected")))
  dev.off()

  # Barplot of Total Number of Cells per Cluster Expressing a Marker - Corrected
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_num_cell_per_cluster_correct.png"), width = 600, height = 400)
  print(markerCellPerCluster(tj, cur_list$mzebra) + ggtitle(paste("Normalized Total Number of Cells per Expressing Expressing", name, "Markers in Cichlid Tooth")))
  dev.off()

  # Heatmap of Expression per Cluster
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_heatmap.png"), width = 1200, height = 1200, res=100)
  print(markerHeatmap(tj, cur_list$mzebra) + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
  dev.off()

  # DotPlot of Expression per Cluster
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_dotplot.png"), width = 800, height = 1600, res=100)
  print(dot_plot + ggtitle(paste0(name, " Markers in Cichlid Tooth")) + coord_flip())
  dev.off()

  # DotPlot BIN per Cluster
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_dotplot_bin.png"), width = 800, height = 400)
  print(dot_bin)
  dev.off()

  # DotPlot BIN per Cluster - Weighted
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_dotplot_bin_weight.png"), width = 800, height = 400)
  print(dot_bin_weight)
  dev.off()

  # DEGs in Markers
  if (! is.null(res[[2]])) {
    # Number of DEGs
    png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_deg.png"), width = 500, height = 400)
    print(res[[2]] + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
    dev.off()

    # LogFC
    png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_deg_logFC.png"), width = 500, height = 400)
    print(res[[3]] + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
    dev.off()

    # % DEGs in Markers
    png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_deg_pct.png"), width = 500, height = 400)
    print(res_pct + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
    dev.off()
  }

  # Write Genes Found
  found_genes = res[[1]]
  colnames(found_genes) = c("mzebra", "cluster", "avg_logFC", "isPos")
  found_genes = left_join(found_genes, cur_list, by = "mzebra")
  found_genes = found_genes[,c("hgnc", "mzebra", "cluster", "avg_logFC", "isPos")]
  write.table(found_genes, paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/hm/tj_", name_ns, "_deg.tsv"), sep="\t", quote = F, row.names = F)
}


for (i in 1:ncol(tal)) {
  name = str_replace(colnames(tal)[i], "/", "_")
  name_ns = str_replace(name, " ", "_")
  print(name)
  
  # cur_list = hgncGood(as.data.frame(tal[which(! is.na(tal[,i])), i])[,1], rownames(jaw), as_df = T)
  # write.table(cur_list, paste0("C:/Users/miles/Downloads/d_tooth/data/markers/igor/im/", name_ns, ".tsv"), sep="\t", quote = F, row.names = F)
  cur_list = read.table(paste0("C:/Users/miles/Downloads/d_tooth/data/markers/igor/im/", name_ns, ".tsv"), sep="\t", header=T)
  
  res = markersInDEGUpDown(jaw_deg_all, cur_list$mzebra)
  res_pct = markersInDEGUpDown(jaw_deg_all, cur_list$mzebra, pct=T)[[2]]
  dot_res = myDotPlot(jaw, cur_list$mzebra)
  dot_plot = dot_res[[1]]
  dot_bin = dot_res[[2]]
  dot_bin_weight = dot_res[[3]]

  # UMAP plot of Expression of Markers per Cell
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_exp_per_cell.png"), width = 500, height = 400)
  print(markerExpPerCell(jaw, cur_list$mzebra) + ggtitle(paste0("Expression of ", name, " Markers per Cell in Cichlid Jaw")))
  dev.off()

  # Boxplot of Genes per Cell per Cluster - Corrected
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_avg_genes_per_cell_per_cluster_correct.png"), width = 600, height = 400)
  print(markerExpPerCellPerCluster(jaw, cur_list$mzebra, correct = T)[[1]] + ggtitle(paste0("Average Number of ", name, " Markers per Cell per Cluster in Cichlid Jaw - Corrected")))
  dev.off()

  # Boxplot of Expression of Genes per Cell per Cluster - Corrected
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_exp_per_cell_per_cluster_correct.png"), width = 600, height = 400)
  print(markerExpPerCellPerCluster(jaw, cur_list$mzebra, correct = T, n_markers = F)[[1]] + ggtitle(paste0("Average Expression of ", name, " Markers per Cell per Cluster in Cichlid Jaw - Corrected")))
  dev.off()

  # Barplot of Total Number of Cells per Cluster Expressing a Marker - Corrected
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_num_cell_per_cluster_correct.png"), width = 600, height = 400)
  print(markerCellPerCluster(jaw, cur_list$mzebra) + ggtitle(paste("Normalized Total Number of Cells per Expressing Expressing", name, "Markers in Cichlid Jaw")))
  dev.off()

  # Heatmap of Expression per Cluster
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_heatmap.png"), width = 1200, height = 1200, res=100)
  print(markerHeatmap(jaw, cur_list$mzebra) + ggtitle(paste0(name, " Markers in Cichlid Jaw")))
  dev.off()

  # DotPlot of Expression per Cluster
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_dotplot.png"), width = 800, height = 1600, res=100)
  print(dot_plot + ggtitle(paste0(name, " Markers in Cichlid Jaw")) + coord_flip())
  dev.off()

  # DotPlot BIN per Cluster
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_dotplot_bin.png"), width = 800, height = 1600, res=100)
  print(dot_bin + ggtitle(paste0(name, " Markers in Cichlid Jaw")))
  dev.off()

  # DotPlot BIN per Cluster - Weighted
  png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_dotplot_bin_weight.png"), width = 800, height = 400)
  print(dot_bin_weight)
  dev.off()

  # DEGs in Markers
  if (! is.null(res[[2]])) {
    # Number of DEGs
    png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_deg.png"), width = 500, height = 400)
    print(res[[2]] + ggtitle(paste0(name, " Markers in Cichlid Jaw")))
    dev.off()

    # LogFC
    png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_deg_logFC.png"), width = 500, height = 400)
    print(res[[3]] + ggtitle(paste0(name, " Markers in Cichlid Jaw")))
    dev.off()

    # % DEGs in Markers
    png(paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_deg_pct.png"), width = 500, height = 400)
    print(res_pct + ggtitle(paste0(name, " Markers in Cichlid Jaw")))
    dev.off()
  }

  # Write Genes Found
  found_genes = res[[1]]
  colnames(found_genes) = c("mzebra", "cluster", "avg_logFC", "isPos")
  found_genes = left_join(found_genes, cur_list, by = "mzebra")
  found_genes = found_genes[,c("hgnc", "mzebra", "cluster", "avg_logFC", "isPos")]
  write.table(found_genes, paste0("C:/Users/miles/Downloads/d_tooth/results/igor_markers/jaw_", name_ns, "_deg.tsv"), sep="\t", quote = F, row.names = F)
}
