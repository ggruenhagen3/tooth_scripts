# Load Packages
library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("biomaRt")
library("stringr")
library("dplyr")
library("ggplot2")
library("qvalue")
library("jaccard")

igor_path = "c:/Users/miles/Downloads/d_tooth/data/igor/"

# Mouse Incisor #########################################################################

igor.mouse.data = read.table(paste0(igor_path, "counts_SS2_mouse_incisor.txt"))
igor.mouse <- CreateSeuratObject(counts = igor.mouse.data, project = "INSR")
igor.mouse$cond = "igor.mouse"
igor.mouse = NormalizeData(igor.mouse, normalization.method = "LogNormalize", scale.factor = 100000)
subset(igor.mouse, subset = nFeature_RNA > 500 & nFeature_RNA < 7500)
igor.mouse = FindVariableFeatures(object = igor.mouse, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
igor.mouse <- ScaleData(object = igor.mouse, vars.to.regress = NULL)
igor.mouse <- RunPCA(igor.mouse, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
igor.mouse <- RunUMAP(igor.mouse, reduction = "pca", dims = 1:20)
igor.mouse <- FindNeighbors(igor.mouse, reduction = "umap", dims = 1:2)
igor.mouse <- FindClusters(igor.mouse, resolution = 0.1)
DimPlot(igor.mouse, reduction="umap", label=T)

# Igor's annotations
igor.mouse.annot = read.table(paste0(igor_path, "annotation_mouse_incisor.txt"), sep="\t", header = F, stringsAsFactors = F)
igor.mouse$annot = igor.mouse.annot[,2]
Idents(igor.mouse) = igor.mouse$annot
DimPlot(igor.mouse, reduction="umap", label=T, pt.size = 1)

# Comparitive Mouse #####################################################################
# Igor's Annotation of Clusters
test = readRDS("C:/Users/miles/Downloads/incisor_molar_data.rds")
comp.annot = read.table(paste0(igor_path, "annotation_comparative_mouse.txt"), header = T, stringsAsFactors = F)
test1 = as.data.frame(test$UMAP)
test1$cluster = factor(comp.annot$id)
colnames(test1) = c("UMAP_1", "UMAP_2", "cluster")
p = ggplot(test1, aes(UMAP_1, UMAP_2, color=cluster)) + geom_point(size=1) + theme_classic() + xlab("UMAP_1") + ylab("UMAP_2")
LabelClusters(p, "cluster", levels(test1$cluster), colour="black", repel=F)

# My own Clustering
igor.comp.data = read.table(paste0(igor_path, "counts_comparative_mouse.txt"))
igor.comp = CreateSeuratObject(counts = igor.comp.data, project = "COMP")
igor.comp$cond = "igor.comp"
igor.comp = NormalizeData(igor.comp, normalization.method = "LogNormalize", scale.factor = 100000)
subset(igor.comp, subset = nFeature_RNA > 500 & nFeature_RNA < 7500)
igor.comp = FindVariableFeatures(object = igor.comp, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
igor.comp <- ScaleData(object = igor.comp, vars.to.regress = NULL)
igor.comp <- RunPCA(igor.comp, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
igor.comp <- RunUMAP(igor.comp, reduction = "pca", dims = 1:20)
igor.comp <- FindNeighbors(igor.comp, reduction = "umap", dims = 1:2)
igor.comp <- FindClusters(igor.comp, resolution = 0.1)
DimPlot(igor.comp, reduction="umap", label=T) + ggtitle("Incisor + Molar: My Clustering")
pulp_cells = as.character(as.vector(names(test$pulp.clusters)))
pulp_cells = gsub("-", ".", gsub("\\.", "", pulp_cells))
pulp_cells[which(! startsWith(pulp_cells, "one_"))] = paste0("X", pulp_cells[which(startsWith(pulp_cells, "5_one"))])
DimPlot(igor.comp, reduction="umap", label=T, cells.highlight = pulp_cells) + ggtitle("Incisor + Molar: My Clustering (Highlight Pulp Cells)")

# Putting Igor's Annotations and Clustering on my Seurat object
igor.comp@reductions$umap@cell.embeddings=as.matrix(test1[,1:2])
rownames(igor.comp@reductions$umap@cell.embeddings) = colnames(igor.comp)
igor.comp$seurat_clusters = comp.annot$id
Idents(igor.comp) = factor(igor.comp$seurat_clusters)
DimPlot(igor.comp, reduction="umap", label=T) + ggtitle("Incisor + Molar: Igor's Clustering")
DimPlot(igor.comp, reduction="umap", label=T,  cells.highlight = pulp_cells) + ggtitle("Incisor + Molar: Igor's Clustering (Highlight Pulp Cells)")
# Human Molar ###########################################################################


################################
# Cluster Igor Mouse w/ Mzebra #
################################
tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
jpool <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/pool_shiny/data/jpool.rds")
jaw <- jpool

tj.mouse  = convertToMouseObj(tj)
tj.mouse[[paste0("NormalizeData.RNA")]] = tj[[paste0("NormalizeData.RNA")]]
tj.mouse = FindVariableFeatures(tj.mouse)

jaw.mouse = convertToMouseObj(jaw)
jaw.mouse[[paste0("NormalizeData.RNA")]] = tj[[paste0("NormalizeData.RNA")]]
jaw.mouse = FindVariableFeatures(jaw.mouse)

# Tooth Merge
igor.mouse.tj = merge(tj.mouse, igor.mouse)
# igor.mouse.tj <- RunCCA(tj.mouse, igor.mouse, renormalize = TRUE, rescale = TRUE)
igor.mouse.tj <- FindVariableFeatures(igor.mouse.tj)
igor.mouse.tj <- ScaleData(object = igor.mouse.tj, vars.to.regress = NULL)
igor.mouse.tj <- RunPCA(igor.mouse.tj, npcs = 30, verbose = FALSE)
igor.mouse.tj <- RunUMAP(igor.mouse.tj, reduction = "pca", dims = 1:20) 
igor.mouse.tj <- FindNeighbors(igor.mouse.tj, reduction = "umap", dims = 1:2)
igor.mouse.tj <- FindClusters(igor.mouse.tj, resolution = 0.1)
DimPlot(igor.mouse.tj, reduction = "umap", split.by = "cond", label = TRUE)

# Tooth Merge - IntegrateData
anchors <- FindIntegrationAnchors(object.list = list(tj.mouse, igor.mouse), reference = 2, dims = 1:30)
igor.mouse.tj <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(igor.mouse.tj) <- "integrated"
igor.mouse.tj <- ScaleData(object = igor.mouse.tj, vars.to.regress = NULL)
igor.mouse.tj <- RunPCA(igor.mouse.tj, npcs = 30, verbose = FALSE)
igor.mouse.tj <- RunUMAP(igor.mouse.tj, reduction = "pca", dims = 1:20) 
igor.mouse.tj <- FindNeighbors(igor.mouse.tj, reduction = "umap", dims = 1:2)
igor.mouse.tj <- FindClusters(igor.mouse.tj, resolution = 0.1)
DimPlot(igor.mouse.tj, reduction = "umap", split.by = "cond", label = TRUE, pt.size = 1)

# Label Tooth Merge by Original Cluster
igor.mouse.tj$annot = as.vector(igor.mouse.tj$seurat_clusters)
igor.mouse.tj$annot[which(names(igor.mouse.tj$annot) %in% colnames(igor.mouse))] = igor.mouse$annot
igor.mouse.tj$annot[which(names(igor.mouse.tj$annot) %in% colnames(tj))] = paste0("tj_", tj$seurat_clusters)
Idents(igor.mouse.tj) = igor.mouse.tj$annot
DimPlot(igor.mouse.tj, reduction = "umap", split.by = "cond", label = TRUE) + ggtitle("Igor Mouse Incisor Clustered with Cichlid Tooth - Simple Merge")

# Jaw Merge
# igor.mouse.jaw = merge(jaw.mouse, igor.mouse)
igor.mouse.jaw <- RunCCA(jaw.mouse, igor.mouse, renormalize = TRUE, rescale = TRUE)
igor.mouse.jaw <- FindVariableFeatures(igor.mouse.jaw)
# igor.mouse.jaw <- ScaleData(object = igor.mouse.jaw, vars.to.regress = NULL)
igor.mouse.jaw <- RunPCA(igor.mouse.jaw, npcs = 30, verbose = FALSE)
igor.mouse.jaw <- RunUMAP(igor.mouse.jaw, reduction = "pca", dims = 1:20) 
igor.mouse.jaw <- FindNeighbors(igor.mouse.jaw, reduction = "umap", dims = 1:2)
igor.mouse.jaw <- FindClusters(igor.mouse.jaw, resolution = 0.1)
DimPlot(igor.mouse.jaw, reduction = "umap", split.by = "cond", label = TRUE)

# Label Jaw Merge by Original Cluster
igor.mouse.jaw$annot = as.vector(igor.mouse.jaw$seurat_clusters)
igor.mouse.jaw$annot[which(names(igor.mouse.jaw$annot) %in% colnames(igor.mouse))] = igor.mouse$annot
igor.mouse.jaw$annot[which(names(igor.mouse.jaw$annot) %in% colnames(jaw))] = paste0("jaw_", jaw$seurat_clusters)
Idents(igor.mouse.jaw) = igor.mouse.jaw$annot
DimPlot(igor.mouse.jaw, reduction = "umap", split.by = "cond", label = TRUE) + ggtitle("Igor Mouse Incisor Clustered with Cichlid Jaw - IntegrateData")

#######################
# Similarity b/w DEGs #
#######################
# Only Positive
Idents(igor.mouse) = igor.mouse$annot
igor.mouse.deg = FindAllMarkers(igor.mouse, only.pos = T)
igor.mouse.deg = igor.mouse.deg[which(igor.mouse.deg$p_val_adj < 0.05),]

Idents(igor.comp) = igor.comp$seurat_clusters
igor.comp.deg = FindAllMarkers(igor.comp, only.pos = T)
igor.comp.deg = igor.comp.deg[which(igor.comp.deg$p_val_adj < 0.05),]

tj.deg = FindAllMarkers(tj, only.pos = T)
tj.deg = tj.deg[which(tj.deg$p_val_adj < 0.05),]
# tj.deg = convertMzebraDFToMouse(tj.deg, 7) # Convert to Mouse genes
tj.deg = hgncMzebraInPlace(tj.deg, 7, rownames(tj))
heatmapComparison(tj.deg, igor.comp.deg, "TJ", "Incisor+Molar", "ig_comp_v_tj", igor_path)

jaw.deg = FindAllMarkers(jaw, only.pos = T)
jaw.deg = jaw.deg[which(jaw.deg$p_val_adj < 0.05),]
# jaw.deg = convertMzebraDFToMouse(jaw.deg, 7) # Convert to Mouse genes
jaw.deg = hgncMzebraInPlace(jaw.deg, 7, rownames(jaw))
heatmapComparison(jaw.deg, igor.comp.deg, "Jaw", "Incisor+Molar", "ig_comp_v_jaw", igor_path)

# A sophisticated point based system. Get points for having same direction DEG (up and down regulated)
# Get more points for lower p_val_adj
Idents(igor.mouse) = igor.mouse$annot
igor.mouse.deg = FindAllMarkers(igor.mouse, only.pos = F)
igor.mouse.deg = igor.mouse.deg[which(igor.mouse.deg$p_val_adj < 0.05),]

tj.deg = FindAllMarkers(tj, only.pos = F)
tj.deg = tj.deg[which(tj.deg$p_val_adj < 0.05),]
tj.deg = convertMzebraDFToMouse(tj.deg, 7) # Convert to Mouse genes
# tj.deg = hgncMzebraInPlace(tj.deg, 7, rownames(tj))
heatmapComparison(tj.deg, igor.mouse.deg, "TJ", "Incisor", "im_vs_tj", igor_path)

jaw.deg = FindAllMarkers(jaw, only.pos = F)
jaw.deg = jaw.deg[which(jaw.deg$p_val_adj < 0.05),]
jaw.deg = convertMzebraDFToMouse(jaw.deg, 7) # Convert to Mouse genes
# jaw.deg = hgncMzebraInPlace(jaw.deg, 7, rownames(jaw))
heatmapComparison(jaw.deg, igor.mouse.deg, "Jaw", "Incisor", "im_vs_jaw", igor_path)

###########################
# Similarity b/w Matrices #
###########################
tj.mouse.avg   = myAverageExpression(tj.mouse, cells = colnames(tj.mouse))
Idents(igor.mouse) = igor.mouse$cond
igor.mouse.avg = myAverageExpression(igor.mouse, cells = colnames(igor.mouse))

igor.mouse.avg = setNames(as.vector(igor.mouse.avg[,1]), rownames(igor.mouse.avg))
tj.mouse.avg   = setNames(as.vector(tj.mouse.avg[,1]), rownames(tj.mouse.avg))
igor.mouse.avg = igor.mouse.avg[which(names(igor.mouse.avg) %in% names(tj.mouse.avg))]
tj.mouse.avg   = tj.mouse.avg  [which(names(tj.mouse.avg)   %in% names(igor.mouse.avg))]

pct.sim = c()
for (gene in names(tj.mouse.avg)) {
  tot = igor.mouse.avg[gene] + tj.mouse.avg[gene]
  diff = abs(igor.mouse.avg[gene] - tj.mouse.avg[gene])
  pct.sim = c( pct.sim, 1 - (diff/tot) )
}
names(pct.sim) = names(tj.mouse.avg)
pct.sim = pct.sim[which(! is.na(pct.sim) )]
mean(pct.sim)

pct.sim.bi = c()
for (gene in names(tj.mouse.avg)) {
  igor_is_zero = igor.mouse.avg[gene] == 0
  tj_is_zero = tj.mouse.avg[gene] == 0
  pct.sim.bi = c( pct.sim.bi, igor_is_zero ==tj_is_zero )
}
names(pct.sim.bi) = names(tj.mouse.avg)
length(pct.sim.bi[which(pct.sim.bi==1)])/length(pct.sim.bi)

##########
# DEGS ##
##########

deg_table = table(deg$gene)
deg$n_gene_appears = deg_table[match(deg$gene, names(deg_table))]

appears = data.frame()
for (gene in unique(deg$gene)) {
  rows = deg[which(deg$gene == gene),]
  appears = rbind(appears, t(c(gene, paste(rows$cluster, collapse = ", ") )))
}
deg$DEG_in = appears$V2[match(deg$gene, appears$V1)]
write.csv(deg, "C:/Users/miles/Downloads/d_tooth/results/hm_pos_degs.csv")
write.csv(deg, "~/scratch/d_tooth/results/hm_pos_degs.csv")

top = data.frame(row.names = 1:100)
top_unique = data.frame(row.names = 1:100)
for (cluster in levels(Idents(hm))) {
  rows = deg[which(deg$cluster == cluster),]
  rows_unique = rows[which(rows$n_gene_appears == 1),]
  top = cbind(top, rows$gene[1:100])
  top_unique = cbind(top_unique, rows_unique$gene[1:100])
}
colnames(top) = levels(Idents(hm))
colnames(top_unique) = levels(Idents(hm))
write.table(top, "~/scratch/d_tooth/data/hm_100.tsv", sep="\t", quote=F, row.names=F)
write.table(top_unique, "~/scratch/d_tooth/data/hm_unique_100.tsv", sep="\t", quote=F, row.names=F)
write.xlsx(top, "C:/Users/miles/Downloads/incsr_100.xlsx", row.names = F)
write.xlsx(top_unique, "C:/Users/miles/Downloads/incsr_unique_100.xlsx", row.names = F)
