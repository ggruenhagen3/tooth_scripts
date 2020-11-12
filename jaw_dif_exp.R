library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("monocle3")
library("EnhancedVolcano")

rna_path <- "C:/Users/miles/Downloads/d_tooth/"

lj.data <- Read10X(data.dir = paste(rna_path, "data/LJ/outs/filtered_feature_bc_matrix/", sep=""))
uj.data <- Read10X(data.dir = paste(rna_path, "data/UJ/outs/filtered_feature_bc_matrix/", sep=""))
jpool.data <- Read10X(data.dir = paste(rna_path, "data/JPOOL/outs/filtered_feature_bc_matrix/", sep=""))

lj <- CreateSeuratObject(counts = lj.data, project = "LJ")
uj <- CreateSeuratObject(counts = uj.data, project = "UJ")
jpool <- CreateSeuratObject(counts = jpool.data, project = "JPOOL")

lj$cond <- "LJ"
lj$sample <- "LJ"
uj$cond <- "UJ"
uj$sample <- "UJ"
jpool$cond <- "JPOOL"
jpool$sample <- "JPOOL"

lj <- NormalizeData(lj, normalization.method = "LogNormalize", scale.factor = 100000)
uj <- NormalizeData(uj, normalization.method = "LogNormalize", scale.factor = 100000)
jpool <- NormalizeData(jpool, normalization.method = "LogNormalize", scale.factor = 100000)

# INDV
combined <- merge(lj, uj, merge.data = TRUE)
combined <- subset(combined, subset = nFeature_RNA > 300)
combined <- subset(combined, subset = nFeature_RNA < 2500)
# There are major signs of doublets here. Removing 85 cells.
combined <- FindVariableFeatures(object = combined, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
combined <- ScaleData(object = combined, vars.to.regress = NULL)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# JackStraw Values: 5.079762e-71  1.019556e-65 3.573636e-113  1.145205e-49  9.412051e-70  4.467056e-35 5.496905e-29  7.158615e-16  5.555104e-20  3.277846e-20  4.046250e-24  1.418974e-39 1.471733e-03  1.804553e-07  8.781465e-17  7.104022e-11  3.826268e-13  8.796275e-12 5.004931e-04  7.970772e-14
# Steep dropoff after PC 12
combined <- RunUMAP(combined, reduction = "pca", dims = 1:12)
combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
combined <- FindClusters(combined, resolution = 0.50)
DimPlot(combined, reduction = "umap", split.by = "cond", label = TRUE)
# saveRDS(combined, "C:/Users/miles/Downloads/d_tooth/tooth_scripts/jaw_shiny/data/jaw.rds")

# POOL
# FeatureScatter(jpool, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
jpool <- subset(jpool, subset = nFeature_RNA > 300)
jpool <- subset(jpool, subset = nFeature_RNA < 2500)
jpool <- FindVariableFeatures(object = jpool, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
jpool <- ScaleData(object = jpool, vars.to.regress = NULL)
jpool <- RunPCA(jpool, npcs = 30, verbose = FALSE)
# JackStraw Values: 1.373999e-74  1.817431e-65 1.758434e-105  1.294784e-51  4.594446e-82  9.763171e-23 1.972570e-26  1.209141e-15  2.042022e-15  9.412616e-20  8.212660e-25  5.087122e-36 1.007234e-04  1.658780e-14  1.484178e-16  1.482937e-11  1.074398e-17  7.816959e-09 1.471733e-03  1.482937e-11
# Steep dropoff after PC 12
jpool <- RunUMAP(jpool, reduction = "pca", dims = 1:12)
jpool <- FindNeighbors(jpool, reduction = "umap", dims = 1:2)
jpool <- FindClusters(jpool, resolution = 0.30)
DimPlot(jpool, reduction = "umap", split.by = "cond", label = TRUE)
# saveRDS(jpool, "C:/Users/miles/Downloads/d_tooth/tooth_scripts/pool_shiny/data/jpool.rds")
# jpool <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/pool_shiny/data/jpool.rds")
# png(paste(rna_path, "/results/jpool_umap.png", sep=""), width = 2000, height = 1500, res = 300, type = "cairo")
# DimPlot(jpool, reduction = "umap", split.by = "cond", label = TRUE)
# dev.off()

# Label Cell Types
Idents(jpool) <- "seurat_clusters"
new.cluster.ids <- c("Epithelial Stem Cells", "Epithelium", "HPSC", "Epithelium", "Epithelium", "Pigmented", "Fibroblast", "Epithelium", "DPSC", "HPSC", "Mature Taste Bud", "Immature Taste Bud")
names(new.cluster.ids) <- levels(jpool)
jpool <- RenameIdents(jpool, new.cluster.ids)
png(filename = paste(rna_path, "results/", "jpool.png", sep=""), width = 1500, height = 900, unit="px", res=200)
DimPlot(jpool, reduction = "umap", cols = c("cadetblue1", "blue", "darkorange", "pink", "firebrick1", "chartreuse", "gold", "purple"), label = TRUE, pt.size=1.5)
dev.off()

# DEG Across Clusters
rna_path <- "C:/Users/miles/Downloads/d_tooth/"
combined.nk.markers <- FindAllMarkers(combined)
combined.nk.markers <- combined.nk.markers[,c(ncol(combined.nk.markers)-1, ncol(combined.nk.markers), 1:(combined.ncol(nk.markers)-2)),]
combined.nk.markers <- combined.nk.markers[which(combined.nk.markers$p_val_adj < 0.05),]
write.table(combined.nk.markers, file = paste(rna_path, "/results/indv_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
jpool.nk.markers <- FindAllMarkers(jpool, only.pos = TRUE)
jpool.nk.markers <- jpool.nk.markers[,c(ncol(jpool.nk.markers)-1, ncol(jpool.nk.markers), 1:(ncol(jpool.nk.markers)-2)),]
jpool.nk.markers <- jpool.nk.markers[which(jpool.nk.markers$p_val_adj < 0.05),]
write.table(jpool.nk.markers, file = paste(rna_path, "/results/pool_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

# DEG Across Conditions
Idents(combined) <- combined$cond
lj.markers <- FindMarkers(combined, ident.1 = "LJ", ident.2 = "UJ", only.pos = TRUE)
uj.markers <- FindMarkers(combined, ident.1 = "UJ", ident.2 = "LJ", only.pos = TRUE)
lj.markers <- lj.markers[which(lj.markers$p_val_adj < .05),]
uj.markers <- uj.markers[which(uj.markers$p_val_adj < .05),]
write.table(lj.markers, file = paste(rna_path, "/results/lj_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = TRUE)
write.table(uj.markers, file = paste(rna_path, "/results/uj_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = TRUE)
# EnhancedVolcano(combined.cond.markers, lab = rownames(combined.cond.markers), x = "avg_logFC", y = "p_val_adj", pCutoff = 0.05, subtitle = NULL, col=c("black", "black", "black", "red3"), transcriptPointSize =3, title="LJ vs UJ")

# Jaw Painting
taste_bud <- c("fgf8a", "ascl1a", "mib1", "calb2a")
pres <- c("pitx2", "wnt10a", "bmp4")
clust_3 <- c("shha", "bmp4", "edar", "pitx2")
cycle <- c("unga", "top2a", "ccnb2", "cdc20")
epi <- c("krt5", "krt15", "krt8", "krt4")
paul <- c("mitfa", "celsr1a", "gad1b", "calb2a")
cur_list <- paul
for (gene in cur_list) {
  print(gene)
  png(paste(rna_path, "/results/painting/jpool/jpool_", gene, "_umap.png", sep=""), width = 1800, height = 1000, unit = "px", res = 200)
  p <- FeaturePlot(jpool, features = gene, reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
  print(p)
  dev.off()
}

# Feature Scatter
obj <- tj
gene1 <- "sox2"
gene2 <- "pitx2"
expr1 <- FetchData(object = obj, vars = gene1)
expr2 <- FetchData(object = obj, vars = gene2)
duo1 <- colnames(obj@assays$RNA@counts[, which(x = expr1 > 0)])
duo2 <- colnames(obj@assays$RNA@counts[, which(x = expr2 > 0)])
duo <- unique(c(duo1, duo2))
cor.test(obj@assays$RNA@counts[gene1,duo], obj@assays$RNA@counts[gene2,duo], method = "pearson")$p.value
png(paste(rna_path, "/results/fs_", gene1, "_", gene2, "_tj.png", sep=""), width = 600, height = 450, unit = "px")
FeatureScatter(obj[,duo], gene1, gene2, slot = "counts", pt.size = 1.5)
dev.off()

# Regress out Cell Cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- hgncGood(s.genes, rownames(jpool@assays$RNA@counts))
g2m.genes <- hgncGood(g2m.genes, rownames(jpool@assays$RNA@counts))
jpool <- CellCycleScoring(jpool, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
jpool <- RunPCA(jpool, features = c(s.genes, g2m.genes))

# Find the Cell Cycle State of Celsr1
expr <- FetchData(object = jpool, vars = "celsr1a", slot = "counts")
celsr1_cells <- colnames(jpool@assays$RNA@counts[, which(x = expr > 0)])
celsr1_cycle <- as.vector(Idents(jpool)[celsr1_cells])
table(celsr1_cycle)/length(celsr1_cells)
table(as.vector(Idents(jpool)))/length(Idents(jpool))

# Recluster celsr1 cells
jaw$cyto <- jaw_cyto$CytoTRACE
jaw_celsr_cells <- colnames(jaw@assays$RNA@counts[,which(jaw@assays$RNA@counts["celsr1a",] > 0)])
matrix2 <- matrix(as.numeric(as.vector(jaw_cyto$CytoTRACE)), nrow = 1, ncol = length(jaw_cyto$CytoTRACE), dimnames = list("cyto", names(jaw_cyto$CytoTRACE)) )
jaw[["cyto"]] <- CreateAssayObject(counts = matrix2)
jaw_celsr1 <- subset(jaw, cells = jaw_celsr_cells)
jaw_celsr1$orig.clust <- jaw$seurat_clusters[jaw_celsr_cells]
jaw_celsr1 <- FindVariableFeatures(object = jaw_celsr1, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
jaw_celsr1 <- ScaleData(object = jaw_celsr1, vars.to.regress = NULL)
jaw_celsr1 <- RunPCA(jaw_celsr1, npcs = 30, verbose = FALSE)
jaw_celsr1 <- RunUMAP(jaw_celsr1, reduction = "pca", dims = 1:12)
jaw_celsr1 <- FindNeighbors(jaw_celsr1, reduction = "umap", dims = 1:2)
jaw_celsr1 <- FindClusters(jaw_celsr1, resolution = 0.30)
Idents(jaw_celsr1) <- "seurat_clusters"
DimPlot(jaw_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(jaw_celsr1) <- jaw_celsr1$orig.clust
DimPlot(jaw_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)

FeaturePlot(jaw_celsr1, features = "cyto", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE, slot = "counts")[[1]] + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1))
FeaturePlot(jaw_celsr1, features = "cdkn1a", reduction = "umap", pt.size = 2, label = TRUE, order = TRUE)

jaw_celsr1_deg <- FindAllMarkers(jaw_celsr1)
sig_3 <- jaw_celsr1_deg[which(jaw_celsr1_deg$p_val_adj < 0.05 & jaw_celsr1_deg$cluster == 3),]
sig_4 <- jaw_celsr1_deg[which(jaw_celsr1_deg$p_val_adj < 0.05 & jaw_celsr1_deg$cluster == 4),]

# Bin by CytoTRACE score
jaw_celsr1$bin <- jaw_celsr1$cyto
jaw_celsr1$bin[which(jaw_celsr1$cyto < quantile(jaw_celsr1$cyto, 0.33))] <- "relative_low"
jaw_celsr1$bin[which(jaw_celsr1$cyto > quantile(jaw_celsr1$cyto, 0.33) & jaw_celsr1$cyto < quantile(jaw_celsr1$cyto, 0.66))] <- "relative_medium"
jaw_celsr1$bin[which(jaw_celsr1$cyto > quantile(jaw_celsr1$cyto, 0.66))] <- "relative_high"

Idents(jaw_celsr1) <- jaw_celsr1$bin
jaw_celsr1_relative_cyto_deg <- FindAllMarkers(jaw_celsr1)
DimPlot(jaw_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)

jaw_celsr1$bin <- factor(jaw_celsr1$bin, levels = c("relative_high", "relative_medium", "relative_low"))
Idents(jaw_celsr1) <- jaw_celsr1$bin
FeatureScatter(jaw_celsr1, "celsr1a", "cyto", slot = "data")
mes_celsr1$bin <- factor(mes_celsr1$bin, levels = c("relative_high", "relative_medium", "relative_low"))
Idents(mes_celsr1) <- mes_celsr1$bin
FeatureScatter(mes_celsr1, "Celsr1", "cyto", slot = "data")

# CYTOTRACE PAPER FIG 3D
test = data.frame(jaw$cyto, jaw_cyto$Counts)
test$celsr = jaw@assays$RNA@counts["celsr1a",]
test$isCelsr = test$celsr > 0
test = test[order(test$jaw.cyto),]
# ggplot(test, aes(jaw.cyto, jaw_cyto.Counts)) + geom_point(alpha = 0.2) + geom_line(color = "red", lwd=1.2, aes(y=rollmean(jaw_cyto.Counts, 50, na.pad=TRUE))) + scale_x_reverse() + xlab("CytoTRACE") + ylab("RNA content per cell") + ggtitle("Cichlid Jaw (CytoTRACE Paper Fig 3D)")
ggplot(test, aes(jaw.cyto, jaw_cyto.Counts, color = celsr, size=celsr)) + geom_point(aes(alpha = celsr)) + geom_line(color = "red", lwd=1.2, aes(y=rollmean(jaw_cyto.Counts, 50, na.pad=TRUE))) + scale_x_reverse() + xlab("CytoTRACE") + ylab("RNA content per cell") + ggtitle("Cichlid Jaw Celsr") + scale_color_viridis()
ggplot(test, aes(jaw.cyto, jaw_cyto.Counts, color = isCelsr)) + geom_point(aes(alpha = isCelsr)) + geom_line(color = "red", lwd=1.2, aes(y=rollmean(jaw_cyto.Counts, 50, na.pad=TRUE))) + scale_x_reverse() + xlab("CytoTRACE") + ylab("RNA content per cell") + ggtitle("Cichlid Jaw Celsr1 as Quiescent Marker") + scale_color_manual(values=c("black", "blue"))

#
data <- read.table(paste(rna_path, "/data/ophir_supp_1.txt", sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE)
for (col in 1:ncol(data)) {
  bio <- colnames(data)[col]
  if ( bio != "X" && ! endsWith(bio, "FC") ) {
    mzebra_genes <- mouseToMzebra(data[,col], rownames(jpool))
    write.table(data.frame(mzebra_genes, rep(bio, length(mzebra_genes))), paste(rna_path, "/data/markers/", bio, ".tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE) 
  }
}

# mouse to mzebra lazy
mouseToMzebra <- function(mouse_genes, mzebra_all) {
  mouse_lower <- tolower(mouse_genes)
  # mzebra_genes <- mzebr
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  mouse  = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = mouse_genes , mart = mouse, attributesL = c("external_gene_name"), martL = mzebra, uniqueRows=T)
  
  return(ensembl_genes[,2])
}