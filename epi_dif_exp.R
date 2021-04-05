library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("monocle3")
library("RColorBrewer")
library("EnhancedVolcano")

rna_path <- "C:/Users/miles/Downloads/d_tooth/"

ctrl.data <- Read10X(data.dir = paste(rna_path, "data/EPI_CTRL/outs/filtered_feature_bc_matrix/", sep=""))
injr.data <- Read10X(data.dir = paste(rna_path, "data/EPI_INJR/outs/filtered_feature_bc_matrix/", sep=""))

ctrl <- CreateSeuratObject(counts = ctrl.data, project = "CTRL")
injr <- CreateSeuratObject(counts = injr.data, project = "INJR")

ctrl$cond <- "CTRL"
injr$sample <- "INJR"
ctrl$cond <- "CTRL"
injr$sample <- "INJR"

ctrl <- NormalizeData(ctrl, normalization.method = "LogNormalize", scale.factor = 100000)
injr <- NormalizeData(injr, normalization.method = "LogNormalize", scale.factor = 100000)

# INDV
epi <- merge(ctrl, injr, merge.data = TRUE)
epi <- subset(epi, subset = nFeature_RNA > 1000)
epi <- subset(epi, subset = nFeature_RNA < 4000)
# There are major signs of doublets here. Removing 85 cells.
epi <- FindVariableFeatures(object = epi, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
epi <- ScaleData(object = epi, vars.to.regress = NULL)
epi <- RunPCA(epi, npcs = 30, verbose = FALSE)
# JackStraw Values: 
# All PC's significant
epi <- RunUMAP(epi, reduction = "pca", dims = 1:30)
epi <- FindNeighbors(epi, reduction = "umap", dims = 1:2)
epi <- FindClusters(epi, resolution = 0.30)
DimPlot(epi, reduction = "umap", split.by = "orig.ident", label = TRUE)
# saveRDS(epi, "C:/Users/miles/Downloads/d_tooth/data/epi_full.rds")
# epi <- readRDS("C:/Users/miles/Downloads/d_tooth/data/epi_full.rds")

# Remove Mesechymal and Immune Clusters
# Actually Todd says to leave them in
# Idents(epi) <- "seurat_clusters"
# mes_and_immune_clusters <- c("14", "1", "16", "17", "18")
# epi_clusters <- levels(epi@active.ident)[which(! levels(epi@active.ident) %in% mes_and_immune_clusters)]
# epi <- subset(epi, cells = WhichCells(epi, idents = epi_clusters))
# saveRDS(epi, "C:/Users/miles/Downloads/d_tooth/data/epi_epi.rds")

# Recluster
# epi <- FindNeighbors(epi, reduction = "umap", dims = 1:2)
# epi <- FindClusters(epi, resolution = 0.30)
# DimPlot(epi, reduction = "umap", split.by = "orig.ident", label = TRUE)
# saveRDS(epi, "C:/Users/miles/Downloads/d_tooth/data/epi.rds")

# DEG Across Clusters
rna_path <- "C:/Users/miles/Downloads/d_tooth/"
epi.nk.markers <- FindAllMarkers(epi)
epi.nk.markers <- epi.nk.markers[,c(ncol(epi.nk.markers)-1, ncol(epi.nk.markers), 1:(epi.ncol(nk.markers)-2)),]
epi.nk.markers <- epi.nk.markers[which(epi.nk.markers$p_val_adj < 0.05),]
write.table(epi.nk.markers, file = paste(rna_path, "/results/indv_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

# DEG Across Conditions
Idents(epi) <- epi$cond
lj.markers <- FindMarkers(epi, ident.1 = "LJ", ident.2 = "UJ", only.pos = TRUE)
uj.markers <- FindMarkers(epi, ident.1 = "UJ", ident.2 = "LJ", only.pos = TRUE)
lj.markers <- lj.markers[which(lj.markers$p_val_adj < .05),]
uj.markers <- uj.markers[which(uj.markers$p_val_adj < .05),]
write.table(lj.markers, file = paste(rna_path, "/results/lj_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = TRUE)
write.table(uj.markers, file = paste(rna_path, "/results/uj_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = TRUE)
# EnhancedVolcano(epi.cond.markers, lab = rownames(epi.cond.markers), x = "avg_logFC", y = "p_val_adj", pCutoff = 0.05, subtitle = NULL, col=c("black", "black", "black", "red3"), transcriptPointSize =3, title="LJ vs UJ")

# Painting
rna_path <- "C:/Users/miles/Downloads/d_tooth/"
adrenergic <- c("Adra2b", "Adra1d", "Adra2c", "Adrbk2", "Adrb3",  "Adra1a", "Adra1b", "Adrb2", "Adrbk1", "Adra2a", "Adrb1")
mes <- c("Serpinf1", "Pcolce", "Msx1", "Col1a1", "Mfap4", "Ctsk", "Nbl1", "Cxcl12", "Prrx2")
immune <- c("Fcer1g", "C1qb", "C1qa", "C1qc", "Tyrobp", "B2m", "H2-D1", "Ftl1", "Lyz2")
g2_m <- c("Ube2c","Nusap1","Cenpf","Prc1","Ccnb1","Top2a","Kpna2","Arl6ip1","Cdc20","Cks2","Cdk1","Spc25","Lockd","Pbk","Birc5","Hmgb2","H2afx","Cenpa","Tuba1c","Hn1")
oee_2 <- c("Tsc22d1","Cxcl14","Prnp","Igfbp5","Gpha2","Gas6","Krt15","Igfbp2","Fst","Fos","Cd59a","Sparc","Dcn","Sostdc1","Lum","Cyr61","Col11a1","Fam132a","Tgfb2","Malat1")
fig6 <- c("Ccnb1", "Birc5", "Igfbpl1", "Ank2", "Sfrp5", "Cldn10")
pro <- c("Prlr")
val_adr <- c("Adrb3", "Src", "Pnpla2", "Abhd5", "Srebf1", "Adcy2", "Adcy7", "Snap23", "Snap29", "Stx3", "Vamp1", "Vamp2", "Vamp3", "Vamp8")
cur_list <- val_adr
for (gene in cur_list) {
  print(gene)
  png(paste(rna_path, "/results/painting/epi/fig6/epi_full_", gene, "_umap.png", sep=""), width = 900, height = 500, unit = "px")
  p <- FeaturePlot(epi, features = gene, reduction = "umap", split.by = "orig.ident", pt.size = 2, label=TRUE, order = TRUE)
  print(p)
  dev.off()
}
# High Quality
for (gene in cur_list) {
  print(gene)
  png(paste(rna_path, "/results/painting/epi/epi_full_", gene, "_umap.png", sep=""), width = 4800, height = 2400, unit="px", res=300)
  p <- FeaturePlot(epi, features = gene, reduction = "umap", split.by = "orig.ident", pt.size = 2, label=TRUE, order = TRUE)
  print(p)
  dev.off()
}

# Feature Scatters
adrenergic <- c("Adra2b", "Adra1d", "Adra2c", "Adrbk2", "Adrb3",  "Adra1a", "Adra1b", "Adrb2", "Adrbk1", "Adra2a", "Adrb1")
val_adr <- c("Adrb3", "Celsr1", "Src", "Pnpla2", "Abhd5", "Srebf1", "Adcy2", "Adcy7", "Snap23", "Snap29", "Stx3", "Vamp1", "Vamp2", "Vamp3", "Vamp8")
results <- data.frame()
Idents(epi) <- "orig.ident"
epi@active.assay <- "RNA"
cur_list <- unique(c(val_adr, adrenergic))
for (i in 1:length(cur_list)) {
  gene = cur_list[i]
  for (j in i:length(cur_list)) {
    gene2 = cur_list[j]
    if(gene != gene2) {
      print(paste(gene, gene2))
      ctrl_cells = WhichCells(epi, idents = "CTRL")
      injr_cells = WhichCells(epi, idents = "INJR")
      
      # Do a Fisher's Exact Test, where you test the ratio of positive gene1 cells
      # in gene2 cells compared to the number of positive cell gene1 cells in all cells.
      all_p  <- twoGeneFisher(gene, gene2, epi)
      ctrl_p <- twoGeneFisher(gene, gene2, epi[, ctrl_cells])
      injr_p <- twoGeneFisher(gene, gene2, epi[, injr_cells])
      
      # Do a Spearman's Correlation between the two genes in the normalized data
      all_cor  = cor.test(epi@assays$RNA@data[gene,], epi@assays$RNA@data[gene2,])
      ctrl_cor = cor.test(epi@assays$RNA@data[gene, ctrl_cells], epi@assays$RNA@data[gene2, ctrl_cells])
      injr_cor = cor.test(epi@assays$RNA@data[gene, injr_cells], epi@assays$RNA@data[gene2, injr_cells])
      
      # Find the jaccard index
      gene_cells = colnames(epi)[which(epi@assays$RNA@counts[gene,] != 0)]
      gene2_cells = colnames(epi)[which(epi@assays$RNA@counts[gene2,] != 0)]
      all_j  = jaccard(colnames(epi) %in% gene_cells, colnames(epi) %in% gene2_cells)
      ctrl_j = jaccard(ctrl_cells %in% gene_cells,    ctrl_cells %in% gene2_cells)
      injr_j = jaccard(injr_cells %in% gene_cells,    injr_cells %in% gene2_cells)
      
      all_j_p  = jaccard.test(colnames(epi) %in% gene_cells, colnames(epi) %in% gene2_cells, method="mca", accuracy=1e-5)
      ctrl_j_p = jaccard.test(ctrl_cells %in% gene_cells,    ctrl_cells %in% gene2_cells, method="mca", accuracy=1e-5)
      injr_j_p = jaccard.test(injr_cells %in% gene_cells,    injr_cells %in% gene2_cells, method="mca", accuracy=1e-5)
      
      results <- rbind(results, t(c(gene, gene2, all_p, ctrl_p, injr_p, all_cor$p.value, ctrl_cor$p.value, injr_cor$p.value, all_j_p$pvalue, ctrl_j_p$pvalue, injr_j_p$pvalue, all_cor$estimate, ctrl_cor$estimate, injr_cor$estimate, all_j, ctrl_j, injr_j)))
      # p_comb <- FeatureScatter(epi, gene, gene2, slot = "counts", col=c("red", "blue"))
      # p_ctrl <- FeatureScatter(epi[,WhichCells(epi, idents = "CTRL")], gene, gene2, slot = "counts", col = "red")
      # p_injr <- FeatureScatter(epi[,WhichCells(epi, idents = "INJR")], gene, gene2, slot = "counts", col = "blue")
      # results <- rbind(results, t(c(gene, gene2, p_comb$plot_env$plot.cor, p_ctrl$plot_env$plot.cor, p_injr$plot_env$plot.cor)))
    }
  }
}
# colnames(results) <- c("gene1", "gene2", "combined_cor", "ctrl_cor", "injr_cor")
# colnames(results) <- c("gene1", "gene2", "all_p", "ctrl_p", "injr_p")
colnames(results) = c("gene1", "gene2", "all_p", "ctrl_p", "injr_p", "all_cor_p", "ctrl_cor_p", "injr_cor_p", "all_j_p", "ctrl_j_p", "injr_j_p", "all_cor", "ctrl_cor", "injr_cor", "all_j", "ctrl_j", "injr_j")
results_rht <- results[which(results$gene1 == "Celsr1"),]
results_rht$all_q  <- p.adjust(as.numeric(as.vector(results_rht$all_p)), method = "hochberg")
results_rht$ctrl_q <- p.adjust(as.numeric(as.vector(results_rht$ctrl_p)), method = "hochberg")
results_rht$injr_q <- p.adjust(as.numeric(as.vector(results_rht$injr_p)), method = "hochberg")
# write.table(results, paste(rna_path, "/results/epi_val_adr_cor.tsv", sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(results_rht, paste(rna_path, "/results/epi_val_adr_fisher_rht.tsv", sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

results$i = 1:nrow(results)
df2 <- results %>% pivot_longer(c("all_p", "all_cor_p", "all_j_p"), names_to = "variable", values_to = "value")
df2$value = as.numeric(as.vector(df2$value))
ggplot(df2, aes(i, value, color=variable)) + geom_point() + geom_line()

results$all_q     = p.adjust(as.vector(results$all_p))
results$all_cor_q = p.adjust(as.vector(results$all_cor_p))
results$all_j_q   = p.adjust(as.vector(results$all_j_p))
results[which(results$all_q < 0.05 & results$all_cor_q < 0.05),]
# Helper Functions
library(biomaRt)
library(rentrez)
library(dplyr)
library(stringr)

# Helper Functions
twoGeneFisher <- function(gene1, gene2, obj) {
  expr <- FetchData(object = obj, vars = gene1, slot = "counts")
  expr2 <- FetchData(object = obj, vars = gene2, slot = "counts")
  gene_1_cells <- colnames(obj@assays$RNA@counts[, which(x = expr > 0)])
  gene_2_cells <- colnames(obj@assays$RNA@counts[, which(x = expr2 > 0)])
  ovlp <- gene_1_cells[which(gene_1_cells %in% gene_2_cells)]
  contig_table <- data.frame(c(length(ovlp), length(gene_2_cells)), c(length(gene_1_cells), ncol(obj)))
  fisher_p <- fisher.test(contig_table)$p.value
  # print(contig_table)
  return(fisher_p)
}

mzebraToMouseInPlace <- function(obj) {
  genes <- rownames(obj@assays$RNA@counts)
  mouse  = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  zfin_genes    <- getLDS(attributes = c("zfin_id_symbol"), filters = "zfin_id_symbol", values = genes , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  hgnc_genes    <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  
  all_mouse <- rbind(ensembl_genes, setNames(zfin_genes, names(ensembl_genes)), setNames(hgnc_genes, names(ensembl_genes)))
  all_mouse <- all_mouse[!duplicated(all_mouse[,2]),]
  all_mouse <- all_mouse[!duplicated(all_mouse[,1]),]
  colnames(all_mouse) <- c("mzebra", "mouse")
  
  new_counts_matrix <- as.matrix(obj@assays$RNA@counts)
  new_data_matrix   <- as.matrix(obj@assays$RNA@counts)
  # new_counts_matrix <- new_counts_matrix[which(rownames(new_counts_matrix) %in% all_mouse[,1]),]
  for (i in 1:nrow(all_mouse)) {
    mzebra_gene <- all_mouse[i,1]
    mouse_gene  <- all_mouse[i,2]
    # print(mzebra_gene)
    # print(mouse_gene)
    # print(ind)
    ind <- which(genes == mzebra_gene)
    rownames(new_counts_matrix)[ind] <- mouse_gene
    rownames(new_data_matrix)[ind]   <- mouse_gene
  }
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  obj_2$seurat_clusters <- obj$seurat_clusters
    
  return(obj_2)
}

# Find Epi Clusters in Tooth Data
# Compare Each TJ cluster to the Epi data
tj <- mzebraToMouseInPlace(tj_backup)
Idents(epi_epi) <- "EPI"
Idents(tj) <- tj$seurat_clusters
all_teeth <- merge(epi_epi, tj, merge.data = TRUE)
num_tj_clusters <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
big_df <- data.frame()
for (clust in 0:num_tj_clusters) {
  clust.markers <- FindMarkers(all_teeth, ident.1 = clust, ident.2 = "EPI", only.pos = TRUE)
  clust.markers <- cbind(clust.markers, rep(clust, nrow(clust.markers)))
  big_df <- rbind(big_df, clust.markers)
}
big_df$gene <- rownames(big_df)
big_df <- big_df[,c(ncol(big_df), ncol(big_df)-1, 1:(ncol(big_df)-2))]
colnames(big_df)[1:2] <- c("gene", "cluster")
write.table(big_df, paste(rna_path, "/results/epi_tj_deg.tsv", sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Find Mes clusters in tooth data
# Compare Each TJ cluster to the Mes data
Idents(mes) <- "seurat_clusters"
mes_mes_clusters <- c("16", "20")
mes_clusters <- levels(mes@active.ident)[which(! levels(mes@active.ident) %in% mes_mes_clusters)]
mes_mes <- subset(mes, cells = WhichCells(mes, idents = mes_clusters))
Idents(mes_mes) <- "MES"
Idents(tj) <- tj$seurat_clusters
all_teeth <- merge(mes_mes, tj, merge.data = TRUE)
num_tj_clusters <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
big_df <- data.frame()
for (clust in 0:num_tj_clusters) {
  clust.markers <- FindMarkers(all_teeth, ident.1 = clust, ident.2 = "MES", only.pos = TRUE)
  clust.markers <- cbind(clust.markers, rep(clust, nrow(clust.markers)))
  big_df <- rbind(big_df, clust.markers)
}
big_df$gene <- rownames(big_df)
big_df <- big_df[,c(ncol(big_df), ncol(big_df)-1, 1:(ncol(big_df)-2))]
colnames(big_df)[1:2] <- c("gene", "tj_cluster")
write.table(big_df, paste(rna_path, "/results/mes_tj_deg.tsv", sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Paint Cell Types
rna_path <- "C:/Users/miles/Downloads/d_tooth/"
Idents(epi) <- "seurat_clusters"
new.cluster.ids <- c("Class 2", "Immune", "Class 1", "Class 3", "Class 1", "Class 3", "Class 1", "Class 3", "Class 3", "Class 3", "Class 2", "Class 3", "Class 1", "Class 3", "Mesenchymal", "Class 3", "Immune", "Immune", "Immune")
names(new.cluster.ids) <- levels(epi)
epi <- RenameIdents(epi, new.cluster.ids)
png(filename = paste(rna_path, "results/", "epi_label.png", sep=""), width = 2400, height = 1800, unit="px", res=300)
DimPlot(epi, reduction = "umap", cols = c("purple", "cadetblue1", "chartreuse", "tomato", "gold"), label = TRUE)
dev.off()

# Recluster with only Celsr1 Cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
epi$cyto <- epi_cyto$CytoTRACE
epi_celsr_cells <- colnames(epi@assays$RNA@counts[,which(epi@assays$RNA@counts["Celsr1",] > 0)])
matrix2 <- matrix(as.numeric(as.vector(epi_cyto$CytoTRACE)), nrow = 1, ncol = length(epi_cyto$CytoTRACE), dimnames = list("cyto", names(epi_cyto$CytoTRACE)) )
epi[["cyto"]] <- CreateAssayObject(counts = matrix2)

epi_celsr1 <- subset(epi, cells = epi_celsr_cells)
epi_celsr1$orig.clust <- epi$seurat_clusters[epi_celsr_cells]
epi_celsr1 <- FindVariableFeatures(object = epi_celsr1, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
epi_celsr1 <- ScaleData(object = epi_celsr1, vars.to.regress = NULL)
epi_celsr1 <- RunPCA(epi_celsr1, npcs = 30, verbose = FALSE)
epi_celsr1 <- RunUMAP(epi_celsr1, reduction = "pca", dims = 1:12)
epi_celsr1 <- FindNeighbors(epi_celsr1, reduction = "umap", dims = 1:2)
epi_celsr1 <- FindClusters(epi_celsr1, resolution = 0.30)
Idents(epi_celsr1) <- "seurat_clusters"
DimPlot(epi_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(epi_celsr1) <- epi_celsr1$orig.clust
DimPlot(epi_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(epi_celsr1) <- epi_celsr1$Phase
DimPlot(epi_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(epi_celsr1) <- epi_celsr1$cyto
DimPlot(epi_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)

png(paste0(rna_path, "/results/epi_recluster_cyto_old_clust.png"), width = 1200, height = 600, res = 150)
Idents(epi_celsr1) <- "orig.clust"
p <- FeaturePlot(epi_celsr1, features = "cyto", reduction = "umap", cols = rev(brewer.pal(11,"Spectral")), split.by = "cond", pt.size = 2, label=TRUE, order = TRUE)
print(p)
dev.off()

# Bin by CytoTRACE score
epi_celsr1$bin <- epi_celsr1$cyto
epi_celsr1$bin[which(epi_celsr1$cyto <= quantile(epi_celsr1$cyto, 0.33))] <- "relative_low"
epi_celsr1$bin[which(epi_celsr1$cyto > quantile(epi_celsr1$cyto, 0.33) & epi_celsr1$cyto <= quantile(epi_celsr1$cyto, 0.66))] <- "relative_medium"
epi_celsr1$bin[which(epi_celsr1$cyto > quantile(epi_celsr1$cyto, 0.66))] <- "relative_high"

Idents(epi_celsr1) <- epi_celsr1$bin
epi_celsr1_relative_cyto_deg <- FindAllMarkers(epi_celsr1)
DimPlot(epi_celsr1, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Cluster mes celsr1 with epi celsr
epi[[paste0("NormalizeData.RNA")]]      <- tj[[paste0("NormalizeData.RNA")]]
combined[[paste0("NormalizeData.RNA")]] <- tj[[paste0("NormalizeData.RNA")]]
matrix4 <- matrix(as.numeric(as.vector(mes_cyto$CytoTRACE)), nrow = 1, ncol = length(mes_cyto$CytoTRACE), dimnames = list("cyto", names(mes_cyto$CytoTRACE)) )
combined[["cyto"]] <- CreateAssayObject(counts = matrix4)
epi$dataset <- "EPI"
combined$dataset <- "MES"
mes_celsr_cells <- colnames(combined@assays$RNA@counts[,which(combined@assays$RNA@counts["Celsr1",] > 0)])
all_celsr1 <- RunCCA(subset(epi, cells = epi_celsr_cells), subset(combined, cells = mes_celsr_cells), assay1 = "RNA", assay2 = "RNA", renormalize = TRUE, rescale = TRUE)
matrix2 <- matrix(as.numeric(as.vector(c(epi_cyto$CytoTRACE[epi_celsr_cells], mes_cyto$CytoTRACE[mes_celsr_cells]))), nrow = 1, ncol = length(c(epi_celsr_cells, mes_celsr_cells)), dimnames = list("cyto", c(epi_celsr_cells, mes_celsr_cells)) )
all_celsr1[["cyto"]] <- CreateAssayObject(counts = matrix2)
# all_celsr1 <- merge(subset(epi, cells = epi_celsr_cells), subset(combined, cells = mes_celsr_cells))
all_celsr1$orig.clust <- c(paste0(epi$seurat_clusters[epi_celsr_cells]), paste0(combined$seurat_clusters[mes_celsr_cells]))
all_celsr1 <- FindVariableFeatures(object = all_celsr1, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
all_celsr1 <- ScaleData(object = all_celsr1, vars.to.regress = NULL)
all_celsr1 <- RunPCA(all_celsr1, npcs = 30, verbose = FALSE)
all_celsr1 <- RunUMAP(all_celsr1, reduction = "pca", dims = 1:12)
all_celsr1 <- FindNeighbors(all_celsr1, reduction = "umap", dims = 1:2)
all_celsr1 <- FindClusters(all_celsr1, resolution = 0.30)
Idents(all_celsr1) <- "seurat_clusters"
DimPlot(all_celsr1, reduction = "umap", split.by = "dataset", label = TRUE)
Idents(all_celsr1) <- all_celsr1$orig.clust
DimPlot(all_celsr1, reduction = "umap", split.by = "dataset", label = TRUE)
p1 <- FeaturePlot(all_celsr1, features = "cyto", split.by = "dataset", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE, slot = "counts")[[1]] + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1))
p2 <- FeaturePlot(all_celsr1, features = "cyto", split.by = "dataset", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE, slot = "counts")[[2]] + scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1))
plot_grid(p1,p2)
# p1 <- FeaturePlot(all_celsr1, features = "cyto", split.by = "dataset", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE, slot = "counts") 
# fix.sc <- scale_color_gradientn( colors = rev(brewer.pal(11,"Spectral")), limits = c(0,1))
# lapply(p1, function (x) x )
# CombinePlots(p2)

# Do cycling vs immune have significantly different CytoTRACE scores?
# Answer: Yes, significant p
cycling_clust <- c("4", "6", "12")
immune_clust <- c("1", "16")
Idents(epi_celsr1) <- "orig.clust"
cycling_clust_cells <- WhichCells(epi_celsr1, idents = cycling_clust)
immune_clust_cells  <- WhichCells(epi_celsr1, idents = immune_clust)
cycling_clust_cyto <- epi_celsr1$cyto[cycling_clust_cells]
immune_clust_cyto <- epi_celsr1$cyto[immune_clust_cells]
t.test(cycling_clust_cyto, immune_clust_cyto)

# Do the cycling Celsr have lower Celsr counts?
# Answer: No, 1.05 v 1.01
cycling_clust_counts <- epi_celsr1@assays$RNA@counts["Celsr1",cycling_clust_cells]
immune_clust_counts <- epi_celsr1@assays$RNA@counts["Celsr1",immune_clust_cells]
mean(cycling_clust_counts)
mean(immune_clust_counts)

# Try to cluster by Celsr1 and Gli1
cycling_clust <- c("2", "4", "6", "12")
pre           <- c("0", "10")
class3        <- c("1", "3", "5", "7", "8", "9", "11", "13", "14", "15", "16", "17", "18")
epi_celsr_cells <- colnames(epi@assays$RNA@counts[,which(epi@assays$RNA@counts["Celsr1",] > 0)])
epi_gli1_cells  <- colnames(epi@assays$RNA@counts[,which(epi@assays$RNA@counts["Gli1",]   > 0)])
epi_trj <- epi[,unique(c(WhichCells(epi, idents = c(cycling_clust, pre)), epi_celsr_cells, epi_gli1_cells))]
epi_trj <- FindVariableFeatures(object = epi_trj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 20)
epi_trj <- ScaleData(object = epi_trj, vars.to.regress = c("Celsr1", "Gli1"))
epi_trj <- RunPCA(epi_trj, npcs = 30, verbose = FALSE, features = c("Celsr1", "Gli1"), umap.method = 'umap-learn', metric = 'correlation')
epi_trj <- RunUMAP(epi_trj, reduction = "pca", features = c("Celsr1", "Gli1"))
epi_trj <- FindNeighbors(epi_trj, reduction = "umap", dims = 1:2)
epi_trj <- FindClusters(epi_trj, resolution = 0.30)
DimPlot(epi_trj, reduction = "umap", split.by = "orig.ident", label = TRUE)
FeaturePlot(epi_trj, "Celsr1", order = TRUE, pt.size = 1.5)
FeaturePlot(epi_trj, "Gli1", order = TRUE, pt.size = 1.5)
results <- CytoTRACE(as.matrix(epi_trj@assays$RNA@counts))
epi_trj$CytoTRACE <- results$CytoTRACE
less_50 <- names(epi_trj$CytoTRACE[which(epi_trj$CytoTRACE < 0.5)])
matrix2 <- matrix(as.numeric(as.vector(epi_trj$CytoTRACE)), nrow = 1, ncol = length(epi_trj$CytoTRACE), dimnames = list("cyto", names(epi_trj$CytoTRACE)) )
matrix3 <- matrix(1 - as.numeric(as.vector(epi_trj$CytoTRACE)), nrow = 1, ncol = length(epi_trj$CytoTRACE), dimnames = list("rev_cyto", names(epi_trj$CytoTRACE)) )
epi_trj[["cyto"]]    <- CreateAssayObject(counts = matrix2)
epi_trj[["rev-cyto"]] <- CreateAssayObject(counts = matrix3)
FeaturePlot(epi_trj, features = "cyto", reduction = "umap", cols = rev(brewer.pal(11,"Spectral")), pt.size = 2, label=TRUE, order = TRUE) + NoLegend()
FeaturePlot(epi_trj, features = "rev-cyto", reduction = "umap", cols = brewer.pal(11,"Spectral"), pt.size = 2, label=TRUE, order = TRUE) + NoLegend()
plotCytoTRACE(results, emb = epi_trj@reductions$umap@cell.embeddings[,1:2], outputDir = "C:/Users/miles/Downloads/d_tooth/results/epi_celsr1_gli1_cyto.png")
# png("/nv/hp10/ggruenhagen3/scratch/d_tooth/results/Celsr1.png")
# p <- FeaturePlot(epi_trj, "Celsr1", order = TRUE, pt.size = 1.5)
# print(p)
# dev.off()

# Recluster with only Gli1 Cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
epi$cyto <- epi_cyto$CytoTRACE
epi_celsr_cells <- colnames(epi@assays$RNA@counts[,which(epi@assays$RNA@counts["Gli1",] > 0)])
matrix2 <- matrix(as.numeric(as.vector(epi_cyto$CytoTRACE)), nrow = 1, ncol = length(epi_cyto$CytoTRACE), dimnames = list("cyto", names(epi_cyto$CytoTRACE)) )
epi[["cyto"]] <- CreateAssayObject(counts = matrix2)

epi_gli1 <- subset(epi, cells = epi_celsr_cells)
epi_gli1$orig.clust <- epi$seurat_clusters[epi_celsr_cells]
epi_gli1 <- FindVariableFeatures(object = epi_gli1, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
epi_gli1 <- ScaleData(object = epi_gli1, vars.to.regress = NULL)
epi_gli1 <- RunPCA(epi_gli1, npcs = 30, verbose = FALSE)
epi_gli1 <- RunUMAP(epi_gli1, reduction = "pca", dims = 1:12)
epi_gli1 <- FindNeighbors(epi_gli1, reduction = "umap", dims = 1:2)
epi_gli1 <- FindClusters(epi_gli1, resolution = 0.30)
Idents(epi_gli1) <- "seurat_clusters"
DimPlot(epi_gli1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(epi_gli1) <- epi_gli1$orig.clust
DimPlot(epi_gli1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(epi_gli1) <- epi_gli1$Phase
DimPlot(epi_gli1, reduction = "umap", split.by = "orig.ident", label = TRUE)
Idents(epi_gli1) <- as.numeric(as.vector(epi_gli1$cyto))
DimPlot(epi_gli1, reduction = "umap", split.by = "orig.ident") + NoLegend()
FeaturePlot(epi_gli1, features = "cyto", reduction = "umap", cols = rev(brewer.pal(11,"Spectral")), split.by = "orig.ident", pt.size = 2, order = TRUE)

cycling_clust <- c("2", "4", "6", "12")
cycling_cells <- WhichCells(epi, idents = c(cycling_clust))
class3_cells  <- WhichCells(epi, idents = c(class3))
class2_cells  <- WhichCells(epi, idents = c(pre))
df <- setNames(data.frame(epi_cyto$CytoTRACE, epi_cyto$CytoTRACErank, epi$orig.ident, names(epi_cyto$CytoTRACE) %in% epi_celsr_cells, names(epi_cyto$CytoTRACE) %in% epi_gli1_cells, names(epi_cyto$CytoTRACE) %in% cycling_cells, names(epi_cyto$CytoTRACE) %in% class2_cells, names(epi_cyto$CytoTRACE) %in% class3_cells), c("CytoTRACE", "CytoTRACErank", "cond", "isCelsr1", "isGli1", "isCycling", "isClass2", "isClass3"))
df$isCelsr1[which(df$isCelsr1)] <- "Celsr1"
df$isGli1[which(df$isGli1)] <- "Gli1"
df$isCycling[which(df$isCycling)] <- "Class1"
df$isClass2[which(df$isClass2)] <- "Class2"
df$isClass3[which(df$isClass3)] <- "Class3"
df2 <- df %>% pivot_longer(c("isCelsr1", "isGli1", "isCycling", "isClass2", "isClass3"), names_to = "variable", values_to = "value")
df2 <- df2[which(df2$value != "FALSE"),]
df2$value <- factor(df2$value, levels=c("Celsr1", "Gli1", "Class1", "Class2", "Class3"))
ggplot(df2, aes(x = value, y = CytoTRACErank, fill = cond)) + geom_boxplot(alpha = 0.6) + geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.3, aes(colour = cond)) + scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Mouse Epithelium: CytoTRACE Rank for Cells Expressing a Gene")
ggplot(df2, aes(x = value, y = CytoTRACE, fill = cond)) + geom_boxplot(alpha = 0.6) + geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.3, aes(colour = cond)) + scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Mouse Epithelium: CytoTRACE Score (0-1) for Cells Expressing a Gene")

# Boxplot of Celsr1 Expression by CytoBin for each Dataset
df <- data.frame(rep("Mouse_Mes", ncol(mes_celsr1)), mes_celsr1$cyto2, mes_celsr1$bin, as.vector(mes_celsr1@assays$RNA@data["Celsr1",]))
colnames(df) <- c("Dataset", "CytoTRACE", "CytoBIN", "Celsr1_Expr")
df <- rbind(df, setNames(data.frame(rep("Mouse_Epi", ncol(epi_celsr1)), epi_celsr1$cyto, epi_celsr1$bin, as.vector(epi_celsr1@assays$RNA@data["Celsr1",])), c("Dataset", "CytoTRACE", "CytoBIN", "Celsr1_Expr")))
df <- rbind(df, setNames(data.frame(rep("Cichlid_Jaw", ncol(jaw_celsr1)), jaw_celsr1$cyto, jaw_celsr1$bin, as.vector(jaw_celsr1@assays$RNA@data["celsr1a",])), c("Dataset", "CytoTRACE", "CytoBIN", "Celsr1_Expr")))
png(paste0(rna_path, "/results/celsr1_cytobin.png"), width = 1800, height = 1200, res = 300)
p <- ggplot(df, aes(Dataset, Celsr1_Expr, colour = CytoBIN, fill = CytoBIN)) + geom_violin(bw = 0.5, position = position_jitterdodge(jitter.width = 0, dodge.width = 1), alpha = 0.6) + geom_jitter(shape=16, position=position_jitterdodge(dodge.width = 1), alpha = 0.5) + stat_summary(fun.y=median, geom="point", size=2, color="red", position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) + ylim(0, 6) + ggtitle("Celsr1 Expression in CytoTRACE Bins per Dataset with only Celsr1+ Cells")
print(p)                                                                                     
dev.off()

# Compare high bin to low bin for random genes
combined$bin <- combined$cyto
combined$bin[which(combined$cyto <= quantile(combined$cyto, 0.33))] <- "low"
combined$bin[which(combined$cyto > quantile(combined$cyto, 0.33) & combined$cyto <= quantile(combined$cyto, 0.66))] <- "medium"
combined$bin[which(combined$cyto > quantile(combined$cyto, 0.66))] <- "high"
epi$bin <- epi$cyto
epi$bin[which(epi$cyto <= quantile(epi$cyto, 0.33))] <- "low"
epi$bin[which(epi$cyto > quantile(epi$cyto, 0.33) & epi$cyto <= quantile(epi$cyto, 0.66))] <- "medium"
epi$bin[which(epi$cyto > quantile(epi$cyto, 0.66))] <- "high"
jaw$bin <- jaw$cyto
jaw$bin[which(jaw$cyto <= quantile(jaw$cyto, 0.33))] <- "low"
jaw$bin[which(jaw$cyto > quantile(jaw$cyto, 0.33) & jaw$cyto <= quantile(jaw$cyto, 0.66))] <- "medium"
jaw$bin[which(jaw$cyto > quantile(jaw$cyto, 0.66))] <- "high"
Idents(combined) <- combined$bin
Idents(epi) <- epi$bin
Idents(jaw) <- jaw$bin
n_genes <- 50
mes_nonzero <- names(which(rowSums(combined@assays$RNA@counts[,]) > 30))
epi_nonzero <- names(which(rowSums(epi@assays$RNA@counts[,]) > 30))
jaw_nonzero <- names(which(rowSums(jaw@assays$RNA@counts[,]) > 30))
ran_genes_epi <- sample(epi_nonzero, n_genes)
ran_genes_mes <- sample(mes_nonzero, n_genes)
ran_genes_jaw <- sample(jaw_nonzero, n_genes)
df <- myAverageExpression(epi, features = ran_genes_epi, cells = epi_celsr_cells)
df <- rbind(df, myAverageExpression(combined, features = ran_genes_mes, cells = mes_celsr_cells))
df <- rbind(df, myAverageExpression(jaw, features = ran_genes_jaw, cells = jaw_celsr_cells))
df$dataset <- c(rep("Epi", n_genes), rep("Mes", n_genes), rep("Jaw", n_genes))
df$highest_bin <- colnames(df)[apply(df,1,which.max)]
df$lowest_bin <- colnames(df)[apply(df,1,which.min)]
ggplot(df, aes(dataset, fill = highest_bin, color = highest_bin)) + geom_bar(alpha = 0.7) + ggtitle(paste("CytoBIN with the highest mean expression for", n_genes, "genes - Celsr1+ Only"))
ggplot(df, aes(dataset, fill = lowest_bin, color = lowest_bin)) + geom_bar(alpha = 0.7) + ggtitle(paste("CytoBIN with the LOWEST mean expression for", n_genes, "genes - Celsr1+ Only"))

all_mouse <- convertMzebraGeneListToMouse(rownames(tj))
success <- 0
mode <- "_pos"
data_slot <- "counts"
df2 <- data.frame()
Idents(jaw) = jaw$bin
Idents(mes) = mes$bin
Idents(combined) = combined$bin
Idents(epi) = epi$bin
for (i in 1:nrow(all_mouse))  {
  tryCatch({
    gene_mzebra <- all_mouse[i,1]
    gene_mouse  <- all_mouse[i,2]
    if (data_slot == "counts") {
      epi_pos <- colnames(epi@assays$RNA@counts[,which(epi@assays$RNA@counts[gene_mouse,] > 0)])
      mes_pos <- colnames(mes@assays$RNA@counts[,which(mes@assays$RNA@counts[gene_mouse,] > 0)])
      jaw_pos <- colnames(jaw@assays$RNA@counts[,which(jaw@assays$RNA@counts[gene_mzebra,] > 0)])
    } else {
      epi_pos <- colnames(epi@assays$RNA@data[,which(epi@assays$RNA@data[gene_mouse,] > 0)])
      mes_pos <- colnames(mes@assays$RNA@data[,which(mes@assays$RNA@data[gene_mouse,] > 0)])
      jaw_pos <- colnames(jaw@assays$RNA@data[,which(jaw@assays$RNA@data[gene_mzebra,] > 0)])
    }
    
    if ( mode == "_pos" ) {
      title = paste0(gene_mzebra, "/", gene_mouse, " Expression in CytoTRACE Bins per Dataset - Positive Only")
      epi_cells <- epi_pos
      mes_cells <- mes_pos
      jaw_cells <- jaw_pos
    } else if ( mode == "_celsr_pos" ) {
      title = paste0(gene_mzebra, "/", gene_mouse, " Expression in CytoTRACE Bins per Dataset - Celsr+ Only")
      epi_cells <- epi_celsr_cells
      mes_cells <- mes_celsr_cells
      jaw_cells <- jaw_celsr_cells
    } else {
      print("Not a correct mode")
      break
    }
    
    df <-           setNames(as.data.frame(epi@assays$RNA@data[gene_mouse,      epi_cells]), c("Expr"))
    df <- rbind(df, setNames(as.data.frame(combined@assays$RNA@data[gene_mouse, mes_cells]), c("Expr")))
    df <- rbind(df, setNames(as.data.frame(jaw@assays$RNA@data[gene_mzebra,     jaw_cells]), c("Expr")))
    df$dataset <- c(rep("Epi", length(epi_cells)), rep("Mes", length(mes_cells)), rep("Jaw", length(jaw_cells)))
    df$CytoBIN <- c(epi$bin[epi_cells], combined$bin[mes_cells], jaw$bin[jaw_cells])
    
    # epi_pos <- colnames(epi@assays$RNA@counts[,which(epi@assays$RNA@counts[gene_mouse,] > 0)])
    # mes_pos <- colnames(mes@assays$RNA@counts[,which(mes@assays$RNA@counts[gene_mouse,] > 0)])
    # jaw_pos <- colnames(jaw@assays$RNA@counts[,which(jaw@assays$RNA@counts[gene_mzebra,] > 0)])
    # epi_cells <- epi_pos
    # mes_cells <- mes_pos
    # jaw_cells <- jaw_pos
    # df <-           setNames(as.data.frame(epi@assays$RNA@counts[gene_mouse,      epi_cells]), c("Expr"))
    # df <- rbind(df, setNames(as.data.frame(combined@assays$RNA@counts[gene_mouse, mes_cells]), c("Expr")))
    # df <- rbind(df, setNames(as.data.frame(jaw@assays$RNA@counts[gene_mzebra,     jaw_cells]), c("Expr")))
    # df$dataset <- c(rep("Epi", length(epi_cells)), rep("Mes", length(mes_cells)), rep("Jaw", length(jaw_cells)))
    # df$CytoBIN <- c(epi$bin[epi_cells], combined$bin[mes_cells], jaw$bin[jaw_cells])
    
    
    if(length(epi_pos[which(epi_pos %in% epi_cells)]) > 10 && length(mes_pos[which(mes_pos %in% mes_cells)]) > 10 && length(jaw_pos[which(jaw_pos %in% jaw_cells)]) > 10) {
      # png(paste0(rna_path, "/results/violin", mode, "/", gene_mzebra, ".png"), width = 2200, height = 1200, res = 300)
      # p <- ggplot(df, aes(dataset, Expr, colour = CytoBIN, fill = CytoBIN)) + geom_violin(bw = 0.5, position = position_jitterdodge(jitter.width = 0, dodge.width = 1), alpha = 0.6) + geom_jitter(shape=16, position=position_jitterdodge(dodge.width = 1), alpha = 0.5) + stat_summary(fun=median, geom="point", size=2, color="red", position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) + ggtitle(title)
      # print(p)
      # dev.off()
      
      epi_all_num_cells <- c(length(epi_cyto$CytoTRACE[which(epi_cyto$CytoTRACE > 0.66)])                                         , length(epi_cyto$CytoTRACE[which(epi_cyto$CytoTRACE < 0.66 & epi_cyto$CytoTRACE > 0.33)])                                         , length(epi_cyto$CytoTRACE[which(epi_cyto$CytoTRACE < 0.33)]))
      epi_pos_num_cells <- c(length(epi_cyto$CytoTRACE[which(epi_cyto$CytoTRACE > 0.66 & names(epi_cyto$CytoTRACE) %in% epi_pos)]), length(epi_cyto$CytoTRACE[which(epi_cyto$CytoTRACE < 0.66 & epi_cyto$CytoTRACE > 0.33 & names(epi_cyto$CytoTRACE) %in% epi_pos)]), length(epi_cyto$CytoTRACE[which(epi_cyto$CytoTRACE < 0.33 & names(epi_cyto$CytoTRACE) %in% epi_pos)]))
      mes_all_num_cells <- c(length(mes_cyto$CytoTRACE[which(mes_cyto$CytoTRACE > 0.66)])                                         , length(mes_cyto$CytoTRACE[which(mes_cyto$CytoTRACE < 0.66 & mes_cyto$CytoTRACE > 0.33)])                                         , length(mes_cyto$CytoTRACE[which(mes_cyto$CytoTRACE < 0.33)]))
      mes_pos_num_cells <- c(length(mes_cyto$CytoTRACE[which(mes_cyto$CytoTRACE > 0.66 & names(mes_cyto$CytoTRACE) %in% mes_pos)]), length(mes_cyto$CytoTRACE[which(mes_cyto$CytoTRACE < 0.66 & mes_cyto$CytoTRACE > 0.33 & names(mes_cyto$CytoTRACE) %in% mes_pos)]), length(mes_cyto$CytoTRACE[which(mes_cyto$CytoTRACE < 0.33 & names(mes_cyto$CytoTRACE) %in% mes_pos)]))
      jaw_all_num_cells <- c(length(jaw_cyto$CytoTRACE[which(jaw_cyto$CytoTRACE > 0.66)])                                         , length(jaw_cyto$CytoTRACE[which(jaw_cyto$CytoTRACE < 0.66 & jaw_cyto$CytoTRACE > 0.33)])                                         , length(jaw_cyto$CytoTRACE[which(jaw_cyto$CytoTRACE < 0.33)]))
      jaw_pos_num_cells <- c(length(jaw_cyto$CytoTRACE[which(jaw_cyto$CytoTRACE > 0.66 & names(jaw_cyto$CytoTRACE) %in% jaw_pos)]), length(jaw_cyto$CytoTRACE[which(jaw_cyto$CytoTRACE < 0.66 & jaw_cyto$CytoTRACE > 0.33 & names(jaw_cyto$CytoTRACE) %in% jaw_pos)]), length(jaw_cyto$CytoTRACE[which(jaw_cyto$CytoTRACE < 0.33 & names(jaw_cyto$CytoTRACE) %in% jaw_pos)]))
      
      df2 <- rbind(df2, cbind(gene = gene_mouse,  myAverageExpression(epi,      slot = data_slot, features = gene_mouse)                   , dataset = "epi", isSubset = "F", cells = t(epi_all_num_cells)))
      df2 <- rbind(df2, cbind(gene = gene_mouse,  myAverageExpression(combined, slot = data_slot, features = gene_mouse)                   , dataset = "mes", isSubset = "F", cells = t(mes_all_num_cells)))
      df2 <- rbind(df2, cbind(gene = gene_mouse,  myAverageExpression(jaw, slot = data_slot, features = gene_mzebra)                       , dataset = "jaw", isSubset = "F", cells = t(jaw_all_num_cells)))
      df2 <- rbind(df2, cbind(gene = gene_mouse,  myAverageExpression(epi, slot = data_slot, features = gene_mouse     , cells = epi_cells), dataset = "epi", isSubset = "T", cells = t(epi_pos_num_cells)))
      df2 <- rbind(df2, cbind(gene = gene_mouse,  myAverageExpression(combined, slot = data_slot, features = gene_mouse, cells = mes_cells), dataset = "mes", isSubset = "T", cells = t(mes_pos_num_cells)))
      df2 <- rbind(df2, cbind(gene = gene_mzebra, myAverageExpression(jaw, slot = data_slot, features = gene_mzebra    , cells = jaw_cells), dataset = "jaw", isSubset = "T", cells = t(jaw_pos_num_cells)))
      success <- success + 1
      print(paste(gene_mzebra, "- Success!"))
      if (success > 50) {
        break
      }
    } else {
      print(paste(gene_mzebra, "(empty)"))
    }
    

  }, error = function(e) {
    print(gene_mouse)
    print(e)
  })
  
}
df2$highest_cells <- colnames(df2)[c("cells.1", "cells.2", "cells.3")][apply(df2,1,which.max)]
df2$highest_bin <- sapply(1:nrow(df2), function(x) c("high", "low", "medium")[which.max(df2[x,c("high", "low", "medium")])])
df2$lowest_bin <- sapply(1:nrow(df2), function(x) c("high", "low", "medium")[which.min(df2[x,c("high", "low", "medium")])])
ggplot(df2[which(df2$isSubset == "T"),], aes(dataset, fill = highest_bin, color = highest_bin)) + geom_bar() + ggtitle(paste("CytoBIN with the HIGHEST mean expression for ViolinPlot/Orthologous genes - Only Positive Cells"))
ggplot(df2[which(df2$isSubset == "T"),], aes(dataset, fill = lowest_bin, color = lowest_bin)) + geom_bar() + ggtitle(paste("CytoBIN with the LOWEST mean expression for ViolinPlot/Orthologous genes - Only Positive Cells"))
ggplot(df2[which(df2$isSubset == "F"),], aes(dataset, fill = highest_bin, color = highest_bin)) + geom_bar() + ggtitle(paste("CytoBIN with the HIGHEST mean expression for ViolinPlot/Orthologous genes - All Cells"))
ggplot(df2[which(df2$isSubset == "F"),], aes(dataset, fill = lowest_bin, color = lowest_bin)) + geom_bar() + ggtitle(paste("CytoBIN with the LOWEST mean expression for ViolinPlot/Orthologous genes - All Cells"))


vp1 <- ggplot(df2[which(df2$isSubset == "T"),], aes(dataset, fill = highest_bin, color = highest_bin)) + geom_bar() + ggtitle(paste("HIGHEST Expression for ViolinPlot Genes - Only Celsr1+ Cells"))
vp2 <- ggplot(df2[which(df2$isSubset == "T"),], aes(dataset, fill = lowest_bin, color = lowest_bin)) + geom_bar() + ggtitle(paste("LOWEST Expression for ViolinPlot genes - Only Celsr1+ Cells"))
vp3 <- ggplot(df2[which(df2$isSubset == "F"),], aes(dataset, fill = highest_bin, color = highest_bin)) + geom_bar() + ggtitle(paste("HIGHEST Expression for ViolinPlot genes - All Cells"))
vp4 <- ggplot(df2[which(df2$isSubset == "F"),], aes(dataset, fill = lowest_bin, color = lowest_bin)) + geom_bar() + ggtitle(paste("LOWEST Expression for ViolinPlot genes - All Cells"))
png(paste0(rna_path, "/results/violin/_high_low.png"), width = 1800, height = 1000, res = 150)
p <- plot_grid(vp1, vp2, vp3, vp4)
print(p)
dev.off()

# View the Avg_Expr, Gene_Count, and Total_Expr for all genes in CytoBINs
# df <- myAverageExpression(epi)
# df <- rbind(df, myAverageExpression(combined))
# df <- rbind(df, myAverageExpression(jaw))
# df$dataset <- c(rep("Epi", nrow(epi)), rep("Mes", nrow(combined)), rep("Jaw", nrow(jaw)))
# df$CytoBIN <- c(epi$bin, combined$bin, jaw$bin)
# df$value <- "Avg_Expr"
# df <- rbind(df, myTotalTrans(epi))
# df <- rbind(df, myTotalTrans(combined))
# df <- rbind(df, myTotalTrans(jaw))
# df$dataset <- c(rep("Epi", nrow(epi)), rep("Mes", nrow(combined)), rep("Jaw", nrow(jaw)), rep("Epi", nrow(epi)), rep("Mes", nrow(combined)), rep("Jaw", nrow(jaw)))
# df$type <- c(rep("Avg_Expr", (nrow(epi)+nrow(combined)+nrow(jaw))), rep("Total_Trans", (nrow(epi)+nrow(combined)+nrow(jaw))))
# df2 <- df %>% pivot_longer(c("medium", "low", "high"), names_to = "variable", values_to = "value")
epi_bin_div <- c(length(rownames(epi)[which(rowSums(epi@assays$RNA@counts[,WhichCells(epi, idents = "medium")]) > 0)]), length(rownames(epi)[which(rowSums(epi@assays$RNA@counts[,WhichCells(epi, idents = "low")]) > 0)]), length(rownames(epi)[which(rowSums(epi@assays$RNA@counts[,WhichCells(epi, idents = "high")]) > 0)]))
mes_bin_div <- c(length(rownames(combined)[which(rowSums(combined@assays$RNA@counts[,WhichCells(combined, idents = "low")]) > 0)]), length(rownames(combined)[which(rowSums(combined@assays$RNA@counts[,WhichCells(combined, idents = "medium")]) > 0)]), length(rownames(combined)[which(rowSums(combined@assays$RNA@counts[,WhichCells(combined, idents = "high")]) > 0)]))
jaw_bin_div <- c(length(rownames(epi)[which(rowSums(epi@assays$RNA@counts[,WhichCells(epi, idents = "low")]) > 0)]), length(rownames(epi)[which(rowSums(epi@assays$RNA@counts[,WhichCells(epi, idents = "medium")]) > 0)]), length(rownames(epi)[which(rowSums(epi@assays$RNA@counts[,WhichCells(epi, idents = "high")]) > 0)]))
epi_df3 <- cbind(as.data.frame(c(colMeans(sapply( myAverageExpression(epi), as.numeric ))[1:3], colSums(sapply( myTotalTrans(epi), as.numeric ))[1:3], epi_bin_div)),           rep(c("medium", "low", "high"),3), rep("Epi", 9), c(rep("Avg_Expr", 3), rep("Total_Trans",3), rep("Trans_Div",3)))
mes_df3 <- cbind(as.data.frame(c(colMeans(sapply( myAverageExpression(combined), as.numeric ))[1:3], colSums(sapply( myTotalTrans(combined), as.numeric ))[1:3], mes_bin_div)), rep(c("low", "medium", "high"),3), rep("Mes", 9), c(rep("Avg_Expr", 3), rep("Total_Trans",3), rep("Trans_Div",3)))
jaw_df3 <- cbind(as.data.frame(c(colMeans(sapply( myAverageExpression(jaw), as.numeric ))[1:3], colSums(sapply( myTotalTrans(jaw), as.numeric ))[1:3], jaw_bin_div)),           rep(c("low", "medium", "high"),3), rep("Jaw", 9), c(rep("Avg_Expr", 3), rep("Total_Trans",3), rep("Trans_Div",3)))
df3 <- rbind(setNames(epi_df3, c("value", "CytoBIN", "dataset", "variable")), setNames(mes_df3, c("value", "CytoBIN", "dataset", "variable")), setNames(jaw_df3, c("value", "CytoBIN", "dataset", "variable")))
df3$CytoBIN <- factor(df3$CytoBIN, levels = c("high", "medium", "low"))
total_p1 <- ggplot(df3[which(df3$variable == "Avg_Expr"),], aes(dataset, value, fill = CytoBIN, color = CytoBIN)) + geom_point(size = 5) + ylab("Average Expression")
total_p2 <- ggplot(df3[which(df3$variable == "Total_Trans"),], aes(dataset, value, fill = CytoBIN, color = CytoBIN)) + geom_point(shape = 17, size = 5) + ylab("Total Transranscripts")
total_p3 <- ggplot(df3[which(df3$variable == "Trans_Div"),], aes(dataset, value, fill = CytoBIN, color = CytoBIN)) + geom_point(shape = 15, size = 5) + ylab("Transcriptional Diversity")
png(paste0(rna_path, "/results/cytobin_stats.png"), width = 2000, height = 500, res = 150)
p <- plot_grid(total_p1, total_p2, total_p3, ncol = 3)
print(p)
dev.off()
# ggplot(df3, aes(dataset, value, shape = variable, color = CytoBIN)) + geom_point(size = 5) + scale_y_continuous("left", sec.axis = sec_axis(trans=~.*1e8, name = "right"))
# ggplot(df2[which(df2$type == "Avg_Expr"),],    aes(x=dataset, y=value, fill = variable, color = variable)) + geom_boxplot(outlier.shape = NA, alpha = 0.5) + ggtitle("Average Expression")
# ggplot(df2[which(df2$type == "Total_Trans"),], aes(x=dataset, y=value, fill = variable, color = variable)) + geom_boxplot(outlier.shape = NA, alpha = 0.5) + ggtitle("Total Transcripts")


# See if DEG pattern for random genes as for Celsr
fitsPattern <- data.frame()
for (i in 1:n_genes) {
  # Epithelium
  ran_gene <- ran_genes_epi[i]
  expr <- FetchData(object = epi, vars = ran_gene, slot = "counts")
  pos_cells <- colnames(epi[, which(x = expr > 0)])
  geneObj <- subset(epi, cells = pos_cells)
  Idents(geneObj) <- geneObj$bin
  geneObjDeg <- FindAllMarkers(geneObj)
  geneObjDeg <- geneObjDeg[which(geneObjDeg$p_val_adj < 0.05),]
  if (nrow(geneObjDeg) > 0) {
    lowLogFC    <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "low")])
    mediumLogFC <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "medium")])
    highLogFC   <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "high")])
    geneFitLow    <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "low" & geneObjDeg$gene == ran_gene)]) == sign(lowLogFC)
    geneFitMedium <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "medium" & geneObjDeg$gene == ran_gene)]) == sign(mediumLogFC)
    geneFitHigh   <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "high" & geneObjDeg$gene == ran_gene)]) == sign(highLogFC)
    geneFitLow <- ifelse(identical(geneFitLow, logical(0)), "NA", as.character(geneFitLow))
    geneFitMedium <- ifelse(identical(geneFitMedium, logical(0)), "NA", as.character(geneFitMedium))
    geneFitHigh <- ifelse(identical(geneFitHigh, logical(0)), "NA", as.character(geneFitHigh))
    newRow <- data.frame(rep("epi", 3), rep(ran_gene, 3), c(geneFitLow, geneFitMedium, geneFitHigh), c("Low", "Medium", "High"))
    colnames(newRow) <- c("dataset", "FitsPattern", "CytoBIN")
    fitsPattern <- rbind(fitsPattern, newRow)
  }

  # Mesenchyme
  ran_gene <- ran_genes_mes[i]
  expr <- FetchData(object = combined, vars = ran_gene, slot = "counts")
  pos_cells <- colnames(combined[, which(x = expr > 0)])
  geneObj <- subset(combined, cells = pos_cells)
  Idents(geneObj) <- geneObj$bin
  geneObjDeg <- FindAllMarkers(geneObj)
  geneObjDeg <- geneObjDeg[which(geneObjDeg$p_val_adj < 0.05),]
  if (nrow(geneObjDeg) > 0) {
    lowLogFC    <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "low")])
    mediumLogFC <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "medium")])
    highLogFC   <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "high")])
    geneFitLow    <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "low" & geneObjDeg$gene == ran_gene)]) == sign(lowLogFC)
    geneFitMedium <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "medium" & geneObjDeg$gene == ran_gene)]) == sign(mediumLogFC)
    geneFitHigh   <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "high" & geneObjDeg$gene == ran_gene)]) == sign(highLogFC)
    geneFitLow <- ifelse(identical(geneFitLow, logical(0)), "NA", as.character(geneFitLow))
    geneFitMedium <- ifelse(identical(geneFitMedium, logical(0)), "NA", as.character(geneFitMedium))
    geneFitHigh <- ifelse(identical(geneFitHigh, logical(0)), "NA", as.character(geneFitHigh))
    newRow <- data.frame(rep("medium", 3), rep(ran_gene, 3), c(geneFitLow, geneFitMedium, geneFitHigh), c("Low", "Medium", "High"))
    colnames(newRow) <- c("dataset", "FitsPattern", "CytoBIN")
    fitsPattern <- rbind(fitsPattern, newRow)
  }
  
  # Jaw
  ran_gene <- ran_genes_jaw[i]
  expr <- FetchData(object = jaw, vars = ran_gene, slot = "counts")
  pos_cells <- colnames(jaw[, which(x = expr > 0)])
  geneObj <- subset(jaw, cells = pos_cells)
  Idents(geneObj) <- geneObj$bin
  geneObjDeg <- FindAllMarkers(geneObj)
  geneObjDeg <- geneObjDeg[which(geneObjDeg$p_val_adj < 0.05),]
  if (nrow(geneObjDeg) > 0) {
    lowLogFC    <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "low")])
    mediumLogFC <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "medium")])
    highLogFC   <- sum(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "high")])
    geneFitLow    <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "low" & geneObjDeg$gene == ran_gene)]) == sign(lowLogFC)
    geneFitMedium <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "medium" & geneObjDeg$gene == ran_gene)]) == sign(mediumLogFC)
    geneFitHigh   <- sign(geneObjDeg$avg_logFC[which(geneObjDeg$cluster == "high" & geneObjDeg$gene == ran_gene)]) == sign(highLogFC)
    geneFitLow <- ifelse(identical(geneFitLow, logical(0)), "NA", as.character(geneFitLow))
    geneFitMedium <- ifelse(identical(geneFitMedium, logical(0)), "NA", as.character(geneFitMedium))
    geneFitHigh <- ifelse(identical(geneFitHigh, logical(0)), "NA", as.character(geneFitHigh))
    newRow <- data.frame(rep("jaw", 3), rep(ran_gene, 3), c(geneFitLow, geneFitMedium, geneFitHigh), c("Low", "Medium", "High"))
    colnames(newRow) <- c("dataset", "FitsPattern", "CytoBIN")
    fitsPattern <- rbind(fitsPattern, newRow)
  }
}
colnames(fitsPattern) <- c("dataset", "gene", "FitsPattern", "CytoBIN")
ggplot(fitsPattern, aes(FitsPattern, fill = CytoBIN)) + geom_bar(stat = "count") + ggtitle("Does the Random Gene Fit the Pattern of the Other Gene in the DEGs?")


#
celsr_p = c()
celsr_j = c()
for ( row in 1:length(gene_names) ) {
  gene2 <- gene_names[row]
  celsr_p = c(celsr_p, jaccard.test(gene_bi[[gene1]], gene_bi[[gene2]], method = "mca", accuracy=1e-5)$pvalue)
  celsr_j = c(celsr_j, jaccard(gene_bi[[gene1]], gene_bi[[gene2]]))
}
names(celsr_j) = gene_names
names(celsr_p) = gene_names
which.min(celsr_p)
celsr_p_0 = celsr_p[which(celsr_p == 0)]
length(celsr_p_0)
head(celsr_j[which(celsr_p == 0)])
celsr_r = jaccard.rahman(celsr_j)
head(sort(celsr_r))

#===================================================================================================
# Epi CTRL v INJR DEGs 03/26/2021 ==================================================================
#===================================================================================================
epi$cond_cluster = paste0(epi$cond, "_", epi$seurat_clusters)
Idents(epi) = epi$cond_cluster
all_cond_cluster_deg = data.frame()
for (cluster in levels(epi$seurat_clusters)) {
  print(cluster)
  if (paste0("CTRL_", cluster) %in% levels(Idents(epi)) && paste0("INJR_", cluster) %in% levels(Idents(epi)) && 
      length(WhichCells(epi, idents = paste0("CTRL_", cluster))) > 10 && length(WhichCells(epi, idents = paste0("INJR_", cluster))) > 10) {
    this_deg = FindMarkers(epi, ident.1 = paste0("CTRL_", cluster), ident.2 = paste0("INJR_", cluster))
    this_deg$gene = rownames(this_deg)
    this_deg$cluster = cluster
    this_deg$num_ctrl_cluster_cells = length(WhichCells(epi, idents = paste0("CTRL_", cluster)))
    this_deg$num_injr_cluster_cells = length(WhichCells(epi, idents = paste0("INJR_", cluster)))
    this_deg$pct_ctrl_cluster_cells = this_deg$num_ctrl_cluster_cells / (this_deg$num_ctrl_cluster_cells + this_deg$num_injr_cluster_cells)
    rownames(this_deg) = paste0(this_deg$gene, "_", cluster)
    all_cond_cluster_deg = rbind(all_cond_cluster_deg, this_deg)
  }
}
all_cond_cluster_deg$upCTRL = all_cond_cluster_deg$avg_logFC > 0
all_cond_cluster_deg = all_cond_cluster_deg[which(all_cond_cluster_deg$p_val_adj < 0.05),]
write.csv(all_cond_cluster_deg, "C:/Users/miles/Downloads/d_tooth/results/paul_epi_ctrl_v_injr_cluster_deg_sig.csv")

#===================================================================================================
# Mes CTRL v CLIPP DEGs 03/26/2021 =================================================================
#===================================================================================================
mes$cond_cluster = paste0(mes$cond, "_", mes$seurat_clusters)
Idents(mes) = mes$cond_cluster
all_cond_cluster_deg = data.frame()
for (cluster in levels(mes$seurat_clusters)) {
  print(cluster)
  if (paste0("CTRL_", cluster) %in% levels(Idents(mes)) && paste0("CLIPP_", cluster) %in% levels(Idents(mes)) && 
      length(WhichCells(mes, idents = paste0("CTRL_", cluster))) > 10 && length(WhichCells(mes, idents = paste0("CLIPP_", cluster))) > 10) {
    this_deg = FindMarkers(mes, ident.1 = paste0("CTRL_", cluster), ident.2 = paste0("CLIPP_", cluster))
    this_deg$gene = rownames(this_deg)
    this_deg$cluster = cluster
    this_deg$num_ctrl_cluster_cells = length(WhichCells(mes, idents = paste0("CTRL_", cluster)))
    this_deg$num_CLIPP_cluster_cells = length(WhichCells(mes, idents = paste0("CLIPP_", cluster)))
    this_deg$pct_ctrl_cluster_cells = this_deg$num_ctrl_cluster_cells / (this_deg$num_ctrl_cluster_cells + this_deg$num_CLIPP_cluster_cells)
    rownames(this_deg) = paste0(this_deg$gene, "_", cluster)
    all_cond_cluster_deg = rbind(all_cond_cluster_deg, this_deg)
  }
}
all_cond_cluster_deg$upCTRL = all_cond_cluster_deg$avg_logFC > 0
all_cond_cluster_deg = all_cond_cluster_deg[which(all_cond_cluster_deg$p_val_adj < 0.05),]
write.csv(all_cond_cluster_deg, "C:/Users/miles/Downloads/d_tooth/results/paul_mes_ctrl_v_CLIPP_cluster_deg_sig.csv")

all_cond_cluster_deg$annot = df$annot[match(all_cond_cluster_deg$cluster, df$cluster)]
test = all_cond_cluster_deg[, c(colnames(all_cond_cluster_deg)[1:7], "annot", colnames(all_cond_cluster_deg)[8:11])]
write.csv(test, "C:/Users/miles/Downloads/d_tooth/results/paul_mes_ctrl_v_CLIPP_cluster_deg_sig.csv")
