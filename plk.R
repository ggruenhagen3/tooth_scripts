#*******************************************************************************
# Load Libraries ===============================================================
#*******************************************************************************
wdstr = substr(getwd(), 1, 12)
switch(wdstr,
       "C:/Users/mil" = { main_path = "C:/Users/miles/Downloads/";        },
       "/home/george" = { main_path = "~/research/"                       },
       "/storage/scr" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/hom" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/cod" = { main_path = "/storage/scratch1/6/ggruenhagen3/" })
switch(wdstr,
       "C:/Users/mil" = { tooth_path = "d_tooth/" },
       "/home/george" = { tooth_path = "tooth/"   },
       "/storage/scr" = { tooth_path = "d_tooth/" },
       "/storage/hom" = { tooth_path = "d_tooth/" },
       "/storage/cod" = { tooth_path = "d_tooth/" })
switch(wdstr,
       "C:/Users/mil" = { gene_info_path = "all_research/" },
       "/home/george" = { gene_info_path = "all_research/" },
       "/storage/scr" = { gene_info_path = "m_zebra_ref/"  },
       "/storage/hom" = { gene_info_path = "m_zebra_ref/"  },
       "/storage/cod" = { gene_info_path = "m_zebra_ref/"  })
full_path = paste0(main_path, tooth_path)
gene_info_path = paste0(main_path, gene_info_path)
brain_dir = paste0(main_path, "brain/")
data_dir  = paste0(full_path, "data/")
out_dir   = paste0(full_path, "results/")
source(paste0(brain_dir, "/brain_scripts/all_f.R"))
setwd(out_dir)

#*******************************************************************************
# Load Objects =================================================================
#*******************************************************************************
gene_info = read.csv(paste0(gene_info_path, "gene_info_3.csv"))
plk = readRDS(paste0(data_dir, "plkall_053023.rds"))
div = readRDS(paste0(data_dir, "div_todd_090723.rds"))
plk_subject = readRDS(paste0(data_dir, "plkall_subject_053023.rds"))

#*******************************************************************************
# Diversity QC =================================================================
#*******************************************************************************
dir_of_sr_dirs = "~/scratch/brain/bs/JTS21/" # Folder where all the individual samples are kept
counts_list = list()
counts_list[["alt1"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS21_ALT_1/outs/filtered_feature_bc_matrix/"))
counts_list[["alt2"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS21_ALT_2/outs/filtered_feature_bc_matrix/"))
counts_list[["alt3"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS21_ALT_3/outs/filtered_feature_bc_matrix/"))

objs = list()
objs[["alt1"]] = CreateSeuratObject(counts_list[["alt1"]])
objs[["alt2"]] = CreateSeuratObject(counts_list[["alt2"]])
objs[["alt3"]] = CreateSeuratObject(counts_list[["alt3"]])

objs[["alt1"]]$sample = "alt1";
objs[["alt2"]]$sample = "alt2";
objs[["alt3"]]$sample = "alt3";

for (s in names(objs)) { objs[[s]] = RenameCells(objs[[s]], paste0(s)) }

# Check the quality
for (i in names(objs)) {
  obj = objs[[i]]
  plot1 = ggplot(data.frame(nCount = obj$nCount_RNA), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
  plot2 = VlnPlot(obj, features = "nCount_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  
  plot3 = ggplot(data.frame(nFeature = obj$nFeature_RNA), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
  plot4 = VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  pdf(paste0("~/scratch/d_tooth/results/div_", i, "_quality.pdf"), width = 8, height = 8)
  print(cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2))
  dev.off()
  print(paste0("# Cells UMI > 300 = ", length(which(obj$nCount_RNA > 300))))
  print(paste0("# Cells Genes > 300 = ", length(which(obj$nFeature_RNA > 300))))
}

# Remove cells
removed.df = data.frame()
for (i in names(objs)) {
  # objs[[i]]$pct.mt = colSums(objs[[i]]@assays$RNA@counts[mito.genes,]) / objs[[i]]$nCount_RNA
  # removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$pct.mt >= 0.05)], reason = "dead (% MT)"))
  if (length(which(objs[[i]]$nFeature_RNA >= 2500)) > 0) { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA >= 2500)], reason = "doublet"))       }
  if (length(which(objs[[i]]$nFeature_RNA <= 250)) > 0)  { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA <= 250)], reason = "dead (# Genes)")) }
  removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(!colnames(objs[[i]]) %in% removed.df$cell)], reason = "kept"))
}
removed.df$sample = toupper(reshape2::colsplit(removed.df$cell, "_", c('1', '2'))[,1])

pdf(paste0("~/scratch/d_tooth/results/div_nuclei_kept.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar() + xlab("Sample") + ylab("Number of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()
pdf(paste0("~/scratch/d_tooth/results/div_nuclei_kept_pct.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar(position = "fill") + xlab("Sample") + ylab("Proportion of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()

div = merge(objs[["alt1"]], objs[["alt3"]])

gtf = read.delim("~/scratch/m_zebra_ref/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = stringr::str_replace(mito.genes, "_", "-")
div$pct.mt = colSums(div@assays$RNA@counts[mito.genes,]) / div$nCount_RNA

div = subset(div, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & pct.mt < 0.05)
print(paste("Number of Cells in div After Filterning:", ncol(div)))
div = NormalizeData(div, normalization.method = "LogNormalize", scale.factor = 10000)
div = SCTransform(div, verbose = TRUE)
div@active.assay = "SCT"
saveRDS(div, "~/research/d_tooth/data/div_filtered_081423.rds")

#*******************************************************************************
# Initial clustering ===========================================================
#*******************************************************************************
dir_of_sr_dirs = "~/scratch/brain/bs/JTS15/" # Folder where all the individual samples are kept
counts_list = list()
counts_list[["p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15_clp12_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["p34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15_clp34_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15_cnt12_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["c34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15_cnt34_bcl_cr6/outs/filtered_feature_bc_matrix/"))

objs = list()
objs[["p12"]] = CreateSeuratObject(counts_list[["p12"]])
objs[["p34"]] = CreateSeuratObject(counts_list[["p34"]])
objs[["c12"]] = CreateSeuratObject(counts_list[["c12"]])
objs[["c34"]] = CreateSeuratObject(counts_list[["c34"]])

objs[["p12"]]$sample = "p12"; objs[["p12"]]$cond = "plk";
objs[["p34"]]$sample = "p34"; objs[["p34"]]$cond = "plk";
objs[["c12"]]$sample = "c12"; objs[["c12"]]$cond = "con";
objs[["c34"]]$sample = "c34"; objs[["c34"]]$cond = "con";

for (s in names(objs)) { objs[[s]] = RenameCells(objs[[s]], paste0(s)) }

# Check the quality
for (i in names(objs)) {
  obj = objs[[i]]
  plot1 = ggplot(data.frame(nCount = obj$nCount_RNA), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
  plot2 = VlnPlot(obj, features = "nCount_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())

  plot3 = ggplot(data.frame(nFeature = obj$nFeature_RNA), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
  plot4 = VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  pdf(paste0("~/scratch/d_tooth/results/plk_", i, "_quality.pdf"), width = 8, height = 8)
  print(cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2))
  dev.off()
  print(paste0("# Cells UMI > 300 = ", length(which(obj$nCount_RNA > 300))))
  print(paste0("# Cells Genes > 300 = ", length(which(obj$nFeature_RNA > 300))))
}

# Remove cells
removed.df = data.frame()
for (i in names(objs)) {
        objs[[i]]$pct.mt = colSums(objs[[i]]@assays$RNA@counts[mito.genes,]) / objs[[i]]$nCount_RNA
        # removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$pct.mt >= 0.05)], reason = "dead (% MT)"))
        if (length(which(objs[[i]]$nFeature_RNA >= 2500)) > 0) { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA >= 2500)], reason = "doublet"))       }
        if (length(which(objs[[i]]$nFeature_RNA <= 250)) > 0)  { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA <= 250)], reason = "dead (# Genes)")) }
        removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(!colnames(objs[[i]]) %in% removed.df$cell)], reason = "kept"))
}
removed.df$sample = toupper(reshape2::colsplit(removed.df$cell, "_", c('1', '2'))[,1])

pdf(paste0("~/scratch/d_tooth/results/plk_nuclei_kept.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar() + xlab("Sample") + ylab("Number of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()
pdf(paste0("~/scratch/d_tooth/results/plk_nuclei_kept_pct.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar(position = "fill") + xlab("Sample") + ylab("Proportion of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()

plk = merge(objs[["p12"]], list(objs[["p34"]],objs[["c12"]],objs[["c34"]]))

gtf = read.delim("~/scratch/m_zebra_ref/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = stringr::str_replace(mito.genes, "_", "-")
plk$pct.mt = colSums(plk@assays$RNA@counts[mito.genes,]) / plk$nCount_RNA

plk = subset(plk, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & pct.mt < 0.05)
print(paste("Number of Cells in plk After Filterning:", ncol(plk)))
plk = NormalizeData(plk, normalization.method = "LogNormalize", scale.factor = 10000)
plk = SCTransform(plk, verbose = TRUE)
plk@active.assay = "SCT"
plk = RunPCA(plk)
plk = RunUMAP(plk, reduction = "pca", dims = 1:50)
plk = FindNeighbors(plk, reduction="umap", dims = 1:2)
plk = FindClusters(plk, resolution = .25)

# Plot clusters and plot them by sample
# pdf(paste0("~/scratch/em/results/em_cluster.pdf"), width = 3.5, height = 3.5)
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_cluster_no_regress.png"), width = 1800, height = 1800, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1) + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_cluster_no_regress_split.png"), width = 3200, height = 1400, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1, split.by = "sample") + coord_fixed() + NoLegend())
dev.off()

# Plot the # of Genes per nuclei to see if any clusters have particularly low quality
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_nFeature_umap.png"), width = 1800, height = 1800, res = 200)
print(FeaturePlot(plk, "nFeature_RNA", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_bad_nFeature_umap.png"), width = 1800, height = 1800, res = 200)
plk$inverse_nFeature_RNA = 1/ plk$nFeature_RNA
print(FeaturePlot(plk, "inverse_nFeature_RNA", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()

# Check for cluster proportion differences
sample.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$sample))
s.clust.sdenom = sweep(sample.cluster.mat, 2, colSums(sample.cluster.mat), "/")
s.clust.cdenom = sample.cluster.mat / rowSums(sample.cluster.mat)
s.clust.sdenom.melt = reshape2::melt(s.clust.sdenom)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.sdenom.melt) = c("Cluster", "Sample", "value")
colnames(s.clust.cdenom.melt) = c("Cluster", "Sample", "value")

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_cluster_prop_by_sample.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Sample)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_sample_prop_by_cluster.png"), width = 1000, height = 1400, res = 200)
s.clust.sdenom.melt$Cluster = factor(s.clust.sdenom.melt$Cluster)
ggplot(s.clust.sdenom.melt, aes(x = Sample, y = value, fill = Cluster, color = Cluster)) + geom_bar(alpha = 0.3, stat = 'identity') + theme_bw() + guides(color = "none")
dev.off()

st.clust.deg$label2 = st.clust.deg$label
st.clust.deg$label2[which(st.clust.deg$label == "LOC101482992")] = "LOC101482992 (ACPT)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101476651")] = "LOC101476651 (KRT24)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101483408")] = "LOC101483408 (SLC5A7)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101485955")] = "LOC101485955 (krt5)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101482992")] = "LOC101482992 (ACPT)"
st.clust.deg$label2 = factor(st.clust.deg$label2, unique(st.clust.deg$label2))
ggplot(st.clust.deg, aes(x = label2, y = -log10(p_val_adj), fill = avg_logFC)) + geom_bar(stat = "identity", width = 0.5) + coord_flip() + theme_bw() + scale_y_continuous(expand = c(0,0)) +  facet_wrap(~ cluster, scales = "free_y") + xlab("") + scale_fill_viridis()

# Check for cluster proportion differences
plk.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$cond))
s.clust.pdenom = sweep(plk.cluster.mat, 2, colSums(plk.cluster.mat), "/")
s.clust.cdenom = plk.cluster.mat / rowSums(plk.cluster.mat)
s.clust.pdenom.melt = reshape2::melt(s.clust.pdenom)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.pdenom.melt) = c("Cluster", "Condition", "value")
colnames(s.clust.cdenom.melt) = c("Cluster", "Condition", "value")

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_cluster_prop_by_plk.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Condition)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk_plk_prop_by_cluster.png"), width = 1000, height = 1400, res = 200)
s.clust.pdenom.melt$Cluster = factor(s.clust.pdenom.melt$Cluster)
ggplot(s.clust.pdenom.melt, aes(x = Condition, y = value, color = Cluster, fill = Cluster)) + geom_bar(alpha = 0.3, stat = 'identity') + theme_bw() + guides(color = "none")
dev.off()

all.deg = data.frame()
Idents(plk) = plk$seurat_clusters
for (i in levels(plk$seurat_clusters)) {
        # markers <- FindMarkers(plk, ident.1 = "plk", ident.2 = "con", group.by = 'cond', subset.ident = i)
        markers <- FindMarkers(plk, ident.1 = "plk", ident.2 = "con", group.by = 'cond', subset.ident = i, min.pct = 1e-100, logfc.threshold = 0, only.pos = F)
        markers$gene = rownames(markers)
        markers$cluster = i
        all.deg = rbind(all.deg, markers)
}
all.deg$upInPlk = all.deg$avg_log2FC > 0
all.deg.sig = all.deg[which(all.deg$p_val_adj < 0.05 & abs(all.deg$avg_log2FC) > 0.25 & all.deg$pct.1 > 0.1),]
write.csv(all.deg.sig, "~/scratch/d_tooth/results/plk_plkvcon_cluster_standard_120922.csv")

# Look at # of PlkvCon DEGs using different thresholds
bulk.deg = read.csv("~/research/tooth/results/plk_plkvcon_plkup_bulk_deg_loose_hgnc.csv")
bulk.deg.num.df = data.frame()
for (this.fc in seq(0, 0.25, by = 0.0025)) {
  for (this.min.pct in seq(0, 0.1, by = 0.001)) {
    this.bulk.deg = bulk.deg[which(bulk.deg$abs.avg_logFC >= this.fc & bulk.deg$pct.1 >= this.min.pct & bulk.deg$p_val_adj < 0.05),]
    bulk.deg.num.df = rbind(bulk.deg.num.df, data.frame(fc = this.fc, min.pct = this.min.pct, n.deg = nrow(this.bulk.deg), avg.neg.log.bon = mean(-log10(this.bulk.deg$p_val_adj))  ))
  }
}
ggplot(bulk.deg.num.df, aes(x = fc, y = min.pct, fill = n.deg)) + geom_raster() + scale_fill_viridis() + theme_bw() + coord_fixed() + scale_x_continuous(name = "Log2FC", expand=c(0,0)) + scale_y_continuous(name = "Percent Expressed", expand=c(0,0))
library(metR)
ggplot(bulk.deg.num.df, aes(fc, min.pct)) + metR::geom_contour_fill(aes(z = n.deg)) + geom_contour(aes(z = n.deg), colour = "black") + metR::geom_text_contour(aes(z = n.deg), stroke = 0.15) + scale_fill_viridis() + theme_bw() + coord_fixed() + scale_x_continuous(name = "Log2FC", expand=c(0,0)) + scale_y_continuous(name = "Percent Expressed", expand=c(0,0))

bulk.deg$isSig = bulk.deg$p_val_adj < 0.05
bulk.deg$neg_log_bon = -log10(bulk.deg$p_val_adj)
bulk.deg$col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)[cut(bulk.deg$avg_log2FC, 100)]
bulk.deg$col[which( abs(bulk.deg$avg_log2FC) < 0.25 | bulk.deg$pct.1 < 0.1 | bulk.deg$p_val_adj >= 0.05)] = "gray80"
bulk.deg$isGray = bulk.deg$col == "gray80"
bulk.deg$label2 = bulk.deg$label
bulk.deg$label2[which(startsWith(bulk.deg$label2, "LOC") & !is.na(bulk.deg$hgnc))] = bulk.deg$hgnc[which(startsWith(bulk.deg$label2, "LOC") & !is.na(bulk.deg$hgnc))]
pdf("~/research/tooth/results/plk_plkvcon_bulk_deg.pdf", width = 6, height = 5)
print(ggplot(bulk.deg, aes(x = avg_logFC, y = neg_log_bon, size = pct.1, color = col, alpha = isGray)) + geom_point() + theme_bw() + scale_x_continuous(name = "Log2FC", expand = c(0.01, 0.01)) + scale_y_continuous(name = "-Log10(Adjusted p)", expand = c(0.03, 0.01)) + scale_color_identity() + geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray40") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") + scale_alpha_manual(values = c(0.8, 0.4)) + guides(alpha = "none", size = guide_legend(title="% Experessed")) + geom_text_repel(data = bulk.deg[which(!bulk.deg$isGray),], aes(label = label2)))
dev.off()


bulk.deg$isSig = bulk.deg$p_val_adj < 0.05
bulk.deg$neg_log_bon = -log10(bulk.deg$p_val_adj)
bulk.deg$col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)[cut(bulk.deg$avg_log2FC, 100)]
bulk.deg$col[which( abs(bulk.deg$avg_log2FC) < 0.1 | bulk.deg$pct.1 < 0.05 | bulk.deg$p_val_adj >= 0.05)] = "gray80"
bulk.deg$isGray = bulk.deg$col == "gray80"
bulk.deg$label2 = bulk.deg$label
bulk.deg$label2[which(startsWith(bulk.deg$label2, "LOC") & !is.na(bulk.deg$hgnc))] = bulk.deg$hgnc[which(startsWith(bulk.deg$label2, "LOC") & !is.na(bulk.deg$hgnc))]
pdf("~/research/tooth/results/plk_plkvcon_bulk_deg_fc1_pct05.pdf", width = 6, height = 5)
print(ggplot(bulk.deg, aes(x = avg_logFC, y = neg_log_bon, size = pct.1, color = col, alpha = isGray)) + geom_point() + theme_bw() + scale_x_continuous(name = "Log2FC", expand = c(0.01, 0.01)) + scale_y_continuous(name = "-Log10(Adjusted p)", expand = c(0.03, 0.01)) + scale_color_identity() + geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray40") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") + scale_alpha_manual(values = c(0.8, 0.4)) + guides(alpha = "none", size = guide_legend(title="% Experessed")) + geom_text_repel(data = bulk.deg[which(!bulk.deg$isGray),], aes(label = label2)))
dev.off()

#*******************************************************************************
# Initial clustering w plk60  ==================================================
#*******************************************************************************
dir_of_sr_dirs = "~/scratch/brain/bs/" # Folder where all the individual samples are kept
counts_list = list()
counts_list[["plk7_p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_clp12_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["plk7_p34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_clp34_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["plk7_c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_cnt12_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["plk7_c34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_cnt34_bcl_cr6/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_p12_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_p34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_p34_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_c12_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_c34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_c34_bcl/outs/filtered_feature_bc_matrix/"))

gtf = read.delim("~/scratch/m_zebra_ref/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = stringr::str_replace(mito.genes, "_", "-")

objs = list()
for (this.name in names(counts_list)) { 
  objs[[this.name]] = CreateSeuratObject(counts_list[[this.name]])
  objs[[this.name]]$sample = this.name
  objs[[this.name]]$exp = stringr::str_split(this.name, "_")[[1]][1]
  objs[[this.name]]$cond   = "plk"
  if (grepl("_c", this.name)) { objs[[this.name]]$cond = "con"}
  objs[[this.name]]$pct.mt = colSums(objs[[this.name]]@assays$RNA@counts[mito.genes,]) / objs[[this.name]]$nCount_RNA
}

for (s in names(objs)) { objs[[s]] = RenameCells(objs[[s]], paste0(s)) }

# Check the quality
for (i in names(objs)) {
  obj = objs[[i]]
  plot1 = ggplot(data.frame(nCount = obj$nCount_RNA), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
  plot2 = VlnPlot(obj, features = "nCount_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  
  plot3 = ggplot(data.frame(nFeature = obj$nFeature_RNA), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
  plot4 = VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  pdf(paste0("~/scratch/d_tooth/results/plk60_", i, "_quality.pdf"), width = 8, height = 8)
  print(cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2))
  dev.off()
  print(paste0("# Cells UMI > 300 = ", length(which(obj$nCount_RNA > 300))))
  print(paste0("# Cells Genes > 300 = ", length(which(obj$nFeature_RNA > 300))))
}


# Remove cells
removed.df = data.frame()
for (i in names(objs)) {
  # objs[[i]]$pct.mt = colSums(objs[[i]]@assays$RNA@counts[mito.genes,]) / objs[[i]]$nCount_RNA
  # removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$pct.mt >= 0.05)], reason = "dead (% MT)"))
  if (length(which(objs[[i]]$nFeature_RNA >= 2500)) > 0) { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA >= 2500)], reason = "doublet", sample = i))       }
  if (length(which(objs[[i]]$nFeature_RNA <= 250)) > 0)  { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA <= 250)], reason = "dead (# Genes)", sample = i)) }
  removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(!colnames(objs[[i]]) %in% removed.df$cell)], reason = "kept", sample = i))
}

pdf(paste0("~/scratch/d_tooth/results/plk60_nuclei_kept.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar() + xlab("Sample") + ylab("Number of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()
pdf(paste0("~/scratch/d_tooth/results/plk60_nuclei_kept_pct.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar(position = "fill") + xlab("Sample") + ylab("Proportion of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL)) + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()

plk = merge(objs[["plk7_p12"]], objs[-1])
# plk$pct.mt = colSums(plk@assays$RNA@counts[mito.genes,]) / plk$nCount_RNA

plk = subset(plk, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & pct.mt < 0.05)
print(paste("Number of Cells in plk After Filterning:", ncol(plk)))
plk = NormalizeData(plk, normalization.method = "LogNormalize", scale.factor = 10000)
plk = SCTransform(plk, verbose = TRUE)
plk@active.assay = "SCT"
plk = RunPCA(plk)
plk = RunUMAP(plk, reduction = "pca", dims = 1:25)
plk = FindNeighbors(plk, reduction="umap", dims = 1:2)
plk = FindClusters(plk, resolution = .25)

pdf("~/scratch/d_tooth/results/plk60_elbow.pdf", width = 5, height = 4)
print(ElbowPlot(plk, ndim = 50))
dev.off()

# Plot clusters and plot them by sample
# pdf(paste0("~/scratch/em/results/em_cluster.pdf"), width = 3.5, height = 3.5)
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cluster_no_regress.png"), width = 1800, height = 1800, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1) + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cluster_no_regress_split.png"), width = 3200, height = 2000, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1, split.by = "sample", ncol = 4) + coord_fixed() + NoLegend())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cluster_split_exp.png"), width = 3200, height = 1800, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1, split.by = "exp", ncol = 2) + coord_fixed() + NoLegend())
dev.off()

# Plot the # of Genes per nuclei to see if any clusters have particularly low quality
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_nFeature_umap.png"), width = 1800, height = 1800, res = 200)
print(FeaturePlot(plk, "nFeature_RNA", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cyto_umap.png"), width = 1800, height = 1800, res = 200)
print(FeaturePlot(plk, "cyto", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_bad_nFeature_umap.png"), width = 1800, height = 1800, res = 200)
plk$inverse_nFeature_RNA = 1/ plk$nFeature_RNA
print(FeaturePlot(plk, "inverse_nFeature_RNA", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()

# Check for cluster proportion differences
sample.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$sample))
s.clust.sdenom = sweep(sample.cluster.mat, 2, colSums(sample.cluster.mat), "/")
s.clust.cdenom = sample.cluster.mat / rowSums(sample.cluster.mat)
s.clust.sdenom.melt = reshape2::melt(s.clust.sdenom)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.sdenom.melt) = c("Cluster", "Sample", "value")
colnames(s.clust.cdenom.melt) = c("Cluster", "Sample", "value")

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cluster_prop_by_sample.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Sample)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_sample_prop_by_cluster.png"), width = 1000, height = 1400, res = 200)
s.clust.sdenom.melt$Cluster = factor(s.clust.sdenom.melt$Cluster)
ggplot(s.clust.sdenom.melt, aes(x = Sample, y = value, fill = Cluster, color = Cluster)) + geom_bar(alpha = 0.3, stat = 'identity') + theme_bw() + guides(color = "none")
dev.off()

st.clust.deg$label2 = st.clust.deg$label
st.clust.deg$label2[which(st.clust.deg$label == "LOC101482992")] = "LOC101482992 (ACPT)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101476651")] = "LOC101476651 (KRT24)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101483408")] = "LOC101483408 (SLC5A7)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101485955")] = "LOC101485955 (krt5)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101482992")] = "LOC101482992 (ACPT)"
st.clust.deg$label2 = factor(st.clust.deg$label2, unique(st.clust.deg$label2))
ggplot(st.clust.deg, aes(x = label2, y = -log10(p_val_adj), fill = avg_logFC)) + geom_bar(stat = "identity", width = 0.5) + coord_flip() + theme_bw() + scale_y_continuous(expand = c(0,0)) +  facet_wrap(~ cluster, scales = "free_y") + xlab("") + scale_fill_viridis()

# Check for cluster proportion differences
plk.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$cond))
s.clust.pdenom = sweep(plk.cluster.mat, 2, colSums(plk.cluster.mat), "/")
s.clust.cdenom = plk.cluster.mat / rowSums(plk.cluster.mat)
s.clust.pdenom.melt = reshape2::melt(s.clust.pdenom)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.pdenom.melt) = c("Cluster", "Condition", "value")
colnames(s.clust.cdenom.melt) = c("Cluster", "Condition", "value")

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cluster_prop_by_plk.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Condition)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_plk_prop_by_cluster.png"), width = 1000, height = 1400, res = 200)
s.clust.pdenom.melt$Cluster = factor(s.clust.pdenom.melt$Cluster)
ggplot(s.clust.pdenom.melt, aes(x = Condition, y = value, color = Cluster, fill = Cluster)) + geom_bar(alpha = 0.3, stat = 'identity') + theme_bw() + guides(color = "none")
dev.off()

plk.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$exp))
s.clust.cdenom = plk.cluster.mat / rowSums(plk.cluster.mat)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.cdenom.melt) = c("Cluster", "Experiment", "value")
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plk60_cluster_prop_by_exp.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Experiment)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()

#*******************************************************************************
# Plk + Div + TT  ==============================================================
#*******************************************************************************
tt = readRDS(paste0(data_dir, "tt_072722.rds"))
b10 = subset(tt, cells = colnames(tt)[which(tt$nFeature_RNA > 250 & tt$nFeature_RNA < 2500 & tt$pct.mt < 0.05 & tt$sample == "B10")])
merged = merge(plk, c(div, b10))

#*******************************************************************************
# Plk + Div  ========================================================
#*******************************************************************************
# /storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/conda_envs/r6/bin/R
.libPaths("/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/conda_envs/r6/lib/R/library/")
merged = merge(plk, div)
merged = SCTransform(merged, assay = "RNA", vars.to.regress = "sample", verbose = T)
merged = RunPCA(merged, dim=50, verbose=F)
merged = RunUMAP(merged, dims=1:50, min.dist=0.2, n.neighbors=30, n.epochs=1000, metric="euclidean")
merged = FindNeighbors(merged, reduction="umap", k.param=30, dims=1:2, n.trees=500, prune.SNN=0)
merged = FindClusters(merged, resolution=0.6, algorithm=2)
saveRDS(merged, "~/scratch/d_tooth/data/plk_div_081823.rds")
merged$project = colnames(merged) %in% colnames(plk)
merged$project = plyr::revalue(as.character(merged$project), replace = c("TRUE" = "plk", "FALSE" = "div"))
Idents(merged) = merged$project
pdf("~/scratch/d_tooth/results/div_plk_sct_merge.pdf", width = 7, height = 7)
DimPlot(merged) + ggplot2::scale_color_manual(values=c("gold", "skyblue")) + ggplot2::coord_fixed()
dev.off()

# Simple Project
plk2 = readRDS("~/scratch/d_tooth/data/plkall_081823.rds")
merged = ProjectUMAP(query = div, query.reduction="pca", reference = plk2, reference.reduction = "pca", reduction.model="umap")
pdf("~/scratch/d_tooth/results/div_plk_project_0907223.pdf", width = 7, height = 7)
DimPlot(merged, reduction = "ref.umap") + coord_fixed()
dev.off()
div_coords = data.frame(merged@reductions$ref.umap@cell.embeddings); div_coords$project = "div";
plk_coords = data.frame(plk2@reductions$umap@cell.embeddings);       plk_coords$project = "plk";
colnames(div_coords) = colnames(plk_coords)
all_coords = rbind(plk_coords, div_coords)
pdf("~/scratch/d_tooth/results/div_plk_project_both_0907223.pdf", width = 14, height = 7)
p1 = ggplot(div_coords, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 0.1, color = "skyblue") + theme_classic() + coord_fixed()
p2 = ggplot(plk_coords, aes(x = UMAP_1, y = UMAP_2))    + geom_point(size = 0.1, color = "gold") + theme_classic() + coord_fixed()
print(cowplot::plot_grid(plotlist = list(p1, p2)))
dev.off()
pdf("~/scratch/d_tooth/results/div_plk_project_both2_0907223.pdf", width = 7, height = 7)
ggplot(all_coords, aes(x = UMAP_1, y = UMAP_2, color = project)) + geom_point(size = 0.1) + theme_classic() + scale_color_manual(values=c("skyblue", "gold")) + coord_fixed()
dev.off()

# Project Anchors
plk2$plk_clusters = plk2$seurat_clusters
pancreas.anchors <- FindIntegrationAnchors(object.list = list(plk2, div), dims = 1:30)
integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30, normalization.method = "SCT")
pancreas.anchors2 <- FindTransferAnchors(reference = integrated, query = div, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions <- TransferData(anchorset = pancreas.anchors2, refdata = integrated$plk_clusters, dims = 1:30)
div <- AddMetaData(div, metadata = predictions)
div <- MapQuery(anchorset = pancreas.anchors2, reference = integrated, query = div, refdata = list(plk_clusters = "plk_clusters"), reference.reduction = "pca", reduction.model = "umap")

pancreas.anchors2 <- FindTransferAnchors(reference = plk2, query = div, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions <- TransferData(anchorset = pancreas.anchors2, refdata = plk2$plk_clusters, dims = 1:30)
div <- AddMetaData(div, metadata = predictions)
div <- MapQuery(anchorset = pancreas.anchors2, reference = plk2, query = div, refdata = list(plk_clusters = "plk_clusters"), reference.reduction = "pca", reduction.model = "umap")
Idents(div) = factor(div$predicted.id, levels = levels(plk2$plk_clusters))
pdf("~/scratch/d_tooth/results/div_plk_project_std_both_0907223.pdf", width = 20, height = 6)
p1 = DimPlot(plk2, label = T) + coord_fixed()
p2 = DimPlot(div, reduction = "ref.umap", label = T) + coord_fixed()
print(cowplot::plot_grid(plotlist = list(p1, p2)))
dev.off()

ca=subset(div, cells = colnames(div)[which(div$species_abr == "ca")])
mz=subset(div, cells = colnames(div)[which(div$species_abr == "mz")])
pc=subset(div, cells = colnames(div)[which(div$species_abr == "pc")])

pdf("~/scratch/d_tooth/results/div_plk_project_std_both3_0907223.pdf", width = 14, height = 12)
p1 = DimPlot(plk2, label = T) + coord_fixed() + NoLegend() + ggtitle("Plk/Ctrl")
p2 = DimPlot(ca, reduction = "ref.umap", label = T) + coord_fixed() + NoLegend() + ggtitle("CA")
p3 = DimPlot(mz, reduction = "ref.umap", label = T) + coord_fixed() + NoLegend() + ggtitle("MZ")
p4 = DimPlot(pc, reduction = "ref.umap", label = T) + coord_fixed() + NoLegend() + ggtitle("PC")
print(cowplot::plot_grid(plotlist = list(p2, p3, p4, p1), ncol = 2))
dev.off()

nperm=10000
this_obj = pc
ca_list = parallel::mclapply(1:nperm, function(x) table(sample(plk2$plk_clusters, ncol(this_obj))), mc.cores=24)
ca_df = do.call('rbind', ca_list)
ca_real = data.frame(table(factor(this_obj$predicted.id, levels = levels(plk2$plk_clusters))))
ca_n_greater = t(ca_df) > ca_real[,2]
ca_real$perm_mean = colMeans(ca_df)
ca_real$n_greater = rowSums(ca_n_greater)
ca_real$p = 1 - (ca_real$n_greater / nperm)
colnames(ca_real) = c("plk_cluster", "real_pc_num_nuclei", "perm_mean_num_nuclei", "n_perm_greater_than_real", "p_value")
ca_real[which(ca_real$p_value < 0.05),]

#*******************************************************************************
# ChooseR Modification  ========================================================
#*******************************************************************************
meds_backup = meds

meds = meds_backup
meds = meds[which( meds$n_neighbor < 50),]
threshold <- max(meds$low_med)
meds %>% dplyr::filter(med >= threshold) %>% dplyr::arrange(n_clusters) %>% tail(n = 1)
meds$id = paste(meds$min_dist, meds$n_neighbor, meds$res, sep = "_")
scores$id = paste(scores$min_dist, scores$n_neighbor, scores$res, sep = "_")
meds = meds %>% dplyr::arrange(n_clusters)
choice <- as.character(meds %>% dplyr::filter(med >= threshold) %>% dplyr::arrange(n_clusters) %>% tail(n = 1) %>% dplyr::pull(id))
print("The best clustering paramaters are:")
print(choice)

git_dir_path = "/storage/home/hcoda1/6/ggruenhagen3/scratch/st/cichlid_st/"
obj_path = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plkall_filtered_unclustered.rds"
results_path = paste0("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/chooser_plk/examples/")
# Load Libraries
source(paste0(git_dir_path, "clustering/chooseR_code/pipeline.R"))
library(Seurat)
library(ggplot2)
library(parallel)
`%>%` <- magrittr::`%>%`
obj = readRDS(obj_path)
# Define the number of PCs to use, and which assay and reduction to use.
npcs <- 2
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
assay <- "SCT"
reduction <- "umap"
batch_num = 2
min_dist = 0.2
n_neighbors = 30
res = 0.6
this.obj <- find_clusters(
  obj,
  reduction = reduction,
  assay = assay,
  resolution = res,
  npcs = npcs,
  batch.num = batch_num,
  min.dist = min_dist,
  n.neighbors = n_neighbors
)
pdf("~/scratch/d_tooth/results/div_test1.pdf", width = 7, height = 7)
DimPlot(this.obj, label = T) + coord_fixed()
dev.off()

#*******************************************************************************
# Initial clustering w All Experiments  ========================================
#*******************************************************************************
dir_of_sr_dirs = "~/scratch/brain/bs/" # Folder where all the individual samples are kept
counts_list = list()
counts_list[["plk7_p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_clp12_bcl_cr7/outs/filtered_feature_bc_matrix/"))
counts_list[["plk7_p34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_clp34_bcl_cr7/outs/filtered_feature_bc_matrix/"))
counts_list[["plk7_c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_cnt12_bcl_cr7/outs/filtered_feature_bc_matrix/"))
counts_list[["plk7_c34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS15/JTS15_cnt34_bcl_cr7/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_p12_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_p34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_p34_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_c12_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk60_c34"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS16/JTS16_c34_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk1_p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS17/JTS17_Plk01_1_2_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk1_c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS17/JTS17_Cont01_1_2_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk3_p12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS17/JTS17_Plk03_1_2_bcl/outs/filtered_feature_bc_matrix/"))
counts_list[["plk3_c12"]] = Read10X(paste0(dir_of_sr_dirs, "/JTS17/JTS17_Cont03_1_2_bcl/outs/filtered_feature_bc_matrix/"))

gtf = read.delim("~/scratch/m_zebra_ref/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = stringr::str_replace(mito.genes, "_", "-")

objs = list()
for (this.name in names(counts_list)) { 
  objs[[this.name]] = CreateSeuratObject(counts_list[[this.name]])
  objs[[this.name]]$sample = this.name
  objs[[this.name]]$exp = stringr::str_split(this.name, "_")[[1]][1]
  objs[[this.name]]$cond   = "plk"
  if (grepl("_c", this.name)) { objs[[this.name]]$cond = "con"}
  objs[[this.name]]$pct.mt = colSums(objs[[this.name]]@assays$RNA@counts[mito.genes,]) / objs[[this.name]]$nCount_RNA
}

for (s in names(objs)) { objs[[s]] = RenameCells(objs[[s]], paste0(s)) }

# Check the quality
for (i in names(objs)) {
  obj = objs[[i]]
  plot1 = ggplot(data.frame(nCount = obj$nCount_RNA), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
  plot2 = VlnPlot(obj, features = "nCount_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  
  plot3 = ggplot(data.frame(nFeature = obj$nFeature_RNA), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
  plot4 = VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
  pdf(paste0("~/scratch/d_tooth/results/plkall_", i, "_quality.pdf"), width = 8, height = 8)
  print(cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2))
  dev.off()
  print(paste0("# Cells UMI > 300 = ", length(which(obj$nCount_RNA > 300))))
  print(paste0("# Cells Genes > 300 = ", length(which(obj$nFeature_RNA > 300))))
}


# Remove cells
removed.df = data.frame()
for (i in names(objs)) {
  # objs[[i]]$pct.mt = colSums(objs[[i]]@assays$RNA@counts[mito.genes,]) / objs[[i]]$nCount_RNA
  # removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$pct.mt >= 0.05)], reason = "dead (% MT)"))
  if (length(which(objs[[i]]$nFeature_RNA >= 2500)) > 0) { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA >= 2500)], reason = "doublet", sample = i))       }
  if (length(which(objs[[i]]$nFeature_RNA <= 250)) > 0)  { removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(objs[[i]]$nFeature_RNA <= 250)], reason = "dead (# Genes)", sample = i)) }
  removed.df = rbind(removed.df, data.frame(cell = colnames(objs[[i]])[which(!colnames(objs[[i]]) %in% removed.df$cell)], reason = "kept", sample = i))
}

pdf(paste0("~/scratch/d_tooth/results/plkall_nuclei_kept.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar() + xlab("Sample") + ylab("Number of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()
pdf(paste0("~/scratch/d_tooth/results/plkall_nuclei_kept_pct.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar(position = "fill") + xlab("Sample") + ylab("Proportion of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL)) + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()

ncount.df = lapply(names(objs), function(x) data.frame(ncount = objs[[x]]$nCount_RNA, my.facet = x))
ncount.df = data.frame(data.table::rbindlist(ncount.df))
nfeature.df = lapply(names(objs), function(x) data.frame(ncount = objs[[x]]$nFeature_RNA, my.facet = x))
nfeature.df = data.frame(data.table::rbindlist(nfeature.df))
ncount.df = ncount.df[which(ncount.df$ncount < 3000),]
nfeature.df = nfeature.df[which(nfeature.df$ncount < 3000),]
pdf("~/scratch/d_tooth/results/plkall_ncount_density.pdf", height = 6, width = 8, onefile = F)
print(ggplot(ncount.df, aes(x = ncount, fill = my.facet)) + geom_density(alpha = 0.5, position = "identity", bins = 100) + theme_classic() + ylab("Nuclei Density") + xlab("Number of UMIs per Nuclei") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + facet_wrap(~my.facet, ncol=4))
dev.off()
pdf("~/scratch/d_tooth/results/plkall_nfeature_density.pdf", height = 6, width = 8, onefile = F)
print(ggplot(nfeature.df, aes(x = ncount, fill = my.facet)) + geom_density(alpha = 0.5, position = "identity", bins = 100) + theme_classic() + geom_vline(xintercept = 250, linetype="dashed") + ylab("Nuclei Density") + xlab("Number of Genes per Nucleus") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + facet_wrap(~my.facet, ncol=4))
dev.off()

plk_samples_to_merge = names(objs)[2:length(names(objs))]
plk_samples_to_merge = plk_samples_to_merge[which(!plk_samples_to_merge %in% c("plk7_p34", "plk7_c34"))]
plk = merge(objs[["plk7_p12"]], objs[plk_samples_to_merge])
# plk$pct.mt = colSums(plk@assays$RNA@counts[mito.genes,]) / plk$nCount_RNA

plk = subset(plk, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & pct.mt < 0.05)
print(paste("Number of Cells in plk After Filterning:", ncol(plk)))
plk = NormalizeData(plk, normalization.method = "LogNormalize", scale.factor = 10000)
plk = SCTransform(plk, verbose = TRUE)
plk@active.assay = "SCT"
plk = RunPCA(plk)
plk = RunUMAP(plk, reduction = "pca", dims = 1:25)
plk = FindNeighbors(plk, reduction="umap", dims = 1:2)
plk = FindClusters(plk, resolution = .25)

pdf("~/scratch/d_tooth/results/plkall_elbow.pdf", width = 5, height = 4)
print(ElbowPlot(plk, ndim = 50))
dev.off()

# Plot clusters and plot them by sample
# pdf(paste0("~/scratch/em/results/em_cluster.pdf"), width = 3.5, height = 3.5)
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cluster_no_regress.png"), width = 1800, height = 1800, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1) + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cluster_no_regress_split.png"), width = 3200, height = 3200, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1, split.by = "sample", ncol = 4) + coord_fixed() + NoLegend())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cluster_split_exp.png"), width = 3200, height = 1800, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 1, split.by = "exp", ncol = 2) + coord_fixed() + NoLegend())
dev.off()

# Plot the # of Genes per nuclei to see if any clusters have particularly low quality
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_nFeature_umap.png"), width = 1800, height = 1800, res = 200)
print(FeaturePlot(plk, "nFeature_RNA", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cyto_umap.png"), width = 1800, height = 1800, res = 200)
print(FeaturePlot(plk, "cyto", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_bad_nFeature_umap.png"), width = 1800, height = 1800, res = 200)
plk$inverse_nFeature_RNA = 1/ plk$nFeature_RNA
print(FeaturePlot(plk, "inverse_nFeature_RNA", order = T, label = TRUE, pt.size = 1) + scale_color_viridis() + coord_fixed())
dev.off()

# Old annotations
plk7 = readRDS(paste0("~/scratch/d_tooth/data/", "plk_120922.rds"))
plk7_annot = read.csv("~/scratch/d_tooth/data/Plk7_Combined_Annotations.csv")
plk$plk7_annot = NA
plk$plk7_annot[paste0("plk7_",colnames(plk7))] = plk7_annot$old_label[match(plk7$seurat_clusters, plk7_annot$old_cluster)]
Idents(plk) = plk$plk7_annot
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_by_plk7_annot.png"), width = 2400, height = 2400, res = 200)
print(DimPlot(plk, reduction = "umap", label = TRUE, pt.size = 0.4, order = T) + coord_fixed())
dev.off()

# Check for cluster proportion differences
sample.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$sample))
s.clust.sdenom = sweep(sample.cluster.mat, 2, colSums(sample.cluster.mat), "/")
s.clust.cdenom = sample.cluster.mat / rowSums(sample.cluster.mat)
s.clust.sdenom.melt = reshape2::melt(s.clust.sdenom)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.sdenom.melt) = c("Cluster", "Sample", "value")
colnames(s.clust.cdenom.melt) = c("Cluster", "Sample", "value")

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cluster_prop_by_sample.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Sample)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_sample_prop_by_cluster.png"), width = 1000, height = 1400, res = 200)
s.clust.sdenom.melt$Cluster = factor(s.clust.sdenom.melt$Cluster)
ggplot(s.clust.sdenom.melt, aes(x = Sample, y = value, fill = Cluster, color = Cluster)) + geom_bar(alpha = 0.3, stat = 'identity') + theme_bw() + guides(color = "none")
dev.off()

# Check for cluster proportion differences
plk.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$cond))
s.clust.pdenom = sweep(plk.cluster.mat, 2, colSums(plk.cluster.mat), "/")
s.clust.cdenom = plk.cluster.mat / rowSums(plk.cluster.mat)
s.clust.pdenom.melt = reshape2::melt(s.clust.pdenom)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.pdenom.melt) = c("Cluster", "Condition", "value")
colnames(s.clust.cdenom.melt) = c("Cluster", "Condition", "value")

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cluster_prop_by_plk.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Condition)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_plk_prop_by_cluster.png"), width = 1000, height = 1400, res = 200)
s.clust.pdenom.melt$Cluster = factor(s.clust.pdenom.melt$Cluster)
ggplot(s.clust.pdenom.melt, aes(x = Condition, y = value, color = Cluster, fill = Cluster)) + geom_bar(alpha = 0.3, stat = 'identity') + theme_bw() + guides(color = "none")
dev.off()

# Check for plk vs control cluster proportion differences by experiment
t_order = c (1, 9, 12, 19, 26, 27, 34, 10, 14, 15, 16, 17, 28, 31, 35, 0, 4, 8, 13, 20, 32, 2, 3, 5, 11, 21, 23, 22, 29, 7, 25, 38, 41, 42, 45, 6, 18, 30, 39, 37, 46, 48, 33, 36, 44, 24, 40, 43, 47, 49)
for (this_exp in unique(plk$exp)) {
  this_cells = colnames(plk)[which(plk$exp == this_exp)]
  plk.cluster.mat = as.matrix(table(plk$seurat_clusters[this_cells], plk$cond[this_cells]))
  s.clust.cdenom = plk.cluster.mat / rowSums(plk.cluster.mat)
  s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
  colnames(s.clust.pdenom.melt) = c("Cluster", "Condition", "value")
  colnames(s.clust.cdenom.melt) = c("Cluster", "Condition", "value")
  s.clust.cdenom.melt$Cluster = factor(s.clust.cdenom.melt$Cluster, levels = sort(unique(s.clust.cdenom.melt$Cluster)))
  s.clust.cdenom.melt$num = reshape2::melt(plk.cluster.mat)$value
  s.clust.cdenom.melt$label_y = s.clust.cdenom.melt$value / 2
  s.clust.cdenom.melt$label_y[which(s.clust.cdenom.melt$Condition == "con")] = 1 - (s.clust.cdenom.melt$value[which(s.clust.cdenom.melt$Condition == "con")]/2)
  # s.clust.cdenom.melt2 = s.clust.cdenom.melt %>% arrange(Cluster, Condition) %>% group_by(Cluster) %>% mutate(label_y = cumsum(value))
  s.clust.cdenom.melt$Cluster = factor(s.clust.cdenom.melt$Cluster, levels = t_order)
  
  Cairo::Cairo(paste0("~/research/tooth/results/plkall_", this_exp, "_cluster_prop_by_plk.png"), width = 3000, height = 800, res = 200)
  print(ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Condition)) + geom_bar(stat = 'identity') + theme_bw() + scale_y_continuous(expand = c(0,0)) + geom_text(aes(y = label_y, label = num), color = "white"))
  dev.off()
}

plk.cluster.mat = as.matrix(table(plk$seurat_clusters, plk$exp))
s.clust.cdenom = plk.cluster.mat / rowSums(plk.cluster.mat)
s.clust.cdenom.melt = reshape2::melt(s.clust.cdenom)
colnames(s.clust.cdenom.melt) = c("Cluster", "Experiment", "value")
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkall_cluster_prop_by_exp.png"), width = 2400, height = 800, res = 200)
ggplot(s.clust.cdenom.melt, aes(x = Cluster, y = value, fill = Experiment)) + geom_bar(stat = 'identity') + theme_bw()
dev.off()

st.clust.deg$label2 = st.clust.deg$label
st.clust.deg$label2[which(st.clust.deg$label == "LOC101482992")] = "LOC101482992 (ACPT)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101476651")] = "LOC101476651 (KRT24)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101483408")] = "LOC101483408 (SLC5A7)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101485955")] = "LOC101485955 (krt5)"
st.clust.deg$label2[which(st.clust.deg$label == "LOC101482992")] = "LOC101482992 (ACPT)"
st.clust.deg$label2 = factor(st.clust.deg$label2, unique(st.clust.deg$label2))
ggplot(st.clust.deg, aes(x = label2, y = -log10(p_val_adj), fill = avg_logFC)) + geom_bar(stat = "identity", width = 0.5) + coord_flip() + theme_bw() + scale_y_continuous(expand = c(0,0)) +  facet_wrap(~ cluster, scales = "free_y") + xlab("") + scale_fill_viridis()

#*******************************************************************************
# Soup =========================================================================
#*******************************************************************************
soup_path_df = data.frame(sample = c("plk7_p12", "plk7_c12", "plk60_p12", "plk60_c12", "plk60_p34", "plk60_c34", "plk1_p12", "plk1_c12", "plk3_p12", "plk3_c12"),
                          path = c("JTS15/JTS15_clp12_bcl_cr7", "JTS15/JTS15_cnt12_bcl_cr7", "JTS16/JTS16_p12_bcl", "JTS16/JTS16_c12_bcl", "JTS16/JTS16_p34_bcl", "JTS16/JTS16_c34_bcl", "JTS17/JTS17_Plk01_1_2_bcl", "JTS17/JTS17_Cont01_1_2_bcl", "JTS17/JTS17_Plk03_1_2_bcl", "JTS17/JTS17_Cont03_1_2_bcl"))
for (i in 1:nrow(soup_path_df)) {
  soup2 = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_pred_2.tsv"), data.table = F)
  soup3 = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_pred_3.tsv"), data.table = F)
  soup2_assign = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_dbl_2.tsv"), data.table = F)
  soup3_assign = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_dbl_3.tsv"), data.table = F)
  soup2$norm2 = as.numeric(! soup2$V2)
  soup2$norm3 = soup3$V2
  same_pct = length(which(soup2$V2 == soup2$norm3)) / nrow(soup2)
  if (same_pct > 0.5) { soup2$norm2 = soup2$V2 }
  soup2$status2 = soup2_assign$status
  soup2$status3 = soup3_assign$status
  soup2$agree = soup2$norm2 == soup2$norm3
  soup2$best = NA
  soup2$best[which(soup2$agree)] = soup2$norm2[which(soup2$agree)] 
  soup2$best[which(!soup2$agree & soup2$status2 == "singlet" & soup2$status3 != "singlet")] = soup2$norm2[which(!soup2$agree & soup2$status2 == "singlet" & soup2$status3 != "singlet")]
  soup2$best[which(!soup2$agree & soup2$status3 != "singlet" & soup2$status3 == "singlet")] = soup2$norm2[which(!soup2$agree & soup2$status2 != "singlet" & soup2$status3 == "singlet")]
  print(paste0( soup_path_df$path[i], "==> % of agreement: ", length(which( soup2$agree )) / nrow(soup2) ))
  print(paste0( soup_path_df$path[i], "==> % of unassigned: ", length(which( is.na(soup2$best) )) / nrow(soup2) ))
  system(paste0("ls -lh ~/scratch/brain/bs/", soup_path_df$path[i], "/outs/filtered_final.bam"))
  write.csv(soup2, paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_assign.csv"))
  print ("--------------------------------------------------------------------------------------------------------")
}
# ref = Matrix::readMM(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_ref.mtx"))
# alt = Matrix::readMM(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_alt.mtx"))

for (j in seq(1, nrow(soup_path_df), by = 2)) {
  vcf_plk = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[j],   "/outs/cluster_genotypes_filter.vcf"), skip = 1710, data.table = F)
  vcf_cnt = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[j+1], "/outs/cluster_genotypes_filter.vcf"), skip = 1710, data.table = F)
  vcf_plk$id = paste0(vcf_plk[,1], "_", vcf_plk[,2])
  vcf_cnt$id = paste0(vcf_cnt[,1], "_", vcf_cnt[,2])
  vcf_merge = merge(vcf_plk, vcf_cnt, by = "id", suffixes = c("_plk", "_cnt"))
  vcf_merge$plk_0_gt = substr(vcf_merge[,"0_plk"], 1, 3)
  vcf_merge$plk_1_gt = substr(vcf_merge[,"1_plk"], 1, 3)
  vcf_merge$cnt_0_gt = substr(vcf_merge[,"0_cnt"], 1, 3)
  vcf_merge$cnt_1_gt = substr(vcf_merge[,"1_cnt"], 1, 3)
  same_pct_0 = length(which(vcf_merge$plk_0_gt == vcf_merge$cnt_0_gt)) / nrow(vcf_merge)
  dif_pct_0  = length(which(vcf_merge$plk_0_gt == vcf_merge$cnt_1_gt)) / nrow(vcf_merge)
  if (same_pct_0 > 0.5) { print(paste0("Soup 0 in Plk is Soup 0 in Cnt: ", same_pct_0, "% agreement")); flop_0 = F}
  if (dif_pct_0  > 0.5) { print(paste0("Soup 0 in Plk is Soup 1 in Cnt: ", dif_pct_0,  "% agreement")); flop_0 = T}
  
  same_pct_1 = length(which(vcf_merge$plk_1_gt == vcf_merge$cnt_1_gt)) / nrow(vcf_merge)
  dif_pct_1  = length(which(vcf_merge$plk_1_gt == vcf_merge$cnt_0_gt)) / nrow(vcf_merge)
  if (same_pct_1 > 0.5) { print(paste0("Soup 1 in Plk is Soup 1 in Cnt: ", same_pct_1, "% agreement")); flop_1 = F}
  if (dif_pct_1  > 0.5) { print(paste0("Soup 1 in Plk is Soup 0 in Cnt: ", dif_pct_1,  "% agreement")); flop_1 = T}
  if (flop_0 != flop_1) { print("We're in trouble!!!") }
  if (flop_0 && flop_1) { print("Both in agreement to flop")}
  print ("--------------------------------------------------------------------------------------------------------")
  
  this_pool = stringr::str_split(soup_path_df$sample[j], "_")[[1]][1]
  soup_plk = read.csv(paste0("~/scratch/brain/bs/", soup_path_df$path[j],   "/outs/soup_assign.csv"))
  soup_cnt = read.csv(paste0("~/scratch/brain/bs/", soup_path_df$path[j+1], "/outs/soup_assign.csv"))
  soup_plk$X = NULL; soup_plk$cond = "plk"; soup_plk$pool = this_pool; soup_plk$sample = soup_path_df$sample[j];
  soup_cnt$X = NULL; soup_cnt$cond = "con"; soup_cnt$pool = this_pool; soup_cnt$sample = soup_path_df$sample[j+1];
  soup_plk$subject = as.numeric(soup_plk$best)+1
  soup_cnt$subject = as.numeric(soup_cnt$best)+1
  if (flop_0 && flop_1) { soup_plk$subject = plyr::revalue(as.character(soup_plk$subject), c("1" = "2", "2" = "1")) }
  write.csv(soup_plk, paste0("~/scratch/d_tooth/data/plk_soup/", soup_path_df$sample[j], "_soup.csv"))
  write.csv(soup_cnt, paste0("~/scratch/d_tooth/data/plk_soup/", soup_path_df$sample[j+1], "_soup.csv"))
  plk$subject_num[which(plk$sample == soup_path_df$sample[j])]   = soup_plk$subject[match(colnames(plk)[which(plk$sample == soup_path_df$sample[j])],   paste0(soup_path_df$sample[j],   "_", soup_plk$V1))]
  plk$subject_num[which(plk$sample == soup_path_df$sample[j+1])] = soup_cnt$subject[match(colnames(plk)[which(plk$sample == soup_path_df$sample[j+1])], paste0(soup_path_df$sample[j+1], "_", soup_cnt$V1))]
}
plk$subject = paste0(plk$exp, "_", plk$subject_num)
plk$subject[which(plk$sample %in% c("plk60_c34", "plk60_p34"))] = paste0("plk60_", as.numeric(plk$subject_num[which(plk$sample %in% c("plk60_c34", "plk60_p34"))])+2)

#*******************************************************************************
# Soup Demux ===================================================================
#*******************************************************************************
dbest = data.table::fread("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_out.best", data.table = F)
dsing = data.table::fread("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_out.single", data.table = F)

dsing_df = data.frame(barcode = dbest[,1])
dsing_df$best1_prob = unlist(parallel::mclapply(dsing_df[,1], function(x) { this_v = dsing$POSTPRB[which(dsing[,1] == x)]; this_v[order(this_v, decreasing = T)][1] }))
dsing_df$best2_prob = unlist(parallel::mclapply(dsing_df[,1], function(x) { this_v = dsing$POSTPRB[which(dsing[,1] == x)]; this_v[order(this_v, decreasing = T)][2] }))
dsing_df$best3_prob = unlist(parallel::mclapply(dsing_df[,1], function(x) { this_v = dsing$POSTPRB[which(dsing[,1] == x)]; this_v[order(this_v, decreasing = T)][3] }))
ggplot(dsing_df, aes(x = best1_prob, y = best2_prob)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_post_prob_best_1_2.pdf", width = 5, height = 5)

ggplot(dbest, aes(x = SNG.LLK1, y = SNG.LLK2)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_ll_best_1_2.pdf", width = 5, height = 5)

ggplot(dbest, aes(x = SNG.LLK1, y = PRB.SNG1)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_test_best_1_2.pdf", width = 5, height = 5)

ggplot(dbest, aes(x = SNG.LLK1, y = N.SNP)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_test1_best_1_2.pdf", width = 5, height = 5)

dbest$ll_dif = dbest$SNG.LLK1 - dbest$SNG.LLK2
dbest$ll_pct = dbest$ll_dif / -dbest$SNG.LLK2
ggplot(dbest, aes(x = SNG.LLK1, y = ll_dif)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_test2_best_1_2.pdf", width = 5, height = 5)
ggplot(dbest, aes(x = PRB.SNG1, y = ll_dif)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_test3_best_1_2.pdf", width = 5, height = 5)
ggplot(dbest, aes(x = SNG.LLK1, y = ll_pct)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_test3_best_1_2.pdf", width = 5, height = 5)
ggplot(dbest, aes(x = PRB.SNG1, y = ll_pct)) + geom_point()
ggsave("~/scratch/brain/bs/JTS21/JTS21_ALT_3/outs/demux_test4_best_1_2.pdf", width = 5, height = 5)

length(which(dsing_df$best1_prob > 0.99)) / nrow(dsing_df)
length(which(dbest$PRB.SNG1 > 0.99 & dbest$ll_pct > 0.1 & dbest$N.SNP > 100)) / nrow(dsing_df)
dbest_good = dbest[which(dbest$PRB.SNG1 > 0.99 & dbest$ll_pct > 0.1 & dbest$N.SNP > 100),]

div = readRDS("~/scratch/d_tooth/data/div_optimal_081623.rds")
alt1_soup = read.csv(paste0("~/scratch/brain/bs/JTS21/JTS21_ALT_1/outs/soup_assign.csv"))
div$barcode = reshape2::colsplit(colnames(div), "_", c('1', '2'))[,2]
div$soup = NA
div$soup[which(div$sample == "alt1")] = alt1_soup$best[match(div$barcode[which(div$sample == "alt1")], alt1_soup$barcode)]
div$species_abr = plyr::revalue(as.character(div$soup), c("0" = "pc", "1" = "la", "2" = "ca"))
div$species_abr[which(div$sample == "alt3")] = dbest_good$SNG.1ST[match(div$barcode[which(div$sample == "alt3")], dbest_good$BARCODE)]

merged = readRDS("~/scratch/d_tooth/data/plk_div_081823.rds")
merged$soup = NA
merged$soup[colnames(div)] = paste0(div$sample, "_", div$species_abr)
merged$to_p = colnames(merged) %in% colnames(div)
merged$to_p[which(merged$to_p)] = merged$soup[which(merged$to_p)]
merged$to_p[which(merged$to_p == FALSE)] = "plk/ctrl data"
merged2 = subset(merged, cells = colnames(merged)[which(!merged$soup %in% c("alt1_NA", "alt3_NA") | !colnames(merged) %in% colnames(div))])
Idents(merged2) = factor(merged2$to_p, levels = c("plk/ctrl data", "alt3_ag", "alt3_mz", "alt3_np", "alt1_ca", "alt1_la", "alt1_pc"))
pdf("~/scratch/d_tooth/results/div_species_abr.pdf", width = 9, height = 9)
DimPlot(merged2, split.by = "to_p", ncol = 3) + coord_fixed() + theme_void() + scale_color_manual(values = c("gray60", scales::hue_pal()(6)))
dev.off()

#*******************************************************************************
# Soup Div =====================================================================
#*******************************************************************************
div$barcode = reshape2::colsplit(colnames(div), "_", c('1', '2'))[,2]
div$soup = NA
soup_path_df = data.frame(sample = c("alt1", "alt3"), path = c("JTS21/JTS21_ALT_1", "JTS21/JTS21_ALT_3"))
# soup_path_df = data.frame(sample = c("alt3"), path = c("JTS21/JTS21_ALT_3"))
for (i in 1:nrow(soup_path_df)) {
  # soup2 = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_pred_2.tsv"), data.table = F)
  # soup3 = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_pred_3.tsv"), data.table = F)
  # soup2_assign = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_dbl_2.tsv"), data.table = F)
  # soup3_assign = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_dbl_3.tsv"), data.table = F)
  soup2 = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_dbl_2.tsv"), data.table = F)
  soup3 = data.table::fread(paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_out_dbl_3.tsv"), data.table = F)
  
  soup2$raw_2 = soup2$assignment
  soup2$raw_3 = soup3$assignment
  match_table = unclass(table(soup2$raw_2, soup2$raw_3))
  match_table = match_table[which(nchar(rownames(match_table)) ==1), which(nchar(colnames(match_table)) ==1)]
  raw_2_to_3 = data.frame(raw_2 = 0:2, raw_3 = unname(apply(match_table, 1, which.max))-1)
  
  soup2$norm2 = soup2$raw_2
  soup2$norm3 = raw_2_to_3$raw_2[match(soup2$raw_3, raw_2_to_3$raw_3)]
  
  soup2$status2 = soup2$status
  soup2$status3 = soup3$status
  soup2$agree = soup2$norm2 == soup2$norm3
  
  soup2$best = NA
  soup2$best[which(soup2$agree)] = soup2$norm2[which(soup2$agree)] 
  soup2$best[which(!soup2$agree & soup2$status2 == "singlet" & soup2$status3 != "singlet")] = soup2$norm2[which(!soup2$agree & soup2$status2 == "singlet" & soup2$status3 != "singlet")]
  soup2$best[which(!soup2$agree & soup2$status3 != "singlet" & soup2$status3 == "singlet")] = soup2$norm2[which(!soup2$agree & soup2$status2 != "singlet" & soup2$status3 == "singlet")]
  print(paste0( soup_path_df$path[i], "==> % of agreement: ", length(which( soup2$agree )) / nrow(soup2) ))
  print(paste0( soup_path_df$path[i], "==> % of unassigned: ", length(which( is.na(soup2$best) )) / nrow(soup2) ))
  system(paste0("ls -lh ~/scratch/brain/bs/", soup_path_df$path[i], "/outs/filtered_final.bam"))
  write.csv(soup2, paste0("~/scratch/brain/bs/", soup_path_df$path[i], "/outs/soup_assign.csv"))
  div$soup[which(div$sample == soup_path_df$sample[i])] = soup2$best[match(div$barcode[which(div$sample == soup_path_df$sample[i])], soup2$barcode)]
  print ("--------------------------------------------------------------------------------------------------------")
}

merged$soup = NA
merged$soup[colnames(div)] = paste0(div$sample, "_", div$soup)
merged$to_p = colnames(merged) %in% colnames(div)
merged$to_p[which(merged$to_p)] = merged$soup[which(merged$to_p)]
merged$to_p[which(merged$to_p == FALSE)] = "plk/ctrl data"
Idents(merged) = merged$to_p
merged2 = subset(merged, cells = colnames(merged)[which(!merged$soup %in% c("alt1_NA", "alt3_NA"))])
pdf("~/scratch/d_tooth/results/div_individuals_split.pdf", width = 9, height = 9)
DimPlot(merged2, split.by = "to_p", ncol = 3) + coord_fixed() + theme_void()
dev.off()
pdf("~/scratch/d_tooth/results/div_individuals_split_large.pdf", width = 9, height = 9)
DimPlot(merged2, split.by = "to_p", ncol = 3, pt.size = 1.5) + coord_fixed() + theme_void()
dev.off()

pdf("~/scratch/d_tooth/results/div_individuals.pdf", width = 7, height = 7)
DimPlot(merged2, ncol = 3) + coord_fixed() + theme_void() + scale_color_manual(values = c("gray60", scales::hue_pal()(6)))
dev.off()

# *** Compare vcf from soup individuals to the vcf of species in the pool *** #
library("vcfR")
# Read soup vcf
soup_vcf2 = data.table::fread("~/scratch/brain/bs/JTS21/JTS21_ALT_1/outs/cluster_genotypes_2.vcf", data.table = F, skip = 1714)
soup_vcf3 = data.table::fread("~/scratch/brain/bs/JTS21/JTS21_ALT_1/outs/cluster_genotypes_3.vcf", data.table = F, skip = 1714)
# Read species vcf
species_vcf = data.table::fread("~/scratch/brain/bs/JTS21/JTS21_ALT_1/outs/cc_la_100.vcf", data.table = F, skip = 1714)

# soup2 and soup3 agree
soup_vcf2$id = paste(soup_vcf2[,1], soup_vcf2[,2], soup_vcf2[,4], soup_vcf2[,5], sep="_")
soup_vcf3$id = paste(soup_vcf3[,1], soup_vcf3[,2], soup_vcf3[,4], soup_vcf3[,5], sep="_")
soup_vcf = soup_vcf2
colnames(soup_vcf)[(ncol(soup_vcf)-3):(ncol(soup_vcf)-1)] = paste0("soup2_", colnames(soup_vcf)[(ncol(soup_vcf)-3):(ncol(soup_vcf)-1)])
soup_vcf[,c("soup3_0", "soup3_1", "soup3_2")] = soup_vcf3[match(soup_vcf2$id, soup_vcf3$id), c("2", "1", "0")]
soup_vcf$soup2_0_gt = reshape2::colsplit(soup_vcf$soup2_0, ":", c('1', '2'))[,1]
soup_vcf$soup2_1_gt = reshape2::colsplit(soup_vcf$soup2_1, ":", c('1', '2'))[,1]
soup_vcf$soup2_2_gt = reshape2::colsplit(soup_vcf$soup2_2, ":", c('1', '2'))[,1]
soup_vcf$soup3_0_gt = reshape2::colsplit(soup_vcf$soup3_0, ":", c('1', '2'))[,1]
soup_vcf$soup3_1_gt = reshape2::colsplit(soup_vcf$soup3_1, ":", c('1', '2'))[,1]
soup_vcf$soup3_2_gt = reshape2::colsplit(soup_vcf$soup3_2, ":", c('1', '2'))[,1]
soup_vcf_simple = soup_vcf[, c("id", "soup2_0_gt", "soup2_1_gt", "soup2_2_gt", "soup3_0_gt", "soup3_1_gt", "soup3_2_gt")]
length(which(soup_vcf_simple$soup2_0_gt == soup_vcf_simple$soup3_0_gt)) / nrow(soup_vcf_simple)

species_vcf$id = paste(species_vcf[,1], species_vcf[,2], species_vcf[,4], species_vcf[,5], sep="_")
soup_vcf_simple$cc = species_vcf[match(soup_vcf_simple$id, species_vcf$id), "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/cc"]
soup_vcf_simple$la = species_vcf[match(soup_vcf_simple$id, species_vcf$id), "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/la/la"]
length(which(!is.na(soup_vcf_simple$cc)))
head(soup_vcf_simple[which(!is.na(soup_vcf_simple$cc)),])

soup_vcf_simple$cc_gt = reshape2::colsplit(soup_vcf_simple$cc, ":", c('1', '2'))[,1]
soup_vcf_simple$la_gt = reshape2::colsplit(soup_vcf_simple$la, ":", c('1', '2'))[,1]

soup_vcf_simple2 = soup_vcf_simple[which(!is.na(soup_vcf_simple$cc) & !is.na(soup_vcf_simple$la)),]
soup_vcf_simple2$cc_ref_num = as.numeric(reshape2::colsplit(reshape2::colsplit(soup_vcf_simple2$cc, ",", c('1', '2'))[,1], ":", c('1', '2'))[,2])
soup_vcf_simple2$cc_alt_num = as.numeric(reshape2::colsplit(reshape2::colsplit(soup_vcf_simple2$cc, ",", c('1', '2'))[,2], ":", c('1', '2'))[,1])
soup_vcf_simple2$la_ref_num = as.numeric(reshape2::colsplit(reshape2::colsplit(soup_vcf_simple2$la, ",", c('1', '2'))[,1], ":", c('1', '2'))[,2])
soup_vcf_simple2$la_alt_num = as.numeric(reshape2::colsplit(reshape2::colsplit(soup_vcf_simple2$la, ",", c('1', '2'))[,2], ":", c('1', '2'))[,1])

sum_df = data.frame()
for (this_soup in c("soup2_0_gt", "soup2_1_gt", "soup2_2_gt")) {
  for (this_species in c("cc", "la")) {
    # alts
    num_alt = length(which(soup_vcf_simple2[, paste0(this_species, "_alt_num")] >= 10 & soup_vcf_simple2[, paste0(this_species, "_ref_num")] == 0))
    v_alt_correct = length(which(soup_vcf_simple2[, paste0(this_species, "_alt_num")] >= 10 & soup_vcf_simple2[, paste0(this_species, "_ref_num")] == 0 & soup_vcf_simple2[, this_soup] == "1/1"))
    v_alt_incorrect = length(which(soup_vcf_simple2[, paste0(this_species, "_alt_num")] >= 10 & soup_vcf_simple2[, paste0(this_species, "_ref_num")] == 0 & soup_vcf_simple2[, this_soup] == "0/0"))
    # refs
    num_ref = length(which(soup_vcf_simple2[, paste0(this_species, "_alt_num")] == 0 & soup_vcf_simple2[, paste0(this_species, "_ref_num")] >= 10))
    v_ref_correct = length(which(soup_vcf_simple2[, paste0(this_species, "_alt_num")] == 0 & soup_vcf_simple2[, paste0(this_species, "_ref_num")] >= 10 & soup_vcf_simple2[, this_soup] == "0/0"))
    v_ref_incorrect = length(which(soup_vcf_simple2[, paste0(this_species, "_alt_num")] == 0 & soup_vcf_simple2[, paste0(this_species, "_ref_num")] >= 10 & soup_vcf_simple2[, this_soup] == "1/1"))
    v_correct_total = v_alt_correct + v_ref_correct
    v_incorrect_total = v_alt_incorrect + v_ref_incorrect
    v_score = v_correct_total - v_incorrect_total
    sum_df = rbind(sum_df, data.frame(this_soup, this_species, v_correct_total, v_incorrect_total, v_score, num_alt, v_alt_correct, v_alt_incorrect, num_ref, v_ref_correct, v_ref_incorrect))
  }
}
sum_df = sum_df[order(sum_df$this_species),]
write.csv(sum_df, "~/scratch/brain/bs/JTS21/JTS21_ALT_1/outs/alt1_matching_summary.csv")

length(which(soup_vcf_simple2$cc_alt_num >= 10 & soup_vcf_simple2$cc_ref_num == 0))
length(which(soup_vcf_simple2$cc_alt_num >= 10 & soup_vcf_simple2$cc_ref_num == 0 & soup_vcf_simple2$soup3_2_gt == "1/1"))
length(which(soup_vcf_simple2$cc_alt_num >= 10 & soup_vcf_simple2$cc_ref_num == 0 & soup_vcf_simple2$soup3_2_gt == "0/0"))

length(which(soup_vcf_simple2$cc_alt_num == 0 & soup_vcf_simple2$cc_ref_num >= 10))
length(which(soup_vcf_simple2$cc_alt_num == 0 & soup_vcf_simple2$cc_ref_num >= 10 & soup_vcf_simple2$soup3_2_gt == "0/0"))
length(which(soup_vcf_simple2$cc_alt_num == 0 & soup_vcf_simple2$cc_ref_num >= 10 & soup_vcf_simple2$soup3_2_gt == "1/1"))

#*******************************************************************************
# CytoTRACE ====================================================================
#*******************************************************************************
library(RColorBrewer)
temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
pal = colorRampPalette(temp)

plk$exp_cond = paste0(plk$exp, ": ", plk$cond)
p_list = FeaturePlot(plk, "cyto", split.by = "exp_cond", combine = F, keep.scale = "all")
p_list2 = list()
for (i in 1:length(p_list)) { p_list2[[i]] = p_list[[i]] + coord_fixed() + scale_color_gradientn(colors = pal(100)) + NoLegend() + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), plot.background = element_rect(fill="#FFFFFF", color = NA)) + ggtitle("") }
Cairo::Cairo(paste0("~/research/tooth/results/plk_cyto_by_time_cond.png"), width = 2400, height = 4800, res = 200)
print(cowplot::plot_grid(plotlist = p_list2, ncol = 2))
dev.off()

#*******************************************************************************
# CellChat =====================================================================
#*******************************************************************************
# Code Snippet for running cellchat on all the experiments
exps = c("plk60", "plk1", "plk3", "plk7")
conds = c("plk", "con")
for (this_exp in exps) {
  for (this_cond in conds) {
    cur_cells = colnames(combined)[which(combined$exp == this_exp & combined$cond == this_cond)]
    data.input = data.input.backup[,cur_cells]
    meta.label = paste0(combined$exp[cur_cells], "_", combined$cond[cur_cells], "_", combined$seurat_clusters[cur_cells])
    this_res = CellChatWeights(x)
    write.csv(this_res[[1]], paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_", this_cond, "_lr_one_050423.csv"))
    write.csv(this_res[[2]], paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_", this_cond, "_one_050423.csv"))
    print(paste0("DONE: ", this_exp, ", ", this_cond))
  }
}

# Plot the celltype interactions
t_order = c (1, 9, 12, 19, 26, 27, 34, 10, 14, 15, 16, 17, 28, 31, 35, 0, 4, 8, 13, 20, 32, 2, 3, 5, 11, 21, 23, 22, 29, 7, 25, 38, 41, 42, 45, 6, 18, 30, 39, 37, 46, 48, 33, 36, 44, 24, 40, 43, 47, 49)
exps = c("plk60", "plk1", "plk3", "plk7")
conds = c("plk", "con")
for (this_exp in exps) {
  for (this_cond in conds) {
    cc_df = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_", this_cond, "_one_050423.csv"))
    cc_df$Sender   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
    cc_df$Receiver = reshape2::colsplit(reshape2::colsplit(cc_df$Receiver, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
    # cc_df$Sender   = factor(cc_df$Sender,   levels = 0:max(as.numeric(cc_df$Sender)))
    # cc_df$Receiver = factor(cc_df$Receiver, levels = 0:max(as.numeric(cc_df$Receiver)))
    cc_df$Sender   = factor(cc_df$Sender,   levels = t_order)
    cc_df$Receiver = factor(cc_df$Receiver, levels = t_order)
    
    fname = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_", this_cond, "_one_062623.pdf")
    ggplot(cc_df, aes(x = Sender, y = Receiver, fill = value)) + geom_raster() + scale_fill_viridis() + ggtitle(paste0("CellChat on Experiment: ", this_exp, " and Condition: ", this_cond)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1)) + coord_fixed() + ggh4x::force_panelsizes(cols = unit(length(unique(cc_df$Sender))/7, "in"), rows = unit(length(unique(cc_df$Sender))/7, "in"))
    ggsave(fname, width = 9, height = 9)
    system(paste0("rclone copy ", fname, " dropbox:BioSci-Streelman/George/Tooth/plk/analysis/cellchat"))
  }
}

# Pluck vs control celltype interaction differences
exps = c("plk60", "plk1", "plk3", "plk7")
for (this_exp in exps) {
  cc_df = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plk_one_050423.csv"))
  cc_df$exp = reshape2::colsplit(cc_df$Sender, "_", c('1', '2'))[,1]
  cc_df$cond   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,1]
  cc_df$Sender   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  cc_df$Receiver = reshape2::colsplit(reshape2::colsplit(cc_df$Receiver, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  cc_df$Sender   = factor(cc_df$Sender,   levels = t_order)
  cc_df$Receiver = factor(cc_df$Receiver, levels = t_order)
  cc_df$id = paste0(cc_df$exp, "_", cc_df$Sender, "_", cc_df$Receiver)
  cc_df_plk = cc_df
  
  cc_df = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_con_one_050423.csv"))
  cc_df$exp = reshape2::colsplit(cc_df$Sender, "_", c('1', '2'))[,1]
  cc_df$cond   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,1]
  cc_df$Sender   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  cc_df$Receiver = reshape2::colsplit(reshape2::colsplit(cc_df$Receiver, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  cc_df$Sender   = factor(cc_df$Sender,   levels = t_order)
  cc_df$Receiver = factor(cc_df$Receiver, levels = t_order)
  cc_df$id = paste0(cc_df$exp, "_", cc_df$Sender, "_", cc_df$Receiver)
  cc_df_con = cc_df
  
  cc_df = cc_df_plk
  cc_df$plk = cc_df$value
  cc_df$con = cc_df_con$value[match(cc_df$id, cc_df_con$id)]
  cc_df$dif = cc_df$plk - cc_df$con
  cc_df = cc_df[which(!is.na(cc_df$dif)),]
  fname = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plkvcon_one_062623.pdf")
  ggplot(cc_df, aes(x = Sender, y = Receiver, fill = dif)) + geom_raster() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(cc_df$dif)), max(abs(cc_df$dif))), oob = scales::squish) + ggtitle(paste0("CellChat Plucked - Control Weights for Experiment: ", this_exp)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1)) + coord_fixed() + ggh4x::force_panelsizes(cols = unit(length(unique(cc_df$Sender))/7, "in"), rows = unit(length(unique(cc_df$Sender))/7, "in"))
  ggsave(fname, width = 9, height = 9)
  system(paste0("rclone copy ", fname, " dropbox:BioSci-Streelman/George/Tooth/plk/analysis/cellchat"))
  # fname2 = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plkvcon_one_050423.csv")
  # write.csv(cc_df, fname2)
  # system(paste0("rclone copy ", fname2, " dropbox:BioSci-Streelman/George/Tooth/plk/analysis/cellchat"))
}

# Pluck vs control LR differences
exps = c("plk60", "plk1", "plk3", "plk7")
for (this_exp in exps) {
  lr_df_plk = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plk_lr_one_050423.csv"))
  lr_df_plk$Sender   = reshape2::colsplit(reshape2::colsplit(lr_df_plk$source,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  lr_df_plk$Receiver = reshape2::colsplit(reshape2::colsplit(lr_df_plk$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  lr_df_plk$id = paste0(lr_df_plk$Sender, "_", lr_df_plk$Receiver, "_", lr_df_plk$ligand, "_", lr_df_plk$receptor)
  
  lr_df_con = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_con_lr_one_050423.csv"))
  lr_df_con$Sender   = reshape2::colsplit(reshape2::colsplit(lr_df_con$source,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  lr_df_con$Receiver = reshape2::colsplit(reshape2::colsplit(lr_df_con$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  lr_df_con$id = paste0(lr_df_con$Sender, "_", lr_df_con$Receiver, "_", lr_df_con$ligand, "_", lr_df_con$receptor)
  
  lr_df_plk[,c("plk_prob", "plk_pval")] = lr_df_plk[,c("prob", "pval")]
  lr_df_plk[,c("con_prob", "con_pval")] = lr_df_con[match(lr_df_plk$id, lr_df_con$id), c("prob", "pval")]
  lr_df_plk = lr_df_plk[which(!is.na(lr_df_plk$con_prob) & lr_df_plk$Sender != lr_df_plk$Receiver),]
  lr_df_plk$prob_all = lr_df_plk$plk_prob + lr_df_plk$con_prob
  lr_df_plk$plk.minus.con = lr_df_plk$plk_prob - lr_df_plk$con_prob
  
  # fname = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plkvcon_lr_one_050423.pdf")
  # ggplot(lr_df_plk, aes(x = prob_all, y = plk.minus.con, color = factor(Sender))) + geom_point(alpha = 0.8) + geom_text_repel(data = lr_df_plk[which(lr_df_plk$prob_all > 0.2 & abs(lr_df_plk$plk.minus.con) > 0.1),], aes(label = id)) + theme_bw() + NoLegend()
  # ggsave(fname, width = 11, height = 7)
  # system(paste0("rclone copy ", fname, " dropbox:BioSci-Streelman/George/Tooth/plk/analysis/cellchat"))
  fname2 = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plkvcon_lr_one_050423.csv")
  write.csv(lr_df_plk, fname2)
  system(paste0("rclone copy ", fname2, " dropbox:BioSci-Streelman/George/Tooth/plk/analysis/cellchat"))
}

#*******************************************************************************
# CC ===========================================================================
#*******************************************************************************
cc_df = data.frame( hgnc = c(cc.genes$s.genes, cc.genes$g2m.genes), class = c(rep("S", length(cc.genes$s.genes)), rep("G2M", length(cc.genes$g2m.genes))) )
cc_df$mz = gene_info$mzebra[match(cc_df$hgnc, gene_info$one_to_one_human)]
plkcc = CellCycleScoring(plk, cc_df$mz[which(cc_df$class == "S" & !is.na(cc_df$mz))], cc_df$mz[which(cc_df$class == "G2M" & !is.na(cc_df$mz))])

library(RColorBrewer)
temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
pal = colorRampPalette(temp)
Cairo::Cairo(paste0("~/scratch/d_tooth/results/near_optimal_2_s.png"), width = 2400, height = 2400, res = 200)
print(FeaturePlot(plkcc, "S.Score", label = TRUE, order = T, pt.size = 0.4) + coord_fixed() + NoLegend() + scale_color_gradientn(colors = pal(100)))
dev.off()

#*******************************************************************************
# Merge w/ TT ==================================================================
#*******************************************************************************
tt = readRDS("~/scratch/d_tooth/data/tt_072722.rds")
plk$dataset = "plk"
tt$dataset  = "tt"
combined = merge(plk, tt)
combined = subset(combined, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & pct.mt < 0.05)
print(paste("Number of Cells in combined After Filterning:", ncol(combined)))
combined = NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined = SCTransform(combined, verbose = TRUE)
combined@active.assay = "SCT"
combined = RunPCA(combined)
combined = RunUMAP(combined, reduction = "pca", dims = 1:50)
combined = FindNeighbors(combined, reduction="umap", dims = 1:2)
combined = FindClusters(combined, resolution = .25)

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkandtt_cluster_split_sample.png"), width = 3200, height = 2800, res = 200)
print(DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 1, split.by = "sample", ncol = 4) + coord_fixed() + NoLegend())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkandtt_cluster_split_datatset.png"), width = 2800, height = 1400, res = 200)
print(DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 1, split.by = "dataset") + coord_fixed() + NoLegend())
dev.off()

tt_annot = read.csv("~/scratch/d_tooth/results/tt_annot_first_pass.csv")
tt_annot$Annotations = paste0(tt_annot$Annotations, " (", tt_annot$Cluster, ")")
combined$plk_cluster = NA
combined$plk_cluster[which(colnames(combined) %in% colnames(plk))] = as.numeric(as.vector(plk$seurat_clusters))[match(colnames(combined)[which(colnames(combined) %in% colnames(plk))], colnames(plk))]
combined$tt_cluster = NA
combined$tt_cluster[which(colnames(combined) %in% colnames(tt))] =  as.numeric(as.vector(tt$seurat_clusters))[match(colnames(combined)[which(colnames(combined) %in% colnames(tt))], colnames(tt))]
combined$tt_annot = combined$tt_cluster
combined$tt_annot = tt_annot$Annotations[match(as.numeric(as.vector(combined$tt_annot)), tt_annot$Cluster)]

Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkandtt_plk_cluster.png"), width = 1400, height = 1400, res = 200)
Idents(combined) = factor(combined$plk_cluster, levels = sort(unique(combined$plk_cluster)))
print(DimPlot(combined, order = T, label = TRUE, pt.size = 1) + coord_fixed())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkandtt_plk_cluster_split_dataset.png"), width = 3200, height = 1400, res = 200)
Idents(combined) = factor(combined$plk_cluster, levels = sort(unique(combined$plk_cluster)))
print(DimPlot(combined, order = T, label = TRUE, pt.size = 1, split.by = "dataset", ncol = 3) + coord_fixed() + NoLegend())
dev.off()
Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkandtt_tt_annot.png"), width = 2400, height = 2400, res = 200)
# Cairo::Cairo(paste0("~/scratch/d_tooth/results/plkandtt_tt_annot.png"), width = 1400, height = 1400, res = 200)
Idents(combined) = factor(combined$tt_annot, levels = sort(unique(combined$tt_annot)))
print(DimPlot(combined, order = T, label = T, pt.size = 1) + coord_fixed())
dev.off()

#*******************************************************************************
# Merge w/ Paul ================================================================
#*******************************************************************************
paul = readRDS("~/research/tooth/data/paul_full.rds")
paul$annot = plyr::revalue(paul$seurat_clusters, replace = c("0" = "Mesenchymal", "1" = "Mesenchymal", "2" = "Mesenchymal", "3" = "Mesenchymal", "4" = "Mesenchymal", "5" = "Mesenchymal", "6" = "Perivascular", "7" = "Mesenchymal", "8" = "Odontoblast", "9" = "Endothelial", "10" = "Macrophages", "11" = "TAC", "12" = "Endothelial", "13" = "Immune", "14" = "Immune", "15" = "Immune", "16" = "Neural Crest", "17" = "?", "18" = "?", "19" = "?", "20" = "?"))
plk.counts = sort(rowSums(plk@assays$RNA@counts), decreasing = T)
plk.counts = data.frame(plk.counts)
plk.counts$mouse = stringr::str_to_title( gene_info$human[match(rownames(plk.counts), gene_info$seurat_name)] )
plk.counts = plk.counts[which(!is.na(plk.counts$mouse) & plk.counts$mouse != "" & plk.counts$mouse %in% rownames(paul@assays$RNA@counts)),]
plk.counts = plk.counts[!duplicated(plk.counts$mouse),]
plk.counts.mat =  plk@assays$RNA@counts[rownames(plk.counts),]
rownames(plk.counts.mat) = plk.counts$mouse
plk.mouse = CreateSeuratObject(counts = plk.counts.mat, meta.data = plk@meta.data)
plk.mouse = NormalizeData(plk.mouse)
plk.mouse = FindVariableFeatures(plk.mouse)
plk.mouse$plk.cluster = plk.mouse$seurat_clusters
plk.mouse$lab = "Talha"
paul.small = CreateSeuratObject(counts = paul@assays$RNA@counts[plk.counts$mouse,], meta.data = paul@meta.data)
paul.small = NormalizeData(paul.small)
paul.small = FindVariableFeatures(paul.small)
paul.small$paul.cluster = paul.small$seurat_clusters
paul.small$lab = "Paul"
paul.small$annot = plyr::revalue(paul.small$seurat_clusters, replace = c("0" = "Mesenchymal", "1" = "Mesenchymal", "2" = "Mesenchymal", "3" = "Mesenchymal", "4" = "Mesenchymal", "5" = "Mesenchymal", "6" = "Perivascular", "7" = "Mesenchymal", "8" = "Odontoblast", "9" = "Endothelial", "10" = "Macrophages", "11" = "TAC", "12" = "Endothelial", "13" = "Immune", "14" = "Immune", "15" = "Immune", "16" = "Neural Crest", "17" = "?", "18" = "?", "19" = "?", "20" = "?"))

paul.plk = RunCCA(paul, plk.mouse, assay1 = "RNA", assay2 = "RNA", renormalize = F, rescale = T)
paul.plk = SCTransform(paul.plk, verbose = TRUE)
paul.plk@active.assay = "SCT"
paul.plk = RunPCA(paul.plk)
paul.plk = RunUMAP(paul.plk, reduction = "pca", dims = 1:50)
paul.plk = FindNeighbors(paul.plk, reduction="umap", dims = 1:2)
paul.plk = FindClusters(paul.plk, resolution = .25)
DimPlot(paul.plk, split.by="lab", label = T)

paul.deg$hgnc = toupper(paul.deg$gene)
plk.deg = read.csv("~/research/tooth/results/plk_cluster_markers_standard_120922.csv")
plk.deg$hgnc = gene_info$human[match(plk.deg$gene, gene_info$seurat_name)]
this.list = list(paul.deg, plk.deg)
names(this.list) = c("paul", "plk")
tmp = bigHeatmap(this.list, pdf.name = "~/research/tooth/results/paul_plk_similarities.pdf", rm.self = T)

#*******************************************************************************
# Colquitt Correlation =========================================================
#*******************************************************************************
incsr = readRDS("~/scratch/d_tooth/data/igor_incsr.rds")
incsr.mes = subset(incsr, cells = colnames(incsr)[which(incsr$annot %in% c("Maturing pulp", "Distal pulp", "Pre-odontoblasts", "Alveolar osteo.", "Dental follicle 1", "Dental follicle 2"))])
incsr.epi = subset(incsr, cells = colnames(incsr)[which(incsr$annot %in% c("Ameloblasts", "SI + SR", "OEE"))])
incsr.immune = subset(incsr, cells = colnames(incsr)[which(incsr$annot %in% c("Lyve1 Macrophages", "Macrophages", "Lymphocytes", "Innate leukocytes"))])


plk.mes =subset(plk, cells = colnames(plk)[which(plk$seurat_clusters %in% c(0,4,5,8,22))])
plk.epi = subset(plk, cells = colnames(plk)[which(plk$seurat_clusters %in% c(1,2,3,7,10,11,12,14,15,16,19,20,26,30,33))])
plk.immune = subset(plk, cells = colnames(plk)[which(plk$seurat_clusters %in% c(6, 17, 24, 27, 28, 31))])

#*******************************************************************************
# Time Point ===================================================================
#*******************************************************************************
Idents(plk) = paste0(plk$cond, " ", plk$exp)
avg_exp = AverageExpression(plk, assays = "SCT")[[1]]
exp_mean = rowMeans(avg_exp)
n_cells = rowSums(plk@assays$RNA@counts[names(exp_mean),])
avg_exp_over_mean = avg_exp / exp_mean
# avg_exp_cool = avg_exp[which( rowSums(avg_exp_over_mean >= 5) > 0 & n_cells > 50 ),]
deg = read.csv("~/scratch/d_tooth/results/plk_master_deg_standard_041223.csv")
bulk_deg_genes = unique(deg$gene[which(deg$level == "bulk" & deg$p_val_adj < 1e-20)])
bulk_deg_genes = bulk_deg_genes[which(bulk_deg_genes %in% rownames(plk))]
avg_exp_cool = avg_exp[bulk_deg_genes,]

exp_df = reshape2::melt(avg_exp_cool)
exp_df$over_mean = reshape2::melt(avg_exp_over_mean[rownames(avg_exp_cool),])[,3]
exp_df[,c("cond", "exp")] = reshape2::colsplit(exp_df[,2], " ", c('1', '2'))
exp_df$exp = factor(as.vector(exp_df$exp), levels = c("plk60", "plk1", "plk3", "plk7"))
ggplot(exp_df, aes(x = exp, y = value, color = Var1, group = Var1)) + geom_point() + geom_line() + facet_wrap(~ cond, ncol = 1) + theme_classic()
ggsave(paste0(out_dir, "deg_genes_over_time_points.pdf"), width = 12, height = 8)
ggplot(exp_df, aes(x = exp, y = over_mean, color = Var1, group = Var1)) + geom_point() + geom_line() + facet_wrap(~ cond, ncol = 1) + theme_classic()
ggsave(paste0(out_dir, "deg_genes_over_time_points_scale.pdf"), width = 12, height = 8)

exp_df2 = exp_df %>% group_by(Var1, exp) %>% summarise(value = mean(over_mean)) %>% mutate(variable = "mean")
exp_df3 = exp_df %>% group_by(Var1, exp) %>% summarise(value = over_mean[1] - over_mean[2]) %>% mutate(variable = "dif")
exp_df4 = data.frame(rbind(exp_df2, exp_df3))
ggplot(exp_df4, aes(x = exp, y = value, color = Var1, group = Var1)) + geom_point() + geom_line() + facet_wrap(~ variable, ncol = 1) + theme_classic()
ggsave(paste0(out_dir, "deg_genes_over_time_points_scale2.pdf"), width = 12, height = 8)

#*******************************************************************************
# Dif Exp ======================================================================
#*******************************************************************************
clusterDeg = function(obj, cluster) {
  this_ident_1 = paste0("plk", "_", cluster)
  this_ident_2 = paste0("con", "_", cluster)
  if (length(which(Idents(obj) == this_ident_1)) >= 3 & length(which(Idents(obj) == this_ident_2)) >= 3) {
    deg = FindMarkers(obj, ident.1=this_ident_1, ident.2=this_ident_2)
    deg$cluster = cluster
    deg$gene = rownames(deg)
    return(deg)
  }
  return(data.frame())
}
Idents(plk1) = paste0(plk1$cond, "_", plk1$seurat_clusters)
plk1_deg = mclapply(as.character(unique(as.vector(plk1$seurat_clusters))), function(x) clusterDeg(plk1, x), mc.cores = 20)
plk1_deg = do.call('rbind', plk1_deg)
plk1_deg_sig = plk1_deg[which(plk1_deg$p_val_adj < 0.05),]

Idents(plk3) = paste0(plk3$cond, "_", plk3$seurat_clusters)
plk3_deg = mclapply(as.character(unique(as.vector(plk3$seurat_clusters))), function(x) clusterDeg(plk3, x), mc.cores = 20)
plk3_deg = do.call('rbind', plk3_deg)
plk3_deg_sig = plk3_deg[which(plk3_deg$p_val_adj < 0.05),]

Idents(plk60) = paste0(plk60$cond, "_", plk60$seurat_clusters)
plk60_deg = mclapply(as.character(unique(as.vector(plk60$seurat_clusters))), function(x) clusterDeg(plk60, x), mc.cores = 20)
plk60_deg = do.call('rbind', plk60_deg)
plk60_deg_sig = plk60_deg[which(plk60_deg$p_val_adj < 0.05),]

Idents(plk7) = paste0(plk7$cond, "_", plk7$seurat_clusters)
plk7_deg = mclapply(as.character(unique(as.vector(plk7$seurat_clusters))), function(x) clusterDeg(plk7, x), mc.cores = 20)
plk7_deg = do.call('rbind', plk7_deg)
plk7_deg_sig = plk7_deg[which(plk7_deg$p_val_adj < 0.05),]

Idents(plk1) = paste0(plk1$cond)
plk1_bulk_deg = FindAllMarkers(plk1, only.pos = T)
plk1_bulk_deg_sig = plk1_bulk_deg[which(plk1_bulk_deg$p_val_adj < 0.05),]
nrow(plk1_bulk_deg_sig)
Idents(plk3) = paste0(plk3$cond)
plk3_bulk_deg = FindAllMarkers(plk3, only.pos = T)
plk3_bulk_deg_sig = plk3_bulk_deg[which(plk3_bulk_deg$p_val_adj < 0.05),]
nrow(plk3_bulk_deg_sig)
Idents(plk60) = paste0(plk60$cond)
plk60_bulk_deg = FindAllMarkers(plk60, only.pos = T)
plk60_bulk_deg_sig = plk60_bulk_deg[which(plk60_bulk_deg$p_val_adj < 0.05),]
nrow(plk60_bulk_deg_sig)
Idents(plk7) = paste0(plk7$cond)
plk7_bulk_deg = FindAllMarkers(plk7, only.pos = T)
plk7_bulk_deg_sig = plk7_bulk_deg[which(plk7_bulk_deg$p_val_adj < 0.05),]
nrow(plk7_bulk_deg_sig)

plk1_deg_sig$exp  = "plk1";  plk1_bulk_deg_sig$exp  = "plk1";
plk3_deg_sig$exp  = "plk3";  plk3_bulk_deg_sig$exp  = "plk3";
plk60_deg_sig$exp = "plk60"; plk60_bulk_deg_sig$exp = "plk60";
plk7_deg_sig$exp  = "plk7";  plk7_bulk_deg_sig$exp  = "plk7";

plk1_deg_sig$level  = "cluster"; plk1_bulk_deg_sig$level  = "bulk";
plk3_deg_sig$level  = "cluster"; plk3_bulk_deg_sig$level  = "bulk";
plk60_deg_sig$level = "cluster"; plk60_bulk_deg_sig$level = "bulk";
plk7_deg_sig$level  = "cluster"; plk7_bulk_deg_sig$level  = "bulk";
plk_deg_master = rbind(plk1_deg_sig,      plk3_deg_sig,      plk60_deg_sig,      plk7_deg_sig,
                       plk1_bulk_deg_sig, plk3_bulk_deg_sig, plk60_bulk_deg_sig, plk7_bulk_deg_sig)


names <- colnames(plk@assays$RNA@counts)
group <- c(rep("CTRL", length(grep("CTRL*", names))), rep("CLIPP", length(grep("CLIPP*", names))))
y <- DGEList(counts = plk@assays$RNA@counts, group = plk$cond)
y <- calcNormFactors(y)
# design <- model.matrix(~group)
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
# y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

gene.means = rowMeans(plk@assays$SCT@counts)
disp = unlist(mclapply(rownames(plk@assays$SCT@counts), function(x) { (var(plk@assays$SCT@counts[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)

library("glmmSeq")
# plk.norm = plk@assays$SCT@data
# disp = setNames(edgeR::estimateDisp()$tagwise.dispersion, rownames())
# 
# dds = DESeq2::DESeqDataSetFromMatrix(plk@assays$RNA@counts, colData = plk@meta.data, design = ~ 1)
# dds = DESeq2::DESeq(dds)
# dispersions = setNames(DESeq2::dispersions(dds), rownames(plk@assays$RNA@counts))
# sizeFactors = DESeq2::estimateSizeFactorsForMatrix(plk@assays$RNA@counts)
# non_zero_genes = rowSums(plk@assays$RNA@counts)
# non_zero_genes = names(non_zero_genes)[which(non_zero_genes >= 5)]
tpm = plk@assays$SCT@data
gene.means = rowMeans(as.matrix(tpm))
disp = unlist(mclapply(rownames(tpm), function(x) { (var(tpm[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)
sizeFactors = colSums(tpm) / mean(tpm)
results = glmmSeq::glmmSeq(~ cond + (1|sample), countdata = tpm, metadata = plk@meta.data, dispersion = disp, progress = T)
stats[order(stats[,'P_cond']),]
write.csv(stats, "~/scratch/d_tooth/results/plk60_7_r12_cond_sample_deg.csv")
length(which(stats[, 'P_cond'] < 0.05))
stats[which(stats[,'P_cond'] < 0.05),]
saveRDS(results, "~/scratch/d_tooth/results/plk60_7_r12_cond_sample_deg.rds")

tpm = plk@assays$SCT@data[,which(plk@meta.data$exp == "plk60")]
gene.means = rowMeans(as.matrix(tpm))
disp = unlist(mclapply(rownames(tpm), function(x) { (var(tpm[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)
sizeFactors = colSums(tpm) / mean(tpm)
results = glmmSeq::glmmSeq(~ cond + (1|sample), countdata = tpm, metadata = plk@meta.data[which(plk@meta.data$exp == "plk60"),], dispersion = disp, progress = T)
length(which(stats[, 'P_cond'] < 0.05))
results <- glmmQvals(results)
saveRDS(results, "~/scratch/d_tooth/results/plk60_r12_cond_sample_deg.rds")

tpm = plk@assays$SCT@data[,which(plk@meta.data$exp == "plk60" & plk@meta.data$sample %in% c("plk60_c12", "plk60_p12"))]
gene.means = rowMeans(as.matrix(tpm))
disp = unlist(mclapply(rownames(tpm), function(x) { (var(tpm[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)
sizeFactors = colSums(tpm) / mean(tpm)
results = glmmSeq::glmmSeq(~ cond + (1|sample), countdata = tpm, metadata = plk@meta.data[which(plk@meta.data$exp == "plk60" & plk@meta.data$sample %in% c("plk60_c12", "plk60_p12")),], dispersion = disp, progress = T)
length(which(stats[, 'P_cond'] < 0.05))
results <- glmmQvals(results)
saveRDS(results, "~/scratch/d_tooth/results/plk60_r1_cond_sample_deg.rds")

tpm = plk@assays$SCT@data[,which(plk@meta.data$exp == "plk7")]
gene.means = rowMeans(as.matrix(tpm))
disp = unlist(mclapply(rownames(tpm), function(x) { (var(tpm[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)
sizeFactors = colSums(tpm) / mean(tpm)
results = glmmSeq::glmmSeq(~ cond + (1|sample), countdata = tpm, metadata = plk@meta.data[which(plk@meta.data$exp == "plk7"),], dispersion = disp, progress = T)
length(which(stats[, 'P_cond'] < 0.05))
results <- glmmQvals(results)
saveRDS(results, "~/scratch/d_tooth/results/plk7_r12_cond_sample_deg.rds")

tpm = plk@assays$SCT@data[,which(plk@meta.data$exp == "plk7" & plk@meta.data$sample %in% c("plk7_c12", "plk7_p12"))]
gene.means = rowMeans(as.matrix(tpm))
disp = unlist(mclapply(rownames(tpm), function(x) { (var(tpm[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)
sizeFactors = colSums(tpm) / mean(tpm)
results = glmmSeq::glmmSeq(~ cond + (1|sample), countdata = tpm, metadata = plk@meta.data[which(plk@meta.data$exp == "plk7" & plk@meta.data$sample %in% c("plk7_c12", "plk7_p12")),], dispersion = disp, progress = T)
length(which(stats[, 'P_cond'] < 0.05))
results <- glmmQvals(results)
saveRDS(results, "~/scratch/d_tooth/results/plk7_r1_cond_sample_deg.rds")

tpm = plk@assays$SCT@data[,which(plk@meta.data$sample %in% c("plk7_c12", "plk7_p12", "plk60_c12", "plk60_p12"))]
gene.means = rowMeans(as.matrix(tpm))
disp = unlist(mclapply(rownames(tpm), function(x) { (var(tpm[x,], na.rm = T)-gene.means[x])/(gene.means[x]**2) }, mc.cores = 24))
names(disp) = names(gene.means)
sizeFactors = colSums(tpm) / mean(tpm)
results = glmmSeq::glmmSeq(~ cond + (1|sample), countdata = tpm, metadata = plk@meta.data[which(plk@meta.data$sample %in% c("plk7_c12", "plk7_p12", "plk60_c12", "plk60_p12")),], dispersion = disp, progress = T)
length(which(stats[, 'P_cond'] < 0.05))
results <- glmmQvals(results)
saveRDS(results, "~/scratch/d_tooth/results/plk60_7_r1_cond_sample_deg.rds")

.libPaths("/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/conda_envs/glmmseq/lib/R/library")
library(glmmSeq)
library(parallel)
set.seed(1234)
plk = readRDS("~/scratch/d_tooth/data/plk60_7_022223.rds")
# data(PEAC_minimal_load)
# disp <- apply(tpm, 1, function(x){
#   (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
# })
# 
# head(disp)
# results <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#                    countdata = tpm,
#                    metadata = metadata,
#                    dispersion = disp,
#                    progress = TRUE)

#*******************************************************************************
# DEG Heatmaps =================================================================
#*******************************************************************************
talha_genes = list()
talha_genes[["60"]] <- c("FOSB", "COL1A2", "MAFB", "DIP2A", "PITX2", "JUP", "PTCHD4", "CELSR1", "FOXO1", "TNFRSF6B", "DSP")
talha_genes[["1"]]  <- c("CD97", "NOTCH1", "TCF7L2", "KRT1")
talha_genes[["3"]]  <- c("SEMA6A", "SATB2", "PLXNA4", "PLEC", "MAP2K4", "AHNAK", "PLXNA2", "SLC11A2")
talha_genes[["7"]]  <- c("PLXNB2", "SCUBE2", "LRP4", "TSPAN6", "SEMA3E", "CHD2", "ST6GALNAC3", "LGR6", "DSC1", "CD74", "SLC9A1", "FOSB")
all_deg = data.frame()
for (this_time in c("60", "1", "3", "7")) {
  this_deg = read.csv(paste0("~/research/tooth/results/plk", this_time, "_plk_v_ctrl_cluster_level_sig_pct.csv"))
  this_deg$time = this_time
  this_deg$y_dif = abs(this_deg$y_plk - this_deg$y_con)
  this_deg$label = this_deg$hgnc
  this_deg$label[which(!this_deg$hgnc %in% talha_genes[[this_time]])] = ""
  this_deg = this_deg[order(this_deg$y_dif, decreasing = T),]
  this_deg$label[which(duplicated(this_deg$label))] = ""
  all_deg = rbind(all_deg, this_deg)
}
all_deg$time = factor(all_deg$time, levels = c("60", "1", "3", "7"))
all_deg$y_dif2 = all_deg$y_dif
all_deg$y_dif2[which(all_deg$y_dif2 > 8)] = 8
all_deg$neg_log_bh = -log10(all_deg$bh)
all_deg$neg_log_bh[which(all_deg$neg_log_bh > 100)] = 100
pdf("~/research/tooth/results/deg_summary.pdf", width = 9, height = 7)
# ggplot(all_deg, aes(x = y_plk, y = y_con, color = y_dif2, size = neg_log_bh, alpha = y_dif2)) + geom_point() + geom_abline(color="gray40", linetype="dashed") + xlab("Mean Expression in Plucked") + ylab("Mean Expression in Control") + scale_color_viridis(option = "plasma") + theme_bw() + geom_text_repel(data = all_deg[which(all_deg$y_dif2 > 4.5),], aes(label=hgnc), size = 5) + guides(size=guide_legend(title=expression("-"*Log["10"]*"(adjusted-p)")), alpha = FALSE, color = guide_legend(title="Dif. in Exp.")) + facet_wrap(~time, scales = "free") 
ggplot(all_deg, aes(x = y_plk, y = y_con, color = y_dif2, size = neg_log_bh, alpha = y_dif2)) + geom_point() + geom_abline(color="gray40", linetype="dashed") + xlab("Mean Expression in Plucked") + ylab("Mean Expression in Control") + scale_color_viridis(option = "plasma") + theme_bw() + geom_text_repel(aes(label=label), size = 5, min.segment.length = unit(0, 'lines')) + guides(size=guide_legend(title=expression("-"*Log["10"]*"(adjusted-p)")), alpha = FALSE, color = guide_legend(title="Dif. in Exp.")) + facet_wrap(~time, scales = "free") 
dev.off()

#*******************************************************************************
# Glmmseq ======================================================================
#*******************************************************************************

# Pluck vs Control
message('Starting glmmseq analysis')
obj = subset(plk_subject, cells = colnames(plk_subject)[which(plk_subject$exp == "plk7")])
obj$pair = obj$subject
n_pairs = length(unique(obj$pair))
n_cond  = length(unique(obj$cond))
glmm_out_dir = "~/scratch/d_tooth/results/plk_glmmseq_plk7_clusters50/"
big_res = data.frame()
for (this_clust in sort(unique(obj$seurat_clusters))) {
  print(this_clust)
  this_file = paste0(glmm_out_dir, "cluster_", this_clust, ".csv")
  if ( file.exists(this_file) ) {
    res = read.csv(this_file)
    res$cluster = this_clust
    # res$q = qvalue::qvalue(res$P_cond)$qvalues
    this_counts = obj@assays$RNA@counts[res$X, which(obj$seurat_clusters == this_clust)]
    this_meta = obj@meta.data[which(obj$seurat_clusters == this_clust),]
    plk_num = rowSums(this_counts[,which(this_meta$cond == "plk")] > 0)
    con_num = rowSums(this_counts[,which(this_meta$cond == "con")] > 0)
    res$num_plk = plk_num
    res$num_con = con_num
    res$pct_plk = plk_num / length(which(this_meta$cond == "plk"))
    res$pct_con = con_num / length(which(this_meta$cond == "con"))
    big_res = rbind(res, big_res)
  }
}
big_res$bh = p.adjust(big_res$P_cond, method = "BH")
big_res$up_pct = big_res$pct_plk
big_res$up_pct[which( big_res$condplk < 0 )] = big_res$pct_con[which( big_res$condplk < 0 )]
big_res$up_num = big_res$num_plk
big_res$up_num[which( big_res$condplk < 0 )] = big_res$num_con[which( big_res$condplk < 0 )]
deg_sig = big_res[which(big_res$bh < 0.05 & big_res$up_pct > 0.1 & big_res$up_num > 5),]
head(deg_sig[order(deg_sig$bh, decreasing = F),])
deg_sig$hgnc = gene_info$human[match(deg_sig$X, gene_info$seurat_name)]
write.csv(big_res, paste0(glmm_out_dir, "all.csv"))
# write.csv(deg_sig, paste0(glmm_out_dir, "all_sig.csv"))
write.csv(deg_sig, paste0(glmm_out_dir, "all_sig_pct.csv"))

bulk = read.csv("plk7_plk_v_ctrl_bulk_sig.csv")
this_counts = obj@assays$RNA@counts[bulk$X, ]
this_meta = obj@meta.data
plk_num = rowSums(this_counts[,which(this_meta$cond == "plk")] > 0)
con_num = rowSums(this_counts[,which(this_meta$cond == "con")] > 0)
bulk$num_plk = plk_num
bulk$con_num = con_num
bulk$pct_plk = plk_num / length(which(this_meta$cond == "plk"))
bulk$pct_con = con_num / length(which(this_meta$cond == "con"))
bulk$up_pct = bulk$pct_plk
bulk$up_pct[which( bulk$condplk < 0 )] = bulk$pct_con[which( bulk$condplk < 0 )]
bulk = bulk[which(bulk$up_pct > 0.1),]
write.csv(bulk, "~/scratch/d_tooth/results/plk7_plk_v_ctrl_bulk_sig_pct.csv")


# Experiment vs Experiment
message('Starting glmmseq analysis')
obj = subset(plk_subject, cells = colnames(plk_subject)[which(plk_subject$cond == "con")])
obj$pair = obj$subject
n_pairs = length(unique(obj$pair))
n_cond  = length(unique(obj$cond))
glmm_out_dir = "~/scratch/d_tooth/results/plk_glmmseq_con_clusters50/"
big_res = data.frame()
for (this_clust in sort(unique(obj$seurat_clusters))) {
  print(this_clust)
  this_file = paste0(glmm_out_dir, "cluster_", this_clust, ".csv")
  if ( file.exists(this_file) ) {
    res = read.csv(this_file)
    res$cluster = this_clust
    # res$q = qvalue::qvalue(res$P_cond)$qvalues
    this_counts = obj@assays$RNA@counts[res$X, which(obj$seurat_clusters == this_clust)]
    this_meta = obj@meta.data[which(obj$seurat_clusters == this_clust),]
    count_list = lapply(unique(this_meta$exp), function(x) data.frame(gene = rownames(this_counts), count = rowSums(this_counts[,which(this_meta$exp == x)] > 0), exp = x))
    count_pct  = lapply(unique(this_meta$exp), function(x) data.frame(gene = rownames(this_counts), pct   = rowSums(this_counts[,which(this_meta$exp == x)] > 0) / length(which(this_meta$exp == x)), exp = x))
    count_df = data.frame(data.table::rbindlist(count_list))
    pct_df   = data.frame(data.table::rbindlist(count_pct))
    count_df = reshape2::dcast(count_df, gene ~ exp, value.var = "count")
    pct_df   = reshape2::dcast(pct_df,   gene ~ exp, value.var = "pct")
    res[, paste0("ncells_",   colnames(count_df)[2:ncol(count_df)])] = count_df[match(res$gene, count_df$gene), 2:ncol(count_df)]
    res[, paste0("pctcells_", colnames(pct_df)[2:ncol(pct_df)])]     = pct_df[match(res$gene, pct_df$gene),     2:ncol(pct_df)]
    big_res = rbind(res, big_res)
  }
}
p_cols = colnames(big_res)[which(startsWith(colnames(big_res), "p."))]
bh_cols = gsub("p\\.", "bh.", p_cols)
for (i in 1:length(p_cols)) { big_res[,bh_cols[i]] = p.adjust(big_res[,p_cols[i]], method = "BH") }

sig_df = data.frame()
for (this_bh in bh_cols) {
  this_est_col = gsub("bh", "estimate", this_bh)
  this_p_col = gsub("bh", "p", this_bh)
  this_split = strsplit(this_bh, "\\.\\.\\.")[[1]]
  exp1 = this_split[2]
  exp2 = this_split[3]
  pos_sig_idx = which(big_res[,this_bh] < 0.05 & big_res[,paste0("pctcells_", exp1)] > 0.1 & big_res[,this_est_col] > 0)
  neg_sig_idx = which(big_res[,this_bh] < 0.05 & big_res[,paste0("pctcells_", exp2)] > 0.1 & big_res[,this_est_col] < 0)
  all_sig_idx = c(pos, neg_sig_idx)
  if (length(neg_sig_idx) > 0) {
    this_df = big_res[neg_sig_idx, c("gene", "cluster", this_p_col, this_bh, this_est_col, paste0("pctcells_", exp1), paste0("pctcells_", exp2), paste0("y_", exp1), paste0("y_", exp2))]
    colnames(this_df) = c("gene", "cluster", "p", "bh", "estimate", "pct.1", "pct.2", "mean.1", "mean.2")
    this_df$exp1Up = F
    this_df$exp1 = exp1
    this_df$exp2 = exp2
    sig_df = rbind(sig_df, this_df)
  }
  if (length(pos_sig_idx) > 0) {
    this_df = big_res[pos_sig_idx, c("gene", "cluster", this_p_col, this_bh, this_est_col, paste0("pctcells_", exp1), paste0("pctcells_", exp2), paste0("y_", exp1), paste0("y_", exp2))]
    colnames(this_df) = c("gene", "cluster", "p", "bh", "estimate", "pct.1", "pct.2", "mean.1", "mean.2")
    this_df$exp1Up = T
    this_df$exp1 = exp1
    this_df$exp2 = exp2
    sig_df = rbind(sig_df, this_df)
  }
}

sig_df = sig_df[order(sig_df$exp1, sig_df$bh, decreasing = F),]
sig_df$hgnc = gene_info$human[match(sig_df$gene, gene_info$seurat_name)]
big_res$hgnc = gene_info$human[match(big_res$gene, gene_info$seurat_name)]
sig_df = sig_df[,c(ncol(sig_df), 1:(ncol(sig_df)-1))]
write.csv(big_res, paste0(glmm_out_dir, "all.csv"))
write.csv(sig_df,  paste0(glmm_out_dir, "all_sig_pct.csv"))

# # Calculate pct.1 and pct.2
# tot = data.frame(table(obj$cond, obj$seurat_clusters))
# colnames(tot) = c("cond", "seurat_clusters", "num_cells")
# count_by_gene = data.frame(gene = rownames(this_counts))
# for (i in 1:nrow(tot)) {
#   this_sum = data.frame(rowSums(this_counts[,which(this_meta$cond == tot$cond[i] & this_meta$seurat_clusters == tot$seurat_clusters[i])]))
#   colnames(this_sum) = i
#   # this_sum = this_sum / length(which(this_meta$seurat_clusters == tot$seurat_clusters[i]))
#   count_by_gene = cbind(this_sum, count_by_gene)
# }

#*******************************************************************************
# Ophir ========================================================================
#*******************************************************************************
# mtx = Matrix::readMM("~/scratch/msc/GSE131204_raw_counts_8594x27998.mtx.gz")
# # mtx  = ReadMtx(mtx = "~/scratch/msc/GSE131204_raw_counts_8594x27998.mtx.gz")
# meta = data.table::fread("~/scratch/msc/GSE131204_cell_info_8594x25.tsv")
# 
# ophir = CreateSeuratObject(counts = mtx, meta.data = meta)

meta = data.table::fread("~/scratch/msc/GSE131204_cell_info_8594x25.tsv", data.table = F)
meta$barcode2 = paste0(meta$barcode, "_1")
ophir_inj = Read10X("~/scratch/msc/ffm_injury/")
ophir_con = Read10X("~/scratch/msc/ffm_control/")
ophir_inj = CreateSeuratObject(ophir_inj)
ophir_con = CreateSeuratObject(ophir_con)
ophir_inj$barcode = colnames(ophir_inj)
ophir_con$barcode = colnames(ophir_con)
ophir_inj$cond = "injured"
ophir_con$cond = "control"
ophir = merge(ophir_inj, ophir_con)
ophir$my_barcode_id = paste0(ophir$cond, "_", ophir$barcode)
meta$my_barcode_id = paste0(meta$condition, "_", meta$barcode)
ophir@meta.data[,colnames(meta)] = meta[match(ophir$my_barcode_id, meta$my_barcode_id), colnames(meta)]
ophir = subset(ophir, cells = colnames(ophir)[which(ophir$pass_quality_filters)])
saveRDS(ophir, "~/scratch/d_tooth/data/ophir.rds")

#*******************************************************************************
# GO ===========================================================================
#*******************************************************************************

# Venn Diagram of Overlapping Cluster Level DEGs
all_hgnc = data.frame()
all_hgnc_list = list()
exps = c("plk60", "plk1", "plk3", "plk7")
for (this_exp in exps) {
  glmm_out_dir = paste0("~/scratch/d_tooth/results/plk_glmmseq_", this_exp, "_clusters50/")
  deg_sig = read.csv(paste0(glmm_out_dir, "all_sig_pct.csv"))
  deg_sig_hgnc = sort(unique(deg_sig$hgnc[which(deg_sig$hgnc != "" & !is.na(deg_sig$hgnc))]))
  all_hgnc = rbind(all_hgnc, data.frame(hgnc = deg_sig_hgnc, exp = this_exp))
  all_hgnc_list[[this_exp]] = deg_sig_hgnc
}
ggvenn::ggvenn(all_hgnc_list, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 0.5, set_name_size = 4)
ggsave("~/scratch/d_tooth/results/plk_v_con_cluster_level_hgnc_venn.pdf", width = 6, height = 6)

# Venn Diagram of Overlapping Bulk DEGs
all_hgnc = data.frame()
all_hgnc_list = list()
exps = c("plk60", "plk1", "plk3", "plk7")
for (this_exp in exps) {
  this_file = paste0("~/scratch/d_tooth/results/", this_exp, "_plk_v_ctrl_bulk_sig_pct.csv")
  deg_sig = read.csv(this_file)
  deg_sig$hgnc = gene_info$human[match(deg_sig$X, gene_info$seurat_name)]
  deg_sig_hgnc = sort(unique(deg_sig$hgnc[which(deg_sig$hgnc != "" & !is.na(deg_sig$hgnc))]))
  all_hgnc = rbind(all_hgnc, data.frame(hgnc = deg_sig_hgnc, exp = this_exp))
  all_hgnc_list[[this_exp]] = deg_sig_hgnc
}
ggvenn::ggvenn(all_hgnc_list, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 0.5, set_name_size = 4)
ggsave("~/scratch/d_tooth/results/plk_v_con_bulk_hgnc_venn.pdf", width = 6, height = 6)

# Cluster Level GO
library("GOfuncR")
all_enrich = data.frame()
exps = c("plk60", "plk1", "plk3", "plk7")
for (this_exp in exps) {
  print(this_exp)
  glmm_out_dir = paste0("~/scratch/d_tooth/results/plk_glmmseq_", this_exp, "_clusters50/")
  deg_sig = read.csv(paste0(glmm_out_dir, "all_sig_pct.csv"))
  deg_sig_hgnc = sort(unique(deg_sig$hgnc[which(deg_sig$hgnc != "" & !is.na(deg_sig$hgnc))]))
  this_enrich = GOfuncR::go_enrich(data.frame(gene = deg_sig_hgnc, test = 1), silent = T)$results
  # this_enrich$bon = p.adjust(this_enrich$results$raw_p_overrep, method = "bonferroni")
  if (this_exp == "plk60") { all_enrich = this_enrich[,c("ontology", "node_id", "node_name", "FWER_overrep")]                           }
  else                     { all_enrich = cbind(all_enrich, this_enrich[match(all_enrich$node_id, this_enrich$node_id), "FWER_overrep"]) }
}

# Plotting
colnames(all_enrich) = c("Ontology", "ID", "Name", paste0(exps, "_FWER"))
all_enrich_sig = all_enrich[which( rowSums(all_enrich[,paste0(exps, "_FWER")] < 0.05) > 0 ),]
write.csv(all_enrich,     "~/scratch/d_tooth/results/plkall_cluster_hdeg_go_unfiltered.csv")
write.csv(all_enrich_sig, "~/scratch/d_tooth/results/plkall_cluster_hdeg_go_sig_in_any.csv")
all_enrich_sig_melt = reshape2::melt(all_enrich_sig, id.vars = c("Ontology", "ID", "Name"))
all_enrich_sig_melt$value2 = (1-all_enrich_sig_melt$value) * 1000
ggplot(all_enrich_sig_melt, aes(x = value2, y = Name, fill = variable)) + geom_bar(stat='identity', position = position_dodge()) + theme_classic() + scale_x_continuous(expand = c(0,0))
ggsave("~/scratch/d_tooth/results/plk_v_con_cluster_level_go.pdf", width = 8, height = 10)

tmp = all_enrich_sig[,paste0(exps, "_FWER")] < 0.05
rownames(tmp) = all_enrich_sig$Name
tmp = hclust(dist(tmp))
all_enrich_sig_melt$sig = all_enrich_sig_melt$value < 0.05
all_enrich_sig_melt$Name = factor(all_enrich_sig_melt$Name, levels = tmp$labels[tmp$order])
ggplot(all_enrich_sig_melt, aes(x = variable, y = Name, size = value2, color = sig)) + geom_point() + theme_classic() + scale_color_manual(values = c("gray60", "goldenrod1")) + scale_size_continuous(range=c(1,4)) + xlab("") + ylab("")
ggsave("~/scratch/d_tooth/results/plk_v_con_cluster_level_go_dot.pdf", width = 8, height = 10)
tm_sel = read.csv("~/scratch/d_tooth/results/plkall_cluster_hdeg_go_sig_in_any_TM_Selections.csv")
all_enrich_sig_melt_tm = all_enrich_sig_melt[which(all_enrich_sig_melt$ID %in% tm_sel$ID),]
ggplot(all_enrich_sig_melt_tm, aes(x = variable, y = Name, size = value2, color = sig)) + geom_point() + theme_classic() + scale_color_manual(values = c("gray60", "goldenrod1")) + scale_size_continuous(range=c(1,4)) + xlab("") + ylab("")
ggsave("~/scratch/d_tooth/results/plk_v_con_cluster_level_go_dot_tm_sel.pdf", width = 8, height = 5)

# Bulk DEGs GO
all_enrich_bulk = data.frame()
for (this_exp in c("plk60", "plk1", "plk3", "plk7")) {
  print(this_exp)
  this_file = paste0("~/scratch/d_tooth/results/", this_exp, "_plk_v_ctrl_bulk_sig_pct.csv")
  deg_sig = read.csv(this_file)
  deg_sig$hgnc = gene_info$human[match(deg_sig$X, gene_info$seurat_name)]
  deg_sig_hgnc = sort(unique(deg_sig$hgnc[which(deg_sig$hgnc != "" & !is.na(deg_sig$hgnc))]))
  this_enrich = GOfuncR::go_enrich(data.frame(gene = deg_sig_hgnc, test = 1), silent = T)$results
  # this_enrich$bon = p.adjust(this_enrich$results$raw_p_overrep, method = "bonferroni")
  if (this_exp == "plk60") { all_enrich_bulk = this_enrich[,c("ontology", "node_id", "node_name", "FWER_overrep")] }
  else                     { all_enrich_bulk = cbind(all_enrich_bulk, this_enrich[match(all_enrich_bulk$node_id, this_enrich$node_id), "FWER_overrep"]) }
}

# Plotting
colnames(all_enrich_bulk) = c("Ontology", "ID", "Name", paste0(exps, "_FWER"))
all_enrich_bulk_sig = all_enrich_bulk[which( rowSums(all_enrich_bulk[,paste0(exps, "_FWER")] < 0.05) > 0 ),]
write.csv(all_enrich,     "~/scratch/d_tooth/results/plkall_bulk_hdeg_go_unfiltered.csv")
write.csv(all_enrich_sig, "~/scratch/d_tooth/results/plkall_bulk_hdeg_go_sig_in_any.csv")
all_enrich_bulk_sig_melt = reshape2::melt(all_enrich_bulk_sig, id.vars = c("Ontology", "ID", "Name"))
all_enrich_bulk_sig_melt$value2 = (1-all_enrich_bulk_sig_melt$value) * 1000
tmp = all_enrich_bulk_sig[,paste0(exps, "_FWER")] < 0.05
rownames(tmp) = all_enrich_bulk_sig$Name
tmp = hclust(dist(tmp))
all_enrich_bulk_sig_melt$sig = all_enrich_bulk_sig_melt$value < 0.05
all_enrich_bulk_sig_melt$Name = factor(all_enrich_bulk_sig_melt$Name, levels = tmp$labels[tmp$order])
ggplot(all_enrich_bulk_sig_melt, aes(x = variable, y = Name, size = value2, color = sig)) + geom_point() + theme_classic() + scale_color_manual(values = c("gray60", "goldenrod1")) + scale_size_continuous(range=c(1,4))
ggsave("~/scratch/d_tooth/results/plk_v_con_bulk_go_dot.pdf", width = 8, height = 20)
tm_sel = read.csv("~/scratch/d_tooth/results/plkall_bulk_hdeg_go_sig_in_any_TM_Selections.csv")
all_enrich_sig_melt_tm = all_enrich_sig_melt[which(all_enrich_sig_melt$ID %in% tm_sel$ID),]
ggplot(all_enrich_sig_melt_tm, aes(x = variable, y = Name, size = value2, color = sig)) + geom_point() + theme_classic() + scale_color_manual(values = c("gray60", "goldenrod1")) + scale_size_continuous(range=c(1,4)) + xlab("") + ylab("")
ggsave("~/scratch/d_tooth/results/plk_v_con_bulk_level_go_dot_tm_sel.pdf", width = 8, height = 5)


# KEGG
# pathways.list <- keggList("pathway", "hsa")
# myKeggTest = function(x) { }
# hsa = search_kegg_organism('hsa', by='kegg_code')
library("diffEnrich")
kegg_hsa = get_kegg('hsa')
kegg_hsa$pathway_to_species$V1 = paste0("path:", kegg_hsa$pathway_to_species$V1)
human_ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_enrich = data.frame()
exps = c("plk60", "plk1", "plk3", "plk7")
for (this_exp in exps) {
  print(this_exp)
  glmm_out_dir = paste0("~/scratch/d_tooth/results/plk_glmmseq_", this_exp, "_clusters50/")
  deg_sig = read.csv(paste0(glmm_out_dir, "all_sig_pct.csv"))
  deg_sig_hgnc = sort(unique(deg_sig$hgnc[which(deg_sig$hgnc != "" & !is.na(deg_sig$hgnc))]))
  mapping = biomaRt::getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), filters = 'hgnc_symbol', values = deg_sig_hgnc, mart = human_ensembl)
  this_enrich = diffEnrich::pathEnrich(gk_obj = kegg_hsa, gene_list = as.character(mapping$entrezgene_id), cutoff = 0.05, N=0)$enrich_table
  this_enrich$bh = p.adjust(this_enrich$enrich_p, method = "BH")
  if (this_exp == "plk60") { all_enrich = this_enrich[,c("KEGG_PATHWAY_ID", "KEGG_PATHWAY_description", "bh")] }
  else                     { all_enrich = cbind(all_enrich, this_enrich[match(all_enrich$KEGG_PATHWAY_ID, this_enrich$KEGG_PATHWAY_ID), "bh"]) }
}
colnames(all_enrich) = c("ID", "Description", paste0(exps, "_bh"))
for (this_exp in exps) { all_enrich[,paste0(this_exp, "_neg_log_bh")] = -log10(all_enrich[,paste0(this_exp, "_bh")]) }
all_enrich_sig = all_enrich[which( rowSums(all_enrich[,paste0(exps, "_bh")] < 0.05) > 0 ),]
all_enrich_sig_melt = reshape2::melt(all_enrich_sig[,c("ID", "Description", paste0(exps, "_neg_log_bh"))], id.vars = c("ID", "Description"))
all_enrich_sig_melt$value2 = all_enrich_sig_melt$value
tmp = all_enrich_sig[,paste0(exps, "_bh")] < 0.05
rownames(tmp) = all_enrich_sig$Description
tmp = hclust(dist(tmp))
all_enrich_sig_melt$sig = all_enrich_sig_melt$value > -log10(0.05)
all_enrich_sig_melt$Description = factor(all_enrich_sig_melt$Description, levels = tmp$labels[tmp$order])
all_enrich_sig_melt$Description = reshape2::colsplit(all_enrich_sig_melt$Description, " - Homo sapiens", c('1', '2'))[,1]
ggplot(all_enrich_sig_melt, aes(x = variable, y = Description, size = value2, color = sig)) + geom_point() + theme_classic() + scale_color_manual(values = c("gray60", "goldenrod1")) + scale_size_continuous(range=c(1,4))
ggsave("~/scratch/d_tooth/results/plk_v_con_cluster_kegg_dot.pdf", width = 10, height = 4)


all_enrich_bulk = data.frame()
for (this_exp in c("plk60", "plk1", "plk3", "plk7")) {
  print(this_exp)
  this_file = paste0("~/scratch/d_tooth/results/", this_exp, "_plk_v_ctrl_bulk_sig_pct.csv")
  deg_sig = read.csv(this_file)
  deg_sig$hgnc = gene_info$human[match(deg_sig$X, gene_info$seurat_name)]
  deg_sig_hgnc = sort(unique(deg_sig$hgnc[which(deg_sig$hgnc != "" & !is.na(deg_sig$hgnc))]))
  mapping = biomaRt::getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), filters = 'hgnc_symbol', values = deg_sig_hgnc, mart = human_ensembl)
  this_enrich = diffEnrich::pathEnrich(gk_obj = kegg_hsa, gene_list = as.character(mapping$entrezgene_id), cutoff = 0.05, N=0)$enrich_table
  this_enrich$bh = p.adjust(this_enrich$enrich_p, method = "BH")
  if (this_exp == "plk60") { all_enrich_bulk = this_enrich[,c("KEGG_PATHWAY_ID", "KEGG_PATHWAY_description", "bh")] }
  else                     { all_enrich_bulk = cbind(all_enrich_bulk, this_enrich[match(all_enrich_bulk$KEGG_PATHWAY_ID, this_enrich$KEGG_PATHWAY_ID), "bh"]) }
}
colnames(all_enrich_bulk) = c("ID", "Description", paste0(exps, "_bh"))
for (this_exp in exps) { all_enrich_bulk[,paste0(this_exp, "_neg_log_bh")] = -log10(all_enrich_bulk[,paste0(this_exp, "_bh")]) }

# Plotting
all_enrich_bulk_sig = all_enrich_bulk[which( rowSums(all_enrich_bulk[,paste0(exps, "_bh")] < 0.05) > 0 ),]
all_enrich_bulk_sig_melt = reshape2::melt(all_enrich_bulk_sig[,c("ID", "Description", paste0(exps, "_neg_log_bh"))], id.vars = c("ID", "Description"))
all_enrich_bulk_sig_melt$value2 = all_enrich_bulk_sig_melt$value
tmp = all_enrich_bulk_sig[,paste0(exps, "_bh")] < 0.05
rownames(tmp) = all_enrich_bulk_sig$Description
tmp = hclust(dist(tmp))
all_enrich_bulk_sig_melt$sig = all_enrich_bulk_sig_melt$value > -log10(0.05)
all_enrich_bulk_sig_melt$Description = factor(all_enrich_bulk_sig_melt$Description, levels = tmp$labels[tmp$order])
all_enrich_bulk_sig_melt$Description = reshape2::colsplit(all_enrich_bulk_sig_melt$Description, " - Homo sapiens", c('1', '2'))[,1]
ggplot(all_enrich_bulk_sig_melt, aes(x = variable, y = Description, size = value2, color = sig)) + geom_point() + theme_classic() + scale_color_manual(values = c("gray60", "goldenrod1")) + scale_size_continuous(range=c(1,4))
ggsave("~/scratch/d_tooth/results/plk_v_con_bulk_kegg_dot.pdf", width = 10, height = 4)

# END KEGG


#*******************************************************************************
# WGCNA ========================================================================
#*******************************************************************************
library("Seurat")
library("SeuratObject")
good_libPaths = .libPaths()
.libPaths(c(.libPaths(), "/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/conda_envs/r4/lib/R/library"))
library("cluster")
library("WGCNA")
.libPaths(good_libPaths)
obj = plk
# non_low_genes = rownames(obj)[which( rowSums(obj@assays$RNA@counts[rownames(obj),]) > 10 )]
# data_mat_c = t(obj@assays$RNA@data[non_low_genes,])
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# # sft_c = pickSoftThreshold(data_mat_c, powerVector = powers, verbose = 5)
# adjacency = adjacency(data_mat_c, type = "signed", power = 1)
# TOM = adjacency
# dissTOM = 1-TOM

# Load Correlation Matrix computed in python
library("rhdf5")
h5f = H5Fopen("~/scratch/d_tooth/results/py_ns/plk_real_cor.h5")
cor_mat = h5f$name
h5closeAll()

# Use WGCNA's adjacency.fromSimilarity function, which will use the power 6
rownames(cor_mat) = colnames(cor_mat) = rownames(obj@assays$SCT@data)
adjacency_6 = adjacency.fromSimilarity(cor_mat, type = "signed")
TOM_6 = adjacency_6
dissTOM_6 = 1-TOM_6
distDissTOM_6 = as.dist(dissTOM_6)
geneTree_6 = hclust(distDissTOM_6, method = "average")

# Regular hierarchical clustering of correlation matrix
adjacency_raw = cor_mat
TOM_raw = adjacency_raw
dissTOM_raw = 1 - TOM_raw
distDissTOM_raw = as.dist(dissTOM_raw)
geneTree_raw = hclust(distDissTOM_raw)

# Cluster
dynamicMods_6     = cutreeDynamic(dendro = geneTree_6, distM = dissTOM_6, pamRespectsDendro = FALSE, minClusterSize = 5)
dynamicMods30_6   = cutreeDynamic(dendro = geneTree_6, distM = dissTOM_6, pamRespectsDendro = FALSE, minClusterSize = 30)
dynamicMods_raw   = cutreeDynamic(dendro = geneTree_raw, distM = dissTOM_raw, pamRespectsDendro = FALSE, minClusterSize = 5)
dynamicMods30_raw = cutreeDynamic(dendro = geneTree_raw, distM = dissTOM_raw, pamRespectsDendro = FALSE, minClusterSize = 30)
df = data.frame(gene = colnames(cor_mat), hgnc = gene_info$human[match(colnames(cor_mat), gene_info$seurat_name)], 
                module_5_6 = dynamicMods_6, module_30_6 = dynamicMods30_6, 
                module_5_raw = dynamicMods_raw, module_30_raw = dynamicMods30_raw, row.names = colnames(cor_mat))
fossil::rand.index(df$module_5_6,  df$module_5_raw)
fossil::rand.index(df$module_30_6, df$module_30_raw)

# DBSCAN
library("dbscan")
cl5_6 = hdbscan(distDissTOM_6, minPts = 5)
cl5_raw = hdbscan(distDissTOM_raw, minPts = 5)
df$dbscan_5_6   = cl5_6$cluster
df$dbscan_5_raw = cl5_raw$cluster
# df$dbscan_membership = cl5$membership_prob
# test = dbscan(distDissTOM_raw, minPts = 5, eps = 0.60)
# km_model <- ClusterR::MiniBatchKmeans(cor_mat, clusters = 50, batch_size = 1000, num_init = 10, max_iters = 100, 
#                                       init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
#                                       verbose = T)
# pred <- ClusterR::predict_MBatchKMeans(cor_mat, km_model$centroids)

# Evaluate Quality of Clustering
sum_df = evaluateClustering(df[,3], distDissTOM_6, cor_mat)
sum_df = rbind(sum_df, evaluateClustering(df[,4], distDissTOM_6,   cor_mat))
sum_df = rbind(sum_df, evaluateClustering(df[,5], distDissTOM_raw, cor_mat))
sum_df = rbind(sum_df, evaluateClustering(df[,6], distDissTOM_raw, cor_mat))
sum_df = rbind(sum_df, evaluateClustering(df[,7], distDissTOM_6,   cor_mat))
sum_df = rbind(sum_df, evaluateClustering(df[,8], distDissTOM_raw, cor_mat))

# Iterative cutree
geneTree = geneTree_raw
distDissTOM = distDissTOM_raw
kdf = data.frame(gene = colnames(cor_mat), hgnc = gene_info$human[match(colnames(cor_mat), gene_info$seurat_name)], row.names = colnames(cor_mat))
my_ks = seq(50, 750, by = 25)
for (i in my_ks) {
  print(i)
  this_clustering = cutree(geneTree, k = i)
  kdf = cbind(kdf, this_clustering)
}
res = parallel::mclapply(3:ncol(kdf), function(x) evaluateClustering(kdf[,x], distDissTOM, cor_mat), mc.cores = 1)
sum_kdf = do.call('rbind', res)

evaluateClustering = function(this_clustering, distDissTOM, cor_mat) {
  print(".")
  cluster_sum = summary(silhouette(this_clustering, distDissTOM))
  this_avg_clus_avg_width = mean(cluster_sum$clus.avg.widths)
  this_median_clus_avg_width = median(cluster_sum$clus.avg.widths)
  this_avg_width = cluster_sum$avg.width
  this_cors = lapply(unique(this_clustering), function(x) as.vector(cor_mat[which(this_clustering==x), which(this_clustering==x)]))
  this_avg_cor = mean(unlist(this_cors))
  this_avg_clus_avg_cor    = mean(  unlist(lapply(1:length(this_cors), function(x) mean(this_cors[[x]]))))
  this_median_clus_avg_cor = median(unlist(lapply(1:length(this_cors), function(x) mean(this_cors[[x]]))))
  return(data.frame(avg_width = this_avg_width, avg_clus_avg_width = this_avg_clus_avg_width, median_clus_avg_width = this_median_clus_avg_width, avg_cor = this_avg_cor, avg_clus_avg_cor = this_avg_clus_avg_cor, median_clus_avg_cor = this_median_clus_avg_cor))
}

# Iterative kmeans clustering
distDissTOM = as.dist(dissTOM)
kdf = data.frame(k = seq(10, 150, by = 5), avg.sil.width = 0)
for (i in 1:nrow(kdf)) {
  print(i)
  kdf$avg.sil.width[i] = pam(distDissTOM, kdf$k[i], diss = T)[["silinfo"]]$avg.width
}
kdfk = data.frame(gene = colnames(cor_mat), hgnc = gene_info$human[match(colnames(cor_mat), gene_info$seurat_name)], row.names = colnames(cor_mat))
my_ks = seq(10, 150, by = 5)
for (i in my_ks) {
  print(i)
  this_clustering = cutree(geneTree, k = i)
  kdfk = cbind(kdfk, this_clustering)
}
res = parallel::mclapply(3:ncol(kdfk), function(x) evaluateClustering(kdfk[,x], distDissTOM, cor_mat), mc.cores = 1)
sum_kdfk = do.call('rbind', res)

#*******************************************************************************
# LR Pathways ==================================================================
#*******************************************************************************
# this_timepoint = "60"
for (this_timepoint in c("60", "1", "3", "7")) {
  print(this_timepoint)
  lr_df = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_plkvcon_lr_one_050423.csv"))
  # lr_df = lr_df[which(lr_df$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7")),]
  # lr_df$pathway_name[which(startsWith(lr_df$pathway_name, "SEMA"))] = "SEMA"
  # lr_df_sum = lr_df %>% group_by(pathway_name, source, target) %>% summarise(mean = mean(plk.minus.con))
  lr_df_sum_source = lr_df %>% group_by(pathway_name, source) %>% summarise(mean = mean(plk.minus.con))
  colnames(lr_df_sum_source)[2] = "cluster_long"
  lr_df_sum_source$node = "source"
  lr_df_sum_target = lr_df %>% group_by(pathway_name, target) %>% summarise(mean = mean(plk.minus.con))
  colnames(lr_df_sum_target)[2] = "cluster_long"
  lr_df_sum_target$node = "target"
  lr_df_sum = rbind(lr_df_sum_source, lr_df_sum_target)
  lr_df_sum$cluster = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$cluster_long, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster = factor(lr_df_sum$cluster, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  # ggplot(lr_df_sum, aes(x = cluster, y = pathway_name, fill = mean)) + geom_raster() + theme_classic() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(lr_df_sum$mean)), max(abs(lr_df_sum$mean))), oob = scales::squish) + facet_wrap(~ node, ncol = 1) + coord_fixed() + ggtitle(paste0("Timepoint: ", this_timepoint))
  # ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_lr_pathways.pdf"), width = 12, height = 12)
  ggplot(lr_df_sum, aes(x = cluster, y = pathway_name, fill = mean)) + geom_raster() + theme_classic() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(lr_df_sum$mean)), max(abs(lr_df_sum$mean))), oob = scales::squish) + facet_wrap(~ node, ncol = 2) + coord_fixed() + ggtitle(paste0("Timepoint: ", this_timepoint))
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_lr_pathways_all.pdf"), width = 25, height = 15)
}

this_timepoint = "60"
lr_df_plk = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_plk_lr_one_050423.csv"))
lr_df_con = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_con_lr_one_050423.csv"))
# lr_df = lr_df[which(lr_df$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7")),]
# lr_df$pathway_name[which(startsWith(lr_df$pathway_name, "SEMA"))] = "SEMA"
# lr_df_sum = lr_df %>% group_by(pathway_name, source, target) %>% summarise(mean = mean(plk.minus.con))
lr_df = lr_df_plk
# lr_df_sum = lr_df %>% group_by(source, target) %>% summarise(mean = n())
lr_df_sum = lr_df %>% group_by(source, target) %>% summarise(mean = mean(prob))
lr_df_sum$cluster_source = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$source, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
lr_df_sum$cluster_source = factor(lr_df_sum$cluster_source, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
lr_df_sum$cluster_target = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
lr_df_sum$cluster_target = factor(lr_df_sum$cluster_target, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
lr_df_sum$id = paste0(lr_df_sum$cluster_source, "_", lr_df_sum$cluster_target)
lr_df_sum_plk = lr_df_sum

lr_df = lr_df_con
lr_df_sum = lr_df %>% group_by(source, target) %>% summarise(mean = n())
lr_df_sum = lr_df %>% group_by(source, target) %>% summarise(mean = mean(prob))
lr_df_sum$cluster_source = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$source, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
lr_df_sum$cluster_source = factor(lr_df_sum$cluster_source, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
lr_df_sum$cluster_target = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
lr_df_sum$cluster_target = factor(lr_df_sum$cluster_target, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
lr_df_sum$id = paste0(lr_df_sum$cluster_source, "_", lr_df_sum$cluster_target)
lr_df_sum_con = lr_df_sum

lr_df_sum = merge(lr_df_sum_plk, lr_df_sum_con, all = T, by = "id", suffixes = c("_plk", "_con"))
lr_df_sum$mean_plk[which(is.na(lr_df_sum$mean_plk))] = 0
lr_df_sum$mean_con[which(is.na(lr_df_sum$mean_con))] = 0
lr_df_sum$mean = lr_df_sum$mean_plk - lr_df_sum$mean_con
lr_df_sum[,c("cluster_source", "cluster_target")] = reshape2::colsplit(lr_df_sum$id, "_", c('1', '2'))
lr_df_sum$cluster_source = factor(lr_df_sum$cluster_source, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
lr_df_sum$cluster_target = factor(lr_df_sum$cluster_target, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
ggplot(lr_df_sum, aes(x = cluster_source, y = cluster_target, fill = mean)) + geom_raster() + theme_classic() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(lr_df_sum$mean)), max(abs(lr_df_sum$mean))), oob = scales::squish) + coord_fixed() + ggtitle(paste0("Timepoint: ", this_timepoint))

for (this_timepoint in c("60", "1", "3", "7")) {
  print(this_timepoint)
  lr_df_plk = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_plk_lr_one_050423.csv"))
  lr_df_con = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_con_lr_one_050423.csv"))
  
  lr_df = lr_df_plk
  lr_df = lr_df[which(lr_df$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7")),]
  lr_df$pathway_name[which(startsWith(lr_df$pathway_name, "SEMA"))] = "SEMA"
  lr_df_sum_source = lr_df %>% group_by(pathway_name, source) %>% summarise(mean = n())
  colnames(lr_df_sum_source)[2] = "cluster_long"
  lr_df_sum_source$node = "source"
  lr_df_sum_target = lr_df %>% group_by(pathway_name, target) %>% summarise(mean = n())
  colnames(lr_df_sum_target)[2] = "cluster_long"
  lr_df_sum_target$node = "target"
  lr_df_sum = rbind(lr_df_sum_source, lr_df_sum_target)
  lr_df_sum$cluster = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$cluster_long, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster = factor(lr_df_sum$cluster, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$id = paste0(lr_df_sum$pathway_name, "_", lr_df_sum$cluster, "_", lr_df_sum$node)
  lr_df_sum_plk = lr_df_sum
  
  lr_df = lr_df_con
  lr_df = lr_df[which(lr_df$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7")),]
  lr_df$pathway_name[which(startsWith(lr_df$pathway_name, "SEMA"))] = "SEMA"
  lr_df_sum_source = lr_df %>% group_by(pathway_name, source) %>% summarise(mean = n())
  colnames(lr_df_sum_source)[2] = "cluster_long"
  lr_df_sum_source$node = "source"
  lr_df_sum_target = lr_df %>% group_by(pathway_name, target) %>% summarise(mean = n())
  colnames(lr_df_sum_target)[2] = "cluster_long"
  lr_df_sum_target$node = "target"
  lr_df_sum = rbind(lr_df_sum_source, lr_df_sum_target)
  lr_df_sum$cluster = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$cluster_long, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster = factor(lr_df_sum$cluster, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$id = paste0(lr_df_sum$pathway_name, "_", lr_df_sum$cluster, "_", lr_df_sum$node)
  lr_df_sum_con = lr_df_sum
  
  lr_df_sum = merge(lr_df_sum_plk, lr_df_sum_con, all = T, by = "id", suffixes = c("_plk", "_con"))
  lr_df_sum$mean_plk[which(is.na(lr_df_sum$mean_plk))] = 0
  lr_df_sum$mean_con[which(is.na(lr_df_sum$mean_con))] = 0
  lr_df_sum$mean = lr_df_sum$mean_plk - lr_df_sum$mean_con
  lr_df_sum[, c("pathway_name", "tmp")] = reshape2::colsplit(lr_df_sum$id, "_", c('1', '2'))
  lr_df_sum[, c("cluster", "node")] = reshape2::colsplit(lr_df_sum$tmp, "_", c('1', '2'))
  lr_df_sum$cluster = factor(as.numeric(lr_df_sum$cluster), levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  
  # ggplot(lr_df_sum, aes(x = cluster, y = pathway_name, fill = mean)) + geom_raster() + theme_classic() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(lr_df_sum$mean)), max(abs(lr_df_sum$mean))), oob = scales::squish) + facet_wrap(~ node, ncol = 2) + coord_fixed() + ggtitle(paste0("Timepoint: ", this_timepoint))
  # ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_lr_pathways_all.pdf"), width = 20, height = 15)
  ggplot(lr_df_sum, aes(x = cluster, y = pathway_name, fill = mean)) + geom_raster() + theme_classic() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(lr_df_sum$mean)), max(abs(lr_df_sum$mean))), oob = scales::squish) + facet_wrap(~ node, ncol = 1) + coord_fixed() + ggtitle(paste0("Timepoint: ", this_timepoint))
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_lr_pathways.pdf"), width = 12, height = 12)
}

# Individual connections
celltype_vector = c ("epi1", "epi1", "epi1", "epi1", "epi1", "epi1", "epi1", 
                     "epi2", "epi2", "epi2", "epi2", "epi2", "epi2", "epi2", "epi2", 
                     "epi3", "epi3", "epi3", "epi3", "epi3", 
                     "pdl", 
                     "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", 
                     "endothelium", "endothelium", 
                     "macrophage", "macrophage", "macrophage", "macrophage", "macrophage",
                     "bone", "bone",
                     "bone macrophage", 
                     "glia", 
                     "taste bud", 
                     "myoblasts", "myoblasts",
                     "nk and t cells", 
                     "unknown", 
                     "pericytes", "pericytes",
                     "lympho-epithelial", "lympho-epithelial", "lympho-epithelial", "lympho-epithelial", "lympho-epithelial")
ct_df = data.frame(ct = celltype_vector, cluster = as.numeric(c("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
for (this_timepoint in c("60", "1")) {
  print(this_timepoint)
  lr_df_plk = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_plk_lr_one_050423.csv"))
  lr_df_con = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_con_lr_one_050423.csv"))
  this_good = readxl::read_xlsx("~/Downloads/cc_good_connections.xlsx", sheet = paste0("Day", this_timepoint))
  this_good$clust_id = paste0(this_good$Sender, "_", this_good$Reciever)
  
  lr_df_sum = lr_df_plk
  lr_df_sum$pathway_name[which(! lr_df_sum$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7"))] = "Other"
  # lr_df_sum = lr_df_sum[which(lr_df_sum$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7")),]
  lr_df_sum$pathway_name[which(startsWith(lr_df_sum$pathway_name, "SEMA"))] = "SEMA"
  lr_df_sum$cluster_source = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$source, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster_source = factor(lr_df_sum$cluster_source, levels = as.numeric(c("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$cluster_target = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster_target = factor(lr_df_sum$cluster_target, levels = as.numeric(c("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$clust_id = paste0(lr_df_sum$cluster_source, "_", lr_df_sum$cluster_target)
  lr_df_sum = lr_df_sum[which(lr_df_sum$clust_id %in% this_good$clust_id),]
  lr_df_sum$ct_source = ct_df$ct[match(lr_df_sum$cluster_source, ct_df$cluster)]
  lr_df_sum$ct_target = ct_df$ct[match(lr_df_sum$cluster_target, ct_df$cluster)]
  lr_df_sum$ct_id = paste0(lr_df_sum$ct_source, " -> ", lr_df_sum$ct_target)
  lr_df_sum_num_rows = lr_df_sum %>% group_by(pathway_name, clust_id, ct_id) %>% summarise(num_rows = n())
  lr_df_sum2 = lr_df_sum_num_rows %>% group_by(pathway_name, ct_id) %>% summarise(mean = mean(num_rows))
  # lr_df_sum2 = lr_df_sum %>% group_by(pathway_name, ct_id) %>% summarise(num = n())
  # lr_df_sum$id = paste0(lr_df_sum$cluster_source, "@", lr_df_sum$cluster_target, "@", lr_df_sum$ligand, "@", lr_df_sum$receptor)
  lr_df_sum2$id = paste0(lr_df_sum2$pathway_name, "@", lr_df_sum2$ct_id)
  lr_df_sum_plk = lr_df_sum2
  
  lr_df_sum = lr_df_con
  lr_df_sum$pathway_name[which(! lr_df_sum$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7"))] = "Other"
  # lr_df_sum = lr_df_sum[which(lr_df_sum$pathway_name %in% c("BMP", "WNT", "ncWNT", "HH", "SHH", "NOTCH", "GAS", "SPP1", "CXCL", "CSPG4", "EGF", "FGF", "IGF", "PDGF", "VEGF", "TGFb", "SEMA3", "SEMA4", "SEMA5", "SEMA6", "SEMA7")),]
  lr_df_sum$pathway_name[which(startsWith(lr_df_sum$pathway_name, "SEMA"))] = "SEMA"
  lr_df_sum$cluster_source = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$source, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster_source = factor(lr_df_sum$cluster_source, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$cluster_target = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster_target = factor(lr_df_sum$cluster_target, levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$clust_id = paste0(lr_df_sum$cluster_source, "_", lr_df_sum$cluster_target)
  lr_df_sum = lr_df_sum[which(lr_df_sum$clust_id %in% this_good$clust_id),]
  lr_df_sum$ct_source = ct_df$ct[match(lr_df_sum$cluster_source, ct_df$cluster)]
  lr_df_sum$ct_target = ct_df$ct[match(lr_df_sum$cluster_target, ct_df$cluster)]
  lr_df_sum$ct_id = paste0(lr_df_sum$ct_source, " -> ", lr_df_sum$ct_target)
  lr_df_sum_num_rows = lr_df_sum %>% group_by(pathway_name, clust_id, ct_id) %>% summarise(num_rows = n())
  lr_df_sum2 = lr_df_sum_num_rows %>% group_by(pathway_name, ct_id) %>% summarise(mean = mean(num_rows))
  # lr_df_sum2 = lr_df_sum %>% group_by(pathway_name, ct_id) %>% summarise(num = n())
  # lr_df_sum$id = paste0(lr_df_sum$cluster_source, "@", lr_df_sum$cluster_target, "@", lr_df_sum$ligand, "@", lr_df_sum$receptor)
  lr_df_sum2$id = paste0(lr_df_sum2$pathway_name, "@", lr_df_sum2$ct_id)
  lr_df_sum_con = lr_df_sum2
  
  lr_df_sum = merge(lr_df_sum_plk, lr_df_sum_con, all = T, by = "id", suffixes = c("_plk", "_con"))
  lr_df_sum[, c("pathway_name", "ct_id")] = reshape2::colsplit(lr_df_sum$id, "@", c('1', '2'))
  lr_df_sum$mean_plk[which(is.na(lr_df_sum$mean_plk))] = 0
  lr_df_sum$mean_con[which(is.na(lr_df_sum$mean_con))] = 0
  # lr_df_sum$prob_plk[which(is.na(lr_df_sum$prob_plk))] = 0
  # lr_df_sum$prob_con[which(is.na(lr_df_sum$prob_con))] = 0
  # lr_df_sum$pathway_name_plk[which(is.na(lr_df_sum$pathway_name_plk))] = lr_df_sum$pathway_name_con[which(is.na(lr_df_sum$pathway_name_plk))] 
  # lr_df_sum$plk.minus.con = lr_df_sum$prob_plk - lr_df_sum$prob_con
  # lr_df_sum[, c("cluster_source", "tmp")] = reshape2::colsplit(lr_df_sum$id, "@", c('1', '2'))
  # lr_df_sum[, c("cluster_target", "tmp2")] = reshape2::colsplit(lr_df_sum$tmp, "@", c('1', '2'))
  # lr_df_sum[, c("ligand", "receptor")] = reshape2::colsplit(lr_df_sum$tmp2, "@", c('1', '2'))
  # lr_df_sum$cluster_source = factor(as.numeric(lr_df_sum$cluster_source), levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  # lr_df_sum$cluster_target = factor(as.numeric(lr_df_sum$cluster_target), levels = as.numeric(c ("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  # lr_df_sum$total_prob = lr_df_sum$prob_plk + lr_df_sum$prob_con
  
  # ggplot(lr_df_sum, aes(x = total_prob, y = plk.minus.con, color = pathway_name_plk)) + geom_point() + theme_classic() + coord_fixed() + ggtitle(paste0("Timepoint: ", this_timepoint)) + geom_label_repel(data = lr_df_sum[which(abs(lr_df_sum$plk.minus.con) > 0.1),], aes(label = id))
  # ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_lr_specific_connections.pdf"), width = 12, height = 8) 
  
  # lr_df_sum_melt = reshape2::melt(lr_df_sum)
  lr_df_sum$plk.minus.con = lr_df_sum$mean_plk-lr_df_sum$mean_con
  my_col_pal = RColorBrewer::brewer.pal(11, "RdBu")
  my_col_pal[6] = "white"
  ggplot(lr_df_sum, aes(x = pathway_name, y = ct_id, fill = plk.minus.con)) + geom_raster() + coord_fixed() + scale_fill_gradientn(colors = rev(my_col_pal), limits = c(-max(abs(lr_df_sum$plk.minus.con)), max(abs(lr_df_sum$plk.minus.con))), oob = scales::squish, na.value = my_col_pal[6]) + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint))
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_pathway_specific_connections.pdf"), width = 10, height = 10)
  
  lr_df_sum = lr_df_sum %>% group_by(ct_id) %>% mutate(pct_plk =  100 *mean_plk/sum(mean_plk), pct_con=  100 *mean_con/sum(mean_con))
  lr_df_sum_melt = reshape2::melt(lr_df_sum)
  lr_df_sum_melt = lr_df_sum_melt[which(lr_df_sum_melt$variable %in% c("pct_plk", "pct_con")),]
  # ggplot(lr_df_sum_melt, aes(x = ct_id, y = value, fill = pathway_name, group=variable)) + geom_bar(stat = "identity", position = position_dodge()) + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint))
  lr_df_sum_melt$variable = plyr::revalue(lr_df_sum_melt$variable, c("pct_plk" = "plk", "pct_con" = "con"))
  all_pathways = sort(unique(lr_df_sum_melt$pathway_name))
  lr_df_sum_melt$pathway_name = factor(lr_df_sum_melt$pathway_name, levels = c(all_pathways[which(all_pathways != "Other")], "Other"))
  lr_df_sum_melt = lr_df_sum_melt %>% arrange(pathway_name)
  lr_df_sum_melt <- lr_df_sum_melt %>% group_by(variable, ct_id) %>% mutate(label_y = 100 - (cumsum(value) - 0.5 * value))
  # lr_df_sum_melt = lr_df_sum_melt %>% arrange(pathway_name)
  # my_cols = c(scales::hue_pal()(length(unique(lr_df_sum_melt$pathway_name))-1), "gray60")
  # my_breaks = c(all_pathways[which(all_pathways != "Other")], "Other")
  # ggplot(lr_df_sum_melt, aes(x = variable, y = value, fill = pathway_name)) + geom_bar(stat = "identity") + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint)) + scale_fill_manual(values = my_cols, breaks = my_breaks) + geom_text(data = lr_df_sum_melt[which(lr_df_sum_melt$value > 20),], aes(x = variable, y = label_y, label = pathway_name), angle = 90, vjust = 0.5) + facet_wrap(~ ct_id,  ncol = nrow(lr_df_sum), labeller = label_wrap_gen(width = 10, multi_line = TRUE))
  ggplot(lr_df_sum_melt, aes(x = variable, y = value, fill = pathway_name)) + geom_bar(stat = "identity") + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint)) + scale_fill_manual(values = c(scales::hue_pal()(length(levels(lr_df_sum_melt$pathway_name))-1), "gray60")) + geom_text(data = lr_df_sum_melt[which(lr_df_sum_melt$value > 15 & lr_df_sum_melt$pathway_name != "Other"),], aes(x = variable, y = label_y, label = pathway_name), angle = 90, vjust = 0.5) + facet_wrap(~ ct_id,  ncol = nrow(lr_df_sum), labeller = label_wrap_gen(width = 10, multi_line = TRUE))
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_pathway_specific_connections_bar.pdf"), width = 12, height = 5)
}


# Individual connections by condition
celltype_vector = c ("epi1", "epi1", "epi1", "epi1", "epi1", "epi1", "epi1", 
                     "epi2", "epi2", "epi2", "epi2", "epi2", "epi2", "epi2", "epi2", 
                     "epi3", "epi3", "epi3", "epi3", "epi3", 
                     "pdl", 
                     "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", 
                     "endothelium", "endothelium", 
                     "macrophage", "macrophage", "macrophage", "macrophage", "macrophage",
                     "bone", "bone",
                     "bone macrophage", 
                     "glia", 
                     "taste bud", 
                     "myoblasts", "myoblasts",
                     "nk and t cells", 
                     "unknown", 
                     "pericytes", "pericytes",
                     "lympho-epithelial", "lympho-epithelial", "lympho-epithelial", "lympho-epithelial", "lympho-epithelial")
ct_df = data.frame(ct = celltype_vector, cluster = as.numeric(c("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
process_lr_df_sum = function(lr_df_sum, this_this_good) {
  lr_df_sum$pathway_name[which(startsWith(lr_df_sum$pathway_name, "SEMA"))] = "SEMA"
  lr_df_sum$cluster_source = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$source, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster_source = factor(lr_df_sum$cluster_source, levels = as.numeric(c("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$cluster_target = as.numeric(reshape2::colsplit(reshape2::colsplit(lr_df_sum$target, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2])
  lr_df_sum$cluster_target = factor(lr_df_sum$cluster_target, levels = as.numeric(c("1", "9", "12", "19", "26", "27", "34", "10", "14", "15", "16", "17", "28", "31", "35", "0", "4", "8", "13", "20", "32", "2", "3", "5", "11", "21", "23", "22", "29", "7", "25", "38", "41", "42", "6", "18", "30", "39", "37", "46", "48", "45", "33", "36", "44", "24", "40", "43", "47", "49")))
  lr_df_sum$clust_id = paste0(lr_df_sum$cluster_source, "_", lr_df_sum$cluster_target)
  lr_df_sum = lr_df_sum[which(lr_df_sum$clust_id %in% this_this_good$clust_id),]
  lr_df_sum$ct_source = ct_df$ct[match(lr_df_sum$cluster_source, ct_df$cluster)]
  lr_df_sum$ct_target = ct_df$ct[match(lr_df_sum$cluster_target, ct_df$cluster)]
  lr_df_sum$ct_id = paste0(lr_df_sum$ct_source, " -> ", lr_df_sum$ct_target)
  lr_df_sum_num_rows = lr_df_sum %>% group_by(pathway_name, clust_id, ct_id) %>% summarise(num_rows = n())
  lr_df_sum2 = lr_df_sum_num_rows %>% group_by(pathway_name, ct_id) %>% summarise(mean = mean(num_rows))
  lr_df_sum2$id = paste0(lr_df_sum2$pathway_name, "@", lr_df_sum2$ct_id)
  return(lr_df_sum2)
}
for (this_timepoint in c("60", "1", "3", "7")) {
  print(this_timepoint)
  lr_df_plk = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_plk_lr_one_050423.csv"))
  lr_df_con = read.csv(paste0("~/research/tooth/results/plkall_plk", this_timepoint, "_con_lr_one_050423.csv"))
  this_good = readxl::read_xlsx("~/Downloads/cc_good_connections_both_cond_talha.xlsx", sheet = paste0("Day", this_timepoint))
  this_good = this_good[, c("Sender", "Reciever", "Plk - High")]
  colnames(this_good)[3] = "plkUp"
  this_good$clust_id = paste0(this_good$Sender, "_", this_good$Reciever)
  this_good$Sender = ct_df$ct[match(this_good$Sender, ct_df$cluster)]
  this_good$Reciever = ct_df$ct[match(this_good$Reciever, ct_df$cluster)]
  this_good$ct_id = paste0(this_good$Sender, " -> ", this_good$Reciever)
  
  lr_df_sum_plk = process_lr_df_sum(lr_df_plk, this_good)
  lr_df_sum_plk$upPlk = lr_df_sum_plk$ct_id %in% this_good$ct_id[which(this_good$plkUp == 1)]
  
  lr_df_sum_con = process_lr_df_sum(lr_df_con, this_good)
  lr_df_sum_con$upPlk = lr_df_sum_con$ct_id %in% this_good$ct_id[which(this_good$plkUp == 1)]
  
  lr_df_sum = merge(lr_df_sum_plk, lr_df_sum_con, all = T, by = "id", suffixes = c("_plk", "_con"))
  lr_df_sum[, c("pathway_name", "ct_id")] = reshape2::colsplit(lr_df_sum$id, "@", c('1', '2'))
  lr_df_sum$mean_plk[which(is.na(lr_df_sum$mean_plk))] = 0
  lr_df_sum$mean_con[which(is.na(lr_df_sum$mean_con))] = 0

  lr_df_sum$plk.minus.con = lr_df_sum$mean_plk-lr_df_sum$mean_con
  lr_df_sum$upPlk_con[which(is.na(lr_df_sum$upPlk_con))] = lr_df_sum$upPlk_plk[which(is.na(lr_df_sum$upPlk_con))]
  lr_df_sum$pathway_name = factor(lr_df_sum$pathway_name)
  my_col_pal = RColorBrewer::brewer.pal(11, "RdBu")
  my_col_pal[6] = "white"
  con_ncol = length(unique(lr_df_sum$ct_id[which(lr_df_sum$upPlk_con == 0)]))
  plk_ncol = length(unique(lr_df_sum$ct_id[which(lr_df_sum$upPlk_con == 1)]))
  all_nrow = length(unique(lr_df_sum$pathway_name))
  # ggplot(lr_df_sum, aes(x = pathway_name, y = ct_id, fill = plk.minus.con)) + geom_tile(width = 1, height = 1) + scale_fill_gradientn(colors = rev(my_col_pal), limits = c(-max(abs(lr_df_sum$plk.minus.con)), max(abs(lr_df_sum$plk.minus.con))), oob = scales::squish, na.value = my_col_pal[6]) + theme_classic() + theme(aspect.ratio = 0.2) + ggtitle(paste0("Timepoint: ", this_timepoint)) + facet_wrap(~ upPlk_con, ncol = 1, scales = "free_y")
  p1 = ggplot(lr_df_sum[which(lr_df_sum$upPlk_con == 0),], aes(x = pathway_name, y = ct_id, fill = plk.minus.con)) + geom_raster() + scale_x_discrete(drop = F) + scale_fill_gradientn(colors = rev(my_col_pal), limits = c(-max(abs(lr_df_sum$plk.minus.con)), max(abs(lr_df_sum$plk.minus.con))), oob = scales::squish, na.value = my_col_pal[6]) + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint)) + theme(axis.line.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x.bottom = element_blank(), axis.title.x = element_blank()) + force_panelsizes(cols = unit(all_nrow/6, "in"), rows = unit(con_ncol/6, "in"))
  p1
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_pathway_specific_connections_con.pdf"), width = 13, height = 7)
  p2 = ggplot(lr_df_sum[which(lr_df_sum$upPlk_con == 1),], aes(x = pathway_name, y = ct_id, fill = plk.minus.con)) + geom_raster() + scale_x_discrete(drop = F) + scale_fill_gradientn(colors = rev(my_col_pal), limits = c(-max(abs(lr_df_sum$plk.minus.con)), max(abs(lr_df_sum$plk.minus.con))), oob = scales::squish, na.value = my_col_pal[6]) + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint)) + force_panelsizes(cols = unit(all_nrow/6, "in"), rows = unit(plk_ncol/6, "in"))  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  p2
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_pathway_specific_connections_plk.pdf"), width = 13, height = 7)
  # cowplot::plot_grid(plotlist=list(p1, p2), ncol = 1)
  
  lr_df_sum = lr_df_sum %>% group_by(ct_id) %>% mutate(pct_plk =  100 *mean_plk/sum(mean_plk), pct_con=  100 *mean_con/sum(mean_con))
  lr_df_sum_melt = reshape2::melt(lr_df_sum)
  lr_df_sum_melt = lr_df_sum_melt[which(lr_df_sum_melt$variable %in% c("pct_plk", "pct_con")),]
  lr_df_sum_melt$variable = plyr::revalue(lr_df_sum_melt$variable, c("pct_plk" = "plk", "pct_con" = "con"))
  all_pathways = sort(unique(lr_df_sum_melt$pathway_name))
  lr_df_sum_melt$pathway_name = factor(lr_df_sum_melt$pathway_name, levels = c(all_pathways[which(all_pathways != "Other")], "Other"))
  lr_df_sum_melt = lr_df_sum_melt %>% arrange(pathway_name)
  lr_df_sum_melt <- lr_df_sum_melt %>% group_by(variable, ct_id) %>% mutate(label_y = 100 - (cumsum(value) - 0.5 * value))
  ggplot(lr_df_sum_melt, aes(x = variable, y = value, fill = pathway_name)) + geom_bar(stat = "identity") + theme_classic() + ggtitle(paste0("Timepoint: ", this_timepoint)) + scale_fill_manual(values = c(scales::hue_pal()(length(levels(lr_df_sum_melt$pathway_name))-1), "gray60")) + geom_text(data = lr_df_sum_melt[which(lr_df_sum_melt$value > 15 & lr_df_sum_melt$pathway_name != "Other"),], aes(x = variable, y = label_y, label = pathway_name), angle = 90, vjust = 0.5) + facet_wrap(~ ct_id,  ncol = nrow(lr_df_sum), labeller = label_wrap_gen(width = 10, multi_line = TRUE))
  ggsave(paste0("~/research/tooth/results/plk", this_timepoint, "_pathway_specific_connections_bar.pdf"), width = 12, height = 5)
}

#*******************************************************************************
# Module Scores at Timepoints ==================================================
#*******************************************************************************
dbscan_all = read.csv("~/research/tooth/results/plkall_modules_hgnc.csv")
this_mod = 16
dbmod = dbscan_all[which(dbscan_all$dbscan == this_mod),]
plk_subject$mod_score = colSums(plk_subject@assays$RNA@counts[dbmod$gene,] > 0)

#*******************************************************************************
# Time * Cond ==================================================================
#*******************************************************************************
# Load DEGs
deg = read.csv(paste0(out_dir, "/sig_082923.csv"))
deg$hgnc = gene_info$human[match(deg$gene, gene_info$seurat_name)]

# Fold Change
deg$y_plk60_dif = log2(deg$y_plk_plk60 / deg$y_con_plk60)
deg$y_plk1_dif  = log2(deg$y_plk_plk1  / deg$y_con_plk1)
deg$y_plk3_dif  = log2(deg$y_plk_plk3  / deg$y_con_plk3)
deg$y_plk7_dif  = log2(deg$y_plk_plk7  / deg$y_con_plk7)
# deg$y_plk60_dif = deg$y_plk_plk60 - deg$y_con_plk60
# deg$y_plk1_dif  = deg$y_plk_plk1  - deg$y_con_plk1
# deg$y_plk3_dif  = deg$y_plk_plk3  - deg$y_con_plk3
# deg$y_plk7_dif  = deg$y_plk_plk7  - deg$y_con_plk7

# Plot
deg$id = paste0(deg$gene, "_", deg$cluster)
rownames(deg) = deg$id
# mean_cols = c("y_con_plk1", "y_plk_plk1", "y_con_plk3", "y_plk_plk3", "y_con_plk60", "y_plk_plk60", "y_con_plk7", "y_plk_plk7")
mean_cols = c("y_plk60_dif", "y_plk1_dif", "y_plk3_dif", "y_plk7_dif")
deg_p = deg[, c("id", "gene", "hgnc", "cluster", "P_cond.exp", "bh", mean_cols)]
# deg_p[,mean_cols] = t(scale(t(deg_p[,mean_cols])))
col_order = hclust(dist(deg_p[,mean_cols]), method = "complete")
deg_p = reshape2::melt(deg_p, id.var = c("id", "gene", "hgnc", "cluster", "P_cond.exp", "bh"))
deg_p$id = factor(deg_p$id, levels = col_order$labels[col_order$order])
deg_p$variable = plyr::revalue(deg_p$variable, c("y_plk60_dif" = "0", "y_plk1_dif" = "1", "y_plk3_dif" = "3", "y_plk7_dif" = "7"))
# ggplot(deg_p, aes(x = id, y = variable, fill = value)) + geom_raster() + coord_fixed() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(deg_p$value)), max(abs(deg_p$value))), oob = scales::squish) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot(deg_p, aes(x = id, y = variable, fill = value)) + geom_raster() + coord_fixed() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(deg_p$value)), max(abs(deg_p$value))), oob = scales::squish, name = expression(""*Log["2"]*"FC")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ylab("Timepoint") + xlab("")
ggsave("~/research/tooth/results/time_interaction_cond_deg_heatmap.pdf", width = 12 , height = 6)

# Add in difference in percent expression
deg[, c("plk60_pct", "plk1_pct", "plk3_pct", "plk7_pct")] = NA
Idents(plk_subject) = paste0(plk_subject$seurat_clusters, plk_subject$exp, plk_subject$cond)
for (i in 1:nrow(deg)) {
  for (j in c("plk60", "plk1", "plk3", "plk7")) {
    plk_pct_exp = FoldChange(plk_subject, features = deg_p$gene[i], ident.1 = paste0(deg$cluster[i], j, "plk"))
    con_pct_exp = FoldChange(plk_subject, features = deg_p$gene[i], ident.1 = paste0(deg$cluster[i], j, "con"))
    pct_dif = plk_pct_exp$pct.1 - con_pct_exp$pct.1
    deg[i, paste0(j, "_pct")] = pct_dif
  }
}
deg_pct = reshape2::melt(deg[, c('id', "plk60_pct", "plk1_pct", "plk3_pct", "plk7_pct")], id.var = c("id"))
deg_p$pct_dif = deg_pct$value
ggplot(deg_p, aes(x = id, y = variable, color = value, size = abs(pct_dif))) + geom_point() + geom_point(data=deg_p[which(deg_p$pct_dif < 0),], color = RColorBrewer::brewer.pal(11, "RdBu")[11], shape=1) + geom_point(data=deg_p[which(deg_p$pct_dif > 0),], color = RColorBrewer::brewer.pal(11, "RdBu")[1], shape=1) + coord_fixed() + scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(deg_p$value)), max(abs(deg_p$value))), oob = scales::squish, name = expression(""*Log["2"]*"FC")) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))  + ylab("Timepoint") + xlab("")
ggsave("~/research/tooth/results/time_interaction_cond_deg_dot.pdf", width = 12 , height = 6)

# Todd's line
ggplot(deg_p, aes(x = variable, y = value, color = id, group = id)) + geom_point() + geom_line()

deg_p2 = deg_p
deg_p2$up = "none"
deg_p2$up[which(deg_p2$id %in% deg$id[which(deg$y_plk60_dif > 0)])] = "Up in 60"
deg_p2$up[which(deg_p2$id %in% deg$id[which(deg$y_plk1_dif > 0)])] = "Up in 1"
deg_p2$up[which(deg_p2$id %in% deg$id[which(deg$y_plk3_dif > 0)])] = "Up in 3"
deg_p2$up[which(deg_p2$id %in% deg$id[which(deg$y_plk7_dif > 0)])] = "Up in 7"
ggplot(deg_p2, aes(x = variable, y = value, color = id, group = id)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_point() + geom_line() + facet_wrap(~ up)  + theme_bw() + ylab(expression(""*Log["2"]*"FC")) + xlab("")
ggsave("~/research/tooth/results/time_interaction_cond_deg_line.pdf", width = 14 , height = 7)

#*******************************************************************************
# Talha: Overlap of Modules and DEGs ===========================================
#*******************************************************************************
mod = read.csv("~/scratch/d_tooth/results/plkall_modules_hgnc.csv")
deg = read.csv("~/scratch/d_tooth/results/deg/plk60_plk_v_ctrl_cluster_level_sig_pct.csv")
deg[,1:2] = NULL
mod_deg = merge(mod[,c("gene", "dbscan")], deg, by = "gene")
mod_deg = mod_deg[,c(1, ncol(mod_deg), 2:(ncol(mod_deg)-1))]
write.csv(mod_deg, "~/scratch/d_tooth/results/deg/mod_ovlp_w_plk60_hdeg.csv")

deg = read.csv("~/scratch/d_tooth/results/plk60_plk_v_ctrl_bulk_sig_pct.csv")
deg[,1:2] = NULL
mod_deg = merge(mod[,c("gene", "dbscan")], deg, by.x = "gene", by.y = "X")
mod_deg = mod_deg[,c(1, ncol(mod_deg), 2:(ncol(mod_deg)-1))]
write.csv(mod_deg, "~/scratch/d_tooth/results/deg/mod_ovlp_w_plk60_bulk_hdeg.csv")

#*******************************************************************************
# Diff Cor =====================================================================
#*******************************************************************************
plk_subject$my_sub = paste0(plk_subject$subject, "_", plk_subject$cond)
df2 = data.table::fread("~/scratch/d_tooth/data/all_ind_cond_cor.csv", data.table = F)
ind_counts = lapply(colnames(df2)[2:ncol(df2)], function(x) rowSums(plk_subject@assays$RNA@counts[,which(plk_subject$my_sub == x)] > 0))
ind_counts = do.call("cbind", ind_counts)
colnames(ind_counts) = colnames(df2)[2:ncol(df2)]
df2[,c("gene1", "gene2")] = reshape2::colsplit(df2$id, "_", c('1', '2'))
my_sub_tot_counts = data.frame(table(plk_subject$my_sub))

all_dif_cor = data.frame()
for (this_exp in exps) {
  print(this_exp)
  this_my_sub = colnames(df2)[which(grepl(this_exp, colnames(df2)))]
  df_exp = df2[,c('id', 'gene1', 'gene2', this_my_sub)]
  
  # Find genes expressed in at least 5 cells in all subjects
  ind_exp = ind_counts[,this_my_sub]
  big_genes = rowSums(ind_exp > 5)
  big_genes = names(big_genes)[which(big_genes == ncol(df_exp)-3)]
  df_exp = df_exp[which(df_exp$gene1 %in% big_genes & df_exp$gene2 %in% big_genes),]
  df_exp = df_exp[which(df_exp$gene1 != df_exp$gene2),]
  df_exp$min_gene1_pos = apply(ind_exp[match(df_exp$gene1, rownames(ind_exp)),], 1, min)
  df_exp$min_gene2_pos = apply(ind_exp[match(df_exp$gene2, rownames(ind_exp)),], 1, min)
  df_exp$min_gene_pos = apply(df_exp[,c("min_gene1_pos", "min_gene2_pos")], 1, min)
  
  # Find gene-combos where all subjects are up in one condition
  plk_max = apply(df_exp[,which(endsWith(colnames(df_exp), "plk"))], 1, max)
  plk_min = apply(df_exp[,which(endsWith(colnames(df_exp), "plk"))], 1, min)
  con_max = apply(df_exp[,which(endsWith(colnames(df_exp), "con"))], 1, max)
  con_min = apply(df_exp[,which(endsWith(colnames(df_exp), "con"))], 1, min)
  con_up_idx = which(con_min > plk_max)
  plk_up_idx = which(plk_min > con_max)
  df_exp = df_exp[c(con_up_idx, plk_up_idx),]
  df_exp$up_cond = c(rep("con", length(con_up_idx)), rep("plk", length(plk_up_idx)))
  
  # Find the minimum difference in correlations
  df_exp$up_dist = con_min[c(con_up_idx, plk_up_idx)] - plk_max[c(con_up_idx, plk_up_idx)]
  df_exp$up_dist[which(df_exp$up_cond == "plk")] = plk_min[plk_up_idx] - con_max[plk_up_idx]
  
  # Find the mean difference in correlations
  df_exp$mean_dif = rowMeans(df_exp[,which(endsWith(colnames(df_exp), "con"))]) - rowMeans(df_exp[,which(endsWith(colnames(df_exp), "plk"))]) 
  df_exp$mean_dif[which(df_exp$up_cond == "plk")] = rowMeans(df_exp[which(df_exp$up_cond == "plk"),which(endsWith(colnames(df_exp), "plk"))]) - rowMeans(df_exp[which(df_exp$up_cond == "plk"),which(endsWith(colnames(df_exp), "con"))]) 
  
  # Add in the sample size
  for (j in this_my_sub) { df_exp[,paste0(j, "_counts")] = my_sub_tot_counts$Freq[which(my_sub_tot_counts$Var1 == j)] }
  
  # Perform Correlation Significance Test on all combinations of subjects
  plk_cols = this_my_sub[which(endsWith(this_my_sub, "plk"))]
  con_cols = this_my_sub[which(endsWith(this_my_sub, "con"))]
  for (plk_sub in plk_cols) {
    for (con_sub in con_cols) {
      this_p = unlist(parallel::mclapply(1:nrow(df_exp), function(x) r_to_p(df_exp[x, plk_sub], df_exp[x, con_sub], df_exp[x, paste0(plk_sub, "_counts")], df_exp[x, paste0(con_sub, "_counts")]), mc.cores = 20))
      df_exp[,paste0(plk_sub, "_", con_sub, "_p")] = this_p
      df_exp[,paste0(plk_sub, "_", con_sub, "_bh")] = p.adjust(this_p, method = 'BH')
    }
  }
  bh_cols = colnames(df_exp)[which(grepl("bh", colnames(df_exp)))]
  df_exp_sig = df_exp[which( rowSums(df_exp[,bh_cols] < 0.05) == length(bh_cols) ),]
  df_exp_sig$bh_max = apply(df_exp_sig[,bh_cols], 1, max)
  df_exp_sig$mean_plk_cor = rowMeans(df_exp_sig[,plk_cols])
  df_exp_sig$mean_con_cor = rowMeans(df_exp_sig[,con_cols])
  df_exp_sig$exp = this_exp
  all_dif_cor = rbind(all_dif_cor, df_exp_sig[,c("exp", "id", "gene1", "gene2", "min_gene_pos", "mean_plk_cor", "mean_con_cor", "up_cond", "mean_dif", "up_dist", "bh_max")])
}
all_dif_cor$hgnc1 = gene_info$human[match(all_dif_cor$gene1, gene_info$seurat_name)]
all_dif_cor$hgnc2 = gene_info$human[match(all_dif_cor$gene2, gene_info$seurat_name)]
write.csv(all_dif_cor, "~/scratch/d_tooth/results/cor_more_pos_bh.csv")

r_to_p = function(r1, r2, n1, n2) {
  # Compare Two Correlation Values using Fisher's Z Transformation Method.
  #' @param r1 correlation 1
  #' @param r2 correlation 2
  #' @param n1 number of samples used for correlation 1
  #' @param n2 number of samples used for correlation 2
  z1 = .5 * (log(1+r1) - log(1-r1))
  z2 = .5 * (log(1+r2) - log(1-r2))
  z3 =(z1 - z2) / sqrt( (1 / (n1 - 3)) + (1 / (n2 - 3)) )
  p = 2*pnorm(-abs(z3))
  return(p)
}

#*******************************************************************************
# Supplement ===================================================================
#*******************************************************************************
plk$remove = is.na(plk$subject_num)
df_remove = unclass(table(plk$seurat_clusters, plk$remove))
df_remove = df_remove / rowSums(df_remove)
df_remove = reshape2::melt(df_remove)
df_remove$Var1 = factor(df_remove$Var1, levels = sort(unique(df_remove$Var1)))
pdf("~/research/tooth/results/demux_remove_cells_by_cluster.pdf", width = 10, height = 3)
ggplot(df_remove, aes(x = Var1, y = value*100, fill = Var2)) + geom_bar(stat = 'identity') + theme_classic() + xlab("") + scale_y_continuous(expand = c(0,0), name = "Percent of Nuclei") + scale_fill_manual(values = c("goldenrod1", "gray40")) + NoLegend()
dev.off()

df_remove = unclass(table(plk$exp, plk$remove))
df_remove = df_remove / rowSums(df_remove)
df_remove = reshape2::melt(df_remove)
df_remove$Var1 = factor(df_remove$Var1, levels = sort(unique(df_remove$Var1)))
pdf("~/research/tooth/results/demux_remove_cells_by_exp.pdf", width = 3, height = 3)
ggplot(df_remove, aes(x = Var1, y = value*100, fill = Var2)) + geom_bar(stat = 'identity') + theme_classic() + xlab("") + scale_y_continuous(expand = c(0,0), name = "Percent of Nuclei") + scale_fill_manual(values = c("goldenrod1", "gray40")) + NoLegend()
dev.off()

df_remove = unclass(table(plk$cond, plk$remove))
df_remove = df_remove / rowSums(df_remove)
df_remove = reshape2::melt(df_remove)
df_remove$Var1 = factor(df_remove$Var1, levels = sort(unique(df_remove$Var1)))
pdf("~/research/tooth/results/demux_remove_cells_by_cond.pdf", width = 2, height = 3)
ggplot(df_remove, aes(x = Var1, y = value*100, fill = Var2)) + geom_bar(stat = 'identity') + theme_classic() + xlab("") + scale_y_continuous(expand = c(0,0), name = "Percent of Nuclei") + scale_fill_manual(values = c("goldenrod1", "gray40")) + NoLegend()
dev.off()

