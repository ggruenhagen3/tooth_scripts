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
plk_subject = readRDS(paste0(data_dir, "plkall_subject_053023.rds"))

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
exps = c("plk60", "plk1", "plk3", "plk7")
conds = c("plk", "con")
for (this_exp in exps) {
  for (this_cond in conds) {
    cc_df = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_", this_cond, "_one_050423.csv"))
    cc_df$Sender   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
    cc_df$Receiver = reshape2::colsplit(reshape2::colsplit(cc_df$Receiver, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
    cc_df$Sender   = factor(cc_df$Sender,   levels = 0:max(as.numeric(cc_df$Sender)))
    cc_df$Receiver = factor(cc_df$Receiver, levels = 0:max(as.numeric(cc_df$Receiver)))
    
    fname = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_", this_cond, "_one_050423.pdf")
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
  cc_df$Sender   = factor(cc_df$Sender,   levels = 0:max(as.numeric(cc_df$Sender)))
  cc_df$Receiver = factor(cc_df$Receiver, levels = 0:max(as.numeric(cc_df$Receiver)))
  cc_df$id = paste0(cc_df$exp, "_", cc_df$Sender, "_", cc_df$Receiver)
  cc_df_plk = cc_df
  
  cc_df = read.csv(paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_con_one_050423.csv"))
  cc_df$exp = reshape2::colsplit(cc_df$Sender, "_", c('1', '2'))[,1]
  cc_df$cond   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,1]
  cc_df$Sender   = reshape2::colsplit(reshape2::colsplit(cc_df$Sender,   "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  cc_df$Receiver = reshape2::colsplit(reshape2::colsplit(cc_df$Receiver, "_", c('1', '2'))[,2], "_", c('1', '2'))[,2]
  cc_df$Sender   = factor(cc_df$Sender,   levels = 0:max(as.numeric(cc_df$Sender)))
  cc_df$Receiver = factor(cc_df$Receiver, levels = 0:max(as.numeric(cc_df$Receiver)))
  cc_df$id = paste0(cc_df$exp, "_", cc_df$Sender, "_", cc_df$Receiver)
  cc_df_con = cc_df
  
  cc_df = cc_df_plk
  cc_df$plk = cc_df$value
  cc_df$con = cc_df_con$value[match(cc_df$id, cc_df_con$id)]
  cc_df$dif = cc_df$plk - cc_df$con
  cc_df = cc_df[which(!is.na(cc_df$dif)),]
  # fname = paste0("~/scratch/d_tooth/results/plkall_", this_exp, "_plkvcon_one_050423.pdf")
  # ggplot(cc_df, aes(x = Sender, y = Receiver, fill = dif)) + geom_raster() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-max(abs(cc_df$dif)), max(abs(cc_df$dif))), oob = scales::squish) + ggtitle(paste0("CellChat Plucked - Control Weights for Experiment: ", this_exp)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1)) + coord_fixed() + ggh4x::force_panelsizes(cols = unit(length(unique(cc_df$Sender))/7, "in"), rows = unit(length(unique(cc_df$Sender))/7, "in"))
  # ggsave(fname, width = 9, height = 9)
  # system(paste0("rclone copy ", fname, " dropbox:BioSci-Streelman/George/Tooth/plk/analysis/cellchat"))
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
# Glmmseq ======================================================================
#*******************************************************************************
library(glmmSeq)
library(Seurat)
library(lme4)
library(lmerTest)
library(parallel)
library(edgeR)
library(DESeq2)
library(glmGamPoi)
library(scran)
# library(BiocParallel)

data = plk_subject
data$pair = data$subject
obj = data
# this_cells = colnames(data)[which(data$seurat_clusters == 0)]
this_cells = colnames(data)

# this_cells  = colnames(data) # TODO
this_counts = data@assays$RNA@counts[,this_cells]
this_meta   = data.frame(data@meta.data[this_cells,])

gene_not_present_in_pairs = lapply(unique(this_meta$pair), function(x) rowSums(this_counts[,which(this_meta$pair == x)]) == 0)
gene_not_present_in_pairs = Reduce(`+`, gene_not_present_in_pairs)
genes_present_in_all_pairs = names(gene_not_present_in_pairs[which(gene_not_present_in_pairs == 0)])
this_counts = this_counts[genes_present_in_all_pairs,]

this_max_cluster_size = max(table(this_meta$seurat_clusters)) # TODO make clusters as a variable
multicoreParam <- MulticoreParam(workers = 24)
dds = DESeqDataSetFromMatrix(countData = this_counts, colData = this_meta, design = ~ cond)
size_factors = calculateSumFactors(this_counts, BPPARAM = multicoreParam, max.cluster.size = this_max_cluster_size, clusters = NULL, ref.clust = NULL, positive = TRUE, scaling = NULL,  min.mean = NULL, subset.row = NULL)
# size_factors = calculateSumFactors(this_counts, max.cluster.size = this_max_cluster_size, clusters = NULL, ref.clust = NULL, positive = TRUE, scaling = NULL,  min.mean = NULL, subset.row = NULL)
sizeFactors(dds) = size_factors
dds = estimateDispersions(dds, fitType = "glmGamPoi", useCR = TRUE, maxit = 100, weightThreshold = 0.01, quiet = FALSE, modelMatrix = NULL, minmu = 1e-06)
disp = as.matrix(mcols(dds))
disp = disp[,11]
names(disp) = genes_present_in_all_pairs

results <- glmmSeq(~ cond + (1|sample/subject), id = "subject",
                   countdata = this_counts,
                   metadata = this_meta,
                   dispersion = disp,
                   removeSingles=FALSE,
                   progress=TRUE,
                   cores = 24)


#
# options(MulticoreParam=quote(MulticoreParam(workers=4)))
# FUN <- function(x) { round(sqrt(x), 4) }
# FUN1 <- function(x) { print(x); sleep(1); }
# bplapply(1:4, FUN)
# bplapply(1:8, FUN1)

glmm_out_dir = "~/scratch/d_tooth/results/plk_glmmseq_plk1_clusters50/"
big_res = data.frame()
for (this_clust in sort(unique(obj$seurat_clusters))) {
  this_file = paste0(glmm_out_dir, "cluster_", this_clust, ".csv")
  if ( file.exists(this_file) ) {
    res = read.csv(this_file)
    res$cluster = this_clust
    big_res = rbind(res, big_res)
  }
}
big_res$bh = p.adjust(big_res$P_cond, method = "BH")
deg_sig = big_res[which(big_res$bh < 0.05),]
head(deg_sig[order(deg_sig$bh, decreasing = F),])
deg_sig$hgnc = gene_info$human[match(deg_sig$X, gene_info$seurat_name)]
write.csv(deg_sig, paste0(glmm_out_dir, "all_sig.csv"))
