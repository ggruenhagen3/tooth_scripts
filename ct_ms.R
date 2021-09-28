#=============================================================================================
# Figure 1 ===================================================================================
#=============================================================================================
# *** DotPlot for Tooth ***
# Tooth Data
# deg_tj = read.table("~/research/tooth/results/tj_annot_unique_long_hgnc_100.tsv", sep="\t", header = T, stringsAsFactors = F)
# deg_tj = deg_tj[-which(deg_tj$gene == "ENSMZEG00005000005"),]
# deg_tj2 = as.data.frame(deg_tj %>% group_by(cluster) %>% slice(head(row_number(), 2)))
# deg_tj2 = deg_tj2[order(deg_tj2$cluster),]

deg_tj2 = data.frame(cluster = c(rep("Glia", 3), rep("Immune", 4), rep("Epithelial", 7), rep("Mesenchymal", 8), rep("Endothelial", 6)),
                     gene = c("tcf7", "tmtc2b", "slit3", "ENSMZEG00005022086", "wdr73", "ENSMZEG00005022005", "HBE1 (1 of many).1", "pitx1", "pitx2", "shha", "odam", "lbh", "itgb3b", 'sema3fb', "plod2 (1 of many)", "sox5", "smpd3", "MMP16", "bmp6", "SMAD6", "sall1a", "dkk1a", "tns2a", "EBF1 (1 of many).1", "tie1", "PTPRB (1 of many).1", "cdh5", "robo4"))
deg_tj2$hgnc = toupper(deg_tj2$gene)
deg_tj2$hgnc = str_replace_all(deg_tj2$hgnc, " \\(1 OF MANY\\)\\.1", "")
deg_tj2$hgnc = str_replace_all(deg_tj2$hgnc, " \\(1 OF MANY\\)", "")
deg_tj2$hgnc[which(deg_tj2$gene == "ENSMZEG00005022086")] = "HBB"
deg_tj2$hgnc[which(deg_tj2$gene == "ENSMZEG00005022005")] = "HBE1A"
deg_tj2$hgnc[which(deg_tj2$gene == "HBE1 (1 of many).1")] = "HBE1B"
deg_tj2 = deg_tj2[order(deg_tj2$cluster),]

# Cluster Info
convert_tj = data.frame(old = c(0, 1, 2, 3, 4, 5, 6, 7),
                        new = c("Immune", "Glia", "Epithelial", "Epithelial", "Epithelial", "Mesenchymal", "Endothelial", "Mesenchymal"),
                        col = c("#E9C46A", "#F6B179", "#2A9D8F", "#2A9D8F", "#2A9D8F", "#E24D28", "#264653", "#E24D28"))

# Set colors for plot
deg_tj2$col = convert_tj$col[match(deg_tj2$cluster, convert_tj$new)]
xtext_unique_cols = rev(unique(convert_tj$col[order(convert_tj$new)]))
xtext_cols = c( as.vector(deg_tj2$col) )
pal = rev(brewer.pal(11, "RdYlBu"))
pal = colorRampPalette(pal)

# Labels for genes
deg_tj2$label = tolower(deg_tj2$hgnc)
deg_tj2$label[which(is.na(deg_tj2$label))] = deg_tj2$gene[which(is.na(deg_tj2$label))]
deg_tj2$label[which(deg_tj2$gene == "EBF1 (1 of many).1")] = "ebf1"
all_labels = deg_tj2$label

# Convert cluster IDs
Idents(tj) = factor(Idents(tj), levels = rev(sort(levels(Idents(tj)))))
dotp = DotPlot(tj, features = deg_tj2$gene) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols, face = "italic"), axis.text.y = element_text(colour = xtext_unique_cols)) + scale_x_discrete(labels = deg_tj2$gene) + scale_x_discrete(labels = all_labels)
print(dotp)

# DotPlot
pdf("~/research/tooth/results/figure_1_dotplot_tj_1.pdf", width = 9, height = 5)
print(dotp + xlab("") + ylab(""))
dev.off()
pdf("~/research/tooth/results/figure_1_dotplot_tj_1_big.pdf", width = 14, height = 5)
print(dotp + xlab("") + ylab(""))
dev.off()

# *** Dotplot for Jaw ***
# Jaw Data
# deg_jaw = read.table("~/research/tooth/results/jaw_annot_unique_long_hgnc_100.tsv", sep="\t", header = T, stringsAsFactors = F)
# deg_jaw = deg_jaw[-which(deg_jaw$gene == "zgc:163121"),]
# deg_jaw = deg_jaw[-which(deg_jaw$gene == "ENSMZEG00005011359"),]
# deg_jaw2 = as.data.frame(deg_jaw %>% group_by(cluster) %>% slice(head(row_number(), 2)))
# jaw_order = rev(c("Mature TB", "Immature TB", "Epithelial", "Pigmented", "Immune", "Mesenchymal"))
# deg_jaw2 = deg_jaw2[order(deg_jaw2$cluster),]
# deg_jaw2 = deg_jaw2[order(match(deg_jaw2$cluster, rev(jaw_order))),]

deg_jaw2 = data.frame(cluster = c(rep("Epithelial", 8), rep("Pigmented", 4), rep("Mesenchymal", 12), rep("Mature TB", 5), rep("Immature TB", 8), rep("Immune", 5)),
                      gene = c("krt5", "krt15", "itga6b", "fgfr2", "epcam (1 of many).1", "tp63 (1 of many).2", "klf6a", "pitx1", "mitfa", "mrc1b (1 of many).3", "slc43a2a", "csf1ra", "COL6A3", "col5a1", "epha3", "prrx1b", "itga1", "sfrp2", "smoc1", "col6a2", "PIEZO2", "col6a1", "synm", "bgnb", "calb2a", "trpm5", "avil", "hcn1", "ENSMZEG00005005666", "sox2", "spire2", "dpp6b", "samd10b", "chl1b (1 of many)", "eps8a", "kif26ba", "CALCB", "tox2", "sh2d3ca", "PRKD3 (1 of many)", "skap1", "jak3"))
deg_jaw2$hgnc = toupper(deg_jaw2$gene)
deg_jaw2$hgnc = str_replace_all(deg_jaw2$hgnc, " \\(1 OF MANY\\)\\.\\d", "")
deg_jaw2$hgnc = str_replace_all(deg_jaw2$hgnc, " \\(1 OF MANY\\)", "")
deg_jaw2$hgnc[which(deg_jaw2$gene == "ENSMZEG00005005666")] = "COL2A1"
jaw_order = rev(c("Mature TB", "Immature TB", "Epithelial", "Pigmented", "Immune", "Mesenchymal"))
deg_jaw2 = deg_jaw2[order(match(deg_jaw2$cluster, rev(jaw_order))),]

# Cluster Info
convert_jaw = data.frame(old = 0:11,
                         new = c("Epithelial", "Epithelial", "Immune", "Epithelial", "Epithelial", "Pigmented", "Mesenchymal", "Epithelial", "Immune", "Immune", "Mature TB", "Immature TB"),
                         col = c("#2A9D8F", "#2A9D8F", "#E9C46A", "#2A9D8F", "#2A9D8F", "#F6B179", "#E24D28", "#2A9D8F", "#E9C46A", "#E9C46A", "#264653", "#264653"),
                         order = c(3, 3, 5, 3, 3, 4, 6, 3, 5, 5, 1, 2))

# Set colors for plot
deg_jaw2$col = convert_jaw$col[match(deg_jaw2$cluster, convert_jaw$new)]
xtext_unique_cols = convert_jaw$col[match(jaw_order, convert_jaw$new)]
xtext_cols = c( as.vector(deg_jaw2$col) )
pal = rev(brewer.pal(11, "RdYlBu"))
pal = colorRampPalette(pal)

# Labels for genes
deg_jaw2$label = tolower(deg_jaw2$hgnc)
deg_jaw2$label[which(is.na(deg_jaw2$label))] = deg_jaw2$gene[which(is.na(deg_jaw2$label))]
deg_jaw2$label[which(deg_jaw2$gene == "mrc1b (1 of many).3")] = "mrc1"
all_labels = deg_jaw2$label

# Convert cluster IDs
Idents(jaw) = factor(Idents(jaw), levels = jaw_order)
dotp = DotPlot(jaw, features = deg_jaw2$gene) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols, face = "italic"), axis.text.y = element_text(colour = xtext_unique_cols)) + scale_x_discrete(labels = deg_jaw2$gene) + scale_x_discrete(labels = all_labels)
print(dotp)

# DotPlot
pdf("~/research/tooth/results/figure_1_dotplot_jaw_1.pdf", width = 8, height = 5)
print(dotp + xlab("") + ylab(""))
dev.off()
pdf("~/research/tooth/results/figure_1_dotplot_jaw_1_big.pdf", width = 14, height = 5)
print(dotp + xlab("") + ylab(""))
dev.off()

# *** Cell Types ***
# Tj + Jaw
# tj_jaw$seurat_clusters = factor(tj_jaw$seurat_clusters, levels = 0:10)
# Idents(tj_jaw) = tj_jaw$seurat_clusters
# # tj_jaw = RenameIdents(tj_jaw, '0' = "Epithelial", '1' = "Epithelial", '2' = "Mesenchymal", '3' = "Immune", '4' = "Immune", "5" = "Epithelial", "6" = "Epithelial", "7" = "Shared", "8"= "Pigmented", "9" = "Mature TB", "10" = "Immature TB")
# # tj_jaw = RenameIdents(tj_jaw, '0' = "Epithelial.J1", '1' = "Epithelial.J2", '2' = "Mesenchymal", '3' = "Immune.T1", '4' = "Immune.J1", "5" = "Epithelial.T1", "6" = "Epithelial.J3", "7" = "Shared", "8"= "Pigmented", "9" = "Mature TB", "10" = "Immature TB")
# svg("~/research/tooth/results/tj_jaw_annot_2.svg")
# pal = colorRampPalette(c("#0a3255ff", "#e6d035ff"))(11)
# pal = colorRampPalette(c("#b7094c", "#0091ad"))(11)
# # pal = c("#9d0208", "#d00000", "#dc2f02", "#e85d04", "#f48c06", "#F2D444", "black", "#432371", "#023e8a", "#0077b6", "#0096c7")
# pal = c("#05668D", "#028090", "#00A896", "#01B698", "#02C39A", "#79DBAC", "#F0F3BD", "#432371", "#023e8a", "#0077b6", "#4ea8de")
# names(pal) = c(1,10,9,0,6,8,3,7,2,4,5)
# pal = pal[order(as.numeric(names(pal)))]
# DimPlot(tj_jaw, pt.size = 1.5, label = F, cols = pal) + coord_fixed() + NoLegend()
# dev.off()

# TJ
tj$annot = tj$seurat_clusters
Idents(tj) = tj$annot
tj = RenameIdents(tj, '0' = "Immune", '1' = "Glia", '2' = "Epithelial", '3' = "Epithelial", '4' = "Epithelial", "5" = "Mesenchymal", "6" = "Endothelial", "7" = "Mesenchymal")
svg("~/research/tooth/results/tj_annot.svg")
pal = convert_tj$col[match(levels(Idents(tj)), convert_tj$new)]
# pal = c("#E9C46A", "#F6B179", "#2A9D8F", "#E24D28", "#264653")
print(DimPlot(tj, order = T, label = F, pt.size = 1.5, cols = pal) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
dev.off()

# Jaw
jaw$annot = jaw$seurat_clusters
Idents(jaw) = jaw$annot
jaw = RenameIdents(jaw, '0' = "Epithelial", '1' = "Epithelial", '2' = "Immune", '3' = "Epithelial", '4' = "Epithelial", "5" = "Pigmented", "6" = "Mesenchymal", "7" = "Epithelial", "8" = "Immune", "9" = "Immune", "10" = "Mature TB", "11" = "Immature TB")
svg("~/research/tooth/results/jaw_annot.svg")
pal = c("#2A9D8F", "#E9C46A", "#F6B179", "#E24D28", "#264653", "#264653")
print(DimPlot(jaw, order = T, label = F, pt.size = 1.5, cols = pal) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
dev.off()

# *** krt ***
p1 = print(FeaturePlot(tj, "krt15", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
p2 = print(FeaturePlot(tj, "krt5", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
p3 = print(FeaturePlot(jaw, "krt15", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
p4 = print(FeaturePlot(jaw, "krt5", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
p_list = list(p1, p2, p3, p4)

svg("~/research/tooth/results/tj_and_jaw_krt.svg")
print(plot_grid(plotlist=p_list, ncol = 2))
dev.off()


svg("~/research/tooth/results/jaw_krt15.svg")
print(FeaturePlot(jaw, "krt15", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
dev.off()
svg("~/research/tooth/results/jaw_krt5.svg")
print(FeaturePlot(jaw, "krt5", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
dev.off()
svg("~/research/tooth/results/tj_krt15.svg")
print(FeaturePlot(tj, "krt15", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
dev.off()
svg("~/research/tooth/results/tj_krt5.svg")
print(FeaturePlot(tj, "krt5", order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
dev.off()

# krt_df = data.frame(sample=c("RT", "RT", "Jaw", "Jaw"), gene=c("pos", "neg", "pos", "neg"), value = c(length(tj_krt15)/ncol(tj), (ncol(tj)- length(tj_krt15))/ncol(tj), length(jaw_krt15)/ncol(jaw), (ncol(jaw) - length(jaw_krt15))/ncol(jaw)))
# krt_df$value = krt_df$value * 100
# ggplot(krt_df, aes(sample, value, fill = gene)) + geom_bar(stat="identity") + xlab("") + ylab("% Cells Expressing krt15") + scale_fill_manual(values = c("lightgray", "blue")) + NoLegend() + theme_bw()
# 
# jaw_krt15 = colnames(jaw)[which(jaw@assays$RNA@counts["krt15",] > 0)]
# tj_krt15 = colnames(tj)[which(tj@assays$RNA@counts["krt15",] > 0)]
# jaw_krt5 = colnames(jaw)[which(jaw@assays$RNA@counts["krt5",] > 0)]
# tj_krt5 = colnames(tj)[which(tj@assays$RNA@counts["krt5",] > 0)]
# svg("~/research/tooth/results/jaw_krt15_all.svg")
# print(DimPlot(jaw, cells.highlight = list(), order = T, label = F) + theme_void() + ggtitle("") + coord_fixed() + NoLegend())
# dev.off()

# *************************************************************************************************************
# Figure 2 ====================================================================================================
# *************************************************************************************************************
incsr_unique_100 = readxl::read_excel("~/Downloads/incsr_unique_100.xlsx")
im_unique_100 = readxl::read_excel("~/Downloads/im_unique_100.xlsx")
hm_unique_100 = readxl::read_excel("~/Downloads/hm_unique_100.xlsx")
incsr_unique_100_mz = list()
other_query_genes = c()
for (query_clust in colnames(incsr_unique_100)) {
  print(paste("INCSR:", query_clust))
  this_degs = incsr_unique_100[, query_clust]
  this_degs_mz = convertHgncDataFrameToMzebra(data.frame(this_degs), gene_column = 1, gene_names = rownames(tj), na.rm = T, return_vect = T)
  this_degs_mz = this_degs_mz[which(! this_degs_mz %in% other_query_genes )]
  incsr_unique_100_mz[[query_clust]] = this_degs_mz
  other_query_genes = c(other_query_genes, this_degs_mz)
}
im_unique_100_mz = list()
other_query_genes = c()
for (query_clust in colnames(im_unique_100)) {
  print(paste("IM:", query_clust))
  this_degs = im_unique_100[, query_clust]
  this_degs_mz = convertHgncDataFrameToMzebra(data.frame(this_degs), gene_column = 1, gene_names = rownames(tj), na.rm = T, return_vect = T)
  this_degs_mz = this_degs_mz[which(! this_degs_mz %in% other_query_genes )]
  im_unique_100_mz[[query_clust]] = this_degs_mz
  other_query_genes = c(other_query_genes, this_degs_mz)
}
hm_unique_100_mz = list()
other_query_genes = c()
for (query_clust in colnames(hm_unique_100)) {
  print(paste("HM:", query_clust))
  this_degs = hm_unique_100[, query_clust]
  this_degs_mz = convertHgncDataFrameToMzebra(data.frame(this_degs), gene_column = 1, gene_names = rownames(tj), na.rm = T, return_vect = T)
  this_degs_mz = this_degs_mz[which(! this_degs_mz %in% other_query_genes )]
  hm_unique_100_mz[[query_clust]] = this_degs_mz
  other_query_genes = c(other_query_genes, this_degs_mz)
}

# Make a heatmap for every query cluster for every target cluster
tj_mat = tj@assays$RNA@counts
tj_mat[which(tj_mat > 1)] = 1
jaw_mat = jaw@assays$RNA@counts
jaw_mat[which(jaw_mat > 1)] = 1
tj_num_df = data.frame()
exp_df = data.frame()
other_query_genes = c()
annot_row = data.frame(query_dataset = c(), row.names = c())
tj = ScaleData(tj, features = unique(c(unlist(incsr_unique_100_mz), unlist(im_unique_100_mz), unlist(hm_unique_100_mz))))
for (query_clust in colnames(incsr_unique_100)[order(colnames(incsr_unique_100), decreasing = F)]) {
  this_degs_mz = incsr_unique_100_mz[[query_clust]]
  query_clust = paste0("(MI) ", query_clust)
  for (i in levels(tj$annot)[order(levels(tj$annot), decreasing = F)]) {
    this_sum = sum(rowSums(tj@assays$RNA@scale.data[this_degs_mz, which(tj$annot == i)]))
    this_pct = this_sum / length(this_degs_mz) / length(which(tj$annot == i))
    exp_df = rbind(exp_df, data.frame(query_dataset = "INCSR", query_clust = query_clust, tj_clust = i, value = this_sum, pct = this_pct, num_markers = length(this_degs_mz), num_cells = length(which(tj$annot == i))))
    annot_row = rbind(annot_row, data.frame(query_dataset = "INCSR", row.names = query_clust))
    # exp_df$pct[which(exp_df$query_clust == query_clust & exp_df$tj_clust == i)] = exp_df$value[which(exp_df$query_clust == query_clust & exp_df$tj_clust == i)] / length(which(tj$annot == i)) / length(this_degs_mz)
    
    num_sum = length(which( rowSums(tj_mat[this_degs_mz, which(tj$annot == i)]) > 0 ))
    tj_num_df = rbind(tj_num_df, data.frame(query_dataset = "INCSR", query_clust = query_clust, tj_clust = i, value = num_sum))
  }
}
for (query_clust in colnames(hm_unique_100)[order(colnames(hm_unique_100), decreasing = F)]) {
  this_degs_mz = hm_unique_100_mz[[query_clust]]
  query_clust = paste0("(HM) ", query_clust)
  for (i in levels(tj$annot)[order(levels(tj$annot), decreasing = F)]) {
    this_sum = sum(rowSums(tj@assays$RNA@scale.data[this_degs_mz, which(tj$annot == i)]))
    this_pct = this_sum / length(this_degs_mz) / length(which(tj$annot == i))
    exp_df = rbind(exp_df, data.frame(query_dataset = "HM", query_clust = query_clust, tj_clust = i, value = this_sum, pct = this_pct, num_markers = length(this_degs_mz), num_cells = length(which(tj$annot == i))))
    annot_row = rbind(annot_row, data.frame(query_dataset = "HM", row.names = query_clust))
    # exp_df$pct[which(exp_df$query_clust == query_clust & exp_df$tj_clust == i)] = exp_df$value[which(exp_df$query_clust == query_clust & exp_df$tj_clust == i)] / length(which(tj$annot == i)) / length(this_degs_mz)
    
    num_sum = length(which( rowSums(tj_mat[this_degs_mz, which(tj$annot == i)]) > 0 ))
    tj_num_df = rbind(tj_num_df, data.frame(query_dataset = "HM", query_clust = query_clust, tj_clust = i, value = num_sum))
  }
}
for (query_clust in colnames(im_unique_100)[order(colnames(im_unique_100), decreasing = F)]) {
  this_degs_mz = im_unique_100_mz[[query_clust]]
  query_clust = paste0("(MIM) ", query_clust)
  for (i in levels(tj$annot)[order(levels(tj$annot), decreasing = F)]) {
    this_sum = sum(rowSums(tj@assays$RNA@scale.data[this_degs_mz, which(tj$annot == i)]))
    this_pct = this_sum / length(this_degs_mz) / length(which(tj$annot == i))
    exp_df = rbind(exp_df, data.frame(query_dataset = "IM", query_clust = query_clust, tj_clust = i, value = this_sum, pct = this_pct, num_markers = length(this_degs_mz), num_cells = length(which(tj$annot == i))))
    annot_row = rbind(annot_row, data.frame(query_dataset = "IM", row.names = query_clust))
    # exp_df$pct[which(exp_df$query_clust == query_clust & exp_df$tj_clust == i)] = exp_df$value[which(exp_df$query_clust == query_clust & exp_df$tj_clust == i)] / length(which(tj$annot == i)) / length(this_degs_mz)
    
    num_sum = length(which( rowSums(tj_mat[this_degs_mz, which(tj$annot == i)]) > 0 ))
    tj_num_df = rbind(tj_num_df, data.frame(query_dataset = "IM", query_clust = query_clust, tj_clust = i, value = num_sum))
  }
}
tj_exp_mat = acast(exp_df, tj_clust ~ query_clust, value.var = "pct")
tj_exp_mat = t(tj_exp_mat)
# pheatmap::pheatmap(exp_mat, cellwidth = 25, cellheight = 25, annotation_names_row = F, show_rownames = T, show_colnames = T, angle_col = "0", cluster_cols = F, cluster_rows = F, scale = "column", colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100), filename = "~/research/tooth/results/tj_exp_scale_row.pdf")
# dev.off()

# Cichlid Jaw Heatmap
jaw_num_df = data.frame()
exp_df = data.frame()
other_query_genes = c()
annot_row = data.frame(query_dataset = c(), row.names = c())
jaw = ScaleData(jaw, features = unique(c(unlist(incsr_unique_100_mz), unlist(im_unique_100_mz), unlist(hm_unique_100_mz))))
for (query_clust in colnames(incsr_unique_100)[order(colnames(incsr_unique_100), decreasing = F)]) {
  this_degs_mz = incsr_unique_100_mz[[query_clust]]
  query_clust = paste0("(MI) ", query_clust)
  for (i in levels(jaw$annot)[order(levels(jaw$annot), decreasing = F)]) {
    this_sum = sum(rowSums(jaw@assays$RNA@scale.data[this_degs_mz, which(jaw$annot == i)]))
    this_pct = this_sum / length(this_degs_mz) / length(which(jaw$annot == i))
    exp_df = rbind(exp_df, data.frame(query_dataset = "INCSR", query_clust = query_clust, jaw_clust = i, value = this_sum, pct = this_pct, num_markers = length(this_degs_mz), num_cells = length(which(jaw$annot == i))))
    annot_row = rbind(annot_row, data.frame(query_dataset = "INCSR", row.names = query_clust))
    # exp_df$pct[which(exp_df$query_clust == query_clust & exp_df$jaw_clust == i)] = exp_df$value[which(exp_df$query_clust == query_clust & exp_df$jaw_clust == i)] / length(which(jaw$annot == i)) / length(this_degs_mz)
    
    num_sum = length(which( rowSums(jaw_mat[this_degs_mz, which(jaw$annot == i)]) > 0 ))
    jaw_num_df = rbind(jaw_num_df, data.frame(query_dataset = "INCSR", query_clust = query_clust, jaw_clust = i, value = num_sum))
  }
}
for (query_clust in colnames(hm_unique_100)[order(colnames(hm_unique_100), decreasing = F)]) {
  this_degs_mz = hm_unique_100_mz[[query_clust]]
  query_clust = paste0("(HM) ", query_clust)
  for (i in levels(jaw$annot)[order(levels(jaw$annot), decreasing = F)]) {
    this_sum = sum(rowSums(jaw@assays$RNA@scale.data[this_degs_mz, which(jaw$annot == i)]))
    this_pct = this_sum / length(this_degs_mz) / length(which(jaw$annot == i))
    exp_df = rbind(exp_df, data.frame(query_dataset = "HM", query_clust = query_clust, jaw_clust = i, value = this_sum, pct = this_pct, num_markers = length(this_degs_mz), num_cells = length(which(jaw$annot == i))))
    annot_row = rbind(annot_row, data.frame(query_dataset = "HM", row.names = query_clust))
    # exp_df$pct[which(exp_df$query_clust == query_clust & exp_df$jaw_clust == i)] = exp_df$value[which(exp_df$query_clust == query_clust & exp_df$jaw_clust == i)] / length(which(jaw$annot == i)) / length(this_degs_mz)
    
    num_sum = length(which( rowSums(jaw_mat[this_degs_mz, which(jaw$annot == i)]) > 0 ))
    jaw_num_df = rbind(jaw_num_df, data.frame(query_dataset = "HM", query_clust = query_clust, jaw_clust = i, value = num_sum))
  }
}
for (query_clust in colnames(im_unique_100)[order(colnames(im_unique_100), decreasing = F)]) {
  this_degs_mz = im_unique_100_mz[[query_clust]]
  query_clust = paste0("(MIM) ", query_clust)
  for (i in levels(jaw$annot)[order(levels(jaw$annot), decreasing = F)]) {
    this_sum = sum(rowSums(jaw@assays$RNA@scale.data[this_degs_mz, which(jaw$annot == i)]))
    this_pct = this_sum / length(this_degs_mz) / length(which(jaw$annot == i))
    exp_df = rbind(exp_df, data.frame(query_dataset = "IM", query_clust = query_clust, jaw_clust = i, value = this_sum, pct = this_pct, num_markers = length(this_degs_mz), num_cells = length(which(jaw$annot == i))))
    annot_row = rbind(annot_row, data.frame(query_dataset = "IM", row.names = query_clust))
    # exp_df$pct[which(exp_df$query_clust == query_clust & exp_df$jaw_clust == i)] = exp_df$value[which(exp_df$query_clust == query_clust & exp_df$jaw_clust == i)] / length(which(jaw$annot == i)) / length(this_degs_mz)
    
    num_sum = length(which( rowSums(jaw_mat[this_degs_mz, which(jaw$annot == i)]) > 0 ))
    jaw_num_df = rbind(jaw_num_df, data.frame(query_dataset = "IM", query_clust = query_clust, jaw_clust = i, value = num_sum))
  }
}
jaw_exp_mat = acast(exp_df, jaw_clust ~ query_clust, value.var = "pct")
jaw_exp_mat = t(jaw_exp_mat)
# pheatmap::pheatmap(exp_mat, cellwidth = 25, cellheight = 25, annotation_names_row = F, show_rownames = T, show_colnames = T, angle_col = "0", cluster_cols = F, cluster_rows = F, scale = "column", colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100), filename = "~/research/tooth/results/jaw_exp_scale_row.pdf")
jaw_exp_mat = jaw_exp_mat[,which(! colnames(jaw_exp_mat) %in% c("Mature TB", "Immature TB", "Pigmented"))]
colnames(jaw_exp_mat) = paste("(CJ)", colnames(jaw_exp_mat))
colnames(tj_exp_mat) = paste("(CT)", colnames(tj_exp_mat))

# Reorder clusters
jaw_rownames = rownames(jaw_exp_mat)[which(startsWith(rownames(jaw_exp_mat), "(HM)"))][order(rownames(jaw_exp_mat)[which(startsWith(rownames(jaw_exp_mat), "(HM)"))], decreasing = T)]
jaw_rownames = c(jaw_rownames, rownames(jaw_exp_mat)[which(startsWith(rownames(jaw_exp_mat), "(MIM)"))][order(rownames(jaw_exp_mat)[which(startsWith(rownames(jaw_exp_mat), "(MIM)"))], decreasing = T)])
jaw_rownames = c(jaw_rownames, rownames(jaw_exp_mat)[which(startsWith(rownames(jaw_exp_mat), "(MI)"))][order(rownames(jaw_exp_mat)[which(startsWith(rownames(jaw_exp_mat), "(MI)"))], decreasing = T)])
jaw_exp_mat = jaw_exp_mat[jaw_rownames, order(colnames(jaw_exp_mat), decreasing = T)]
tj_rownames = rownames(tj_exp_mat)[which(startsWith(rownames(tj_exp_mat), "(HM)"))][order(rownames(tj_exp_mat)[which(startsWith(rownames(tj_exp_mat), "(HM)"))], decreasing = T)]
tj_rownames = c(tj_rownames, rownames(tj_exp_mat)[which(startsWith(rownames(tj_exp_mat), "(MIM)"))][order(rownames(tj_exp_mat)[which(startsWith(rownames(tj_exp_mat), "(MIM)"))], decreasing = T)])
tj_rownames = c(tj_rownames, rownames(tj_exp_mat)[which(startsWith(rownames(tj_exp_mat), "(MI)"))][order(rownames(tj_exp_mat)[which(startsWith(rownames(tj_exp_mat), "(MI)"))], decreasing = T)])
tj_exp_mat = tj_exp_mat[tj_rownames, order(colnames(tj_exp_mat), decreasing = T)]

tj_jaw_exp_mat = cbind(tj_exp_mat, jaw_exp_mat)
pheatmap::pheatmap(tj_jaw_exp_mat, border_color = NA, cellwidth = 25, cellheight = 25, annotation_names_row = F, show_rownames = T, show_colnames = T, angle_col = "45", cluster_cols = F, cluster_rows = F, scale = "column", colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100), filename = "~/research/tooth/results/tj_jaw_annot_exp_scale_col.pdf")

p = pheatmap::pheatmap(t(tj_jaw_exp_mat), border_color = NA, annotation_names_row = F, show_rownames = T, show_colnames = T, angle_col = "45", cluster_cols = F, cluster_rows = F, scale = "row", colorRampPalette(viridis(100))(100))
p_df = melt(p$gtable$grobs[[1]]$children[[1]]$gp$fill)
p_df$Var1 = factor(p_df$Var1, levels = unique(p_df$Var1)[order(unique(p_df$Var1), decreasing = T)])
p2 = ggplot(p_df, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pdf("~/research/tooth/results/tj_jaw_annot_exp_scale_col_viridis.pdf", width = 10, height = 10)
print(p2)
dev.off()

# Genes that are present in cichlids that are top 100 unique DEGs in mammals
na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}

test = unlist(c(im_unique_100, hm_unique_100, incsr_unique_100))
test = test[which( !is.na(test) )]
ens_converter = convertHgncDataFrameToMzebra(data.frame(gene = test), 1, rownames(tj))

tj_hits = list()
test2 = hm_unique_100_mz[["PDL"]][which( rowSums(tj_mat[hm_unique_100_mz[["PDL"]], which(tj$annot == "Mesenchymal")]) > 0 )]
tj_hits[["Mesenchymal"]] = ens_converter$gene[match(test2, ens_converter$mz)]
test2 = hm_unique_100_mz[["Immune cells"]][which( rowSums(tj_mat[hm_unique_100_mz[["Immune cells"]], which(tj$annot == "Immune")]) > 0 )]
tj_hits[["Immune"]] = ens_converter$gene[match(test2, ens_converter$mz)]
test2 = hm_unique_100_mz[["Immune cells"]][which( rowSums(tj_mat[hm_unique_100_mz[["Immune cells"]], which(tj$annot == "Glia")]) > 0 )]
tj_hits[["Glia"]] = ens_converter$gene[match(test2, ens_converter$mz)]
test2 = im_unique_100_mz[["Epithelial cells"]][which( rowSums(tj_mat[im_unique_100_mz[["Epithelial cells"]], which(tj$annot == "Epithelial")]) > 0 )]
tj_hits[["Epithelial"]] = ens_converter$gene[match(test2, ens_converter$mz)]
test2 = im_unique_100_mz[["Endothelial cells"]][which( rowSums(tj_mat[im_unique_100_mz[["Endothelial cells"]], which(tj$annot == "Endothelial")]) > 0 )]
tj_hits[["Endothelial"]] = ens_converter$gene[match(test2, ens_converter$mz)]

tj_hits_hgnc = makePaddedDataFrame(tj_hits)
write.csv(tj_hits_hgnc, "~/research/tooth/results/tj_hgnc_mammal_ovlp.csv")

jaw_hits = list()
test2 = im_unique_100_mz[["Pulp cells"]][which( rowSums(jaw_mat[hm_unique_100_mz[["Pulp cells"]], which(jaw$annot == "Mesenchymal")]) > 0 )]
jaw_hits[["Mesenchymal"]] = ens_converter$gene[match(test2, ens_converter$mz)]
test2 = hm_unique_100_mz[["Immune cells"]][which( rowSums(jaw_mat[hm_unique_100_mz[["Immune cells"]], which(jaw$annot == "Immune")]) > 0 )]
jaw_hits[["Immune"]] = ens_converter$gene[match(test2, ens_converter$mz)]
test2 = im_unique_100_mz[["Epithelial cells"]][which( rowSums(jaw_mat[im_unique_100_mz[["Epithelial cells"]], which(jaw$annot == "Epithelial")]) > 0 )]
jaw_hits[["Epithelial"]] = ens_converter$gene[match(test2, ens_converter$mz)]

jaw_hits_hgnc = makePaddedDataFrame(jaw_hits)
write.csv(jaw_hits_hgnc, "~/research/tooth/results/jaw_hgnc_mammal_ovlp.csv")



# *************************************************************************************************************
# Figure 3 ====================================================================================================
# *************************************************************************************************************
pdf("~/research/tooth/results/tj_jaw_cyto.pdf", width = 5.5, height = 5)
FeaturePlot(tj_jaw, "cyto", order = T) + scale_color_gradientn(colors = pal(100)) + ggtitle("") + theme_void()
dev.off()
pdf("~/research/tooth/results/incsr_cyto.pdf", width = 5.5, height = 5)
FeaturePlot(igor_incsr, "cyto", order = T) + scale_color_gradientn(colors = pal(100)) + ggtitle("") + theme_void()
dev.off()
pdf("~/research/tooth/results/incsr_epi_cyto.pdf", width = 5.5, height = 5)
FeaturePlot(igor_incsr_epi, "cyto", order = T, pt.size = 1.5) + scale_color_gradientn(colors = pal(100)) + ggtitle("") + theme_void()
dev.off()
pdf("~/research/tooth/results/im_cyto.pdf", width = 5.5, height = 5)
FeaturePlot(im, "cyto", order = T) + scale_color_gradientn(colors = pal(100)) + ggtitle("") + theme_void()
dev.off()

pdf("~/research/tooth/results/tj_jaw_cyto_celsr1.pdf", width = 5.5, height = 5)
Idents(tj_jaw) = ""
myFeaturePlot(tj_jaw, "cyto2", na.blank = T, my.col.pal = pal) + ggtitle("") + theme_void()
dev.off()
tj$cyto2 = NA
tj$cyto2[which(tj@assays$RNA@counts['celsr1a',] > 0)] = tj$cyto[which(tj@assays$RNA@counts['celsr1a',] > 0)]
pdf("~/research/tooth/results/tj_cyto_celsr1.pdf", width = 5.5, height = 5)
Idents(tj) = ""
myFeaturePlot(tj, "cyto2", na.blank = T, my.col.pal = pal, my.pt.size = 1.5) + ggtitle("") + theme_void()
dev.off()
jaw$cyto2 = NA
jaw$cyto2[which(jaw@assays$RNA@counts['celsr1a',] > 0)] = jaw$cyto[which(jaw@assays$RNA@counts['celsr1a',] > 0)]
pdf("~/research/tooth/results/jaw_cyto_celsr1.pdf", width = 5.5, height = 5)
Idents(jaw) = ""
myFeaturePlot(jaw, "cyto2", na.blank = T, my.col.pal = pal, my.pt.size = 1.5) + ggtitle("") + theme_void()
dev.off()
pdf("~/research/tooth/results/incsr_cyto_celsr1.pdf", width = 5.5, height = 5)
Idents(igor_incsr) = ""
igor_incsr$cyto2 = NA
igor_incsr$cyto2[which(igor_incsr@assays$RNA@counts['Celsr1',] > 0)] = igor_incsr$cyto[which(igor_incsr@assays$RNA@counts['Celsr1',] > 0)]
myFeaturePlot(igor_incsr, "cyto2", na.blank = T, my.col.pal = pal) + ggtitle("") + theme_void()
dev.off()
pdf("~/research/tooth/results/incsr_epi_cyto_celsr1.pdf", width = 5.5, height = 5)
Idents(igor_incsr_epi) = ""
igor_incsr_epi$cyto2 = NA
igor_incsr_epi$cyto2[which(igor_incsr_epi@assays$RNA@counts['Celsr1',] > 0)] = igor_incsr_epi$cyto[which(igor_incsr_epi@assays$RNA@counts['Celsr1',] > 0)]
myFeaturePlot(igor_incsr_epi, "cyto2", na.blank = T, my.col.pal = pal, my.pt.size = 1.5) + ggtitle("") + theme_void()
dev.off()
pdf("~/research/tooth/results/im_cyto_celsr1.pdf", width = 5.5, height = 5)
Idents(im) = ""
im$cyto2 = NA
im$cyto2[which(im@assays$RNA@counts['CELSR1',] > 0)] = im$cyto[which(im@assays$RNA@counts['CELSR1',] > 0)]
myFeaturePlot(im, "cyto2", na.blank = T, my.col.pal = pal) + ggtitle("") + theme_void()
dev.off()

igor_incsr = readRDS("~/research/tooth/data/igor_incsr.rds")
im = readRDS("~/research/tooth/data/igor_incsr_molar.rds")
igor_incsr_epi = readRDS("~/research/tooth/data/igor_incsr_epi.rds")
igor_incsr$cyto = CytoTRACE::CytoTRACE(as.matrix(igor_incsr@assays$RNA@counts))
im$cyto = CytoTRACE::CytoTRACE(as.matrix(im@assays$RNA@counts))

# Calculate CytoTRACE for Celsr1- cells
incsr_no_celsr1 = subset(igor_incsr, cells = colnames(igor_incsr)[which(igor_incsr@assays$RNA@counts["Celsr1",] == 0)])
incsr_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(incsr_no_celsr1@assays$RNA@counts))
im_no_celsr1 = subset(im, cells = colnames(im)[which(im@assays$RNA@counts["CELSR1",] == 0)])
im_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(im_no_celsr1@assays$RNA@counts))
tj_no_celsr1 = subset(tj, cells = colnames(tj)[which(tj@assays$RNA@counts["celsr1a",] == 0)])
tj_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(tj_no_celsr1@assays$RNA@counts))
jaw_no_celsr1 = subset(jaw, cells = colnames(jaw)[which(jaw@assays$RNA@counts["celsr1a",] == 0)])
jaw_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(jaw_no_celsr1@assays$RNA@counts))
# tj_jaw = merge(tj, jaw)

# Calculate CytoTRACE for Celsr1+ cells
incsr_celsr1 = subset(igor_incsr, cells = colnames(igor_incsr)[which(igor_incsr@assays$RNA@counts["Celsr1",] > 0)])
incsr_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(incsr_celsr1@assays$RNA@counts))
im_celsr1 = subset(im, cells = colnames(im)[which(im@assays$RNA@counts["CELSR1",] > 0)])
im_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(im_celsr1@assays$RNA@counts))
tj_celsr1 = subset(tj, cells = colnames(tj)[which(tj@assays$RNA@counts["celsr1a",] > 0)])
tj_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(tj_celsr1@assays$RNA@counts))
jaw_celsr1 = subset(jaw, cells = colnames(jaw)[which(jaw@assays$RNA@counts["celsr1a",] > 0)])
jaw_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(jaw_celsr1@assays$RNA@counts))
# tj_jaw_celsr1 = subset(tj_jaw, cells = colnames(tj_jaw)[which(tj_jaw@assays$RNA@counts["celsr1a",] > 0)])
# tj_jaw_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(tj_jaw_celsr1@assays$RNA@counts))

# Calculate Number of cells
incsr_mat = incsr_no_celsr1@assays$RNA@counts
incsr_mat[which(incsr_mat > 1)] = 1
incsr_celsr1_mat = incsr_celsr1@assays$RNA@counts
incsr_celsr1_mat[which(incsr_celsr1_mat > 1)] = 1

im_mat = im_no_celsr1@assays$RNA@counts
im_mat[which(im_mat > 1)] = 1
im_celsr1_mat = im_celsr1@assays$RNA@counts
im_celsr1_mat[which(im_celsr1_mat > 1)] = 1

cg_df = data.frame(gene = names(incsr_celsr1_cyto$cytoGenes)[which( names(incsr_celsr1_cyto$cytoGenes) %in% names(incsr_no_celsr1_cyto$cytoGenes) )])
cg_df$all = incsr_no_celsr1_cyto$cytoGenes[match(cg_df$gene, names(incsr_no_celsr1_cyto$cytoGenes))]
cg_df$celsr1 = incsr_celsr1_cyto$cytoGenes[match(cg_df$gene, names(incsr_celsr1_cyto$cytoGenes))]
cg_df$dif = cg_df$celsr1 - cg_df$all
cg_df$abs_dif = abs(cg_df$dif)
cg_df$p = r_to_p(cg_df$all, cg_df$celsr1, ncol(incsr_mat), ncol(incsr_celsr1_mat))
cg_df$q = p.adjust(cg_df$p, method = "BH")
cg_df$all_num = rowSums(incsr_mat)[match(cg_df$gene, rownames(incsr_mat))]
cg_df$celsr1_num = rowSums(incsr_celsr1_mat)[match(cg_df$gene, rownames(incsr_celsr1_mat))]
# ggplot(cg_df, aes(x = all, y = celsr1, color = abs(all-celsr1))) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( abs(cg_df$all - cg_df$celsr1) > 0.7 ),], aes(label = gene), color = "black") + ggtitle("Igor Mouse Incisor")
# ggplot(cg_df, aes(x = celsr1_num, y = abs_dif, color = abs_dif)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$abs_dif > 0.6 ),], aes(label = gene), color = "black") + ggtitle("Igor Mouse Incisor")
ggplot(cg_df, aes(x = celsr1_num, y = dif, color = abs_dif)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$dif > 0.55 ),], aes(label = gene), color = "black") + ggtitle("Igor Mouse Incisor")
cg_df_incsr = cg_df

cg_df = data.frame(gene = names(im_celsr1_cyto$cytoGenes)[which( names(im_celsr1_cyto$cytoGenes) %in% names(im_no_celsr1_cyto$cytoGenes) )])
cg_df$all = im_no_celsr1_cyto$cytoGenes[match(cg_df$gene, names(im_no_celsr1_cyto$cytoGenes))]
cg_df$celsr1 = im_celsr1_cyto$cytoGenes[match(cg_df$gene, names(im_celsr1_cyto$cytoGenes))]
cg_df$dif = cg_df$celsr1 - cg_df$all
cg_df$abs_dif = abs(cg_df$dif)
cg_df$p = r_to_p(cg_df$all, cg_df$celsr1, ncol(im_mat), ncol(im_celsr1_mat))
cg_df$q = p.adjust(cg_df$p, method = "BH")
cg_df$all_num = rowSums(im_mat)[match(cg_df$gene, rownames(im_mat))]
cg_df$celsr1_num = rowSums(im_celsr1_mat)[match(cg_df$gene, rownames(im_celsr1_mat))]
# ggplot(cg_df, aes(x = all, y = celsr1, color = abs(all-celsr1))) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( abs(cg_df$all - cg_df$celsr1) > 0.7 ),], aes(label = gene), color = "black") + ggtitle("Igor Mouse Incisor + Molar")
# ggplot(cg_df, aes(x = celsr1_num, y = abs_dif, color = abs_dif)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$abs_dif > 0.6 ),], aes(label = gene), color = "black") + ggtitle("Igor Mouse Incisor + Molar")
ggplot(cg_df, aes(x = celsr1_num, y = dif, color = abs_dif)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$dif > 0.50 ),], aes(label = gene), color = "black") + ggtitle("Igor Mouse Incisor + Molar")
cg_df$incsr_dif = cg_df_incsr$dif[match(str_to_title(cg_df$gene), cg_df_incsr$gene)]
cg_df$tmp =  cg_df$dif + cg_df$incsr_dif
ggplot(cg_df, aes(x = dif, y = incsr_dif, color = tmp)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$tmp > 0.6 ),], aes(label = gene), color = "black") + xlab("Mouse Incisor & Molar Celsr1+ minus Celsr1- R") + ylab("Mouse Incisor Celsr1+ minus Celsr1- R")
cg_df_im = cg_df

common_celsr1_genes = rownames(tj_celsr1_cyto$exprMatrix)[which(rownames(tj_celsr1_cyto$exprMatrix) %in% rownames(jaw_celsr1_cyto$exprMatrix))]
common_celsr1_exprMatrix = cbind(tj_celsr1_cyto$exprMatrix[common_celsr1_genes,], jaw_celsr1_cyto$exprMatrix[common_celsr1_genes,])
common_celsr1_cyto = c(tj_celsr1_cyto$CytoTRACE, jaw_celsr1_cyto$CytoTRACE)
c_celsr1_cytoGenes = unlist(lapply(common_celsr1_genes, function(x) cor(common_celsr1_exprMatrix[x,], common_celsr1_cyto) ))
names(c_celsr1_cytoGenes) = common_celsr1_genes

common_all_genes = rownames(tj_no_celsr1_cyto$exprMatrix)[which(rownames(tj_no_celsr1_cyto$exprMatrix) %in% rownames(jaw_no_celsr1_cyto$exprMatrix))]
common_all_exprMatrix = cbind(tj_no_celsr1_cyto$exprMatrix[common_all_genes,], jaw_no_celsr1_cyto$exprMatrix[common_all_genes,])
common_all_cyto = c(tj_no_celsr1_cyto$CytoTRACE, jaw_no_celsr1_cyto$CytoTRACE)
c_all_cytoGenes = unlist(lapply(common_all_genes, function(x) cor(common_all_exprMatrix[x,], common_all_cyto) ))
names(c_all_cytoGenes) = common_all_genes

common_common_genes = common_celsr1_genes[which(common_celsr1_genes %in% common_all_genes)]
cg_df = data.frame(gene = common_common_genes, all = c_all_cytoGenes[match(common_common_genes, names(c_all_cytoGenes))], celsr1 = c_celsr1_cytoGenes[match(common_common_genes, names(c_celsr1_cytoGenes))])
cg_df$dif = cg_df$celsr1 - cg_df$all
cg_df$abs_dif = abs(cg_df$dif)
cg_df$p = r_to_p(cg_df$all, cg_df$celsr1, length(common_celsr1_cyto), length(common_all_cyto))
cg_df$q = p.adjust(cg_df$p, method = "BH")
ggplot(cg_df, aes(x = all, y = celsr1, color = abs_dif)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$abs_dif > 0.6 ),], aes(label = gene), color = "black") + ggtitle("Cichlid Tooth and Jaw")
cg_df$hgnc = ens_gene_info$hgnc[match(cg_df$gene, ens_gene_info$gene_raw)]
cg_df$incsr_dif = cg_df_incsr$dif[match(cg_df$hgnc, toupper(cg_df_incsr$gene))]
cg_df$im_dif = cg_df_im$dif[match(cg_df$hgnc, cg_df_im$gene)]
cg_df$tmp =  cg_df$dif + cg_df$incsr_dif
cg_df$tmp_im =  cg_df$dif + cg_df$im_dif
ggplot(cg_df, aes(x = dif, y = incsr_dif, color = tmp)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(data = cg_df[which( cg_df$tmp > 0.6 ),], aes(label = hgnc), color = "black") + xlab("Cichlid Celsr1+ minus Celsr1- R") + ylab("Mouse Incisor Celsr1+ minus Celsr1- R")

cg_df$incsr_q = cg_df_incsr$q[match(cg_df$hgnc, toupper(cg_df_incsr$gene))]
cg_df$im_q = cg_df_im$q[match(cg_df$hgnc, toupper(cg_df_im$gene))]
cdh3_df = data.frame(dataset = c("Cichlid", "MI", "MIM"), gene = "CDH3", q = c(cg_df$q[which(cg_df$hgnc == "CDH3")], cg_df$incsr_q[which(cg_df$hgnc == "CDH3")], cg_df$im_q[which(cg_df$hgnc == "CDH3")]))
cdh3_df$neg_log_10 = -log10(cdh3_df$q)
ggplot(cdh3_df, aes(x = gene, y =neg_log_10, fill = dataset)) + geom_bar(stat = "identity", position = position_dodge2()) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + xlab("") + ylab("-Log10 Adjusted P")

up_celsr_genes = cg_df$hgnc[which(cg_df$q < 0.05 & cg_df$incsr_q < 0.05 & cg_df$im_q < 0.05 & cg_df$dif > 0 & cg_df$incsr_dif > 0 & cg_df$im_dif > 0)]
up_celsr_df = data.frame(dataset = c(rep("Cichlid", length(up_celsr_genes)), rep("MI", length(up_celsr_genes)), rep("MIM", length(up_celsr_genes))),
                         gene = rep(up_celsr_genes, 3),
                         isPos = "Celsr1+",
                         data = c(cg_df$celsr1[match(up_celsr_genes, cg_df$hgnc)], cg_df_incsr$celsr1[match(str_to_title(up_celsr_genes), cg_df_incsr$gene)], cg_df_im$celsr1[match(up_celsr_genes, cg_df_im$gene)]))
up_celsr_df = rbind(up_celsr_df, data.frame(dataset = c(rep("Cichlid", length(up_celsr_genes)), rep("MI", length(up_celsr_genes)), rep("MIM", length(up_celsr_genes))),
                                            gene = rep(up_celsr_genes, 3),
                                            isPos = "Celsr1-",
                                            data = c(cg_df$all[match(up_celsr_genes, cg_df$hgnc)], cg_df_incsr$all[match(str_to_title(up_celsr_genes), cg_df_incsr$gene)], cg_df_im$all[match(up_celsr_genes, cg_df_im$gene)])))
up_celsr_df$dataset_gene = paste0(up_celsr_df$dataset, "_", up_celsr_df$gene)
up_celsr_df %>% ggplot(aes(x = gene, y = data, color = dataset, group = dataset_gene, shape = isPos)) + geom_point(stat = "identity", position = position_dodge(width = 0.5), size = 2.5) + geom_line(position = position_dodge(width = 0.5), arrow = arrow(length=unit(0.30,"cm"), ends = "first", angle = 20)) + xlab("") + ylab("Correlation w/ CytoTRACE") + ggtitle("Correlation of Genes w/ CytoTRACE in Celsr1+ Cells (Top Point) vs Celsr1- Cells (Bottom Point)") + theme(plot.subtitle = element_text("Celsr1+ Cells = Top Point, Celsr1- Cells = Bottom Point")) + scale_shape_manual(values = c(15,19)) + theme(legend.title=element_blank())
up_celsr_df %>% ggplot(aes(x = dataset, y = data, color = dataset, group = dataset_gene, shape = isPos)) + geom_point(stat = "identity", position = position_dodge(width = 0.5), size = 2.5) + geom_line(position = position_dodge(width = 0.5), arrow = arrow(length=unit(0.30,"cm"), ends = "first", angle = 20)) + xlab("") + ylab("Correlation w/ CytoTRACE") + ggtitle("Correlation of Genes w/ CytoTRACE in Celsr1+ Cells (Top Point) vs Celsr1- Cells (Bottom Point)") + theme(plot.subtitle = element_text("Celsr1+ Cells = Top Point, Celsr1- Cells = Bottom Point")) + scale_shape_manual(values = c(15,19)) + theme(legend.title=element_blank()) + facet_wrap(~gene, ncol = length(unique(up_celsr_df$gene))) + theme_bw()
# ggplot(cg_df2, aes(x = dif, y = incsr_dif, color = tmp)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(aes(label = hgnc), color = "black") + xlab("Cichlid Celsr1+ minus Celsr1- R") + ylab("Mouse Incisor Celsr1+ minus Celsr1- R")
# ggplot(cg_df2, aes(x = dif, y = im_dif, color = tmp_im)) + geom_point() + scale_color_gradientn(colors = plasma(50)) + geom_text_repel(aes(label = hgnc), color = "black") + xlab("Cichlid Celsr1+ minus Celsr1- R") + ylab("Mouse Incisor & Molar Celsr1+ minus Celsr1- R")
# + geom_smooth(method = "lm", se  = F, color = "#ad544e") + ggpubr::stat_cor(method = "pearson", hjust=0, vjust=0, label.x = -0.1, label.y = -0.08, color = "#ad544e") + coord_cartesian(clip = 'off')


# CDH3 - Method 1
p_gene =  "Cdh3"
mz_gene = ens_gene_info$gene_raw[match(toupper(p_gene), ens_gene_info$hgnc)]
gene_df = data.frame(dataset = c(rep("Cichlid", length(common_celsr1_cyto)), rep("(MI)", length(incsr_celsr1_cyto$CytoTRACE)), rep("(MIM)", length(im_celsr1_cyto$CytoTRACE))),
                     expr = c(common_celsr1_exprMatrix[mz_gene,], incsr_celsr1_cyto$exprMatrix[p_gene,], im_celsr1_cyto$exprMatrix[toupper(p_gene),]),
                     cyto = c(common_celsr1_cyto, incsr_celsr1_cyto$CytoTRACE, im_celsr1_cyto$CytoTRACE),
                     cor = c(rep(cor(common_celsr1_exprMatrix[mz_gene,], common_celsr1_cyto), length(common_celsr1_cyto)), rep(cor(incsr_celsr1_cyto$exprMatrix[p_gene,], incsr_celsr1_cyto$CytoTRACE), length(incsr_celsr1_cyto$CytoTRACE)), rep(cor(im_celsr1_cyto$exprMatrix[toupper(p_gene),], im_celsr1_cyto$CytoTRACE), length(im_celsr1_cyto$CytoTRACE))),
                     gene = p_gene)
ggplot(gene_df, aes(x = cyto, y = expr, color = dataset)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm", se  = F) + ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", hjust=0, vjust=0, label.x = 1, label.y = c(0.6, 0.3, 0), size = 5, geom = "label") + coord_cartesian(clip = 'off') + xlim(expand = c(0,1)) + ggtitle(paste0(p_gene, " in Celsr1+ Cells"))
gene_df = data.frame(dataset = c(rep("Cichlid", length(common_all_cyto)), rep("(MI)", length(incsr_no_celsr1_cyto$CytoTRACE)), rep("(MIM)", length(im_no_celsr1_cyto$CytoTRACE))),
                     expr = c(common_all_exprMatrix[mz_gene,], incsr_no_celsr1_cyto$exprMatrix[p_gene,], im_no_celsr1_cyto$exprMatrix[toupper(p_gene),]),
                     cyto = c(common_all_cyto, incsr_no_celsr1_cyto$CytoTRACE, im_no_celsr1_cyto$CytoTRACE),
                     cor = c(rep(cor(common_all_exprMatrix[mz_gene,], common_all_cyto), length(common_all_cyto)), rep(cor(incsr_no_celsr1_cyto$exprMatrix[p_gene,], incsr_no_celsr1_cyto$CytoTRACE), length(incsr_no_celsr1_cyto$CytoTRACE)), rep(cor(im_no_celsr1_cyto$exprMatrix[toupper(p_gene),], im_no_celsr1_cyto$CytoTRACE), length(im_no_celsr1_cyto$CytoTRACE))),
                     gene = p_gene)
ggplot(gene_df, aes(x = cyto, y = expr, color = dataset)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm", se  = F) + ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", hjust=0, vjust=0, label.x = 1, label.y = c(0.6, 0.3, 0), size = 5, geom = "label") + coord_cartesian(clip = 'off') + xlim(expand = c(0,1)) + ggtitle(paste0(p_gene, " in Celsr1- Cells")) + scale_color_manual(values = c("#d19a97", "#65a679", "#6d87b3"))

# CDH3 - Method 2
gene_df = data.frame(dataset = c(rep("Cichlid", length(common_celsr1_cyto)+length(common_all_cyto))),
                     isPos = c(rep("Celsr1+", length(common_celsr1_cyto)), rep("Celsr1-", length(common_all_cyto))),
                     expr = c(common_celsr1_exprMatrix[mz_gene,], common_all_exprMatrix[mz_gene,]),
                     cyto = c(common_celsr1_cyto, common_all_cyto),
                     cor = c(rep(cor(common_celsr1_exprMatrix[mz_gene,], common_celsr1_cyto), length(common_celsr1_cyto)), rep(cor(common_all_exprMatrix[mz_gene,], common_all_cyto), length(common_all_cyto))),
                     gene = p_gene)
ggplot(gene_df, aes(x = cyto, y = expr, color = isPos)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm", se  = F) + ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", hjust=0, vjust=0, label.x = 1, label.y = c(0.6, 0.3, 0), size = 5, geom = "label") + coord_cartesian(clip = 'off') + xlim(expand = c(0,1)) + ggtitle(paste0(p_gene, " in Celsr1+ vs Celsr1- Cells in Cichlid")) + scale_color_manual(values = c(hue_pal()(3)[1], "darkgray")) + theme(legend.title=element_blank()) 
gene_df = data.frame(dataset = c(rep("(MI)", length(incsr_celsr1_cyto$CytoTRACE)+length(incsr_no_celsr1_cyto$CytoTRACE))),
                     isPos = c(rep("Celsr1+", length(incsr_celsr1_cyto$CytoTRACE)), rep("Celsr1-", length(incsr_no_celsr1_cyto$CytoTRACE))),
                     expr = c(incsr_celsr1_cyto$exprMatrix[p_gene,], incsr_no_celsr1_cyto$exprMatrix[p_gene,]),
                     cyto = c(incsr_celsr1_cyto$CytoTRACE, incsr_no_celsr1_cyto$CytoTRACE),
                     cor = c(rep(cor(incsr_celsr1_cyto$CytoTRACE, incsr_celsr1_cyto$CytoTRACE), length(incsr_celsr1_cyto$CytoTRACE)), rep(cor(incsr_no_celsr1_cyto$exprMatrix[p_gene,], incsr_no_celsr1_cyto$CytoTRACE), length(incsr_no_celsr1_cyto$CytoTRACE))),
                     gene = p_gene)
ggplot(gene_df, aes(x = cyto, y = expr, color = isPos)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm", se  = F) + ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", hjust=0, vjust=0, label.x = 1, label.y = c(0.6, 0.3, 0), size = 5, geom = "label") + coord_cartesian(clip = 'off') + xlim(expand = c(0,1)) + ggtitle(paste0(p_gene, " in Celsr1+ vs Celsr1- Cells in (MI)")) + scale_color_manual(values = c(hue_pal()(3)[2], "darkgray")) + theme(legend.title=element_blank()) 
gene_df = data.frame(dataset = c(rep("(MIM)", length(im_celsr1_cyto$CytoTRACE)+length(im_no_celsr1_cyto$CytoTRACE))),
                     isPos = c(rep("Celsr1+", length(im_celsr1_cyto$CytoTRACE)), rep("Celsr1-", length(im_no_celsr1_cyto$CytoTRACE))),
                     expr = c(im_celsr1_cyto$exprMatrix[p_gene,], im_no_celsr1_cyto$exprMatrix[p_gene,]),
                     cyto = c(im_celsr1_cyto$CytoTRACE, im_no_celsr1_cyto$CytoTRACE),
                     cor = c(rep(cor(im_celsr1_cyto$CytoTRACE, im_celsr1_cyto$CytoTRACE), length(im_celsr1_cyto$CytoTRACE)), rep(cor(im_no_celsr1_cyto$exprMatrix[p_gene,], im_no_celsr1_cyto$CytoTRACE), length(im_no_celsr1_cyto$CytoTRACE))),
                     gene = p_gene)
ggplot(gene_df, aes(x = cyto, y = expr, color = isPos)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm", se  = F) + ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", hjust=0, vjust=0, label.x = 1, label.y = c(0.6, 0.3, 0), size = 5, geom = "label") + coord_cartesian(clip = 'off') + xlim(expand = c(0,1)) + ggtitle(paste0(p_gene, " in Celsr1+ vs Celsr1- Cells in (MIM)")) + scale_color_manual(values = c(hue_pal()(3)[3], "darkgray")) + theme(legend.title=element_blank()) 


# Tooth and Jaw Merging
tj = readRDS("~/research/tooth/data/tj.rds")
jaw = readRDS("~/research/tooth/data/jpool.rds")
tj_jaw = merge(tj, jaw)
tj_jaw_celsr1_deg = mergeToothJawCytoBIN2(tj, jaw, "celsr1a")
igor_incsr_celsr1_deg = createCytoBINsInGene(igor_incsr, "Celsr1")
im_celsr1_deg = createCytoBINsInGene(im, "CELSR1")
igor_incsr_epi_celsr1_deg = createCytoBINsInGene(igor_incsr_epi, "Gli1")

hm_celsr1_deg_low = hm_celsr1_deg[which(hm_celsr1_deg$cluster == "Low"),]
hm_celsr1_deg_medium = hm_celsr1_deg[which(hm_celsr1_deg$cluster == "Medium"),]
hm_celsr1_deg_high = hm_celsr1_deg[which(hm_celsr1_deg$cluster == "High"),]
im_celsr1_deg_low = im_celsr1_deg[which(im_celsr1_deg$cluster == "Low"),]
im_celsr1_deg_medium = im_celsr1_deg[which(im_celsr1_deg$cluster == "Medium"),]
im_celsr1_deg_high = im_celsr1_deg[which(im_celsr1_deg$cluster == "High"),]
igor_incsr_celsr1_deg_low = igor_incsr_celsr1_deg[which(igor_incsr_celsr1_deg$cluster == "Low"),]
igor_incsr_celsr1_deg_medium = igor_incsr_celsr1_deg[which(igor_incsr_celsr1_deg$cluster == "Medium"),]
igor_incsr_celsr1_deg_high = igor_incsr_celsr1_deg[which(igor_incsr_celsr1_deg$cluster == "High"),]

common_deg = degWithHgncAndDescription(tj_jaw_celsr1_deg, rownames(tj))
common_deg_low = common_deg[which(common_deg$cluster == "Low"),]
common_deg_low$im_cluster = im_celsr1_deg_low$cluster_sign[match(common_deg_low$hgnc, im_celsr1_deg_low$gene)]
common_deg_low$incsr_cluster = igor_incsr_celsr1_deg_low$cluster_sign[match(common_deg_low$hgnc, toupper(igor_incsr_celsr1_deg_low$gene))]
common_deg_low = common_deg_low[which(common_deg_low$cluster_sign == common_deg_low$im_cluster & common_deg_low$cluster_sign == common_deg_low$incsr_cluster),]
common_deg_medium = common_deg[which(common_deg$cluster == "Medium"),]
common_deg_medium$im_cluster = im_celsr1_deg_medium$cluster_sign[match(common_deg_medium$hgnc, im_celsr1_deg_medium$gene)]
common_deg_medium$incsr_cluster = igor_incsr_celsr1_deg_medium$cluster_sign[match(common_deg_medium$hgnc, toupper(igor_incsr_celsr1_deg_medium$gene))]
common_deg_medium = common_deg_medium[which(common_deg_medium$cluster_sign == common_deg_medium$im_cluster & common_deg_medium$cluster_sign == common_deg_medium$incsr_cluster),]
common_deg_high = common_deg[which(common_deg$cluster == "High"),]
common_deg_high$im_cluster = im_celsr1_deg_high$cluster_sign[match(common_deg_high$hgnc, im_celsr1_deg_high$gene)]
common_deg_high$incsr_cluster = igor_incsr_celsr1_deg_high$cluster_sign[match(common_deg_high$hgnc, toupper(igor_incsr_celsr1_deg_high$gene))]
common_deg_high = common_deg_high[which(common_deg_high$cluster_sign == common_deg_high$im_cluster & common_deg_high$cluster_sign == common_deg_high$incsr_cluster),]

common_deg = common_deg[which(common_deg$cluster_sign == common_deg$im_cluster & common_deg$cluster_sign == common_deg$incsr_cluster),]

common_deg_backup = common_deg
common_deg = common_deg[which(common_deg$cluster_sign == common_deg$im_cluster & common_deg$cluster_sign == common_deg$incsr_cluster),]

common_deg = read.csv("~/research/tooth/results/celsr_cytobin_common_deg.csv")
mz_high_degs = common_deg$gene[which(!common_deg$gene %in% c("jag1b", "ENSMZEG00005001584") )]
mouse_high_degs = str_to_title( unique(common_deg$hgnc) )
human_high_degs = unique(common_deg$hgnc)

cytobin_df = data.frame(bin = tj_jaw_celsr1$bin, exp = colSums(tj_jaw_celsr1@assays$RNA@data[mz_high_degs,]), cyto = tj_jaw_celsr1$cyto, dataset = "Cichlid")
cytobin_df = rbind(cytobin_df, data.frame(bin = igor_incsr_celsr1$bin, exp = colSums(igor_incsr_celsr1@assays$RNA@data[mouse_high_degs,]), cyto = igor_incsr_celsr1$cyto, dataset = "MI"))
# cytobin_df = rbind(cytobin_df, data.frame(bin = igor_incsr_epi_celsr1$bin, exp = igor_incsr_epi_celsr1@assays$RNA@data["Ashl",], cyto = igor_incsr_epi_celsr1$cyto, dataset = "MIE"))
cytobin_df = rbind(cytobin_df, data.frame(bin = im_celsr1$bin, exp = colSums(im_celsr1@assays$RNA@data[human_high_degs,]), cyto = im_celsr1$cyto, dataset = "MIM"))
cytobin_df$bin = plyr::revalue(cytobin_df$bin, replace = c("relative_high" = "High", "relative_medium" = "Medium", "relative_low" = "Low"))
cytobin_df$bin = factor(cytobin_df$bin, levels = c("High", "Medium", "Low"))

cdh3_df = data.frame(dataset = c("Cichlid", "MI", "MIM"),
                     gene = "Cdh3",
                     p_val_adj = c(tj_jaw_celsr1_deg$p_val_adj[which(tj_jaw_celsr1_deg$gene == "ENSMZEG00005018585" & tj_jaw_celsr1_deg$avg_log2FC > 0)], igor_incsr_celsr1_deg$p_val_adj[which(igor_incsr_celsr1_deg$gene == "Cdh3" & igor_incsr_celsr1_deg$avg_log2FC > 0)], im_celsr1_deg$p_val_adj[which(im_celsr1_deg$gene == "CDH3" & im_celsr1_deg$avg_log2FC > 0)]))
cdh3_df$neg_log_10 = -log10(cdh3_df$p_val_adj)
ggplot(cdh3_df, aes(x = gene, y =neg_log_10, fill = dataset)) + geom_bar(stat = "identity", position = position_dodge2()) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + xlab("") + ylab("-Log10 Adjusted P")

pdf("~/research/tooth/results/celsr1_cyto.pdf", width = 8, height = 5)
ggplot(cytobin_df, aes(x=dataset, y = cyto, fill = bin, color = bin)) + geom_boxplot(alpha = 0.6) + geom_jitter(position=position_dodge2(width=0.75), alpha=0.5) + scale_color_manual(values = c(temp[10], "darkgoldenrod1", temp[2])) + scale_fill_manual(values = c(temp[10], "darkgoldenrod1", temp[2])) + xlab("") + ylab("CytoTRACE") + theme(panel.grid.major.x = element_blank()) + theme_bw()
dev.off()

pdf("~/research/tooth/results/celsr1_cytobin_deg.pdf", width = 8, height = 5)
ggplot(cytobin_df, aes(x=dataset, y = exp, fill = bin, color = bin)) + geom_boxplot(alpha = 0.6) + geom_jitter(position=position_dodge2(width=0.75), alpha=0.5) + scale_color_manual(values = c(temp[10], "darkgoldenrod1", temp[2])) + scale_fill_manual(values = c(temp[10], "darkgoldenrod1", temp[2])) + xlab("") + ylab("Normalized Expression") + ggtitle("") + theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major.x = element_blank())
dev.off()

igor_incsr_epi_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/igor_incsr_epi_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
incsr_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/incsr_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
im_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/im_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
# tj_jaw_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/tj_jaw_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")
tj_jaw_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/cichlid_celsr1_cytoBIN_deg_pos_w_toppgene.xlsx", sheet = "ToppGene")
igor_incsr_epi_high_toppgene = openxlsx::read.xlsx("~/research/tooth/results/igor_incsr_epi_celsr1_cytoBIN_deg_pos.xlsx", sheet = "high toppgene")

incsr_high_toppgene = incsr_high_toppgene[which(incsr_high_toppgene$`q-value.Bonferroni` < 0.05),]
im_high_toppgene = im_high_toppgene[which(im_high_toppgene$`q-value.Bonferroni` < 0.05),]
tj_jaw_high_toppgene = tj_jaw_high_toppgene[which(tj_jaw_high_toppgene$`q-value.Bonferroni` < 0.05),]

tj_jaw_high_toppgene$inIncsr = tj_jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & tj_jaw_high_toppgene$Name %in% incsr_high_toppgene$Name
tj_jaw_high_toppgene$inIM = tj_jaw_high_toppgene$Category %in% c("GO: Molecular Function", "GO: Biological Process", "GO: Cellular Component") & tj_jaw_high_toppgene$Name %in% im_high_toppgene$Name
common_high = tj_jaw_high_toppgene$ID[which( tj_jaw_high_toppgene$inIM & tj_jaw_high_toppgene$inIncsr )]

# toppgene_df = data.frame(category = tj_jaw_high_toppgene$Category[match(common_high, tj_jaw_high_toppgene$ID)], id = common_high, name = tj_jaw_high_toppgene$Name[match(common_high, tj_jaw_high_toppgene$ID)], cichlid_q = tj_jaw_high_toppgene$`q-value.FDR.B&H`[match(common_high, tj_jaw_high_toppgene$ID)], cichlid_p = tj_jaw_high_toppgene$`p-value`[match(common_high, tj_jaw_high_toppgene$ID)], cichlid_neg_log10_p = - log10(tj_jaw_high_toppgene$`p-value`[match(common_high, tj_jaw_high_toppgene$ID)]),
#                          mi_q = incsr_high_toppgene$`q-value.FDR.B&H`[match(common_high, incsr_high_toppgene$ID)], mi_p = incsr_high_toppgene$`p-value`[match(common_high, incsr_high_toppgene$ID)], mi_neg_log10_p = - log10(incsr_high_toppgene$`p-value`[match(common_high, incsr_high_toppgene$ID)]),
#                          mim_q = im_high_toppgene$`q-value.FDR.B&H`[match(common_high, im_high_toppgene$ID)], mim_p = im_high_toppgene$`p-value`[match(common_high, im_high_toppgene$ID)], mim_neg_log10_p = -log10(im_high_toppgene$`p-value`[match(common_high, im_high_toppgene$ID)]))
# toppgene_df = data.frame(category = tj_jaw_high_toppgene$Category[match(common_high, tj_jaw_high_toppgene$ID)], id = common_high, name = tj_jaw_high_toppgene$Name[match(common_high, tj_jaw_high_toppgene$ID)], cichlid_q = tj_jaw_high_toppgene$`q-value.Bonferroni`[match(common_high, tj_jaw_high_toppgene$ID)], cichlid_p = tj_jaw_high_toppgene$`p-value`[match(common_high, tj_jaw_high_toppgene$ID)], cichlid_neg_log10_p = - log10(tj_jaw_high_toppgene$`p-value`[match(common_high, tj_jaw_high_toppgene$ID)]),
#                          mi_q = incsr_high_toppgene$`q-value.Bonferroni`[match(common_high, incsr_high_toppgene$ID)], mi_p = incsr_high_toppgene$`p-value`[match(common_high, incsr_high_toppgene$ID)], mi_neg_log10_p = - log10(incsr_high_toppgene$`p-value`[match(common_high, incsr_high_toppgene$ID)]),
#                          mim_q = im_high_toppgene$`q-value.Bonferroni`[match(common_high, im_high_toppgene$ID)], mim_p = im_high_toppgene$`p-value`[match(common_high, im_high_toppgene$ID)], mim_neg_log10_p = -log10(im_high_toppgene$`p-value`[match(common_high, im_high_toppgene$ID)]))
# toppgene_df = toppgene_df[which(toppgene_df$cichlid_q < 0.05 & toppgene_df$mi_q < 0.05 & toppgene_df$mim_q < 0.05),]
# toppgene_df = toppgene_df[order(toppgene_df$cichlid_p),]
# toppgene_df = toppgene_df[order(toppgene_df$cichlid_p)[1:100],]
toppgene_df = data.frame(dataset = "Cichlid", id = common_high, name = tj_jaw_high_toppgene$Name[match(common_high, tj_jaw_high_toppgene$ID)], q = tj_jaw_high_toppgene$`q-value.FDR.B&H`[match(common_high, tj_jaw_high_toppgene$ID)], p = tj_jaw_high_toppgene$`p-value`[match(common_high, tj_jaw_high_toppgene$ID)])
toppgene_df = rbind(toppgene_df, data.frame(dataset = "MI", id = common_high, name = tj_jaw_high_toppgene$Name[match(common_high, tj_jaw_high_toppgene$ID)], q = incsr_high_toppgene$`q-value.FDR.B&H`[match(common_high, incsr_high_toppgene$ID)], p = incsr_high_toppgene$`p-value`[match(common_high, incsr_high_toppgene$ID)]))
# toppgene_df = rbind(toppgene_df, data.frame(dataset = "MIE", id = common_high, name = tj_jaw_high_toppgene$Name[match(common_high, tj_jaw_high_toppgene$ID)], q = igor_incsr_epi_high_toppgene$`q-value.FDR.B&H`[match(common_high, igor_incsr_epi_high_toppgene$ID)], p = igor_incsr_epi_high_toppgene$`p-value`[match(common_high, igor_incsr_epi_high_toppgene$ID)]))
toppgene_df = rbind(toppgene_df, data.frame(dataset = "MIM", id = common_high, name = tj_jaw_high_toppgene$Name[match(common_high, tj_jaw_high_toppgene$ID)], q = im_high_toppgene$`q-value.FDR.B&H`[match(common_high, im_high_toppgene$ID)], p = im_high_toppgene$`p-value`[match(common_high, im_high_toppgene$ID)]))
toppgene_df$neg_log10_p = - log10(toppgene_df$p)

toppgene_df$name2 = paste0( sapply(toppgene_df$name, function(x) substr(x, 1, 35)) , "...")
# pal = c("#2A9D8F", "#E9C46A", "#F6B179", "#E24D28", "#264653", "#264653")
pal = c("#2A9D8F", "#E9C46A", "#F6B179", "#E24D28", "#264653")

pdf("~/research/tooth/results/enriched_categories.pdf", width = 12, height = 8)
ggplot(toppgene_df, aes(x=name, y=neg_log10_p, fill = dataset, color = dataset)) + geom_bar(alpha = 0.85, stat = "identity", width = 0.5) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.grid.major.x = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 0.25, colour = "black")) + scale_y_continuous(expand = c(0,0), name = expression(-Log["10"]*" P")) + xlab("") + scale_fill_manual(values = c("darkgoldenrod1", "#20bf55ff", "#01baefff", "#9594edff")) + scale_color_manual(values = c("darkgoldenrod1", "#20bf55ff", "#01baefff", "#9594edff"))
dev.off()


# pdf("~/research/tooth/results/enriched_categories.pdf", width = 10, height = 10)
# go_sums = sapply(1:nrow(toppgene_df), function(x) toppgene_df$cichlid_neg_log10_p[x] + toppgene_df$mi_neg_log10_p[x] + toppgene_df$mim_neg_log10_p[x])
# toppgene_df = toppgene_df[order(go_sums), ]
# circos.par("gap.after" = c(rep(0, 11), 30), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.02, 0.02))
# circos.initialize(sectors = toppgene_df$name, xlim = c(0,1))
# # circos.track(ylim = c(0, max(toppgene_df$cichlid_neg_log10_p)), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
# circos.track(ylim = c(0, ceiling(max(go_sums))), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
#   pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
#   
#   circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index], col = "#00A08A90", border = "#00A08A")
#   circos.rect(CELL_META$cell.xlim[1], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index], CELL_META$cell.xlim[2], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mi_neg_log10_p[CELL_META$sector.numeric.index], col = "#FF000090", border = "#FF0000")
#   circos.rect(CELL_META$cell.xlim[1], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mi_neg_log10_p[CELL_META$sector.numeric.index], CELL_META$cell.xlim[2], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mi_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mim_neg_log10_p[CELL_META$sector.numeric.index], col = "#F2AD0090",  border = "#F2AD00")
#   
#   circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(1),
#               toppgene_df$name[CELL_META$sector.numeric.index], facing = "clockwise", niceFacing = TRUE,
#               adj = c(1, 0.5), cex = 0.6)
#   # circos.segments(CELL_META$cell.xlim[1], max(go_sums), CELL_META$cell.xlim[2], max(go_sums), col = "gray20", lty = 1) # top
#   circos.segments(CELL_META$cell.xlim[1], ceiling(max(go_sums)/2), CELL_META$cell.xlim[2], ceiling(max(go_sums)/2), col = "gray80", lty = 2) # middle
#   circos.segments(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 0, col = "gray20", lty = 1) # bottom
#   if (CELL_META$sector.numeric.index == 1) {
#     circos.segments(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[1], max(go_sums), col = "gray20", lty = 1)
#   } else if  (CELL_META$sector.numeric.index == nrow(toppgene_df)) {
#     # circos.segments(CELL_META$cell.xlim[1], max(go_sums), CELL_META$cell.xlim[2]+CELL_META$xcenter, max(go_sums), col = "gray20", lty = 1)
#     circos.segments(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2]+CELL_META$xcenter, 0, col = "gray20", lty = 1)
#     circos.segments(CELL_META$cell.xlim[1], ceiling(max(go_sums)/2), CELL_META$cell.xlim[2]+CELL_META$xcenter, ceiling(max(go_sums)/2), col = "gray80", lty = 2) # middle
#     # circos.segments(CELL_META$cell.xlim[2]+CELL_META$xcenter, 0, CELL_META$cell.xlim[2]+CELL_META$xcenter, max(go_sums), col = "gray20", lty = 1)
#   }
# })
# circos.yaxis(at = c(0, ceiling(max(go_sums)/2), ceiling(max(go_sums))), labels = F, sector.index = toppgene_df$name[1], track.index = 1, side = "left", labels.niceFacing = F)
# dev.off()

# Cor
tj_cytoCor = cytoCor(tj, "celsr1a")
incsr_cytoCor = cytoCor(incsr, "Celsr1")
im_cytoCor = cytoCor(im, "CELSR1")
ens_gene_info = read.csv("~/research/all_research/ens_hgnc_description.csv")
celsr1_df = data.frame(tj = tj_cytoCor)
celsr1_df$hgnc = ens_gene_info$hgnc[match(rownames(celsr1_df), ens_gene_info$gene_raw)]
celsr1_df = celsr1_df[which( !is.na(celsr1_df$hgnc) ),]
celsr1_df$im = im_cytoCor[celsr1_df$hgnc]
celsr1_df$incsr = incsr_cytoCor[str_to_title(celsr1_df$hgnc)]
celsr1_df = celsr1_df[which(! is.na(celsr1_df$im) & ! is.na(celsr1_df$incsr) ),]
celsr1_df$tj_rank[sort(celsr1_df$tj, decreasing = T, index.return = T)$ix] = 1:nrow(celsr1_df)
celsr1_df$im_rank[sort(celsr1_df$im, decreasing = T, index.return = T)$ix] = 1:nrow(celsr1_df)
celsr1_df$incsr_rank[sort(celsr1_df$incsr, decreasing = T, index.return = T)$ix] = 1:nrow(celsr1_df)
celsr1_df$tj_bot_rank[sort(celsr1_df$tj, decreasing = F, index.return = T)$ix] = 1:nrow(celsr1_df)
celsr1_df$im_bot_rank[sort(celsr1_df$im, decreasing = F, index.return = T)$ix] = 1:nrow(celsr1_df)
celsr1_df$incsr_bot_rank[sort(celsr1_df$incsr, decreasing = F, index.return = T)$ix] = 1:nrow(celsr1_df)

celsr1_df_top = celsr1_df[which( abs(celsr1_df$tj_rank - celsr1_df$im_rank) < 150 & abs(celsr1_df$tj_rank - celsr1_df$incsr_rank) < 150),]
temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
pal = colorRampPalette(temp)

# INCSR
time_df_p = data.frame()
gene = "ptprf"
# gene_cyto_df = data.frame(data = incsr@assays$RNA@data[str_to_title(gene),which(incsr@assays$RNA@data[str_to_title(gene),] > 0)], cyto = incsr$cyto, dataset = "MI", gene = "Celsr1")
gene_cyto_df = data.frame(data = incsr@assays$RNA@data[str_to_title(gene),which(incsr@assays$RNA@data[str_to_title("Celsr1"),] > 0)], cyto = incsr$cyto[which(incsr@assays$RNA@data[str_to_title("Celsr1"),] > 0)], dataset = "MI", gene = "Celsr1")
loessMod75 <- loess(gene_cyto_df[,"data"] ~ gene_cyto_df[,"cyto"], span=0.75)
smoothed75 <- predict(loessMod75)
# this_df = gene_cyto_df[,c("cyto", "dataset", "gene")]
# this_df$smooth = smoothed75
# time_df_p = rbind(time_df_p, this_df)
gene_cyto_df$smooth = smoothed75
time_df_p = rbind(time_df_p, gene_cyto_df)

# IM
# gene_cyto_df = data.frame(data = im@assays$RNA@data[toupper(gene),], cyto = im$cyto, dataset = "MIM", gene = "Celsr1")
gene_cyto_df = data.frame(data = im@assays$RNA@data[toupper(gene),which(im@assays$RNA@data[toupper("CELSR1"),] > 0)], cyto = im$cyto[which(im@assays$RNA@data[toupper("CELSR1"),] > 0)], dataset = "MIM", gene = "Celsr1")
loessMod75 <- loess(gene_cyto_df[,"data"] ~ gene_cyto_df[,"cyto"], span=0.75)
smoothed75 <- predict(loessMod75)
# this_df = gene_cyto_df[,c("cyto", "dataset", "gene")]
# this_df$smooth = smoothed75
# time_df_p = rbind(time_df_p, this_df)
gene_cyto_df$smooth = smoothed75
time_df_p = rbind(time_df_p, gene_cyto_df)

#TJ
mz_gene = paste0(gene, "b")
# gene_cyto_df = data.frame(data = tj@assays$RNA@data[paste0(gene, ".1"),], cyto = tj$cyto, dataset = "CT", gene = "Celsr1")
gene_cyto_df = data.frame(data = tj@assays$RNA@data[mz_gene, which(tj@assays$RNA@data["celsr1a",] > 0)], cyto = tj$cyto[which(tj@assays$RNA@data["celsr1a",] > 0)], dataset = "CT", gene = "Celsr1")
loessMod75 <- loess(gene_cyto_df[,"data"] ~ gene_cyto_df[,"cyto"], span=0.75)
smoothed75 <- predict(loessMod75)
# this_df = gene_cyto_df[,c("cyto", "dataset", "gene")]
# this_df$smooth = smoothed75
# time_df_p = rbind(time_df_p, this_df)
gene_cyto_df$smooth = smoothed75
time_df_p = rbind(time_df_p, gene_cyto_df)

time_df_p = time_df_p[order(time_df_p$cyto),]
time_df_p$data[which(time_df_p$data == 0)] = NA
ggplot(time_df_p, aes(color = dataset)) + geom_line(data = time_df_p, aes(cyto, smooth), size = 1.5) + geom_point(data = time_df_p, aes(cyto, data), alpha = 0.1) + ggtitle(gene) + theme_bw()

gene_cyto_df = data.frame(data = incsr@assays$RNA@data["Celsr1",], cyto = incsr$cyto, dataset = "MI", gene = "Celsr1")
gene_cyto_df = rbind(gene_cyto_df, data.frame(data = im@assays$RNA@data["CELSR1",], cyto = im$cyto, dataset = "MIM", gene = "Celsr1"))
gene_cyto_df = rbind(gene_cyto_df, data.frame(data = tj@assays$RNA@data["celsr1a",], cyto = tj$cyto, dataset = "CT", gene = "Celsr1"))
ggplot(gene_cyto_df, aes(cyto, data)) + geom_point() + geom_smooth(formula = y ~ x, method='loess') + theme_bw()

#***************************************************************************************************************************
# All Combined =============================================================================================================
#***************************************************************************************************************************
# Merge Cichlid Tooth, Cichlid Jaw, Mouse Incisor, Mouse Incisor + Molar
tj_hgnc = readRDS("~/research/tooth/data/tj_hgnc.rds")
jaw_hgnc = readRDS("~/research/tooth/data/jaw_hgnc.rds")
incsr_hgnc = readRDS("~/research/tooth/data/incsr_hgnc.rds")
im = readRDS("~/research/tooth/data/igor_incsr_molar.rds")
allt = merge(tj_hgnc, jaw_hgnc, incsr_hgnc, im)
allt = merge(tj_hgnc, list(jaw_hgnc, incsr_hgnc, im), add.cell.ids = c("ct", "cj", "mi", 'mim'))
saveRDS(allt, "~/research/tooth/data/allt.rds")

# Read Merged Data
allt = readRDS("~/research/tooth/data/allt.rds")
# allt_cyto = CytoTRACE::CytoTRACE(as.matrix(allt@assays$RNA@counts))
# allt$cyto = allt_cyto$CytoTRACE

# Find DEGs b/w High Bins
res = createCytoBINsInGene(allt, "CELSR1")
allt_celsr1_ctyo = res[[1]]
allt_celsr1 = res[[2]]
Idents(allt_celsr1) = allt_celsr1$bin
celsr1_deg = FindAllMarkers(allt_celsr1, only.pos = F, logfc.threshold = 0)
celsr1_deg$isSig = celsr1_deg$p_val_adj < 0.05
celsr1_deg$cluster_isSig = paste0(celsr1_deg$cluster, "_", celsr1_deg$isSig)
celsr1_deg_sig = celsr1_deg[which(celsr1_deg$p_val_adj < 0.05),]
celsr1_deg_sig_pos = celsr1_deg_sig[which(celsr1_deg_sig$avg_log2FC > 0),]

# Volcano Plot
celsr1_deg$col = plyr::revalue(celsr1_deg$cluster, replace = c("Low" = rev(brewer.pal(11,"Spectral"))[2], "Medium" = "gold", "High" = rev(brewer.pal(11,"Spectral"))[10]))
celsr1_deg$col[which(celsr1_deg$p_val_adj >= 0.05)] = "gray"
celsr1_deg$abs_pct_dif = abs(celsr1_deg$pct.1 - celsr1_deg$pct.2)
celsr1_deg$col[which( (celsr1_deg$abs_pct_dif < 0.05 | abs(celsr1_deg$avg_log2FC) < 0.25) )] = "gray"
pdf("~/scratch/d_tooth/results/igor/allt/cm_celsr1_sig_thresh.pdf", width = 6, height = 6)
ggplot(celsr1_deg, aes(x = avg_log2FC, y = -log10(p_val), color = col)) + geom_point(alpha = 0.3) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_y_sqrt() + theme_light() + scale_color_identity()
dev.off()

# Expression of CytoBIN DEGs vs CytoTRACE
library(ggpmisc)
celsr1_dif_genes = unique(celsr1_deg_sig$gene)
cytobin_df = data.frame(bin = allt_celsr1$bin, exp = colSums(allt_celsr1@assays$RNA@data[celsr1_deg_sig$gene[which(celsr1_deg_sig$cluster == "Low" & celsr1_deg_sig$avg_log2FC > 0)],]), bin_deg = "Low", cyto = allt_celsr1$cyto, dataset = "All")
cytobin_df = rbind(cytobin_df, data.frame(bin = allt_celsr1$bin, exp = colSums(allt_celsr1@assays$RNA@data[celsr1_deg_sig$gene[which(celsr1_deg_sig$cluster == "Medium" & celsr1_deg_sig$avg_log2FC > 0)],]), bin_deg = "Medium", cyto = allt_celsr1$cyto, dataset = "All"))
cytobin_df = rbind(cytobin_df, data.frame(bin = allt_celsr1$bin, exp = colSums(allt_celsr1@assays$RNA@data[celsr1_deg_sig$gene[which(celsr1_deg_sig$cluster == "High" & celsr1_deg_sig$avg_log2FC > 0)],]), bin_deg = "High", cyto = allt_celsr1$cyto, dataset = "All"))
cytobin_df = rbind(cytobin_df, data.frame(bin = allt_celsr1$bin, exp = colSums(allt_celsr1@assays$RNA@data[ unique(celsr1_deg_sig$gene[which(celsr1_deg_sig$cluster == "High" & celsr1_deg_sig$avg_log2FC > 0)], celsr1_deg_sig$gene[which(celsr1_deg_sig$cluster == "Low" & celsr1_deg_sig$avg_log2FC < 0)]),]), bin_deg = "Test", cyto = allt_celsr1$cyto, dataset = "All"))
cytobin_df$bin = factor(cytobin_df$bin, levels = c("Low", "Medium", "High"))
cytobin_df$bin_deg = factor(cytobin_df$bin_deg, levels = c("Low", "Medium", "High", "Test"))
png("allt_celsr1_boxplot.png", width = 900, height = 600, res = 100)
ggplot(cytobin_df, aes(x=bin_deg, y = exp, fill = bin, color = bin)) + geom_boxplot(alpha = 0.6) + geom_point(stat = "identity", position = position_jitterdodge(), alpha = 0.2) + scale_color_manual(values = c(temp[2], temp[6], temp[10]), guide = F) + scale_fill_manual(values = c(temp[2], temp[6], temp[10])) + xlab("CytoBIN DEGs") + ylab("Normalized Expression") + ggtitle("") + theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major.x = element_blank())
dev.off()
png("allt_celsr1_trendline_test.png", width = 900, height = 600, res = 100)
ggplot(cytobin_df[which(cytobin_df$bin_deg == "Test"),], aes(x=cyto, y = exp, color = cyto)) + geom_point(alpha = 0.6) + scale_color_gradientn(colors = pal(50), guide = F) + theme_light() + geom_smooth(method = "lm", color = "gray40") + stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label.., sep = "~~~")), parse = T) + ylab("Normalized Expression") + xlab("CytoTRACE")
dev.off()

# All in Celsr1-
allt_no_celsr1 = subset(allt, cells = colnames(allt)[which(! colnames(allt) %in% colnames(allt_celsr1))])
allt_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(allt_no_celsr1@assays$RNA@counts))
allt_no_celsr1$cyto = allt_no_celsr1_cyto$CytoTRACE

# Correlations
# allt_no_celsr1_cor = newCytoCor(allt_no_celsr1)
# write.csv(allt_no_celsr1_cor, "~/research/tooth/results/allt_no_celsr1_cor.csv")
# allt_celsr1_cor = newCytoCor(allt_celsr1)
# write.csv(allt_celsr1_cor, "~/research/tooth/results/allt_celsr1_cor.csv")
allt_no_celsr1_cor = read.csv("~/research/tooth/results/allt_no_celsr1_cor.csv")[,2]
names(allt_no_celsr1_cor) = read.csv("~/research/tooth/results/allt_no_celsr1_cor.csv")[,1]
allt_celsr1_cor = read.csv("~/research/tooth/results/allt_celsr1_cor.csv")[,2]
names(allt_celsr1_cor) = read.csv("~/research/tooth/results/allt_celsr1_cor.csv")[,1]

# Correlation in Celsr1+ vs Celsr1-
common_allt_gene = names(allt_celsr1_cor)[which(names(allt_celsr1_cor) %in% names(allt_no_celsr1_cor))]
allt_cor_df = data.frame(gene = common_allt_gene, celsr1 = allt_celsr1_cor[match(common_allt_gene, names(allt_celsr1_cor))], no_celsr1 = allt_no_celsr1_cor[match(common_allt_gene, names(allt_no_celsr1_cor))])
allt_cor_df$dif = allt_cor_df$celsr1 - allt_cor_df$no_celsr1
ggplot(allt_cor_df, aes(x = no_celsr1, y = celsr1, color = abs(dif))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + xlab("Correlation of Gene w/ CytoTRACE in Celsr1- Cells") + ylab("Correlation of Gene w/ CytoTRACE in Celsr1+ Cells") + geom_text_repel(data = allt_cor_df[which(allt_cor_df$dif > 0.55),], aes(label = gene), color = plasma(100)[10])

mat_allt_celsr1 = allt_celsr1@assays$RNA@counts
mat_allt_celsr1[which(mat_allt_celsr1 > 1)] = 1
allt_celsr1_cells = rowSums(mat_allt_celsr1)
mat_allt_no_celsr1 = allt_no_celsr1@assays$RNA@counts
mat_allt_no_celsr1[which(mat_allt_no_celsr1 > 1)] = 1
allt_no_celsr1_cells = rowSums(mat_allt_no_celsr1)
allt_cor_df$num_in_celsr1 = ncol(allt_celsr1)
allt_cor_df$num_in_no_celsr1 = ncol(allt_no_celsr1)
# allt_cor_df$num_in_celsr1 = allt_celsr1_cells[match(allt_cor_df$gene, names(allt_celsr1_cells))]
# allt_cor_df$num_in_no_celsr1 = allt_no_celsr1_cells[match(allt_cor_df$gene, names(allt_no_celsr1_cells))]
allt_cor_df$p = r_to_p(allt_cor_df$celsr1, allt_cor_df$no_celsr1, allt_cor_df$num_in_celsr1, allt_cor_df$num_in_no_celsr1)
allt_cor_df$bon = p.adjust(allt_cor_df$p, method = "bonferroni")
allt_cor_df$neg_log_bon = -log10(allt_cor_df$bon)
allt_cor_df$isSig = allt_cor_df$bon < 0.05
pal = colorRamp(plasma(10))
allt_cor_df$col = rgb(pal( allt_cor_df$neg_log_bon/max(allt_cor_df$neg_log_bon) ), maxColorValue = 255)
allt_cor_df$col[which(allt_cor_df$neg_log_bon < 30 | abs(allt_cor_df$dif) < 0.3)] = "darkgray"
png("allt_yc_vs_nc_cor_p.png", width = 800, height = 600, res = 100)
# ggplot(allt_cor_df, aes(x = dif, y = neg_log_bon, color = neg_log_bon)) + geom_point() + geom_hline(yintercept = -log10(0.05), linetype="dashed") + xlab("R Celsr1+ - R Celsr1-") + ylab(expression(-Log["10"]*" Adjusted P")) + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + geom_text_repel(data = allt_cor_df[which(allt_cor_df$neg_log_bon > 60),], aes(label = gene))
ggplot(allt_cor_df, aes(x = dif, y = neg_log_bon, color = col)) + geom_point() + geom_hline(yintercept = 30, linetype="dashed") + geom_vline(xintercept = c(-0.3, 0.3), linetype="dashed") + xlab("R Celsr1+ - R Celsr1-") + ylab(expression(-Log["10"]*" Adjusted P")) + scale_color_identity() + theme_bw() + geom_text_repel(data = allt_cor_df[which(allt_cor_df$neg_log_bon > 60),], aes(label = gene))
dev.off()


#
allt_no_celsr1$bin <- allt_no_celsr1$cyto
allt_no_celsr1$bin[which(allt_no_celsr1$cyto <= quantile(allt_no_celsr1$cyto, 0.33))] <- "Low"
allt_no_celsr1$bin[which(allt_no_celsr1$cyto > quantile(allt_no_celsr1$cyto, 0.33) & allt_no_celsr1$cyto <= quantile(allt_no_celsr1$cyto, 0.66))] <- "Medium"
allt_no_celsr1$bin[which(allt_no_celsr1$cyto > quantile(allt_no_celsr1$cyto, 0.66))] <- "High"
Idents(allt_no_celsr1) = allt_no_celsr1$bin
no_celsr1_deg = FindAllMarkers(allt_no_celsr1, only.pos = F, logfc.threshold = 0)
no_celsr1_deg$isSig = no_celsr1_deg$p_val_adj < 0.05
no_celsr1_deg$cluster_isSig = paste0(no_celsr1_deg$cluster, "_", no_celsr1_deg$isSig)
no_celsr1_deg_sig = no_celsr1_deg[which(no_celsr1_deg$p_val_adj < 0.05),]
no_celsr1_deg_sig_pos = no_celsr1_deg_sig[which(no_celsr1_deg_sig$avg_log2FC > 0),]

# Compare allt Celsr1+ vs Celsr1- CytoBIN DEGs
celsr1_deg = read.csv("allt_celsr1_deg.csv")
rownames(celsr1_deg) = celsr1_deg$X
celsr1_deg$X = NULL
celsr1_deg_sig = celsr1_deg[which(celsr1_deg$p_val_adj < 0.05),]
no_celsr1_deg = read.csv("allt_no_celsr1_deg.csv")
rownames(no_celsr1_deg) = no_celsr1_deg$X
no_celsr1_deg$X = NULL
no_celsr1_deg_sig = no_celsr1_deg[which(no_celsr1_deg$p_val_adj < 0.05),]

celsr1_deg_sig_low    = celsr1_deg_sig[which(celsr1_deg_sig$cluster == "Low"),]
celsr1_deg_sig_medium = celsr1_deg_sig[which(celsr1_deg_sig$cluster == "Medium"),]
celsr1_deg_sig_high   = celsr1_deg_sig[which(celsr1_deg_sig$cluster == "High"),]
no_celsr1_deg_sig_low    = no_celsr1_deg_sig[which(no_celsr1_deg_sig$cluster == "Low"),]
no_celsr1_deg_sig_medium = no_celsr1_deg_sig[which(no_celsr1_deg_sig$cluster == "Medium"),]
no_celsr1_deg_sig_high   = no_celsr1_deg_sig[which(no_celsr1_deg_sig$cluster == "High"),]

celsr1_y_n_not_common = celsr1_deg_sig_low[which(! celsr1_deg_sig_low$gene %in% no_celsr1_deg_sig_low$gene ),]
celsr1_y_n_not_common = rbind(celsr1_y_n_not_common, celsr1_deg_sig_medium[which(! celsr1_deg_sig_medium$gene %in% no_celsr1_deg_sig_medium$gene ),])
celsr1_y_n_not_common = rbind(celsr1_y_n_not_common, celsr1_deg_sig_high[which(! celsr1_deg_sig_high$gene %in% no_celsr1_deg_sig_high$gene ),])
write.csv(celsr1_y_n_not_common, "~/scratch/d_tooth/results/igor/allt/allt_celsr1_deg_sig_and_not_a_no_celsr1_deg.csv")

num_celsr1_deg_sig = data.frame(aggregate(gene ~ cluster, data = celsr1_deg_sig, length))
colnames(num_celsr1_deg_sig)[2] = "All Celsr1+ DEG"
num_celsr1_deg_sig$not_common = aggregate(gene ~ cluster, data = celsr1_y_n_not_common, length)[,2]
colnames(num_celsr1_deg_sig)[3] = "Celsr1+ DEG and Not Celsr1-"
num_celsr1_deg_sig = reshape2::melt(num_celsr1_deg_sig)
png("~/scratch/d_tooth/results/igor/allt/allt_celsr1_deg_sig_vs_no_celsr1_big.png", width = 800, height = 600, res = 100)
print(ggplot(num_celsr1_deg_sig, aes(x=cluster, y = value, color = variable, fill = variable)) + geom_bar(stat = "identity", alpha = 0.8, position=position_dodge2()) + scale_color_manual(values = c("#2dcc84", "gold")) + scale_fill_manual(values = c("#2dcc84", "gold")) + xlab("") + ylab("Number of DEGs"))
dev.off()

allt$celsr1_bin = "None"
allt$celsr1_bin[colnames(allt_celsr1)] = paste0("y_", allt_celsr1$bin)
allt$celsr1_bin[colnames(allt_no_celsr1)] = paste0("n_", allt_no_celsr1$bin)
Idents(allt) = allt$celsr1_bin
hvh_deg = FindMarkers(allt, ident.1="y_High", ident.2="n_High")
mvm_deg = FindMarkers(allt, ident.1="y_Medium", ident.2="n_Medium")
lvl_deg = FindMarkers(allt, ident.1="y_Low", ident.2="n_Low")

hvh_deg$cluster = "High"
mvm_deg$cluster = "Medium"
lvl_deg$cluster = "Low"
hvh_deg$gene = rownames(hvh_deg)
mvm_deg$gene = rownames(mvm_deg)
lvl_deg$gene = rownames(lvl_deg)
bin_deg_dif = rbind(hvh_deg, mvm_deg, lvl_deg)
bin_deg_dif_sig = bin_deg_dif[which(bin_deg_dif$p_val_adj < 0.05),]
write.csv(bin_deg_dif, "~/scratch/d_tooth/results/igor/allt/allt_celsr1_vs_no_celsr1_deg.csv")
p_df = as.data.frame(table(bin_deg_dif_sig$cluster))
png("~/scratch/d_tooth/results/igor/allt/allt_y_vs_n_celsr1_deg_sig.png", width = 800, height = 600, res = 100)
print(ggplot(p_df, aes(x=Var1, y = Freq)) + geom_bar(stat = "identity", alpha = 0.8, position=position_dodge2()) + xlab("") + ylab("Number of DEGs"))
dev.off()

# Human CytoBIN DEG
hm = readRDS("~/scratch/d_tooth/data/hm.rds")
res = createCytoBINsInGene(hm, "CELSR1")
hm_celsr1_ctyo = res[[1]]
hm_celsr1 = res[[2]]
Idents(hm_celsr1) = hm_celsr1$bin
hm_celsr1_deg = FindAllMarkers(hm_celsr1, only.pos = F, logfc.threshold = 0)
hm_celsr1_deg$isSig = hm_celsr1_deg$p_val_adj < 0.05
hm_celsr1_deg$cluster_isSig = paste0(hm_celsr1_deg$cluster, "_", hm_celsr1_deg$isSig)
hm_celsr1_deg_sig = hm_celsr1_deg[which(hm_celsr1_deg$isSig),]

# Find Human CytoBIN DEGs that are also Cichlid+Mouse CytoBIN DEGs
hm_celsr1_deg_sig_low    = hm_celsr1_deg_sig[which(hm_celsr1_deg_sig$cluster == "Low"),]
hm_celsr1_deg_sig_medium = hm_celsr1_deg_sig[which(hm_celsr1_deg_sig$cluster == "Medium"),]
hm_celsr1_deg_sig_high   = hm_celsr1_deg_sig[which(hm_celsr1_deg_sig$cluster == "High"),]

celsr1_allt_hm_common = celsr1_deg_sig_low[which(! celsr1_deg_sig_low$gene %in% hm_celsr1_deg_sig_low$gene ),]
celsr1_allt_hm_common = rbind(celsr1_allt_hm_common, celsr1_deg_sig_medium[which(! celsr1_deg_sig_medium$gene %in% hm_celsr1_deg_sig_medium$gene ),])
celsr1_allt_hm_common = rbind(celsr1_allt_hm_common, celsr1_deg_sig_high[which(! celsr1_deg_sig_high$gene %in% hm_celsr1_deg_sig_high$gene ),])
write.csv(celsr1_allt_hm_common, "~/scratch/d_tooth/results/igor/allt/allt_celsr1_deg_sig_vs_hm_celsr1_deg_sig_not_common.csv")

num_celsr1_deg_sig = data.frame(aggregate(gene ~ cluster, data = celsr1_deg_sig, length))
colnames(num_celsr1_deg_sig)[2] = "all"
num_celsr1_deg_sig$not_common = aggregate(gene ~ cluster, data = celsr1_allt_hm_common, length)[,2]
num_celsr1_deg_sig = reshape2::melt(num_celsr1_deg_sig)
png("~/scratch/d_tooth/results/igor/allt//allt_celsr1_deg_sig_vs_hm_celsr1_deg_sig_not_common.png", width = 800, height = 600, res = 100)
print(ggplot(num_celsr1_deg_sig, aes(x=cluster, y = value, color = variable, fill = variable)) + geom_bar(stat = "identity", alpha = 0.8, position=position_dodge2()) + scale_color_manual(values = c("#2dcc84", "orange")) + scale_fill_manual(values = c("#2dcc84", "orange")) + xlab("") + ylab("Number of DEGs"))
dev.off()

# Human Celsr1+ Correlation vs Cichlid+Mouse Celsr1+ Correlation
# hm_celsr1_cor = newCytoCor(hm_celsr1)
allt_celsr1_cor = read.csv("allt_celsr1_cor.csv")[,2]
names(allt_celsr1_cor) = read.csv("allt_celsr1_cor.csv")[,1]
hm_celsr1_cor = read.csv("hm_celsr1_cor.csv")[,2]
names(hm_celsr1_cor) = read.csv("hm_celsr1_cor.csv")[,1]

common_gene = names(allt_celsr1_cor)[which(names(allt_celsr1_cor) %in% names(hm_celsr1_cor))]
cor_df = data.frame(gene = common_gene, c_m = allt_celsr1_cor[match(common_gene, names(allt_celsr1_cor))], hm = hm_celsr1_cor[match(common_gene, names(hm_celsr1_cor))])
cor_df$dif = cor_df$c_m - cor_df$hm
png("~/scratch/d_tooth/results/igor/allt/cm_vs_hm_cyto_cor.png", width = 800, height = 600, res = 100)
ggplot(cor_df, aes(x = hm, y = c_m, color = abs(dif))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + xlab("R of Gene w/ CytoTRACE in Human Celsr1+ Cells") + ylab("R of Gene w/ CytoTRACE in Cichlid+Mouse Celsr1+ Cells") + geom_text_repel(data = cor_df[which(cor_df$dif > 0.55),], aes(label = gene), color = plasma(100)[10])
dev.off()

common_gene = names(hm_celsr1_cor)[which(names(hm_celsr1_cor) %in% names(hm_no_celsr1_cor))]
cor_df = data.frame(gene = common_gene, hm_celsr1 = hm_celsr1_cor[match(common_gene, names(hm_celsr1_cor))], hm_no_celsr1 = hm_no_celsr1_cor[match(common_gene, names(hm_no_celsr1_cor))])
cor_df$dif = cor_df$hm_celsr1 - cor_df$hm_no_celsr1
cor_df$num_in_celsr1 = ncol(hm_celsr1)
cor_df$num_in_no_celsr1 = ncol(hm_no_celsr1)
# allt_cor_df$num_in_celsr1 = allt_celsr1_cells[match(allt_cor_df$gene, names(allt_celsr1_cells))]
# allt_cor_df$num_in_no_celsr1 = allt_no_celsr1_cells[match(allt_cor_df$gene, names(allt_no_celsr1_cells))]
cor_df$p = r_to_p(cor_df$hm_celsr1, cor_df$hm_no_celsr1, cor_df$num_in_celsr1, cor_df$num_in_no_celsr1)
cor_df$bon = p.adjust(cor_df$p, method = "bonferroni")
cor_df$neg_log_bon = -log10(cor_df$bon)
cor_df$isSig = cor_df$bon < 0.05
pal = colorRamp(plasma(10))
cor_df$col = rgb(pal( cor_df$neg_log_bon/max(cor_df$neg_log_bon) ), maxColorValue = 255)
cor_df$col[which(cor_df$neg_log_bon < 30 | abs(cor_df$dif) < 0.3)] = "darkgray"
png("hm_yc_vs_nc_cor_p.png", width = 800, height = 600, res = 100)
# ggplot(allt_cor_df, aes(x = dif, y = neg_log_bon, color = neg_log_bon)) + geom_point() + geom_hline(yintercept = -log10(0.05), linetype="dashed") + xlab("R Celsr1+ - R Celsr1-") + ylab(expression(-Log["10"]*" Adjusted P")) + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + geom_text_repel(data = allt_cor_df[which(allt_cor_df$neg_log_bon > 60),], aes(label = gene))
ggplot(cor_df, aes(x = dif, y = neg_log_bon, color = col)) + geom_point() + geom_hline(yintercept = 30, linetype="dashed") + geom_vline(xintercept = c(-0.3, 0.3), linetype="dashed") + xlab("R Celsr1+ - R Celsr1-") + ylab(expression(-Log["10"]*" Adjusted P")) + scale_color_identity() + theme_bw() + geom_text_repel(data = cor_df[which(cor_df$neg_log_bon > 60),], aes(label = gene))
dev.off()
allt_cor_df[,c("hm_bon")] = cor_df[match(allt_cor_df$gene, cor_df$gene), c("bon")]
allt_cor_df$hm_bon[which( is.na(allt_cor_df$hm_bon) )] = 1
allt_cor_df$hm_neg_log_bon = -log10(allt_cor_df$hm_bon)
png("allt_hm_yc_vs_nc_cor_p.png", width = 800, height = 600, res = 100)
# ggplot(allt_cor_df, aes(x = dif, y = neg_log_bon, color = neg_log_bon)) + geom_point() + geom_hline(yintercept = -log10(0.05), linetype="dashed") + xlab("R Celsr1+ - R Celsr1-") + ylab(expression(-Log["10"]*" Adjusted P")) + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + geom_text_repel(data = allt_cor_df[which(allt_cor_df$neg_log_bon > 60),], aes(label = gene))
ggplot(allt_cor_df, aes(x = hm_neg_log_bon, y = neg_log_bon, color = col)) + geom_point() + geom_hline(yintercept = 30, linetype="dashed") + geom_vline(xintercept = 30, linetype="dashed") + xlab(expression(-Log["10"]*" Adjusted P of Celsr1+ vs Celsr1- R Dif in Human")) + ylab(expression(-Log["10"]*" Adjusted P of Celsr1+ vs Celsr1- R Dif in C+M")) + scale_color_identity() + theme_bw() + geom_text_repel(data = allt_cor_df[which(allt_cor_df$neg_log_bon > 60),], aes(label = gene))
dev.off()

# CM Celsr1+ Cor and Number of Cells
mat_allt_celsr1 = allt_celsr1@assays$RNA@counts
mat_allt_celsr1[which(mat_allt_celsr1 > 1)] = 1
allt_celsr1_cells = rowSums(mat_allt_celsr1)
allt_cor_num_df = data.frame(gene = names(allt_celsr1_cor), cor = unname(allt_celsr1_cor), num_cells = unname(allt_celsr1_cells[match(names(allt_celsr1_cor), names(allt_celsr1_cells))]))
png("allt_cor_num_cells.png", width = 800, height = 600, res = 100)
ggplot(allt_cor_num_df, aes(x = num_cells, y = cor, color = cor)) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + xlab("Number of Cells") + ylab("Correlation w/ CytoTRACE")
dev.off()

# Human Celsr1+ CytoBINs vs Cichlid+Mouse Celsr1+ CytoBINs
cm_hm_celsr1 = merge(allt_celsr1, hm_celsr1)
cm_hm_celsr1$celsr1_bin = "None"
cm_hm_celsr1$celsr1_bin[colnames(allt_celsr1)] = paste0("cm_", allt_celsr1$bin)
cm_hm_celsr1$celsr1_bin[colnames(hm_celsr1)] = paste0("hm_", hm_celsr1$bin)
Idents(cm_hm_celsr1) = cm_hm_celsr1$celsr1_bin
cm_hm_hvh_deg = FindMarkers(cm_hm_celsr1, ident.1="cm_High", ident.2="hm_High")
cm_hm_mvm_deg = FindMarkers(cm_hm_celsr1, ident.1="cm_Medium", ident.2="hm_Medium")
cm_hm_lvl_deg = FindMarkers(cm_hm_celsr1, ident.1="cm_Low", ident.2="hm_Low")

cm_hm_hvh_deg$cluster = "High"
cm_hm_mvm_deg$cluster = "Medium"
cm_hm_lvl_deg$cluster = "Low"
cm_hm_hvh_deg$gene = rownames(cm_hm_hvh_deg)
cm_hm_mvm_deg$gene = rownames(cm_hm_mvm_deg)
cm_hm_lvl_deg$gene = rownames(cm_hm_lvl_deg)
bin_deg_dif2 = rbind(cm_hm_hvh_deg, cm_hm_mvm_deg, cm_hm_lvl_deg)
bin_deg_dif2_sig = bin_deg_dif2[which(bin_deg_dif2$p_val_adj < 0.05),]
write.csv(bin_deg_dif2, "~/scratch/d_tooth/results/igor/allt/allt_celsr1_vs_hm_celsr1_deg.csv")
p_df = as.data.frame(table(bin_deg_dif2_sig$cluster))
png("~/scratch/d_tooth/results/igor/allt/allt_celsr1_vs_hm_celsr1_deg_sig.png", width = 800, height = 600, res = 100)
print(ggplot(p_df, aes(x=Var1, y = Freq)) + geom_bar(stat = "identity", alpha = 0.8, position=position_dodge2()) + xlab("") + ylab("Number of DEGs"))
dev.off()

# Human vs Cichlid+Mouse Celsr1+ Volcano Plot
bin_deg_dif2$col = plyr::revalue(bin_deg_dif2$cluster, replace = c("High" = temp[10], "Medium" = temp[6], "Low" = temp[2]))
bin_deg_dif2$col[which(bin_deg_dif2$p_val_adj >= 0.05)] = "darkgray"
bin_deg_dif2$isSig = bin_deg_dif2$p_val_adj >= 0.05
png("~/scratch/d_tooth/results/igor/allt/allt_celsr1_vs_hm_celsr1_deg_sig_volcano.png", width = 800, height = 600, res = 100)
ggplot(bin_deg_dif2, aes(x = avg_log2FC, y = -log10(p_val), color = col, alpha = isSig)) + geom_point() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_y_sqrt() + theme_light() + scale_color_identity() + scale_alpha_manual(values = c(0.6, 0.1), guide = F)
dev.off()

# HM Celsr1-
hm_no_celsr1 = subset(hm, cells = colnames(hm)[which(! colnames(hm) %in% colnames(hm_celsr1))])
hm_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(hm_no_celsr1@assays$RNA@counts))
hm_no_celsr1$cyto = hm_no_celsr1_cyto$CytoTRACE

# R diff in C+M (Celsr1+ - Celsr1-) vs R diff in (Human Celsr1+ - C+M Celsr1+)
common_gene = names(allt_celsr1_cor)[which(names(allt_celsr1_cor) %in% names(hm_celsr1_cor))]
cor_df = data.frame(gene = common_gene, c_m = allt_celsr1_cor[match(common_gene, names(allt_celsr1_cor))], hm = hm_celsr1_cor[match(common_gene, names(hm_celsr1_cor))])
cor_df$dif = cor_df$c_m - cor_df$hm
common_allt_gene = names(allt_celsr1_cor)[which(names(allt_celsr1_cor) %in% names(allt_no_celsr1_cor))]
allt_cor_df = data.frame(gene = common_allt_gene, celsr1 = allt_celsr1_cor[match(common_allt_gene, names(allt_celsr1_cor))], no_celsr1 = allt_no_celsr1_cor[match(common_allt_gene, names(allt_no_celsr1_cor))])
allt_cor_df$dif = allt_cor_df$celsr1 - allt_cor_df$no_celsr1
allt_cor_df$hm_celsr1 = cor_df$hm[match(allt_cor_df$gene, cor_df$gene)]
allt_cor_df = allt_cor_df[which( !is.na(allt_cor_df$hm_celsr1) ),]
allt_cor_df$cm_hm_celsr1_dif = allt_cor_df$celsr1 - allt_cor_df$hm_celsr1
colnames(allt_cor_df)[4] = "celsr1_no_celsr1_dif"
png("~/scratch/d_tooth/results/igor/allt/r_dif.png", width = 800, height = 600, res = 100)
ggplot(allt_cor_df, aes(x = celsr1_no_celsr1_dif, y = cm_hm_celsr1_dif, color = celsr1_no_celsr1_dif + cm_hm_celsr1_dif)) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + xlab("R Celsr1+ - R Celsr1- in C+M") + ylab("(R Celsr1+ C+M) - (R Celsr1+ Human)") + geom_text_repel(data=allt_cor_df[which(allt_cor_df$celsr1_no_celsr1_dif + allt_cor_df$cm_hm_celsr1_dif > 0.8),], aes(label = gene))
dev.off()

# R diff in C+M (Celsr1+ - Celsr1-) vs R diff in Human (Celsr1+ - Celsr1-)
common_gene = names(hm_celsr1_cor)[which(names(hm_celsr1_cor) %in% names(hm_no_celsr1_cor))]
cor_df = data.frame(gene = common_gene, hm_celsr1 = hm_celsr1_cor[match(common_gene, names(hm_celsr1_cor))], hm_no_celsr1 = hm_no_celsr1_cor[match(common_gene, names(hm_no_celsr1_cor))])
cor_df$dif = cor_df$hm_celsr1 - cor_df$hm_no_celsr1
common_allt_gene = names(allt_celsr1_cor)[which(names(allt_celsr1_cor) %in% names(allt_no_celsr1_cor))]
allt_cor_df = data.frame(gene = common_allt_gene, celsr1 = allt_celsr1_cor[match(common_allt_gene, names(allt_celsr1_cor))], no_celsr1 = allt_no_celsr1_cor[match(common_allt_gene, names(allt_no_celsr1_cor))])
allt_cor_df$cm_dif = allt_cor_df$celsr1 - allt_cor_df$no_celsr1
allt_cor_df[,c("hm_celsr1", "hm_no_celsr1")] = cor_df[match(allt_cor_df$gene, cor_df$gene), c("hm_celsr1", "hm_no_celsr1")]
allt_cor_df = allt_cor_df[which( !is.na(allt_cor_df$hm_celsr1) ),]
allt_cor_df$hm_dif = allt_cor_df$hm_celsr1 - allt_cor_df$hm_no_celsr1
allt_cor_df$dif_dif = allt_cor_df$cm_dif - allt_cor_df$hm_dif
png("~/scratch/d_tooth/results/igor/allt/hm_celsr1_vs_no_celsr1_cor.png", width = 800, height = 600, res = 100)
ggplot(cor_df, aes(x = hm_no_celsr1, y = hm_celsr1, color = abs(hm_celsr1 - hm_no_celsr1))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + xlab("R Celsr1- Human") + ylab("R Celsr1+ Human") + geom_text_repel(data=cor_df[which(cor_df$hm_celsr1 - cor_df$hm_no_celsr1 > 0.5),], aes(label = gene))
dev.off()
png("~/scratch/d_tooth/results/igor/allt/allt_hm_celsr1_vs_no_celsr1_cor.png", width = 800, height = 600, res = 100)
ggplot(allt_cor_df, aes(x = hm_dif, y = cm_dif, color = abs(cm_dif - hm_dif))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + theme_bw() + xlab("R Celsr1+ - R Celsr1- in Human") + ylab("R Celsr1+ - R Celsr1- Cichlid & Mouse") + geom_text_repel(data=allt_cor_df[which(allt_cor_df$cm_dif - allt_cor_df$hm_dif > 0.5),], aes(label = gene))
dev.off()

# Talha's Genius Idea: (Celsr1+ vs Celsr1- in C+M) vs (Celsr1+ vs Celsr1- in Human)
hm_no_celsr1$bin <- hm_no_celsr1$cyto
hm_no_celsr1$bin[which(hm_no_celsr1$cyto <= quantile(hm_no_celsr1$cyto, 0.33))] <- "Low"
hm_no_celsr1$bin[which(hm_no_celsr1$cyto > quantile(hm_no_celsr1$cyto, 0.33) & hm_no_celsr1$cyto <= quantile(hm_no_celsr1$cyto, 0.66))] <- "Medium"
hm_no_celsr1$bin[which(hm_no_celsr1$cyto > quantile(hm_no_celsr1$cyto, 0.66))] <- "High"
Idents(hm_no_celsr1) = hm_no_celsr1$bin
hm_no_celsr1_deg = FindAllMarkers(hm_no_celsr1, only.pos = F, logfc.threshold = 0)
hm_no_celsr1_deg$isSig = hm_no_celsr1_deg$p_val_adj < 0.05
hm_no_celsr1_deg$cluster_isSig = paste0(hm_no_celsr1_deg$cluster, "_", hm_no_celsr1_deg$isSig)
hm_no_celsr1_deg_sig = hm_no_celsr1_deg[which(hm_no_celsr1_deg$p_val_adj < 0.05),]
write.csv(hm_no_celsr1_deg, "hm_no_celsr1_deg.csv")

hm_celsr1_deg = read.csv("hm_celsr1_deg.csv")
hm_celsr1_deg_sig = hm_celsr1_deg[which(hm_celsr1_deg$p_val_adj < 0.05),]
hm_celsr1_deg_sig_low    = hm_celsr1_deg_sig[which(hm_celsr1_deg_sig$cluster == "Low"),]
hm_celsr1_deg_sig_medium = hm_celsr1_deg_sig[which(hm_celsr1_deg_sig$cluster == "Medium"),]
hm_celsr1_deg_sig_high   = hm_celsr1_deg_sig[which(hm_celsr1_deg_sig$cluster == "High"),]
hm_no_celsr1_deg_sig_low    = hm_no_celsr1_deg_sig[which(hm_no_celsr1_deg_sig$cluster == "Low"),]
hm_no_celsr1_deg_sig_medium = hm_no_celsr1_deg_sig[which(hm_no_celsr1_deg_sig$cluster == "Medium"),]
hm_no_celsr1_deg_sig_high   = hm_no_celsr1_deg_sig[which(hm_no_celsr1_deg_sig$cluster == "High"),]

hm_celsr1_y_n_not_common = hm_celsr1_deg_sig_low[which(! hm_celsr1_deg_sig_low$gene %in% hm_no_celsr1_deg_sig_low$gene ),]
hm_celsr1_y_n_not_common = rbind(hm_celsr1_y_n_not_common, hm_celsr1_deg_sig_medium[which(! hm_celsr1_deg_sig_medium$gene %in% hm_no_celsr1_deg_sig_medium$gene ),])
hm_celsr1_y_n_not_common = rbind(hm_celsr1_y_n_not_common, hm_celsr1_deg_sig_high[which(! hm_celsr1_deg_sig_high$gene %in% hm_no_celsr1_deg_sig_high$gene ),])
write.csv(hm_celsr1_y_n_not_common, "hm_celsr1_y_n_not_common_big.csv")

hm_celsr1_y_n_not_common$cluster = factor(hm_celsr1_y_n_not_common$cluster, levels = c("Low", "Medium", "High"))
num_celsr1_deg_sig = data.frame(aggregate(gene ~ cluster, data = hm_celsr1_deg_sig, length))
colnames(num_celsr1_deg_sig)[2] = "All Celsr1+ DEG"
aggr_tmp = aggregate(gene ~ cluster, data = hm_celsr1_y_n_not_common, length)
num_celsr1_deg_sig$not_common = aggr_tmp$gene[match(num_celsr1_deg_sig$cluster, aggr_tmp$cluster)]
num_celsr1_deg_sig$not_common[which(is.na(num_celsr1_deg_sig$not_common))] = 0
colnames(num_celsr1_deg_sig)[3] = "Celsr1+ DEG and Not Celsr1-"
num_celsr1_deg_sig = reshape2::melt(num_celsr1_deg_sig)
png("~/scratch/d_tooth/results/igor/allt/hm_celsr1_deg_sig_vs_no_celsr1_big.png", width = 800, height = 600, res = 100)
print(ggplot(num_celsr1_deg_sig, aes(x=cluster, y = value, color = variable, fill = variable)) + geom_bar(stat = "identity", alpha = 0.8, position=position_dodge2()) + scale_color_manual(values = c("#2dcc84", "gold")) + scale_fill_manual(values = c("#2dcc84", "gold")) + xlab("") + ylab("Number of DEGs"))
dev.off()

allt_celsr1_y_n_not_common =  read.csv("allt_celsr1_deg_sig_and_not_a_no_celsr1_deg.csv")
allt_celsr1_y_n_not_common_low = allt_celsr1_y_n_not_common[which(allt_celsr1_y_n_not_common$cluster == "Low"),]
allt_celsr1_y_n_not_common_med = allt_celsr1_y_n_not_common[which(allt_celsr1_y_n_not_common$cluster == "Medium"),]
allt_celsr1_y_n_not_common_high = allt_celsr1_y_n_not_common[which(allt_celsr1_y_n_not_common$cluster == "High"),]
hm_celsr1_y_n_not_common_low = hm_celsr1_y_n_not_common[which(hm_celsr1_y_n_not_common$cluster == "Low"),]
hm_celsr1_y_n_not_common_med = hm_celsr1_y_n_not_common[which(hm_celsr1_y_n_not_common$cluster == "Medium"),]
hm_celsr1_y_n_not_common_high = hm_celsr1_y_n_not_common[which(hm_celsr1_y_n_not_common$cluster == "High"),]

cm_celsr1_spe = allt_celsr1_y_n_not_common_low[which(! allt_celsr1_y_n_not_common_low$gene %in% hm_celsr1_y_n_not_common_low$gene ),]
cm_celsr1_spe = rbind(cm_celsr1_spe, allt_celsr1_y_n_not_common_med[which(! allt_celsr1_y_n_not_common_med$gene %in% hm_celsr1_y_n_not_common_med$gene ),])
cm_celsr1_spe = rbind(cm_celsr1_spe, allt_celsr1_y_n_not_common_high[which(! allt_celsr1_y_n_not_common_high$gene %in% hm_celsr1_y_n_not_common_high$gene ),])
write.csv(cm_celsr1_spe, "cm_celsr1_special.csv")

# Avg LogFC of Human Celsr1+ vs C+M CElsr1+ in CM Celsr1+ Special Genes
hm_celsr1_pct_fc_low = pct_dif_avg_logFC(hm_celsr1, cells.1 = colnames(hm_celsr1)[which(hm_celsr1$bin == "Low")], cells.2 = colnames(hm_celsr1)[which(hm_celsr1$bin != "Low")])
hm_celsr1_pct_fc_med = pct_dif_avg_logFC(hm_celsr1, cells.1 = colnames(hm_celsr1)[which(hm_celsr1$bin == "Medium")], cells.2 = colnames(hm_celsr1)[which(hm_celsr1$bin != "Medium")])
hm_celsr1_pct_fc_high = pct_dif_avg_logFC(hm_celsr1, cells.1 = colnames(hm_celsr1)[which(hm_celsr1$bin == "High")], cells.2 = colnames(hm_celsr1)[which(hm_celsr1$bin != "High")])
allt_celsr1_pct_fc_low = pct_dif_avg_logFC(allt_celsr1, cells.1 = colnames(allt_celsr1)[which(allt_celsr1$bin == "Low")], cells.2 = colnames(allt_celsr1)[which(allt_celsr1$bin != "Low")])
allt_celsr1_pct_fc_med = pct_dif_avg_logFC(allt_celsr1, cells.1 = colnames(allt_celsr1)[which(allt_celsr1$bin == "Medium")], cells.2 = colnames(allt_celsr1)[which(allt_celsr1$bin != "Medium")])
allt_celsr1_pct_fc_high = pct_dif_avg_logFC(allt_celsr1, cells.1 = colnames(allt_celsr1)[which(allt_celsr1$bin == "High")], cells.2 = colnames(allt_celsr1)[which(allt_celsr1$bin != "High")])

hm_celsr1_pct_fc_low$cluster = "Low"
hm_celsr1_pct_fc_med$cluster = "Medium"
hm_celsr1_pct_fc_high$cluster = "High"
hm_celsr1_pct_fc_low$dataset = "HM"
hm_celsr1_pct_fc_med$dataset = "HM"
hm_celsr1_pct_fc_high$dataset = "HM"

allt_celsr1_pct_fc_low$cluster = "Low"
allt_celsr1_pct_fc_med$cluster = "Medium"
allt_celsr1_pct_fc_high$cluster = "High"
allt_celsr1_pct_fc_low$dataset = "CM"
allt_celsr1_pct_fc_med$dataset = "CM"
allt_celsr1_pct_fc_high$dataset = "CM"

cm_celsr1_spe_pct_fc = cbind(hm_celsr1_pct_fc_low[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Low")], hm_celsr1_pct_fc_low$gene),], allt_celsr1_pct_fc_low[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Low")], allt_celsr1_pct_fc_low$gene),])
cm_celsr1_spe_pct_fc = rbind(cm_celsr1_spe_pct_fc, cbind(hm_celsr1_pct_fc_med[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Medium")], hm_celsr1_pct_fc_med$gene),], allt_celsr1_pct_fc_med[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Medium")], allt_celsr1_pct_fc_med$gene),]))
cm_celsr1_spe_pct_fc = rbind(cm_celsr1_spe_pct_fc, cbind(hm_celsr1_pct_fc_high[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "High")], hm_celsr1_pct_fc_high$gene),], allt_celsr1_pct_fc_high[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "High")], allt_celsr1_pct_fc_high$gene),]))
colnames(cm_celsr1_spe_pct_fc) = c(paste0("hm.", colnames(hm_celsr1_pct_fc_low)), paste0("cm.", colnames(allt_celsr1_pct_fc_low)))
cm_celsr1_spe_pct_fc$hm.avg_logFC[which( is.na(cm_celsr1_spe_pct_fc$hm.avg_logFC))] = 0
cm_celsr1_spe_pct_fc$cm.cluster = factor(cm_celsr1_spe_pct_fc$cm.cluster, levels = c("Low", "Medium", "High"))

png("~/scratch/d_tooth/results/igor/allt/allt_and_hm_fc_at_cm_celsr1_special_big.png", width = 800, height = 600, res = 100)
ggplot(cm_celsr1_spe_pct_fc, aes(x = hm.avg_logFC, y = cm.avg_logFC, color = abs(cm.avg_logFC - hm.avg_logFC))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + geom_text_repel(data=cm_celsr1_spe_pct_fc[which( abs(cm_celsr1_spe_pct_fc$cm.avg_logFC - cm_celsr1_spe_pct_fc$hm.avg_logFC) > 4 ),], aes(label = cm.genes)) + xlab("Human Celsr1+ Average LogFC") + ylab("Cichlid and Mouse Celsr1+ Average LogFC") + theme_bw() 
dev.off()
png("~/scratch/d_tooth/results/igor/allt/allt_and_hm_fc_at_cm_celsr1_special_cluster_big.png", width = 800, height = 600, res = 100)
ggplot(cm_celsr1_spe_pct_fc, aes(x = hm.avg_logFC, y = cm.avg_logFC, color = cm.cluster)) + geom_point(alpha = 0.7) + scale_color_manual(values = c(temp[2], temp[6], temp[10]), guide = F) + geom_text_repel(data=cm_celsr1_spe_pct_fc[which( abs(cm_celsr1_spe_pct_fc$cm.avg_logFC - cm_celsr1_spe_pct_fc$hm.avg_logFC) > 4 ),], aes(label = cm.genes)) + xlab("Human Celsr1+ Average LogFC") + ylab("Cichlid and Mouse Celsr1+ Average LogFC") + theme_bw() 
dev.off()

# Avg LogFC of C+M Celsr1- vs C+M CElsr1+ in CM Celsr1+ Special Genes
allt_no_celsr1_pct_fc_low = pct_dif_avg_logFC(allt_no_celsr1, cells.1 = colnames(allt_no_celsr1)[which(allt_no_celsr1$bin == "Low")], cells.2 = colnames(allt_no_celsr1)[which(allt_no_celsr1$bin != "Low")])
allt_no_celsr1_pct_fc_med = pct_dif_avg_logFC(allt_no_celsr1, cells.1 = colnames(allt_no_celsr1)[which(allt_no_celsr1$bin == "Medium")], cells.2 = colnames(allt_no_celsr1)[which(allt_no_celsr1$bin != "Medium")])
allt_no_celsr1_pct_fc_high = pct_dif_avg_logFC(allt_no_celsr1, cells.1 = colnames(allt_no_celsr1)[which(allt_no_celsr1$bin == "High")], cells.2 = colnames(allt_no_celsr1)[which(allt_no_celsr1$bin != "High")])
allt_no_celsr1_pct_fc_low$cluster = "Low"
allt_no_celsr1_pct_fc_med$cluster = "Medium"
allt_no_celsr1_pct_fc_high$cluster = "High"
allt_no_celsr1_pct_fc_low$dataset = "CM-"
allt_no_celsr1_pct_fc_med$dataset = "CM-"
allt_no_celsr1_pct_fc_high$dataset = "CM-"

cm_celsr1_spe_pct_fc = cbind(allt_no_celsr1_pct_fc_low[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Low")], allt_no_celsr1_pct_fc_low$gene),], allt_celsr1_pct_fc_low[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Low")], allt_celsr1_pct_fc_low$gene),])
cm_celsr1_spe_pct_fc = rbind(cm_celsr1_spe_pct_fc, cbind(allt_no_celsr1_pct_fc_med[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Medium")], allt_no_celsr1_pct_fc_med$gene),], allt_celsr1_pct_fc_med[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "Medium")], allt_celsr1_pct_fc_med$gene),]))
cm_celsr1_spe_pct_fc = rbind(cm_celsr1_spe_pct_fc, cbind(allt_no_celsr1_pct_fc_high[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "High")], allt_no_celsr1_pct_fc_high$gene),], allt_celsr1_pct_fc_high[match(cm_celsr1_spe$gene[which(cm_celsr1_spe$cluster == "High")], allt_celsr1_pct_fc_high$gene),]))
colnames(cm_celsr1_spe_pct_fc) = c(paste0("nc.", colnames(allt_no_celsr1_pct_fc_high)), paste0("yc.", colnames(allt_celsr1_pct_fc_low)))
cm_celsr1_spe_pct_fc$nc.avg_logFC[which( is.na(cm_celsr1_spe_pct_fc$nc.avg_logFC) )] = 0
cm_celsr1_spe_pct_fc$yc.cluster = factor(cm_celsr1_spe_pct_fc$yc.cluster, levels = c("Low", "Medium", "High"))
cm_celsr1_spe_pct_fc$avg_logFC_dif = cm_celsr1_spe_pct_fc$yc.avg_logFC - cm_celsr1_spe_pct_fc$nc.avg_logFC
cm_celsr1_spe_pct_fc$pct_dif_dif = cm_celsr1_spe_pct_fc$yc.pct_dif - cm_celsr1_spe_pct_fc$nc.pct_dif

png("~/scratch/d_tooth/results/igor/allt/allt_yc_and_nc_fc_at_cm_celsr1_special_big.png", width = 800, height = 600, res = 100)
ggplot(cm_celsr1_spe_pct_fc, aes(x = nc.avg_logFC, y = yc.avg_logFC, color = abs(yc.avg_logFC - nc.avg_logFC))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + geom_text_repel(data=cm_celsr1_spe_pct_fc[which( abs(cm_celsr1_spe_pct_fc$yc.avg_logFC - cm_celsr1_spe_pct_fc$nc.avg_logFC) > 4 ),], aes(label = yc.genes)) + xlab("Cichlid and Mouse Celsr1- Average LogFC") + ylab("Cichlid and Mouse Celsr1+ Average LogFC") + theme_bw() 
dev.off()
png("~/scratch/d_tooth/results/igor/allt/allt_yc_and_nc_fc_at_cm_celsr1_special_cluster_big.png", width = 800, height = 600, res = 100)
ggplot(cm_celsr1_spe_pct_fc, aes(x = nc.avg_logFC, y = yc.avg_logFC, color = yc.cluster)) + geom_point() + scale_color_manual(values = c(temp[2], temp[6], temp[10]), guide = F) + geom_text_repel(data=cm_celsr1_spe_pct_fc[which( abs(cm_celsr1_spe_pct_fc$yc.avg_logFC - cm_celsr1_spe_pct_fc$nc.avg_logFC) > 4 ),], aes(label = yc.genes)) + xlab("Cichlid and Mouse Celsr1- Average LogFC") + ylab("Cichlid and Mouse Celsr1+ Average LogFC") + theme_bw() 
dev.off()
png("~/scratch/d_tooth/results/igor/allt/allt_yc_and_nc_pct_fc_at_cm_celsr1_special_big.png", width = 800, height = 600, res = 100)
ggplot(cm_celsr1_spe_pct_fc, aes(x = pct_dif_dif, y = avg_logFC_dif, color = abs(avg_logFC_dif - pct_dif_dif/10))) + geom_point() + scale_color_gradientn(colors = plasma(100), guide = F) + geom_text_repel(data=cm_celsr1_spe_pct_fc[which( abs(cm_celsr1_spe_pct_fc$yc.avg_logFC - cm_celsr1_spe_pct_fc$nc.avg_logFC) > 2.8 ),], aes(label = yc.genes)) + xlab("% Difference Dif (Celsr1+ vs Celsr1- in C+M)") + ylab("Average LogFC Dif (Celsr1+ vs Celsr1- in C+M)") + theme_bw() 
dev.off()

# Avg LogFC Diff vs Pct Dif Diff
png("~/scratch/d_tooth/results/igor/allt/test.png", width = 800, height = 600, res = 100)
ggplot(cm_celsr1_spe_pct_fc, aes(x = pct_dif_dif, y = abs(avg_logFC_dif), color = abs(avg_logFC_dif) > 0.25 & abs(pct_dif_dif) > 10)) + geom_point() + geom_text_repel(data=cm_celsr1_spe_pct_fc[which( abs(cm_celsr1_spe_pct_fc$yc.avg_logFC - cm_celsr1_spe_pct_fc$nc.avg_logFC) > 4 ),], aes(label = yc.genes)) + xlab("% Difference Dif (Celsr1+ vs Celsr1- in C+M)") + ylab("Average LogFC Dif (Celsr1+ vs Celsr1- in C+M)") + theme_bw() 
dev.off()

# Combo Approach
combo = cm_celsr1_spe[which( cm_celsr1_spe$gene %in% allt_cor_df$gene[which(allt_cor_df$bon < 0.05 & allt_cor_df$hm_bon >= 0.05)] ),]
combo[,c("r_celsr1", "r_no_celsr1", "r_dif", "num_celsr1", "num_no_celsr1", "p", "cm_bon", "hm_bon")] = allt_cor_df[match(combo$gene, allt_cor_df$gene), c("celsr1", "no_celsr1", "dif", "num_in_celsr1", "num_in_no_celsr1", "p", "bon", "hm_bon")]
combo = combo[order(combo$p_val_adj),]


#*************************************************************************************************
# CM Celsr1+ Special ToppGene Circle Plot ========================================================
#*************************************************************************************************
toppgene_df = read.table("C:/Users/miles/Downloads/cm_celsr1_spe_big_toppgene.txt", sep = "\t", header = T, nrow = 19490)
toppgene_df = toppgene_df[which(! toppgene_df$Category %in% c("Pubmed", "Pathway")),]
toppgene_df = toppgene_df[which(toppgene_df$q.value.Bonferroni < 0.05),]
toppgene_df$neg_log_bon = -log10(toppgene_df$q.value.Bonferroni)

pdf("C:/Users/miles/Downloads/cm_celsr1_spe_big_toppgene.pdf", width = 10, height = 10)
toppgene_df = toppgene_df[order(toppgene_df$q.value.Bonferroni, decreasing = F), ]
toppgene_df = toppgene_df[1:100,]
toppgene_df$col = plyr::revalue(toppgene_df$Category, replace = c("GO: Biological Process" = "#00A08A", "GO: Cellular Component" = "#FF0000", "GO: Molecular Function" = "#F2AD00", "Mouse Phenotype" = "#264653"))
circos.par("gap.after" = c(rep(0, 99), 20), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.02, 0.02))
circos.initialize(sectors = toppgene_df$Name, xlim = c(0,1))
# circos.track(ylim = c(0, max(toppgene_df$cichlid_neg_log10_p)), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
circos.track(ylim = c(0, ceiling(max(toppgene_df$neg_log_bon))), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))

  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], toppgene_df$neg_log_bon[CELL_META$sector.numeric.index], col = paste0(toppgene_df$col[CELL_META$sector.numeric.index], "90"), border = toppgene_df$col[CELL_META$sector.numeric.index])
  # circos.rect(CELL_META$cell.xlim[1], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index], CELL_META$cell.xlim[2], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mi_neg_log10_p[CELL_META$sector.numeric.index], col = "#FF000090", border = "#FF0000")
  # circos.rect(CELL_META$cell.xlim[1], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mi_neg_log10_p[CELL_META$sector.numeric.index], CELL_META$cell.xlim[2], toppgene_df$cichlid_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mi_neg_log10_p[CELL_META$sector.numeric.index]+toppgene_df$mim_neg_log10_p[CELL_META$sector.numeric.index], col = "#F2AD0090",  border = "#F2AD00")

  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(1),
              toppgene_df$Name[CELL_META$sector.numeric.index], facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  circos.segments(CELL_META$cell.xlim[1],ceiling(max(toppgene_df$neg_log_bon)), CELL_META$cell.xlim[2], ceiling(max(toppgene_df$neg_log_bon)), col = "gray20", lty = 1) # top
  circos.segments(CELL_META$cell.xlim[1], ceiling(max(toppgene_df$neg_log_bon)/2), CELL_META$cell.xlim[2], ceiling(max(toppgene_df$neg_log_bon)/2), col = "gray80", lty = 2) # middle
  circos.segments(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 0, col = "gray20", lty = 1) # bottom
})
circos.yaxis(at = c(0, ceiling(max(toppgene_df$neg_log_bon)/2), ceiling(max(toppgene_df$neg_log_bon))), labels = c(0, ceiling(max(toppgene_df$neg_log_bon)/2), ceiling(max(toppgene_df$neg_log_bon))), sector.index = toppgene_df$Name[1], track.index = 1, side = "left", labels.niceFacing = F)
dev.off()

#*************************************************************************************************
# Helper Functions ===============================================================================
#*************************************************************************************************
createCytoBINsInGene = function(obj, gene) {
  gene_pos_cells = colnames(obj)[which(obj@assays$RNA@counts[gene,] > 0)]
  gene_obj = subset(obj, cells = gene_pos_cells)
  gene_obj_cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj@assays$RNA@counts))
  gene_obj$cyto = gene_obj_cyto$CytoTRACE
  
  gene_obj$bin <- gene_obj$cyto
  gene_obj$bin[which(gene_obj$cyto <= quantile(gene_obj$cyto, 0.33))] <- "Low"
  gene_obj$bin[which(gene_obj$cyto > quantile(gene_obj$cyto, 0.33) & gene_obj$cyto <= quantile(gene_obj$cyto, 0.66))] <- "Medium"
  gene_obj$bin[which(gene_obj$cyto > quantile(gene_obj$cyto, 0.66))] <- "High"
  gene_obj$bin <- gene_obj$cyto
  gene_obj$bin[which(gene_obj$cyto <= 0.33)] <- "Low"
  gene_obj$bin[which(gene_obj$cyto > 0.33 & gene_obj$cyto <= 0.66)] <- "Medium"
  gene_obj$bin[which(gene_obj$cyto > 0.66)] <- "High"
  return(list(gene_obj_cyto, gene_obj))
  # Idents(gene_obj) = gene_obj$bin
  # deg = FindAllMarkers(gene_obj, only.pos = F, logfc.threshold = 0)
  # deg = deg[which(deg$p_val_adj < 0.05),]
  # deg$cluster_sign = paste0(sign(deg$avg_log2FC), deg$cluster)
  # return(deg)
}

mergeToothJawCytoBIN = function(obj1, obj2, gene) {
  gene_pos_cells = colnames(obj1)[which(obj1@assays$RNA@counts[gene,] > 0)]
  gene_obj1 = subset(obj1, cells = gene_pos_cells)
  gene_obj1$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj1@assays$RNA@counts))$CytoTRACE
  
  gene_pos_cells = colnames(obj2)[which(obj2@assays$RNA@counts[gene,] > 0)]
  gene_obj2 = subset(obj2, cells = gene_pos_cells)
  gene_obj2$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj2@assays$RNA@counts))$CytoTRACE
  
  merged = merge(gene_obj1, gene_obj2)
  merged$bin <- merged$cyto
  merged$bin[which(merged$cyto <= quantile(merged$cyto, 0.33))] <- "Low"
  merged$bin[which(merged$cyto > quantile(merged$cyto, 0.33) & merged$cyto <= quantile(merged$cyto, 0.66))] <- "Medium"
  merged$bin[which(merged$cyto > quantile(merged$cyto, 0.66))] <- "High"
  print(quantile(merged$cyto, 0.66))
  
  Idents(merged) = merged$bin
  deg = FindAllMarkers(merged, only.pos = F, logfc.threshold = 0)
  deg = deg[which(deg$p_val_adj < 0.05),]
  deg$cluster_sign = paste0(sign(deg$avg_log2FC), deg$cluster)
  return(deg)
}
mergeToothJawCytoBIN2 = function(obj1, obj2, gene) {
  gene_pos_cells = colnames(obj1)[which(obj1@assays$RNA@counts[gene,] > 0)]
  gene_obj1 = subset(obj1, cells = gene_pos_cells)
  gene_obj1$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj1@assays$RNA@counts))$CytoTRACE
  gene_obj1$bin = gene_obj1$cyto
  gene_obj1$bin[which(gene_obj1$cyto <= quantile(gene_obj1$cyto, 0.33))] <- "Low"
  gene_obj1$bin[which(gene_obj1$cyto > quantile(gene_obj1$cyto, 0.33) & gene_obj1$cyto <= quantile(gene_obj1$cyto, 0.66))] <- "Medium"
  gene_obj1$bin[which(gene_obj1$cyto > quantile(gene_obj1$cyto, 0.66))] <- "High"
  print(quantile(gene_obj1$cyto, 0.66))
  
  gene_pos_cells = colnames(obj2)[which(obj2@assays$RNA@counts[gene,] > 0)]
  gene_obj2 = subset(obj2, cells = gene_pos_cells)
  gene_obj2$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj2@assays$RNA@counts))$CytoTRACE
  gene_obj2$bin = gene_obj2$cyto
  gene_obj2$bin[which(gene_obj2$cyto <= quantile(gene_obj2$cyto, 0.33))] <- "Low"
  gene_obj2$bin[which(gene_obj2$cyto > quantile(gene_obj2$cyto, 0.33) & gene_obj2$cyto <= quantile(gene_obj2$cyto, 0.66))] <- "Medium"
  gene_obj2$bin[which(gene_obj2$cyto > quantile(gene_obj2$cyto, 0.66))] <- "High"
  print(quantile(gene_obj2$cyto, 0.66))
  
  merged = merge(gene_obj1, gene_obj2)
  Idents(merged) = merged$bin
  deg = FindAllMarkers(merged, only.pos = F, logfc.threshold = 0)
  deg = deg[which(deg$p_val_adj < 0.05),]
  deg$cluster_sign = paste0(sign(deg$avg_log2FC), deg$cluster)
  
  cytobin_df = data.frame(bin = gene_obj2$bin, exp = colSums(gene_obj2@assays$RNA@data[c("ENSMZEG00005018585", "ppp2r2ca", "mpzl2b"),]), cyto = gene_obj2$cyto, dataset = "Cichlid")
  cytobin_df = rbind(cytobin_df, data.frame(bin = gene_obj2$bin, exp = colSums(gene_obj2@assays$RNA@data[c("ENSMZEG00005018585", "ppp2r2ca", "mpzl2b"),]), cyto = gene_obj2$cyto, dataset = "Cichlid"))
  
  return(deg)
}
mergeToothJawCytoBIN3 = function(obj1, obj2, gene) {
  gene_pos_cells = colnames(obj1)[which(obj1@assays$RNA@counts[gene,] > 0)]
  gene_obj1 = subset(obj1, cells = gene_pos_cells)
  # gene_obj1$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj1@assays$RNA@counts))$CytoTRACE
  gene_obj1$bin = gene_obj1$cyto
  gene_obj1$bin[which(gene_obj1$cyto <= quantile(gene_obj1$cyto, 0.33))] <- "Low"
  gene_obj1$bin[which(gene_obj1$cyto > quantile(gene_obj1$cyto, 0.33) & gene_obj1$cyto <= quantile(gene_obj1$cyto, 0.66))] <- "Medium"
  gene_obj1$bin[which(gene_obj1$cyto > quantile(gene_obj1$cyto, 0.66))] <- "High"
  
  gene_pos_cells = colnames(obj2)[which(obj2@assays$RNA@counts[gene,] > 0)]
  gene_obj2 = subset(obj2, cells = gene_pos_cells)
  # gene_obj2$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj2@assays$RNA@counts))$CytoTRACE
  gene_obj2$bin = gene_obj2$cyto
  gene_obj2$bin[which(gene_obj2$cyto <= quantile(gene_obj2$cyto, 0.33))] <- "Low"
  gene_obj2$bin[which(gene_obj2$cyto > quantile(gene_obj2$cyto, 0.33) & gene_obj2$cyto <= quantile(gene_obj2$cyto, 0.66))] <- "Medium"
  gene_obj2$bin[which(gene_obj2$cyto > quantile(gene_obj2$cyto, 0.66))] <- "High"
  
  merged = merge(gene_obj1, gene_obj2)
  Idents(merged) = merged$bin
  deg = FindAllMarkers(merged, only.pos = F, logfc.threshold = 0)
  deg = deg[which(deg$p_val_adj < 0.05),]
  deg$cluster_sign = paste0(sign(deg$avg_log2FC), deg$cluster)
  return(deg)
}
mergeToothJawCytoBIN4 = function(obj1, obj2, gene) {
  gene_pos_cells = colnames(obj1)[which(obj1@assays$RNA@counts[gene,] > 0)]
  gene_obj1 = subset(obj1, cells = gene_pos_cells)
  # gene_obj1$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj1@assays$RNA@counts))$CytoTRACE
  
  gene_pos_cells = colnames(obj2)[which(obj2@assays$RNA@counts[gene,] > 0)]
  gene_obj2 = subset(obj2, cells = gene_pos_cells)
  # gene_obj2$cyto = CytoTRACE::CytoTRACE(as.matrix(gene_obj2@assays$RNA@counts))$CytoTRACE
  
  merged = merge(gene_obj1, gene_obj2)
  merged$bin <- merged$cyto
  merged$bin[which(merged$cyto <= quantile(merged$cyto, 0.33))] <- "Low"
  merged$bin[which(merged$cyto > quantile(merged$cyto, 0.33) & merged$cyto <= quantile(merged$cyto, 0.66))] <- "Medium"
  merged$bin[which(merged$cyto > quantile(merged$cyto, 0.66))] <- "High"
  
  Idents(merged) = merged$bin
  deg = FindAllMarkers(merged, only.pos = F, logfc.threshold = 0)
  deg = deg[which(deg$p_val_adj < 0.05),]
  deg$cluster_sign = paste0(sign(deg$avg_log2FC), deg$cluster)
  return(deg)
}
cytoCor = function(obj, gene) {
  gene_pos_cells = colnames(obj)[which(obj@assays$RNA@counts[gene,] > 0)]
  gene_obj = subset(obj, cells = gene_pos_cells)
  mat = gene_obj@assays$RNA@counts
  mat[which(mat > 1)] = 1
  mat_rowsums = rowSums(mat)
  non_zero_genes = names(mat_rowsums)[which(mat_rowsums >= 3)]
  gene_cors = unlist(mclapply(non_zero_genes, function(x) cor(gene_obj@assays$RNA@data[x,], gene_obj$cyto), mc.cores = detectCores()))
  names(gene_cors) = non_zero_genes
  return(gene_cors)
}
newCytoCor = function(obj) {
  mat = obj@assays$RNA@counts
  mat[which(mat > 1)] = 1
  mat_rowsums = rowSums(mat)
  non_zero_genes = names(mat_rowsums)[which(mat_rowsums >= 3)]
  gene_cors = unlist(mclapply(non_zero_genes, function(x) cor(obj@assays$RNA@data[x,], obj$cyto), mc.cores = detectCores()))
  names(gene_cors) = non_zero_genes
  return(gene_cors)
}
