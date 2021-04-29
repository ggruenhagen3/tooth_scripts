#=============================================================================================
# Figure 1 ===================================================================================
#=============================================================================================
# *** DotPlot for Tooth ***
# Tooth Data
# deg_tj = read.table("~/research/tooth/results/tj_annot_unique_long_hgnc_100.tsv", sep="\t", header = T, stringsAsFactors = F)
# deg_tj = deg_tj[-which(deg_tj$gene == "ENSMZEG00005000005"),]
# deg_tj2 = as.data.frame(deg_tj %>% group_by(cluster) %>% slice(head(row_number(), 2)))
# deg_tj2 = deg_tj2[order(deg_tj2$cluster),]

deg_tj2 = data.frame(cluster = c(rep("Glia", 3), rep("Immune", 4), rep("Epithelial", 8), rep("Mesenchymal", 8), rep("Endothelial", 6)),
                     gene = c("tcf7", "tmtc2b", "slit3", "ENSMZEG00005022086", "wdr73", "ENSMZEG00005022005", "HBE1 (1 of many).1", "pitx1", "pitx2", "shha", "odam", "lbh", "itgb3b", 'sema3fb', "cldn10a", "plod2 (1 of many)", "sox5", "smpd3", "MMP16", "bmp6", "SMAD6", "sall1a", "dkk1a", "tns2a", "EBF1 (1 of many).1", "tie1", "PTPRB (1 of many).1", "cdh5", "robo4"))
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
pal = c("#E9C46A", "#F6B179", "#2A9D8F", "#E24D28", "#264653")
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