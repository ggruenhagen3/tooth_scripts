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
brain_dir = paste0(main_path, "brain/")
data_dir  = paste0(main_path, "tooth/data/")
out_dir   = paste0(main_path, "tooth/results/")
if (main_path == "/storage/scratch1/6/ggruenhagen3/") { data_dir = "/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/tooth/data/" }
source(paste0(brain_dir, "/brain_scripts/all_f.R"))
setwd(out_dir)

#*******************************************************************************
# Initial clustering ===========================================================
#*******************************************************************************

# Load Data 
b10.counts = Read10X(paste0(data_dir, "JTS13_B10_nuc/outs/filtered_feature_bc_matrix/"))
d2.counts = Read10X(paste0(data_dir, "JTS13_D2_nuc/outs/filtered_feature_bc_matrix/"))
d10.counts = Read10X(paste0(data_dir, "JTS13_D10_nuc/outs/filtered_feature_bc_matrix/"))
g10.counts = Read10X(paste0(data_dir, "JTS13_G10_nuc/outs/filtered_feature_bc_matrix/"))

b10 = CreateSeuratObject(counts = b10.counts, project = "B10")
d2 = CreateSeuratObject(counts = d2.counts, project = "D2")
d10 = CreateSeuratObject(counts = d10.counts, project = "D10")
g10 = CreateSeuratObject(counts = g10.counts, project = "G10")

# Percent mitochondrial
gtf = read.delim("~/research/all_research/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = str_replace(mito.genes, "_", "-")
b10$pct.mt = colSums(b10@assays$RNA@counts[mito.genes,]) / b10$nCount_RNA
d2$pct.mt  = colSums(d2@assays$RNA@counts[mito.genes,])  / d2$nCount_RNA
d10$pct.mt = colSums(d10@assays$RNA@counts[mito.genes,]) / d10$nCount_RNA
g10$pct.mt = colSums(g10@assays$RNA@counts[mito.genes,]) / g10$nCount_RNA

# Add metadata
b10$sample = "B10"; b10$hasRT = T; b10$species.abr = "TM"; b10$species = "Tropheus sp mauve"; 
d2$sample  = "D2";  d2$hasRT  = T; d2$species.abr  = "TM"; d2$species  = "Tropheus sp mauve"; 
d10$sample = "D10"; d10$hasRT = T; d10$species.abr = "TM"; d10$species = "Tropheus sp mauve"; 
g10$sample = "G10"; g10$hasRT = F; g10$species.abr = "MF"; g10$species = "Metriclima fainzilberi"; 

# Load TJ and JPOOL data for comparison purposes
tj.tmp.counts = Read10X(paste0(data_dir, "TJ/outs/filtered_feature_bc_matrix/"))
tj.tmp = CreateSeuratObject(counts = tj.tmp.counts, project = "TJ")
jpool.tmp.counts = Read10X(paste0(data_dir, "JPOOL/outs/filtered_feature_bc_matrix/"))
jpool.tmp = CreateSeuratObject(counts = jpool.tmp.counts, project = "JPOOL")

# Plot the Data Quality: nCount and nFeature
ncount.df = rbind(data.frame(ncount = b10$nCount_RNA, sample = "B10", my.facet = "B10"), data.frame(ncount = tj.tmp$nCount_RNA, sample = "TJ", my.facet = "B10"), data.frame(ncount = jpool.tmp$nCount_RNA, sample = "JPOOL", my.facet = "B10"))
ncount.df = rbind(ncount.df, data.frame(ncount = d2$nCount_RNA, sample = "D2", my.facet = "D2"), data.frame(ncount = tj.tmp$nCount_RNA, sample = "TJ", my.facet = "D2"), data.frame(ncount = jpool.tmp$nCount_RNA, sample = "JPOOL", my.facet = "D2"))
ncount.df = rbind(ncount.df, data.frame(ncount = d10$nCount_RNA, sample = "D10", my.facet = "D10"), data.frame(ncount = tj.tmp$nCount_RNA, sample = "TJ", my.facet = "D10"), data.frame(ncount = jpool.tmp$nCount_RNA, sample = "JPOOL", my.facet = "D10"))
ncount.df = rbind(ncount.df, data.frame(ncount = g10$nCount_RNA, sample = "G10", my.facet = "G10"), data.frame(ncount = tj.tmp$nCount_RNA, sample = "TJ", my.facet = "G10"), data.frame(ncount = jpool.tmp$nCount_RNA, sample = "JPOOL", my.facet = "G10"))
ncount.df$sample = factor(ncount.df$sample, levels = c("TJ", "JPOOL", "B10", "D2", "D10", "G10"))
ncount.df = ncount.df[which(ncount.df$ncount < 5000),]
pdf(paste0(out_dir, "tt_ncount_density.pdf"), height = 8, width = 6, onefile = F)
print(ggplot(ncount.df, aes(x = ncount, fill = sample)) + geom_density(alpha = 0.5, position = "identity", bins = 100) + theme_classic() + ylab("Nuclei Density") + xlab("Number of UMIs per Nuclei") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + facet_wrap(~my.facet, ncol=1))
dev.off()
pdf(paste0(out_dir, "tt_ncount_histogram.pdf"), height = 8, width = 6, onefile = F)
print(ggplot(ncount.df, aes(x = ncount, fill = sample)) + geom_histogram(alpha = 0.5, position = "identity", bins = 100) + theme_classic() + ylab("Number of Nuclei") + xlab("Number of UMIs per Nuclei") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + facet_wrap(~my.facet, ncol=1))
dev.off()

nfeature.df = rbind(data.frame(nfeature = b10$nFeature_RNA, sample = "B10", my.facet = "B10"), data.frame(nfeature = tj.tmp$nFeature_RNA, sample = "TJ", my.facet = "B10"), data.frame(nfeature = jpool.tmp$nFeature_RNA, sample = "JPOOL", my.facet = "B10"))
nfeature.df = rbind(nfeature.df, data.frame(nfeature = d2$nFeature_RNA, sample = "D2", my.facet = "D2"), data.frame(nfeature = tj.tmp$nFeature_RNA, sample = "TJ", my.facet = "D2"), data.frame(nfeature = jpool.tmp$nFeature_RNA, sample = "JPOOL", my.facet = "D2"))
nfeature.df = rbind(nfeature.df, data.frame(nfeature = d10$nFeature_RNA, sample = "D10", my.facet = "D10"), data.frame(nfeature = tj.tmp$nFeature_RNA, sample = "TJ", my.facet = "D10"), data.frame(nfeature = jpool.tmp$nFeature_RNA, sample = "JPOOL", my.facet = "D10"))
nfeature.df = rbind(nfeature.df, data.frame(nfeature = g10$nFeature_RNA, sample = "G10", my.facet = "G10"), data.frame(nfeature = tj.tmp$nFeature_RNA, sample = "TJ", my.facet = "G10"), data.frame(nfeature = jpool.tmp$nFeature_RNA, sample = "JPOOL", my.facet = "G10"))
nfeature.df$sample = factor(nfeature.df$sample, levels = c("TJ", "JPOOL", "B10", "D2", "D10", "G10"))
nfeature.df = nfeature.df[which(nfeature.df$nfeature < 5000),]
pdf(paste0(out_dir, "tt_nfeature_density.pdf"), height = 8, width = 6, onefile = F)
print(ggplot(nfeature.df, aes(x = nfeature, fill = sample)) + geom_density(alpha = 0.5, position = "identity", bins = 100) + theme_classic() + ylab("Nuclei Density") + xlab("Number of Genes per Nuclei") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + facet_wrap(~my.facet, ncol=1))
dev.off()
pdf(paste0(out_dir, "tt_nfeature_histogram.pdf"), height = 8, width = 6, onefile = F)
print(ggplot(nfeature.df, aes(x = nfeature, fill = sample)) + geom_histogram(alpha = 0.5, position = "identity", bins = 100) + theme_classic() + ylab("Number of Nuclei") + xlab("Number of Genes per Nuclei") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + facet_wrap(~my.facet, ncol=1))
dev.off()

# Rename cells: Add the sample name to the cell names
b10 = RenameCells(b10, "b10")
d2  = RenameCells(d2, "d2")
d10 = RenameCells(d10, "d10")
g10 = RenameCells(g10, "g10")

# Remove cells
removed.df = data.frame(cell = colnames(b10)[which(b10$pct.mt > 0.02)], reason = "dead (% MT)")
removed.df = rbind(removed.df, data.frame(cell = colnames(d2)[which(d2$pct.mt   > 0.02)], reason = "dead (% MT)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d10)[which(d10$pct.mt > 0.02)], reason = "dead (% MT)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(g10)[which(g10$pct.mt > 0.02)], reason = "dead (% MT)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(b10)[which(b10$nFeature_RNA > 2500)], reason = "doublet"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d2)[which(d2$nFeature_RNA   > 2500)], reason = "doublet"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d10)[which(d10$nFeature_RNA > 2500)], reason = "doublet"))
removed.df = rbind(removed.df, data.frame(cell = colnames(g10)[which(g10$nFeature_RNA > 2500)], reason = "doublet"))
removed.df = rbind(removed.df, data.frame(cell = colnames(b10)[which(b10$nFeature_RNA < 200)], reason = "dead (# Genes)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d2)[which(d2$nFeature_RNA   < 200)], reason = "dead (# Genes)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d10)[which(d10$nFeature_RNA < 200)], reason = "dead (# Genes)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(g10)[which(g10$nFeature_RNA < 200)], reason = "dead (# Genes)"))
removed.df = rbind(removed.df, data.frame(cell = colnames(b10)[which(!colnames(b10) %in% removed.df$cell)], reason = "kept"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d2)[which(!colnames(d2)   %in% removed.df$cell)], reason = "kept"))
removed.df = rbind(removed.df, data.frame(cell = colnames(d10)[which(!colnames(d10) %in% removed.df$cell)], reason = "kept"))
removed.df = rbind(removed.df, data.frame(cell = colnames(g10)[which(!colnames(g10) %in% removed.df$cell)], reason = "kept"))
removed.df$sample = toupper(reshape2::colsplit(removed.df$cell, "_", c('1', '2'))[,1])

pdf(paste0(out_dir, "tt_nuclei_kept.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar() + xlab("Sample") + ylab("Number of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()
pdf(paste0(out_dir, "tt_nuclei_kept_pct.pdf"), height = 4, width = 4, onefile = F)
print(ggplot(removed.df, aes(x = sample, fill = reason)) + geom_bar(position = "fill") + xlab("Sample") + ylab("Proportion of Nuclei") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("#023047", "#219ebc", "#8ecae6", "#ffb703"), name = NULL))
dev.off()

b10 = subset(b10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & pct.mt < 0.02)
d2  = subset(d2,  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & pct.mt < 0.02)
d10 = subset(d10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & pct.mt < 0.02)
g10 = subset(g10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & pct.mt < 0.02)

tt = merge(b10, c(d2, d10, g10))
tt = SCTransform(tt)
tt = RunPCA(tt, npcs = 30, verbose = FALSE)
tt = RunUMAP(tt, dims = 1:30, verbose = FALSE)
# tt = FindNeighbors(tt, dims = 1:30)
tt = FindNeighbors(tt, reduction = "umap", dims = 1:2)
tt = FindClusters(tt, resolution = 0.3, verbose = FALSE)
DimPlot(tt, label = TRUE) + coord_fixed() + theme_void()
Idents(tt) = tt$sample
DimPlot(tt, label = FALSE) + coord_fixed() + theme_void()

tt.deg=FindAllMarkers(tt, only.pos = T)
tt.deg$hgnc = gene_info$human[match(tt.deg$gene, gene_info$mzebra)]
tt.deg %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(tt, features = top10$gene) + NoLegend()

tt.rt = merge(b10, c(d2, d10))
tt.rt = SCTransform(tt.rt)
tt.rt = RunPCA(tt.rt, npcs = 30, verbose = FALSE)
tt.rt = RunUMAP(tt.rt, dims = 1:30, verbose = FALSE)
# tt.rt = FindNeighbors(tt.rt, dims = 1:30)
tt.rt = FindNeighbors(tt.rt, reduction = "umap", dims = 1:2)
tt.rt = FindClusters(tt.rt, resolution = 0.3, verbose = FALSE)
DimPlot(tt.rt, label = TRUE) + coord_fixed() + theme_void()
Idents(tt.rt) = tt.rt$sample
DimPlot(tt.rt, label = FALSE) + coord_fixed() + theme_void()

tt.ct = g10
tt.ct = SCTransform(tt.ct)
# tt = SCTransform(tt, vars.to.regress = "sample")
tt.ct = RunPCA(tt.ct, npcs = 30, verbose = FALSE)
tt.ct = RunUMAP(tt.ct, dims = 1:30, verbose = FALSE)
# tt.ct = FindNeighbors(tt.ct, dims = 1:30)
tt.ct = FindNeighbors(tt.ct, reduction = "umap", dims = 1:2)
tt.ct = FindClusters(tt.ct, resolution = 0.3, verbose = FALSE)
DimPlot(tt.ct, label = TRUE) + coord_fixed() + theme_void()

#*******************************************************************************
# Integration w/ Old Data ======================================================
#*******************************************************************************
tj.tmp.counts = Read10X(paste0(data_dir, "TJ/outs/filtered_feature_bc_matrix/"))
jpool.tmp.counts = Read10X(paste0(data_dir, "JPOOL/outs/filtered_feature_bc_matrix/"))
good_gene_info = gene_info[which(! is.na(gene_info$ens) & gene_info$ens != "" & gene_info$ens %in% rownames(tj.tmp.counts) ),]
good_gene_info = good_gene_info[which(! duplicated(good_gene_info$mzebra) ),]
good_gene_info = good_gene_info[which(! duplicated(good_gene_info$ens) ),]
good_ncbi = good_gene_info$mzebra
good_ens = good_gene_info$ens

tj.tmp.counts = tj.tmp.counts[good_ens,]
jpool.tmp.counts = jpool.tmp.counts[good_ens,]
rownames(jpool.tmp.counts) = rownames(tj.tmp.counts) = good_gene_info$mzebra

b10.counts.merge = b10.counts[good_ncbi,]
d2.counts.merge  = d2.counts[good_ncbi,]
d10.counts.merge = d10.counts[good_ncbi,]
g10.counts.merge = g10.counts[good_ncbi,]

tj.tmp = CreateSeuratObject(counts = tj.tmp.counts, project = "TJ")
jpool.tmp = CreateSeuratObject(counts = jpool.tmp.counts, project = "JPOOL")

b10.merge = CreateSeuratObject(counts = b10.counts.merge, project = "b10.merge")
d2.merge = CreateSeuratObject(counts = d2.counts.merge, project = "d2.merge")
d10.merge = CreateSeuratObject(counts = d10.counts.merge, project = "d10.merge")
g10.merge = CreateSeuratObject(counts = g10.counts.merge, project = "g10.merge")

b10.merge$sample = "B10"
d2.merge$sample = "D2"
d10.merge$sample = "D10"
g10.merge$sample = "G10"
tj.tmp$sample = "TJ"
jpool.tmp$sample = "JPOOL"

b10.merge = RenameCells(b10.merge, "b10")
d2.merge  = RenameCells(d2.merge, "d2")
d10.merge = RenameCells(d10.merge, "d10")
g10.merge = RenameCells(g10.merge, "g10")
tj.tmp = RenameCells(tj.tmp, "TJ")
jpool.tmp = RenameCells(jpool.tmp, "JPOOL")

b10.merge = subset(b10.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
d2.merge  = subset(d2.merge,  subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
d10.merge = subset(d10.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
g10.merge = subset(g10.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
tj.tmp = subset(tj.tmp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
jpool.tmp = subset(jpool.tmp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

tt.w.old = merge(b10.merge, c(d2.merge, d10.merge, g10.merge, tj.tmp, jpool.tmp))
tt.w.old = SCTransform(tt.w.old)
tt.w.old = RunPCA(tt.w.old, npcs = 30, verbose = FALSE)
tt.w.old = RunUMAP(tt.w.old, dims = 1:30, verbose = FALSE)
# tt.w.old = FindNeighbors(tt.w.old, dims = 1:30)
tt.w.old = FindNeighbors(tt.w.old, reduction = "umap", dims = 1:2)
tt.w.old = FindClusters(tt.w.old, resolution = 0.3, verbose = FALSE)
DimPlot(tt.w.old, label = TRUE) + coord_fixed() + theme_void()
Idents(tt.w.old) = tt.w.old$sample
DimPlot(tt.w.old, label = FALSE) + coord_fixed() + theme_void()

tt.w.old$new_cluster = "none"
tt.w.old$new_cluster[colnames(tt)] = as.vector(tt$seurat_clusters)
tt.w.old$new_cluster = factor(tt.w.old$new_cluster, c("none", levels(tt$seurat_clusters)))
Idents(tt.w.old) = tt.w.old$new_cluster
DimPlot(tt.w.old, label = T) + coord_fixed() + theme_void()

tt.w.old.tj$orig_cluster = "none"
tt.w.old.tj$orig_cluster[paste0(colnames(tj), "-1")] = as.vector(tj$seurat_clusters)
tt.w.old.tj$orig_cluster = factor(tt.w.old.tj$orig_cluster, levels = c("none", as.character(1:10)))
Idents(tt.w.old.tj) = tt.w.old.tj$orig_cluster
DimPlot(tt.w.old.tj, label = TRUE, cols = c("gray40", hue_pal()(10)[2:10])) + coord_fixed() + theme_void()

#*******************************************************************************
# Cluster Annotations ==========================================================
#*******************************************************************************
tal = as.data.frame(readxl::read_excel("~/research/tooth/data/mouse_and_human_100_unique_cell_type_degs.xlsx", skip = 1, sheet = 3))
tal = as.data.frame(readxl::read_excel("~/research/tooth/data/Specific and enriched Dental Genes.xlsx", skip = 0, sheet = 1))
tal.melt = reshape2::melt(as.matrix(tal))
colnames(tal.melt) = c("id", "cluster", "gene")
tal.melt$cluster = as.vector(tal.melt$cluster)
tal.melt$hgnc = toupper(tal.melt$gene)
tt.deg %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> top100
this.list = list(top100, tal.melt)
names(this.list) = c("TT", "MIM")
res = bigHeatmap(this.list, rm.self = T)

amb_markers = rbind(data.frame(hgnc = c("Sox2", "Bex1", "Six1", "Moxd1", "Sox11", "Etv4"), mzebra = c("sox2", "", "six1", "LOC101477797", "sox11", "etv4"), class = "SR"), 
                    data.frame(hgnc = c("Ascl5", "Moxd1", "Sox11", "Etv4"), mzebra = c("", "LOC101477797", "sox11", "etv4"), class = "Pre-Am"),
                    data.frame(hgnc = c("Ascl5", "Ell2", "Satb2", "Sox21", "Mafb"), mzebra = c("", "ell2", "satb2", "sox21", "mafb"), class = "Secr"),
                    data.frame(hgnc = c("Mafb", "Foxq1", "Foxo1", "Cdkn2b", "Runx2"), mzebra = c("mafb", "LOC101479104", "LOC101469585", "LOC101471656", "runx2"), class = "Mat"),
                    data.frame(hgnc = c("Foxq1", "Foxo1", "Cdkn2b", "Runx2", "Klf5", "Prmd1"), mzebra = c("LOC101479104", "LOC101469585", "LOC101471656", "runx2", "LOC101479207", ""), class = "Post"),
                    data.frame(hgnc = c("Enam", "Klk4", "Odam", "Gm17660"), mzebra = c("LOC101482019", "", "LOC101468182", ""), class = "Known"))
FeaturePlot(tt, features = amb_markers$mzebra[which(amb_markers$mzebra != "" & amb_markers$class == "SR")], order = T)

nestin_gene = "LOC101479744"
eng_gene = "LOC101484493"

tal = as.data.frame(readxl::read_excel("~/research/tooth/data/mouse_and_human_100_unique_cell_type_degs.xlsx", skip = 1, sheet = 1))
tj_deg_all = tt.deg
tj_deg_all$cluster = as.vector(tj_deg_all$cluster)
for (i in 1:ncol(tal)) {
        name = str_replace(colnames(tal)[i], "/", "_")
        name_ns = str_replace(name, " ", "_")
        print(name)
        
        cur_list = gene_info[match(toupper(tal[which(tal[,i] != "NA"),i]), gene_info$human),]
        # cur_list = gene_info[match(tal[which(tal[,i] != "NA"),i], gene_info$human),]
        cur_list = cur_list[which( !is.na(cur_list$mzebra) ),]
        cur_list = cur_list[which(!duplicated(cur_list$mzebra)),]
        
        res = markersInDEGUpDown(tj_deg_all, cur_list$mzebra)
        res_pct = markersInDEGUpDown(tj_deg_all, cur_list$mzebra, pct=T)[[2]]
        dot_res = myDotPlot(tt, cur_list$mzebra)
        dot_plot = dot_res[[1]]
        dot_bin = dot_res[[2]]
        dot_bin_weight = dot_res[[3]]
        
        # UMAP plot of Expression of Markers per Cell
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_exp_per_cell.png"), width = 500, height = 400)
        print(markerExpPerCell(tt, cur_list$mzebra) + ggtitle(paste0("Expression of ", name, " Markers per Cell in Cichlid Tooth")))
        dev.off()
        
        # Boxplot of Genes per Cell per Cluster - Corrected
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_avg_genes_per_cell_per_cluster_correct.png"), width = 600, height = 400)
        print(markerExpPerCellPerCluster(tt, cur_list$mzebra, correct = T)[[1]] + ggtitle(paste0("Average Number of ", name, " Markers per Cell per Cluster in Cichlid Tooth - Corrected")))
        dev.off()
        
        # Boxplot of Expression of Genes per Cell per Cluster - Corrected
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_exp_per_cell_per_cluster_correct.png"), width = 600, height = 400)
        print(markerExpPerCellPerCluster(tt, cur_list$mzebra, correct = T, n_markers = F)[[1]] + ggtitle(paste0("Average Expression of ", name, " Markers per Cell per Cluster in Cichlid Tooth - Corrected")))
        dev.off()
        
        # Barplot of Total Number of Cells per Cluster Expressing a Marker - Corrected
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_num_cell_per_cluster_correct.png"), width = 600, height = 400)
        print(markerCellPerCluster(tt, cur_list$mzebra) + ggtitle(paste("Normalized Total Number of Cells per Expressing Expressing", name, "Markers in Cichlid Tooth")))
        dev.off()
        
        # Heatmap of Expression per Cluster
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_heatmap.png"), width = 1200, height = 1200, res=100)
        print(markerHeatmap(tt, cur_list$mzebra, myslot = "data") + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
        dev.off()
        
        # DotPlot of Expression per Cluster
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_dotplot.png"), width = 800, height = 1600, res=100)
        print(dot_plot + ggtitle(paste0(name, " Markers in Cichlid Tooth")) + coord_flip())
        dev.off()
        
        # DotPlot BIN per Cluster
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_dotplot_bin.png"), width = 800, height = 400)
        print(dot_bin)
        dev.off()
        
        # DotPlot BIN per Cluster - Weighted
        Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_dotplot_bin_weight.png"), width = 800, height = 400)
        print(dot_bin_weight)
        dev.off()
        
        # DEGs in Markers
        if (! is.null(res[[2]])) {
                # Number of DEGs
                Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_deg.png"), width = 500, height = 400)
                print(res[[2]] + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
                dev.off()
                
                # LogFC
                Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_deg_logFC.png"), width = 500, height = 400)
                print(res[[3]] + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
                dev.off()
                
                # % DEGs in Markers
                Cairo(paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_deg_pct.png"), width = 500, height = 400)
                print(res_pct + ggtitle(paste0(name, " Markers in Cichlid Tooth")))
                dev.off()
        }
        
        # Write Genes Found
        # found_genes = res[[1]]
        # colnames(found_genes) = c("mzebra", "cluster", "avg_logFC", "isPos")
        # found_genes = left_join(found_genes, cur_list, by = "mzebra")
        # found_genes = found_genes[,c("hgnc", "mzebra", "cluster", "avg_logFC", "isPos")]
        # write.table(found_genes, paste0("~/research/tooth/results/tt_cluster_annot/tt_", name_ns, "_deg.tsv"), sep="\t", quote = F, row.names = F)
}
