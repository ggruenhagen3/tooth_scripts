library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("CytoTRACE")
library("ggplot2")
library("RColorBrewer")
library("tidyverse")

combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/combined.Rds")
epi      <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/epi.rds")
jpool    <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/jpool.rds")
tj       <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")

epi$cond[is.na(epi$cond)] <- "INJR"
combined$dataset <- "mesenchyme"
epi$dataset      <- "epithelium"
jpool$dataset    <- "jaw"
tj$dataset       <- "tooth"
all_obj <- list(combined, epi, jpool, tj)

# Top 50 Genes (Least Differentiated)
# # all_top_genes <- lapply(0:length(all_obj), function(x) c())
# for (i in 1:length(all_obj)) {
#   obj <- all_obj[[i]]
#   obj_str <- unique(obj$dataset)
#   results <- CytoTRACE(as.matrix(obj@assays$RNA@counts))
#   genes <- rownames(obj@assays$RNA@counts)
#   
#   this_result <- data.frame()
#   for (gene in genes) {
#     gene_cells <- names(obj@assays$RNA@counts[gene,which(obj@assays$RNA@counts[gene,] != 0)])
#     if (length(gene_cells) > 30) {
#       gene_p <- "NA"
#       try(gene_p <- t.test(results$CytoTRACE[gene_cells], results$CytoTRACE)$p.value, silent = TRUE)
#       this_result <- rbind(this_result, t(c( obj_str, gene, length(gene_cells), mean(results$CytoTRACE[gene_cells]), gene_p  )))
#     }
#   } # end gene for
#   colnames(this_result) <- c("dataset", "gene", "num_cells", "mean", "p")
#   this_result <- this_result[order(this_result$mean),]
#   this_result <- this_result[1:50,]
#   # write.table(this_result, paste("/nv/hp10/ggruenhagen3/scratch/d_tooth/", obj_str, "_top_50.tsv", sep=""), ssep = "\t", quote = FALSE)
# } # end object for

# Celsr1 by Cluster
all_celsr1 <- c()
all_results <- c()
all_results_cluster <- data.frame()
# vlndata <- data.frame()
for (obj in all_obj) {
  obj_str <- unique(obj$dataset)
  results <- CytoTRACE(as.matrix(obj@assays$RNA@counts))
  # saveRDS(results, paste("/nv/hp10/ggruenhagen3/scratch/d_tooth/results/", obj_str, "_cyto.rds", sep=""))
  
  all_results <- c(all_results, results$CytoTRACE)

  if (obj_str == "mesenchyme" || obj_str == "epithelium") {
    expr <- FetchData(object = obj, vars = "Celsr1", slot = "counts")
  } else {
    expr <- FetchData(object = obj, vars = "celsr1a", slot = "counts")
  }
  celsr1_cells <- colnames(obj@assays$RNA@counts[, which(x = expr > 0)])

  conds <- unique(obj$cond)
  num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
  for (cond in conds) {
    Idents(obj) <- "cond"
    this_cond_cells <- WhichCells(obj, idents = c(cond))
    Idents(obj) <- "seurat_clusters"
    vlndata <- data.frame()
    for (cluster in 0:num_clusters) {
      this_cluster_cells <- WhichCells(obj, idents = c(cluster))
      this_cells <- this_cond_cells[which(this_cond_cells %in% this_cluster_cells)]
      celsr1_cluster <- celsr1_cells[which(celsr1_cells %in% this_cells)]
      celsr1_cluster_p <- "NA"
      try(celsr1_cluster_p <- t.test(results$CytoTRACE[celsr1_cluster], results$CytoTRACE[this_cells])$p.value, silent = TRUE)
      
      celsr1_row <- data.frame( rep(obj_str, length(celsr1_cluster)), rep(cond, length(celsr1_cluster)), rep(cluster, length(celsr1_cluster)), rep("Celsr1", length(celsr1_cluster)), results$CytoTRACE[celsr1_cluster] )
      colnames(celsr1_row) <- c("dataset", "cond", "cluster", "isCelsr1", "cyto")
      vlndata <- rbind( vlndata, celsr1_row )
      all_row <- data.frame( rep(obj_str, length(this_cells)), rep(cond, length(this_cells)), rep(cluster, length(this_cells)), rep("All", length(this_cells)), results$CytoTRACE[this_cells] )
      colnames(all_row) <- c("dataset", "cond", "cluster", "isCelsr1", "cyto")
      vlndata <- rbind( vlndata, all_row )
      
      all_results_cluster <- rbind(all_results_cluster, t(c( obj_str, cond, cluster, mean(results$CytoTRACE[celsr1_cluster]), mean(results$CytoTRACE[this_cells]), celsr1_cluster_p )))
    } # end cluster for
    
    colnames(vlndata) <- c("dataset", "cond", "cluster", "isCelsr1", "cyto")
    vlndata$cluster <- as.factor(vlndata$cluster)
    filename <-  paste("/nv/hp10/ggruenhagen3/scratch/d_tooth/results/vln_", obj_str, "_", cond, ".png", sep="")
    png(filename, width = 5000, height = 1600, res = 200)
    # p <- ggplot(vlndata, aes(x=cluster, y=cyto, fill=isCelsr1)) + geom_violin(position="dodge") + geom_boxplot(width=0.3, position="dodge") + ggtitle(paste(obj_str, "-", cond, ": Celsr1 vs All")) + theme_classic()
    p <- ggplot(vlndata, aes(x=cluster, y=cyto, fill=isCelsr1, color=isCelsr1)) + geom_boxplot(position=position_dodge2(), alpha=0.5) + geom_jitter(shape=16, position=position_jitterdodge(), aes(alpha = isCelsr1)) + ggtitle(paste(obj_str, "-", cond, ": Celsr1 vs All")) + theme_classic() + scale_alpha_discrete(range = c(0.8, 0.3))
    print(p)
    dev.off()
    system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/laptop_backup/d_tooth/results/", sep=""))
    
    
  } # end condition for
  # # plotCytoTRACE(results, emb = obj@reductions$umap@cell.embeddings[,1:2], outputDir = "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/")
}
# colnames(vlndata) <- c("dataset", "cond", "cluster", "isCelsr1", "cyto")
# vlndata$cluster <- as.factor(vlndata$cluster)
# filename <-  "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/test.png"
# png(filename, width = 1800, height = 1200, res = 150)
# p <- ggplot(vlndata[which(vlndata$dataset == "epithelium" & vlndata$cond == "INJR"),], aes(x=cluster, y=cyto, fill=isCelsr1)) + geom_violin(position=position_dodge(1)) + geom_boxplot(width=0.1)
# print(p)
# dev.off()
# system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/laptop_backup/d_tooth/results/", sep=""))
#   
colnames(all_results_cluster) <- c("dataset", "cond", "cluster", "mean_of_celsr1_in_cluster", "mean_of_cluster", "p")
all_results_cluster$p_sig_greater <- as.vector(all_results_cluster$mean_of_celsr1_in_cluster) > as.vector(all_results_cluster$mean_of_cluster) & as.vector(all_results_cluster$p) < 0.05
write.table(all_results_cluster, "/nv/hp10/ggruenhagen3/scratch/d_tooth/results/celsr1_cluster.tsv", sep="\t", quote = FALSE, row.names = FALSE)

#################
# Local PC Code #
#################
mes_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/mesenchyme_cyto.rds")
epi_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/epithelium_cyto.rds")
jaw_cyto <- readRDS("C:/Users/miles/Downloads/d_tooth/data/jaw_cyto.rds")
tj_cyto  <- readRDS("C:/Users/miles/Downloads/d_tooth/data/tooth_cyto.rds")
epi$cond[is.na(epi$cond)] <- "INJR"
combined$cyto <- mes_cyto$CytoTRACE
epi$cyto <- epi_cyto$CytoTRACE
jaw$cyto <- jaw_cyto$CytoTRACE
tj$cyto = tj_cyto$CytoTRACE
all_cyto <- list(mes_cyto, epi_cyto, jaw_cyto, tj_cyto)
all_obj  <- list(mes, epi, jpool, tj)
all_names <- c("mes", "epi", "jaw", "tooth")
# all_cyto <- list(mes_cyto, epi_cyto)
# all_obj  <- list(mes, epi)
# all_names <- c("mes", "epi")

rna_path <- "C:/Users/miles/Downloads/d_tooth/"

all_top_50 <- data.frame()
all_top_200 <- data.frame()
all_bot_50 <- data.frame()
all_bot_200 <- data.frame()
for (i in 1:length(all_cyto)) {
  cyto <- all_cyto[[i]]
  name <- all_names[i]
  obj <- all_obj[[i]]
  
  # # Add the cyto trace data to the object
  # matrix2 <- matrix(as.numeric(as.vector(cyto$CytoTRACE)), nrow = 1, ncol = length(cyto$CytoTRACE), dimnames = list("cyto", names(cyto$CytoTRACE)) )
  # obj[["cyto"]] <- CreateAssayObject(counts = matrix2)
  # 
  # # Save Plot
  # png(paste(rna_path, "/results/", name, "_cyto.png", sep=""), width = 4800, height = 2400, unit="px", res=300)
  # p <-   FeaturePlot(obj, features = "cyto", reduction = "umap", cols = rev(brewer.pal(11,"Spectral")), split.by = "cond", pt.size = 2, label=TRUE, order = TRUE)
  # print(p)
  # dev.off()
  
  top_gene_cyto <- c()
  for (i in 1:50) {
    gene <- names(cyto$cytoGenes[i])
    expr <- FetchData(object = obj, vars = gene, slot = "counts")
    gene_cells <- colnames(obj@assays$RNA@counts[, which(x = expr > 0)])
    top_gene_cyto <- c(top_gene_cyto, mean(cyto$CytoTRACE[gene_cells]))
  }
  
  bot_gene_cyto <- c()
  for (i in (length(cyto$cytoGenes)-200+1):length(cyto$cytoGenes)) {
    gene <- names(cyto$cytoGenes[i])
    expr <- FetchData(object = obj, vars = gene, slot = "counts")
    gene_cells <- colnames(obj@assays$RNA@counts[, which(x = expr > 0)])
    bot_gene_cyto <- c(bot_gene_cyto, mean(cyto$CytoTRACE[gene_cells]))
  }
  
  newRow_top_50 <- data.frame(rep(name, 50), names(cyto$cytoGenes[1:50]), top_gene_cyto[1:50], cyto$cytoGenes[1:50], names(cyto$gcsGenes[1:50]), cyto$gcsGenes[1:50])
  colnames(newRow_top_50) <- c("dataset", "top50_cytoTRACE_genes", "avg_ctyoTRACE_score", "Pearson_Correlation", "top50_gcs_genes", "Pearson_Correlation")
  all_top_50 <- rbind(all_top_50, newRow_top_50)
  
  newRow <- data.frame(rep(name, 200), names(cyto$cytoGenes)[1:200], top_gene_cyto[1:200], cyto$cytoGenes[1:200], names(cyto$gcsGenes)[1:200], cyto$gcsGenes[1:200])
  colnames(newRow) <- c("dataset", "top200_cytoTRACE_genes", "avg_ctyoTRACE_score", "Pearson_Correlation", "top200_gcs_genes", "Pearson_Correlation")
  all_top_200 <- rbind(all_top_200, newRow)
  
  newRow_bot_50 <- data.frame(rep(name, 50), tail(names(cyto$cytoGenes), n=50), tail(bot_gene_cyto, n=50), tail(cyto$cytoGenes, n=50), tail(names(cyto$gcsGenes), n=50), tail(cyto$gcsGenes, n=50))
  colnames(newRow_bot_50) <- c("dataset", "bot50_cytoTRACE_genes", "avg_ctyoTRACE_score", "Pearson_Correlation", "bot50_gcs_genes", "Pearson_Correlation")
  all_bot_50 <- rbind(all_bot_50, newRow_bot_50)
  
  newRow <- data.frame(rep(name, 200), tail(names(cyto$cytoGenes), n=200), tail(bot_gene_cyto, n=200), tail(cyto$cytoGenes, n=200), tail(names(cyto$gcsGenes), n=200), tail(cyto$gcsGenes, n=200))
  colnames(newRow) <- c("dataset", "bot200_cytoTRACE_genes", "avg_ctyoTRACE_score", "Pearson_Correlation", "bot200_gcs_genes", "Pearson_Correlation")
  all_bot_200 <- rbind(all_bot_200, newRow)
  
  # Plot top 10 and bot 10
  colnames(newRow_top_50) <- c("dataset", "gene", "Avg_CytoTRACE_Score", "Abs_Pearson_Correlation", "gcs", "gcs_Pearson_Correlation")
  colnames(newRow_bot_50) <- c("dataset", "gene", "Avg_CytoTRACE_Score", "Abs_Pearson_Correlation", "gcs", "gcs_Pearson_Correlation")
  newRow_bot_50$Abs_Pearson_Correlation <- abs(newRow_bot_50$Abs_Pearson_Correlation)
  df <- rbind(tail(newRow_bot_50, 10), head(newRow_top_50, 10))
  df$isTop <- c(rep(FALSE, 10), rep(TRUE, 10))
  df2 <- df %>% pivot_longer(c("Abs_Pearson_Correlation", "Avg_CytoTRACE_Score"), names_to = "variable", values_to = "value")
  # ggplot(df, aes(x=gene, y=Abs_Pearson_Correlation, fill=isTop)) + geom_bar(stat="identity") + coord_flip() + theme_classic()
  # png(paste(rna_path, "/results/", name, "_cyto_top_bot_10.png", sep=""), width = 900, height = 1800, unit="px", res=150)
  p <- ggplot(df2, aes(x=gene, y=value, color=variable, fill=isTop)) + geom_bar(stat="identity", position=position_dodge2(), size = 0.8, alpha=0.6) + ylim(0,1) + coord_flip() + theme_classic() + scale_fill_manual(values = c("blue", "red")) + scale_color_manual(values = c("gold","black")) + scale_x_discrete(expand=c(0,0))
  print(p)
  # dev.off()
}
write.table(all_top_50, paste(rna_path, "/results/cyto_top_50.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)
write.table(all_top_200, paste(rna_path, "/results/cyto_top_200.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)
write.table(all_bot_50, paste(rna_path, "/results/cyto_bot_50.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)
write.table(all_bot_200, paste(rna_path, "/results/cyto_bot_200.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

# Boxplot of CytoTRACE score per cell type
Idents(mes) <- "seurat_clusters"
new.cluster.ids <- c("Pre-odontoblast", "Pre-odontoblast", "Pre-odontoblast", "MSC", "MSC", "TAC", "Perivascular","MSC", "Perivascular", "Immune", "Immune", "TAC", "Immune", "Immune", "Immune", "Immune", "Neural Crest", "Immune", "Immune", "Immune", "Epithelial")
names(new.cluster.ids) <- levels(mes)
mes <- RenameIdents(mes, new.cluster.ids)
mes_cols <- c("pink", "cadetblue1", "chartreuse", "tomato", "gold", "purple", "orange")
cell_type_df <- cbind(as.numeric(as.vector(mes_cyto$CytoTRACE)), as.vector(Idents(mes)))
colnames(cell_type_df) <- c("CytoTRACE", "Cell_Type")
cell_type_df <- as.data.frame(cell_type_df)
cell_type_df$CytoTRACE <- as.numeric(as.vector(cell_type_df$CytoTRACE))
cell_type_df$Cell_Type <- factor(cell_type_df$Cell_Type, levels=levels(Idents(mes)))
rownames(cell_type_df) <- NULL
ggplot(cell_type_df, aes(x=Cell_Type, y=CytoTRACE, fill = Cell_Type, color=Cell_Type)) + geom_boxplot(alpha=0.8) + ylim(0,1) + scale_x_discrete(limits=rev(levels(Idents(mes)))) + coord_flip() + theme_classic() + scale_fill_manual(values = mes_cols) + scale_color_manual(values = mes_cols)

cell_type_df$isLow   = cell_type_df$CytoTRACE >= 0 & cell_type_df$CytoTRACE < 0.33
cell_type_df$isMed   = cell_type_df$CytoTRACE >= 0.33 & cell_type_df$CytoTRACE < 0.66
cell_type_df$isHigh  = cell_type_df$CytoTRACE >= 0.66

obj=incsr_2
cell_type_df_2 = data.frame()
for (cell_type in levels(Idents(obj))) {
  print(cell_type)
  rows = data.frame(rep(cell_type, 3), c("Low", "Medium", "High"),
                    c(length(cell_type_df$isLow[which(cell_type_df$Cell_Type == cell_type & cell_type_df$isLow == T)]), 
                      length(cell_type_df$isMed[which(cell_type_df$Cell_Type == cell_type & cell_type_df$isMed == T)]), 
                      length(cell_type_df$isHigh[which(cell_type_df$Cell_Type == cell_type & cell_type_df$isHigh == T)])))
  colnames(rows) = c("Cell_Type", "CytoBIN", "value")
  rows$pct = c(rows$value[1]/sum(rows$value), rows$value[2]/sum(rows$value), rows$value[3]/sum(rows$value))
  rows$pct = rows$pct*100
  cell_type_df_2 = rbind(cell_type_df_2, rows)
}
# colnames(cell_type_df_2) = c("Cell_Type", "CytoBIN", "value", "pct")
cell_type_df_2$CytoBIN = factor(cell_type_df_2$CytoBIN, levels = c("High", "Medium", "Low"))
ggplot(cell_type_df_2, aes(x=Cell_Type, y=value, fill=CytoBIN)) + geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=brewer.pal(11,"Spectral")[c(2,6,10)]) + ggtitle("Number of Cells in CytoBINs per Cluster - Mes Reclustered")
ggplot(cell_type_df_2, aes(x=Cell_Type, y=pct, fill=CytoBIN)) + geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=brewer.pal(11,"Spectral")[c(2,6,10)]) + ggtitle("Percent of Cells in CytoBINs per Cluster - Mes Reclustered")

# Find a TF present in tooth celsr1 cells but not jaw celsr 1 cells
expr <- FetchData(object = tj, vars = "celsr1a", slot = "counts")
tj_celsr_cells <- colnames(tj@assays$RNA@counts[, which(x = expr > 0)])
tj_celsr <- tj[,tj_celsr_cells]
expr <- FetchData(object = jpool, vars = "celsr1a", slot = "counts")
jpool_celsr_cells <- colnames(jpool@assays$RNA@counts[, which(x = expr > 0)])
jpool_celsr <- jpool[,jpool_celsr_cells]
genes_in_tj_celsr    <- rownames(tj_celsr)[which(rowSums(tj_celsr@assays$RNA@counts[,]) != 0)]
genes_in_jpool_celsr <- rownames(jpool_celsr)[which(rowSums(jpool_celsr@assays$RNA@counts[,]) != 0)]
genes_just_tj_celsr <- genes_in_tj_celsr[which(! genes_in_tj_celsr %in% genes_in_jpool_celsr)]
genes_just_jpool_celsr <- genes_in_jpool_celsr[which(! genes_in_jpool_celsr %in% genes_in_tj_celsr)]

# ARE THERE ANY GENES EXPRESSED IN ALL (OR ALMOST ALL) CELSR1 CELLS THAT AREN'T 
# EXPRESSED IN OTHER CELLS?
common_celsr <- genes_in_jpool_celsr[which(genes_in_jpool_celsr %in% genes_in_tj_celsr)]
tj_non_celsr    <- tj[,colnames(tj)[which(! colnames(tj) %in% tj_celsr_cells)]]
jpool_non_celsr <- jpool[,colnames(jpool)[which(! colnames(jpool) %in% jpool_celsr_cells)]]
genes_tj_non_celsr    <- rownames(tj_non_celsr)[which(rowSums(tj_non_celsr@assays$RNA@counts[,]) != 0)]
genes_jpool_non_celsr <- rownames(jpool_non_celsr)[which(rowSums(jpool_non_celsr@assays$RNA@counts[,]) != 0)]
unique_celsr <- common_celsr[which( (! common_celsr %in% genes_tj_non_celsr) & (! common_celsr %in% genes_jpool_non_celsr) )]
# unique_celsr == "ENSMZEG00005016016", "nucb2b", "CPA1", "SELE (1 of many).1", "ENSMZEG00005006843", "ENSMZEG00005013205","celsr1a", "fam92a1", "ENSMZEG00005012661", "ENSMZEG00005020998", "ENSMZEG00005027674", "ENSMZEG00005026148"

# ARE THERE GENES OR COMBOS OF GENES WE COULD USE TO DIFFERENTIATE BETWEEN TOOTH AND JAW CELLS?
genes_tj    <- rownames(tj)[which(rowSums(tj@assays$RNA@counts[,]) != 0)]
genes_jpool <- rownames(jpool)[which(rowSums(jpool@assays$RNA@counts[,]) != 0)]
unique_tj    <- genes_tj[which(! genes_tj %in% genes_jpool)]
unique_jpool <- genes_jpool[which(! genes_jpool %in% genes_tj)]

threshold = 0.8
tj_cells <- length(colnames(tj))
jpool_cells <- length(colnames(jpool))
tj.pct.gene <- data.frame()
jpool.pct.gene <- data.frame()
for (gene in unique_tj) {
  # Find out of the unique genes are expressed broadly in the sample
  gene_cells <- names(which(tj@assays$RNA@counts[gene,] != 0))
  tj.pct.gene <- rbind(tj.pct.gene, t(c(gene, length(gene_cells)/tj_cells)))
}
for (gene in unique_jpool) {
  # Find out of the unique genes are expressed broadly in the sample
  gene_cells <- names(which(jpool@assays$RNA@counts[gene,] != 0))
  jpool.pct.gene <- rbind(jpool.pct.gene, t(c(gene, length(gene_cells)/jpool_cells)))
}
tj.pct.gene$V2    <- as.numeric(as.vector(tj.pct.gene$V2))
jpool.pct.gene$V2 <- as.numeric(as.vector(jpool.pct.gene$V2))

tj_jaw_celsr <- merge(tj_celsr, jpool_celsr, merge.data = TRUE)
tj_jaw_celsr <- FindVariableFeatures(object = tj_jaw_celsr, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
tj_jaw_celsr <- ScaleData(object = tj_jaw_celsr, vars.to.regress = NULL)
tj_jaw_celsr <- RunPCA(tj_jaw_celsr, npcs = 30, verbose = FALSE)
tj_jaw_celsr <- RunUMAP(tj_jaw_celsr, reduction = "pca", dims = 1:12)
tj_jaw_celsr <- FindNeighbors(tj_jaw_celsr, reduction = "umap", dims = 1:2)
tj_jaw_celsr <- FindClusters(tj_jaw_celsr, resolution = 0.30)
DimPlot(tj_jaw_celsr, reduction = "umap", split.by = "cond", label = TRUE)
FeaturePlot(tj_jaw_celsr, feature = "krt15", reduction = "umap", split.by = "cond", label = TRUE)

m_prog = c("Cdh6", "Lrp11", "Vat1l", "Fam19a4", "Hey2", "Cpne5", "Mki67")
m_stem = c("Lgr5", "Lrig1", "Sox2", "Sfrp5", "Grp", "Rhoc", "Disc1", "Fez1")
m_matr = c("Gm17660", "Slc5a8", "Ptpn22", "Klk4", "Gpr155", "Slc34a2")

cell_types = list()
cell_types[["prog"]] = m_prog
cell_types[["stem"]] = m_stem
cell_types[["matr"]] = m_matr

cytogene_rank = data.frame()
for (i in 1:length(cell_types)) {
  name = names(cell_types)[i]
  genes = cell_types[[name]]
  newRow = data.frame(rep(name, length(genes)), genes, which(names(results$cytoGenes) %in% genes), results$cytoGenes[which(names(results$cytoGenes) %in% genes)])
  cytogene_rank = rbind(cytogene_rank, newRow)
}
colnames(cytogene_rank) = c("Cell_Type", "Gene", "Cyto_Rank", "Correlation")
ggplot(cytogene_rank, aes(Cell_Type, Cyto_Rank, fill=Cell_Type, color=Cell_Type)) + geom_boxplot(alpha=0.8) + ggtitle("Igor Incisor Epithelium")
ggplot(cytogene_rank, aes(Cell_Type, Correlation, fill=Cell_Type, color=Cell_Type)) + geom_boxplot(alpha=0.8) + ggtitle("Igor Incisor Epithelium")
