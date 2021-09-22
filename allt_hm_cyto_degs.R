# Set Folder
setwd("~/scratch/d_tooth/results/igor/allt/")

# Load Libraries
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))

# Helper Functions
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
}

# CM Celsr1+ DEGs
allt = readRDS("~/research/tooth/data/allt.rds")
res = createCytoBINsInGene(allt, "CELSR1")
allt_celsr1_ctyo = res[[1]]
allt_celsr1 = res[[2]]
Idents(allt_celsr1) = allt_celsr1$bin
celsr1_deg = FindAllMarkers(allt_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
celsr1_deg$isSig = celsr1_deg$p_val_adj < 0.05
celsr1_deg$cluster_isSig = paste0(celsr1_deg$cluster, "_", celsr1_deg$isSig)
write.csv(celsr1_deg, "allt_celsr1_deg_big.csv")
print("Done CM CElsr1+ DEG")

# CM Celsr1- DEGs
allt_no_celsr1 = subset(allt, cells = colnames(allt)[which(! colnames(allt) %in% colnames(allt_celsr1))])
allt_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(allt_no_celsr1@assays$RNA@counts))
allt_no_celsr1$cyto = allt_no_celsr1_cyto$CytoTRACE
allt_no_celsr1$bin <- allt_no_celsr1$cyto
allt_no_celsr1$bin[which(allt_no_celsr1$cyto <= quantile(allt_no_celsr1$cyto, 0.33))] <- "Low"
allt_no_celsr1$bin[which(allt_no_celsr1$cyto > quantile(allt_no_celsr1$cyto, 0.33) & allt_no_celsr1$cyto <= quantile(allt_no_celsr1$cyto, 0.66))] <- "Medium"
allt_no_celsr1$bin[which(allt_no_celsr1$cyto > quantile(allt_no_celsr1$cyto, 0.66))] <- "High"
Idents(allt_no_celsr1) = allt_no_celsr1$bin
no_celsr1_deg = FindAllMarkers(allt_no_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
no_celsr1_deg$isSig = no_celsr1_deg$p_val_adj < 0.05
no_celsr1_deg$cluster_isSig = paste0(no_celsr1_deg$cluster, "_", no_celsr1_deg$isSig)
write.csv(no_celsr1_deg, "allt_no_celsr1_deg_big.csv")
print("Done CM CElsr1- DEG")

# HM Celsr1+ DEG
hm = readRDS("~/scratch/d_tooth/data/hm.rds")
res = createCytoBINsInGene(hm, "CELSR1")
hm_celsr1_ctyo = res[[1]]
hm_celsr1 = res[[2]]
Idents(hm_celsr1) = hm_celsr1$bin
hm_celsr1_deg = FindAllMarkers(hm_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
hm_celsr1_deg$isSig = hm_celsr1_deg$p_val_adj < 0.05
hm_celsr1_deg$cluster_isSig = paste0(hm_celsr1_deg$cluster, "_", hm_celsr1_deg$isSig)
write.csv(hm_celsr1_deg, "hm_celsr1_deg_big.csv")
print("Done HM CElsr1+ DEG")

# HM Celsr1- DEG
hm_no_celsr1 = subset(hm, cells = colnames(hm)[which(! colnames(hm) %in% colnames(hm_celsr1))])
hm_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(hm_no_celsr1@assays$RNA@counts))
hm_no_celsr1$cyto = hm_no_celsr1_cyto$CytoTRACE
hm_no_celsr1$bin <- hm_no_celsr1$cyto
hm_no_celsr1$bin[which(hm_no_celsr1$cyto <= quantile(hm_no_celsr1$cyto, 0.33))] <- "Low"
hm_no_celsr1$bin[which(hm_no_celsr1$cyto > quantile(hm_no_celsr1$cyto, 0.33) & hm_no_celsr1$cyto <= quantile(hm_no_celsr1$cyto, 0.66))] <- "Medium"
hm_no_celsr1$bin[which(hm_no_celsr1$cyto > quantile(hm_no_celsr1$cyto, 0.66))] <- "High"
Idents(hm_no_celsr1) = hm_no_celsr1$bin
hm_no_celsr1_deg = FindAllMarkers(hm_no_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
hm_no_celsr1_deg$isSig = hm_no_celsr1_deg$p_val_adj < 0.05
hm_no_celsr1_deg$cluster_isSig = paste0(hm_no_celsr1_deg$cluster, "_", hm_no_celsr1_deg$isSig)
write.csv(hm_no_celsr1_deg, "hm_no_celsr1_deg_big.csv")
print("Done HM CElsr1- DEG")
print("All Done")