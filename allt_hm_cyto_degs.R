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
createCytoBINsInGeneNeg = function(obj, gene) {
  gene_neg_cells = colnames(obj)[which(obj@assays$RNA@counts[gene,] > 0)]
  gene_obj = subset(obj, cells = gene_neg_cells)
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

#*******************************************************************************
# All Combined 1 ===============================================================
#*******************************************************************************

# # CM Celsr1+ DEGs
# allt = readRDS("allt.rds")
# res = createCytoBINsInGene(allt, "CELSR1")
# allt_celsr1_ctyo = res[[1]]
# allt_celsr1 = res[[2]]
# Idents(allt_celsr1) = allt_celsr1$bin
# celsr1_deg = FindAllMarkers(allt_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
# celsr1_deg$isSig = celsr1_deg$p_val_adj < 0.05
# celsr1_deg$cluster_isSig = paste0(celsr1_deg$cluster, "_", celsr1_deg$isSig)
# write.csv(celsr1_deg, "allt_celsr1_deg_big.csv")
# print("Done CM CElsr1+ DEG")
# 
# # CM Celsr1- DEGs
# allt_no_celsr1 = subset(allt, cells = colnames(allt)[which(! colnames(allt) %in% colnames(allt_celsr1))])
# allt_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(allt_no_celsr1@assays$RNA@counts))
# allt_no_celsr1$cyto = allt_no_celsr1_cyto$CytoTRACE
# allt_no_celsr1$bin <- allt_no_celsr1$cyto
# allt_no_celsr1$bin[which(allt_no_celsr1$cyto <= quantile(allt_no_celsr1$cyto, 0.33))] <- "Low"
# allt_no_celsr1$bin[which(allt_no_celsr1$cyto > quantile(allt_no_celsr1$cyto, 0.33) & allt_no_celsr1$cyto <= quantile(allt_no_celsr1$cyto, 0.66))] <- "Medium"
# allt_no_celsr1$bin[which(allt_no_celsr1$cyto > quantile(allt_no_celsr1$cyto, 0.66))] <- "High"
# Idents(allt_no_celsr1) = allt_no_celsr1$bin
# no_celsr1_deg = FindAllMarkers(allt_no_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
# no_celsr1_deg$isSig = no_celsr1_deg$p_val_adj < 0.05
# no_celsr1_deg$cluster_isSig = paste0(no_celsr1_deg$cluster, "_", no_celsr1_deg$isSig)
# write.csv(no_celsr1_deg, "allt_no_celsr1_deg_big.csv")
# print("Done CM CElsr1- DEG")
# 
# # HM Celsr1+ DEG
# hm = readRDS("~/scratch/d_tooth/data/hm.rds")
# res = createCytoBINsInGene(hm, "CELSR1")
# hm_celsr1_ctyo = res[[1]]
# hm_celsr1 = res[[2]]
# Idents(hm_celsr1) = hm_celsr1$bin
# hm_celsr1_deg = FindAllMarkers(hm_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
# hm_celsr1_deg$isSig = hm_celsr1_deg$p_val_adj < 0.05
# hm_celsr1_deg$cluster_isSig = paste0(hm_celsr1_deg$cluster, "_", hm_celsr1_deg$isSig)
# write.csv(hm_celsr1_deg, "hm_celsr1_deg_big.csv")
# print("Done HM CElsr1+ DEG")
# 
# # HM Celsr1- DEG
# hm_no_celsr1 = subset(hm, cells = colnames(hm)[which(! colnames(hm) %in% colnames(hm_celsr1))])
# hm_no_celsr1_cyto = CytoTRACE::CytoTRACE(as.matrix(hm_no_celsr1@assays$RNA@counts))
# hm_no_celsr1$cyto = hm_no_celsr1_cyto$CytoTRACE
# hm_no_celsr1$bin <- hm_no_celsr1$cyto
# hm_no_celsr1$bin[which(hm_no_celsr1$cyto <= quantile(hm_no_celsr1$cyto, 0.33))] <- "Low"
# hm_no_celsr1$bin[which(hm_no_celsr1$cyto > quantile(hm_no_celsr1$cyto, 0.33) & hm_no_celsr1$cyto <= quantile(hm_no_celsr1$cyto, 0.66))] <- "Medium"
# hm_no_celsr1$bin[which(hm_no_celsr1$cyto > quantile(hm_no_celsr1$cyto, 0.66))] <- "High"
# Idents(hm_no_celsr1) = hm_no_celsr1$bin
# hm_no_celsr1_deg = FindAllMarkers(hm_no_celsr1, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
# hm_no_celsr1_deg$isSig = hm_no_celsr1_deg$p_val_adj < 0.05
# hm_no_celsr1_deg$cluster_isSig = paste0(hm_no_celsr1_deg$cluster, "_", hm_no_celsr1_deg$isSig)
# write.csv(hm_no_celsr1_deg, "hm_no_celsr1_deg_big.csv")
# print("Done HM CElsr1- DEG")
# print("All Done")

#*******************************************************************************
# All Combined 2 ===============================================================
#*******************************************************************************
# Second time through this time calculating CytoTRACE separately.
# Read Data
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
data_folder = "~/scratch/d_tooth/data/"
allt_folder = "~/scratch/d_tooth/results/igor/allt/"
setwd(allt_folder)
tjh    = readRDS(paste0(data_folder, "tj_hgnc.rds"))
jawh   = readRDS(paste0(data_folder, "jaw_hgnc.rds"))
incsrh = readRDS(paste0(data_folder, "incsr_hgnc.rds"))
im     = readRDS(paste0(data_folder, "igor_incsr_molar.rds"))
hm     = readRDS(paste0(data_folder, "hm.rds"))
allt   = readRDS(paste0(allt_folder, "allt.rds"))

# Celsr1+ Cells
tjh_yc_res    = createCytoBINsInGene(tjh, "CELSR1")
jawh_yc_res   = createCytoBINsInGene(jawh, "CELSR1")
incsrh_yc_res = createCytoBINsInGene(incsrh, "CELSR1")
im_yc_res    = createCytoBINsInGene(im, "CELSR1")
hm_yc_res    = createCytoBINsInGene(hm, "CELSR1")

tjh_yc    = tjh_yc_res[[2]]
jawh_yc   = jawh_yc_res[[2]]
incsrh_yc = incsrh_yc_res[[2]]
im_yc     = im_yc_res[[2]]
hm_yc     = hm_yc_res[[2]]

allt_yc = subset(allt, cells = colnames(allt)[which(allt@assays$RNA@counts["CELSR1",] > 0)])
allt_yc$bin = "none"
allt_yc$bin[paste0("ct_", colnames(tjh_yc)[which(tjh_yc$bin       == "Low")])]    = "Low"
allt_yc$bin[paste0("ct_", colnames(tjh_yc)[which(tjh_yc$bin       == "Medium")])] = "Medium"
allt_yc$bin[paste0("ct_", colnames(tjh_yc)[which(tjh_yc$bin       == "High")])]   = "High"
allt_yc$bin[paste0("cj_", colnames(jawh_yc)[which(jawh_yc$bin     == "Low")])]    = "Low"
allt_yc$bin[paste0("cj_", colnames(jawh_yc)[which(jawh_yc$bin     == "Medium")])] = "Medium"
allt_yc$bin[paste0("cj_", colnames(jawh_yc)[which(jawh_yc$bin     == "High")])]   = "High"
allt_yc$bin[paste0("mi_", colnames(incsrh_yc)[which(incsrh_yc$bin == "Low")])]    = "Low"
allt_yc$bin[paste0("mi_", colnames(incsrh_yc)[which(incsrh_yc$bin == "Medium")])] = "Medium"
allt_yc$bin[paste0("mi_", colnames(incsrh_yc)[which(incsrh_yc$bin == "High")])]    = "High"
allt_yc$bin[paste0("mim_", colnames(im_yc)[which(im_yc$bin        == "Low")])]    = "Low"
allt_yc$bin[paste0("mim_", colnames(im_yc)[which(im_yc$bin        == "Medium")])] = "Medium"
allt_yc$bin[paste0("mim_", colnames(im_yc)[which(im_yc$bin        == "High")])]   = "High"
Idents(allt_yc) = allt_yc$bin
saveRDS(allt_yc, "allt_yc.rds")
saveRDS(hm_yc,   "hm_yc.rds")

print("Finding DEGs in Allt Celsr1+")
allt_yc_deg = FindAllMarkers(allt_yc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
allt_yc_deg$pct_dif = allt_yc_deg$pct.1 - allt_yc_deg$pct.2
allt_yc_deg$abs_pct_dif = allt_yc_deg$pct.1 - allt_yc_deg$pct.2
allt_yc_deg$isSig = allt_yc_deg$p_val_adj < 0.05
allt_yc_deg$cluster_isSig = paste0(allt_yc_deg$cluster, "_", allt_yc_deg$isSig)
write.csv(allt_yc_deg, "allt_yc_deg.csv")
print("Done")

print("Finding DEGs in HM Celsr1+")
hm_yc_deg = FindAllMarkers(hm_yc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
hm_yc_deg$pct_dif = hm_yc_deg$pct.1 - hm_yc_deg$pct.2
hm_yc_deg$abs_pct_dif = hm_yc_deg$pct.1 - hm_yc_deg$pct.2
hm_yc_deg$isSig = hm_yc_deg$p_val_adj < 0.05
hm_yc_deg$cluster_isSig = paste0(hm_yc_deg$cluster, "_", hm_yc_deg$isSig)
write.csv(hm_yc_deg, "hm_yc_deg.csv")
print("Done")

# Celsr1- Cells
tjh_nc_res    = createCytoBINsInGeneNeg(tjh, "CELSR1")
jawh_nc_res   = createCytoBINsInGeneNeg(jawh, "CELSR1")
incsrh_nc_res = createCytoBINsInGeneNeg(incsrh, "CELSR1")
imh_nc_res    = createCytoBINsInGeneNeg(im, "CELSR1")
hmh_nc_res    = createCytoBINsInGeneNeg(hm, "CELSR1")

tjh_nc    = tjh_nc_res[[2]]
jawh_nc   = jawh_nc_res[[2]]
incsrh_nc = incsrh_nc_res[[2]]
imh_nc    = imh_nc_res[[2]]
hm_nc     = hmh_nc_res[[2]]

allt_nc = subset(allt, cells = colnames(allt)[which(allt@assays$RNA@counts["CELSR1",] == 0)])
allt_nc$bin = "none"
allt_nc$bin[paste0("ct_", colnames(tjh_nc)[which(tjh_nc$bin       == "Low")])]    = "Low"
allt_nc$bin[paste0("ct_", colnames(tjh_nc)[which(tjh_nc$bin       == "Medium")])] = "Medium"
allt_nc$bin[paste0("ct_", colnames(tjh_nc)[which(tjh_nc$bin       == "High")])]   = "High"
allt_nc$bin[paste0("cj_", colnames(jawh_nc)[which(jawh_nc$bin     == "Low")])]    = "Low"
allt_nc$bin[paste0("cj_", colnames(jawh_nc)[which(jawh_nc$bin     == "Medium")])] = "Medium"
allt_nc$bin[paste0("cj_", colnames(jawh_nc)[which(jawh_nc$bin     == "High")])]   = "High"
allt_nc$bin[paste0("mi_", colnames(incsrh_nc)[which(incsrh_nc$bin == "Low")])]    = "Low"
allt_nc$bin[paste0("mi_", colnames(incsrh_nc)[which(incsrh_nc$bin == "Medium")])] = "Medium"
allt_nc$bin[paste0("mi_", colnames(incsrh_nc)[which(incsrh_nc$bin == "High")])]   = "High"
allt_nc$bin[paste0("mim_", colnames(im_nc)[which(im_nc$bin        == "Low")])]    = "Low"
allt_nc$bin[paste0("mim_", colnames(im_nc)[which(im_nc$bin        == "Medium")])] = "Medium"
allt_nc$bin[paste0("mim_", colnames(im_nc)[which(im_nc$bin        == "High")])]   = "High"
Idents(allt_nc) = allt_nc$bin
saveRDS(allt_nc, "allt_nc.rds")
saveRDS(hm_nc,   "hm_nc.rds")

print("Finding DEGs in Allt Celsr1-")
allt_nc_deg = FindAllMarkers(allt_nc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
allt_nc_deg$pct_dif = allt_nc_deg$pct.1 - allt_nc_deg$pct.2
allt_nc_deg$abs_pct_dif = allt_nc_deg$pct.1 - allt_nc_deg$pct.2
allt_nc_deg$isSig = allt_nc_deg$p_val_adj < 0.05
allt_nc_deg$cluster_isSig = paste0(allt_nc_deg$cluster, "_", allt_nc_deg$isSig)
write.csv(allt_nc_deg, "allt_nc_deg.csv")
print("Done")

print("Finding DEGs in HM Celsr1-")
hm_nc_deg = FindAllMarkers(hm_nc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
hm_nc_deg$pct_dif = hm_nc_deg$pct.1 - hm_nc_deg$pct.2
hm_nc_deg$abs_pct_dif = hm_nc_deg$pct.1 - hm_nc_deg$pct.2
hm_nc_deg$isSig = hm_nc_deg$p_val_adj < 0.05
hm_nc_deg$cluster_isSig = paste0(hm_nc_deg$cluster, "_", hm_nc_deg$isSig)
write.csv(hm_nc_deg, "hm_nc_deg.csv")
print("Done")
print("All Done")