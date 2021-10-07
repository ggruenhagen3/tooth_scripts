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
  gene_neg_cells = colnames(obj)[which(obj@assays$RNA@counts[gene,] == 0)]
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
# Body =========================================================================
#*******************************************************************************
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
data_folder = "~/scratch/d_tooth/data/"
allt_folder = "~/scratch/d_tooth/results/igor/allt/"
setwd(allt_folder)
tjh    = readRDS(paste0(data_folder, "tj_hgnc.rds"))
jawh   = readRDS(paste0(data_folder, "jaw_hgnc.rds"))
incsrh = readRDS(paste0(data_folder, "incsr_hgnc.rds"))
im     = readRDS(paste0(data_folder, "igor_incsr_molar.rds"))

imi = subset(im, cells = colnames(im)[which(im$anatomy == "Incisor")])
imm = subset(im, cells = colnames(im)[which(im$anatomy == "Molar")])

regen = merge(tjh, list(jawh, incsrh, imi), add.cell.ids = c("ct", "cj", "mi", 'imi'))

# Celsr1+ Cells
tjh_yc_res    = createCytoBINsInGene(tjh, "CELSR1")
jawh_yc_res   = createCytoBINsInGene(jawh, "CELSR1")
incsrh_yc_res = createCytoBINsInGene(incsrh, "CELSR1")
imi_yc_res    = createCytoBINsInGene(imi, "CELSR1")
imm_yc_res    = createCytoBINsInGene(imm, "CELSR1")

tjh_yc    = tjh_yc_res[[2]]
jawh_yc   = jawh_yc_res[[2]]
incsrh_yc = incsrh_yc_res[[2]]
imi_yc    = imi_yc_res[[2]]
imm_yc    = imm_yc_res[[2]]

regen_yc = subset(allt, cells = colnames(regen)[which(regen@assays$RNA@counts["CELSR1",] > 0)])
regen_yc$bin = "none"
regen_yc$bin[paste0("ct_", colnames(tjh_yc)[which(tjh_yc$bin       == "Low")])]    = "Low"
regen_yc$bin[paste0("ct_", colnames(tjh_yc)[which(tjh_yc$bin       == "Medium")])] = "Medium"
regen_yc$bin[paste0("ct_", colnames(tjh_yc)[which(tjh_yc$bin       == "High")])]   = "High"
regen_yc$bin[paste0("cj_", colnames(jawh_yc)[which(jawh_yc$bin     == "Low")])]    = "Low"
regen_yc$bin[paste0("cj_", colnames(jawh_yc)[which(jawh_yc$bin     == "Medium")])] = "Medium"
regen_yc$bin[paste0("cj_", colnames(jawh_yc)[which(jawh_yc$bin     == "High")])]   = "High"
regen_yc$bin[paste0("mi_", colnames(incsrh_yc)[which(incsrh_yc$bin == "Low")])]    = "Low"
regen_yc$bin[paste0("mi_", colnames(incsrh_yc)[which(incsrh_yc$bin == "Medium")])] = "Medium"
regen_yc$bin[paste0("mi_", colnames(incsrh_yc)[which(incsrh_yc$bin == "High")])]   = "High"
regen_yc$bin[paste0("imi_", colnames(imi_yc)[which(imi_yc$bin      == "Low")])]    = "Low"
regen_yc$bin[paste0("imi_", colnames(imi_yc)[which(imi_yc$bin      == "Medium")])] = "Medium"
regen_yc$bin[paste0("imi_", colnames(imi_yc)[which(imi_yc$bin      == "High")])]   = "High"
Idents(regen_yc) = regen_yc$bin
Idents(imm_yc) = imm_yc$bin
saveRDS(regen_yc, "regen_yc.rds")
saveRDS(imm_yc,   "imm_yc.rds")

print("Finding DEGs in Regen Celsr1+")
regen_yc_deg = FindAllMarkers(regen_yc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
regen_yc_deg$pct_dif = regen_yc_deg$pct.1 - regen_yc_deg$pct.2
regen_yc_deg$abs_pct_dif = regen_yc_deg$pct.1 - regen_yc_deg$pct.2
regen_yc_deg$isSig = regen_yc_deg$p_val_adj < 0.05
regen_yc_deg$cluster_isSig = paste0(regen_yc_deg$cluster, "_", regen_yc_deg$isSig)
write.csv(regen_yc_deg, "regen_yc_deg.csv")
print("Done")

print("Finding DEGs in IMM Celsr1+")
Idents(imm_yc) = imm_yc$bin
imm_yc_deg = FindAllMarkers(imm_yc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
imm_yc_deg$pct_dif = imm_yc_deg$pct.1 - imm_yc_deg$pct.2
imm_yc_deg$abs_pct_dif = imm_yc_deg$pct.1 - imm_yc_deg$pct.2
imm_yc_deg$isSig = imm_yc_deg$p_val_adj < 0.05
imm_yc_deg$cluster_isSig = paste0(imm_yc_deg$cluster, "_", imm_yc_deg$isSig)
write.csv(imm_yc_deg, "imm_yc_deg.csv")
print("Done")

# Celsr1- Cells
tjh_nc_res    = createCytoBINsInGeneNeg(tjh, "CELSR1")
jawh_nc_res   = createCytoBINsInGeneNeg(jawh, "CELSR1")
incsrh_nc_res = createCytoBINsInGeneNeg(incsrh, "CELSR1")
imi_nc_res    = createCytoBINsInGeneNeg(imi, "CELSR1")
imm_nc_res    = createCytoBINsInGeneNeg(imm, "CELSR1")

tjh_nc    = tjh_nc_res[[2]]
jawh_nc   = jawh_nc_res[[2]]
incsrh_nc = incsrh_nc_res[[2]]
imi_nc    = imi_nc_res[[2]]
imm_nc     = imm_nc_res[[2]]

regen_nc = subset(regen, cells = colnames(regen)[which(regen@assays$RNA@counts["CELSR1",] == 0)])
regen_nc$bin = "none"
regen_nc$bin[paste0("ct_", colnames(tjh_nc)[which(tjh_nc$bin       == "Low")])]    = "Low"
regen_nc$bin[paste0("ct_", colnames(tjh_nc)[which(tjh_nc$bin       == "Medium")])] = "Medium"
regen_nc$bin[paste0("ct_", colnames(tjh_nc)[which(tjh_nc$bin       == "High")])]   = "High"
regen_nc$bin[paste0("cj_", colnames(jawh_nc)[which(jawh_nc$bin     == "Low")])]    = "Low"
regen_nc$bin[paste0("cj_", colnames(jawh_nc)[which(jawh_nc$bin     == "Medium")])] = "Medium"
regen_nc$bin[paste0("cj_", colnames(jawh_nc)[which(jawh_nc$bin     == "High")])]   = "High"
regen_nc$bin[paste0("mi_", colnames(incsrh_nc)[which(incsrh_nc$bin == "Low")])]    = "Low"
regen_nc$bin[paste0("mi_", colnames(incsrh_nc)[which(incsrh_nc$bin == "Medium")])] = "Medium"
regen_nc$bin[paste0("mi_", colnames(incsrh_nc)[which(incsrh_nc$bin == "High")])]   = "High"
regen_nc$bin[paste0("imi_", colnames(imi_nc)[which(imi_nc$bin      == "Low")])]    = "Low"
regen_nc$bin[paste0("imi_", colnames(imi_nc)[which(imi_nc$bin      == "Medium")])] = "Medium"
regen_nc$bin[paste0("imi_", colnames(imi_nc)[which(imi_nc$bin      == "High")])]   = "High"
Idents(regen_nc) = regen_nc$bin
Idents(imm_nc) = imm_nc$bin
saveRDS(regen_nc, "regen_nc.rds")
saveRDS(imm_nc,   "imm_nc.rds")

print("Finding DEGs in Regen Celsr1-")
regen_nc_deg = FindAllMarkers(regen_nc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
regen_nc_deg$pct_dif = regen_nc_deg$pct.1 - regen_nc_deg$pct.2
regen_nc_deg$abs_pct_dif = regen_nc_deg$pct.1 - regen_nc_deg$pct.2
regen_nc_deg$isSig = regen_nc_deg$p_val_adj < 0.05
regen_nc_deg$cluster_isSig = paste0(regen_nc_deg$cluster, "_", regen_nc_deg$isSig)
write.csv(regen_nc_deg, "regen_nc_deg.csv")
print("Done")

print("Finding DEGs in IMM Celsr1-")
imm_nc = readRDS("imm_nc.rds")
Idents(imm_nc) = imm_nc$bin
imm_nc_deg = FindAllMarkers(imm_nc, only.pos = F, logfc.threshold = 0, min.pct = 0.01)
imm_nc_deg$pct_dif = imm_nc_deg$pct.1 - imm_nc_deg$pct.2
imm_nc_deg$abs_pct_dif = imm_nc_deg$pct.1 - imm_nc_deg$pct.2
imm_nc_deg$isSig = imm_nc_deg$p_val_adj < 0.05
imm_nc_deg$cluster_isSig = paste0(imm_nc_deg$cluster, "_", imm_nc_deg$isSig)
write.csv(imm_nc_deg, "imm_nc_deg.csv")
print("Done")
print("All Done")