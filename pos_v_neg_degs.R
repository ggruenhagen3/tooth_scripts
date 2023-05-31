# Helper Functions
createCytoBINsInGene = function(obj) {
  obj$cyto = 0
  obj$bin = 0
  for (s in c("TJ", "JPOOL", "IMI")) {
    this_cyto = CytoTRACE::CytoTRACE(as.matrix(obj@assays$RNA@counts[which(obj$sample == s)]))
    this_bin = this_cyto
    this_bin[which(this_cyto <= quantile(this_cyto, 0.33))] <- "Low"
    this_bin[which(this_cyto > quantile(this_cyto, 0.33) & this_cyto <= quantile(this_cyto, 0.66))] <- "Medium"
    this_bin[which(this_cyto > quantile(this_cyto, 0.66))] <- "High"
    this_bin <- this_cyto
    this_bin[which(this_cyto <= 0.33)] <- "Low"
    this_bin[which(this_cyto > 0.33 & this_cyto <= 0.66)] <- "Medium"
    this_bin[which(this_cyto > 0.66)] <- "High"
    
    obj$cyto[which(obj$sample == s)] = this_cyto
    obj$bin[which(obj$sample == s)] = this_bin
  }
  
  return(obj)
}

# Body
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
data_folder = "~/scratch/d_tooth/data/"
allt_folder = "~/scratch/d_tooth/results/igor/allt/"
setwd(allt_folder)

regen_yc = readRDS("regen_yc2.rds")
regen_nc = readRDS("regen_nc2.rds")

# regen_yc$celsr1 = "yc"
# regen_nc$celsr1 = "nc"
# 
# regen_yc$celsr1_bin = paste0(regen_yc$celsr1, "_", regen_yc$bin)
# regen_nc$celsr1_bin = paste0(regen_nc$celsr1, "_", regen_nc$bin)
# 
# regen_full = merge(regen_yc, regen_nc)
# 
# Idents(regen_full) = regen_full$celsr1_bin
# yc_nc_bin_full = data.frame()
# for (bin in c("High", "Medium", "Low")) {
#   print(bin)
#   bin_res = FindMarkers(regen_full, ident.1 = paste0("yc_", bin ), ident.2 = paste0("nc_", bin ))
#   bin_res$gene = rownames(bin_res)
#   bin_res$cluster = bin
#   bin_res$bin = bin
#   yc_nc_bin_full = rbind(yc_nc_bin_full, bin_res)
# }
# write.csv(yc_nc_bin_full, "yc_nc_bin_simple_deg_all.csv")
# 
# yc_nc_bin_full_sig = yc_nc_bin_full[which(yc_nc_bin_full$p_val_adj < 0.05),]
# write.csv(yc_nc_bin_full_sig, "yc_nc_bin_simple_deg_sig.csv")

regen_yc$sample[which(is.na(regen_yc$sample))] = "IMI"
regen_yc$isCichlid = F
regen_yc$isCichlid[which(regen_yc$sample %in% c("TJ", "JPOOL"))] = T

regen_yc$isMes = F
regen_yc$isMes[which(regen_yc$annot %in% c("Pulp cells", "PDL"))] = T
regen_yc$isMes[which(regen_yc$sample == "TJ" & regen_yc$seurat_clusters %in% c("5", "7"))] = T
regen_yc$isMes[which(regen_yc$sample == "JPOOL" & regen_yc$seurat_clusters %in% c("6"))] = T
regen_yc_mes = subset(regen_yc, cells = colnames(regen_yc)[which(regen_yc$isMes)])

Idents(regen_yc_mes) = regen_yc_mes$bin
regen_yc_mes_deg = FindAllMarkers(regen_yc_mes)
regen_yc_mes_deg_sig = regen_yc_mes_deg[which(regen_yc_mes_deg$p_val_adj < 0.05),]
write.csv(regen_yc_mes_deg, "regen_yc_mes_bin_deg.csv")
write.csv(regen_yc_mes_deg_sig, "regen_yc_mes_bin_deg_sig.csv")

regen_yc$isEpi = F
regen_yc$isEpi[which(regen_yc$annot %in% c("Epithelial cells"))] = T
regen_yc$isEpi[which(regen_yc$sample == "TJ" & regen_yc$seurat_clusters %in% c("2", "3", "4"))] = T
regen_yc$isEpi[which(regen_yc$sample == "JPOOL" & regen_yc$seurat_clusters %in% c("0", "1", "3", "4", "7"))] = T
regen_yc_epi = subset(regen_yc, cells = colnames(regen_yc)[which(regen_yc$isEpi)])

Idents(regen_yc_epi) = regen_yc_epi$bin
regen_yc_epi_deg = FindAllMarkers(regen_yc_epi)
regen_yc_epi_deg_sig = regen_yc_epi_deg[which(regen_yc_epi_deg$p_val_adj < 0.05),]
write.csv(regen_yc_epi_deg, "regen_yc_epi_bin_deg.csv")
write.csv(regen_yc_epi_deg_sig, "regen_yc_epi_bin_deg_sig.csv")

regen_yc$MesEpi = "Other"
regen_yc$MesEpi[which(regen_yc$isMes)] = "Mesenchyme"
regen_yc$MesEpi[which(regen_yc$isEpi)] = "Epithelium"
df = aggregate(nCount_RNA ~ MesEpi + bin, regen_yc@meta.data, length)
bin_df = aggregate(nCount_RNA ~ bin, regen_yc@meta.data, length)
df$prop_of_bin = df$nCount_RNA / bin_df$nCount_RNA[match(df$bin, bin_df$bin)] * 100
df$bin = factor(df$bin, levels = c("High", "Medium", "Low"))

temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
png("~/scratch/d_tooth/results/igor/allt/third_draft/regen_yc_mes_epi.png")
ggplot(df, aes(x = bin, y = prop_of_bin, fill = MesEpi)) + geom_bar(stat = 'identity', position = position_dodge2()) + scale_fill_manual(values = brewer.pal(9,"Blues")[c(3, 5, 8)]) + xlab("") + ylab("Proportion of Cells in Bin") + theme_bw() + scale_y_continuous(expand = c(0,0))
dev.off()