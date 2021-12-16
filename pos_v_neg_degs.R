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

regen_yc$isMes = F
regen_yc$isMes[which(regen_yc$annot %in% c("Pulp cells", "PDL"))] = T
regen_yc$isMes[which(regen_yc$sample == "TJ" & regen_yc$seurat_clusters %in% c("5", "7"))] = T
regen_yc$isMes[which(regen_yc$sample == "JPOOL" & regen_yc$seurat_clusters %in% c("6"))] = T
regen_yc_mes = subset(regen_yc, cells = colnames(regen_yc)[which(im$isMes)])

regen_yc_mes$bin <- regen_yc_mes$cyto
regen_yc_mes$bin[which(regen_yc_mes$cyto <= quantile(regen_yc_mes$cyto, 0.33))] <- "Low"
regen_yc_mes$bin[which(regen_yc_mes$cyto > quantile(regen_yc_mes$cyto, 0.33) & regen_yc_mes$cyto <= quantile(regen_yc_mes$cyto, 0.66))] <- "Medium"
regen_yc_mes$bin[which(regen_yc_mes$cyto > quantile(regen_yc_mes$cyto, 0.66))] <- "High"
regen_yc_mes$bin <- regen_yc_mes$cyto
regen_yc_mes$bin[which(regen_yc_mes$cyto <= 0.33)] <- "Low"
regen_yc_mes$bin[which(regen_yc_mes$cyto > 0.33 & regen_yc_mes$cyto <= 0.66)] <- "Medium"
regen_yc_mes$bin[which(regen_yc_mes$cyto > 0.66)] <- "High"

Idents(regen_yc_mes) = regen_yc_mes$bins
regen_yc_mes_deg = FindAllMarkers(regen_yc_mes)
regen_yc_mes_deg_sig = regen_yc_mes_deg[which(regen_yc_mes_deg$p_val_adj < 0.05),]
write.csv(regen_yc_mes_deg, "regen_yc_mes_bin_deg.csv")
write.csv(regen_yc_mes_deg_sig, "regen_yc_mes_bin_deg_sig.csv")

regen_yc$isEpi = F
regen_yc$isEpi[which(regen_yc$annot %in% c("Epithelial cells"))] = T
regen_yc$isEpi[which(regen_yc$sample == "TJ" & regen_yc$seurat_clusters %in% c("2", "3", "4"))] = T
regen_yc$isEpi[which(regen_yc$sample == "JPOOL" & regen_yc$seurat_clusters %in% c("0", "1", "3", "4", "7"))] = T
regen_yc_epi = subset(regen_yc, cells = colnames(regen_yc)[which(im$isEpi)])

regen_yc_epi$bin <- regen_yc_epi$cyto
regen_yc_epi$bin[which(regen_yc_epi$cyto <= quantile(regen_yc_epi$cyto, 0.33))] <- "Low"
regen_yc_epi$bin[which(regen_yc_epi$cyto > quantile(regen_yc_epi$cyto, 0.33) & regen_yc_epi$cyto <= quantile(regen_yc_epi$cyto, 0.66))] <- "Medium"
regen_yc_epi$bin[which(regen_yc_epi$cyto > quantile(regen_yc_epi$cyto, 0.66))] <- "High"
regen_yc_epi$bin <- regen_yc_epi$cyto
regen_yc_epi$bin[which(regen_yc_epi$cyto <= 0.33)] <- "Low"
regen_yc_epi$bin[which(regen_yc_epi$cyto > 0.33 & regen_yc_epi$cyto <= 0.66)] <- "Medium"
regen_yc_epi$bin[which(regen_yc_epi$cyto > 0.66)] <- "High"

Idents(regen_yc_epi) = regen_yc_epi$bins
regen_yc_epi_deg = FindAllMarkers(regen_yc_epi)
regen_yc_epi_deg_sig = regen_yc_epi_deg[which(regen_yc_epi_deg$p_val_adj < 0.05),]
write.csv(regen_yc_epi_deg, "regen_yc_epi_bin_deg.csv")
write.csv(regen_yc_epi_deg_sig, "regen_yc_epi_bin_deg_sig.csv")
