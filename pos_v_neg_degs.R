rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
data_folder = "~/scratch/d_tooth/data/"
allt_folder = "~/scratch/d_tooth/results/igor/allt/"
setwd(allt_folder)

regen_yc = readRDS("regen_yc2.rds")
regen_nc = readRDS("regen_nc2.rds")

regen_yc$celsr1 = "yc"
regen_nc$celsr1 = "nc"

regen_yc$celsr1_bin = paste0(regen_yc$celsr1, "_", regen_yc$bin)
regen_nc$celsr1_bin = paste0(regen_nc$celsr1, "_", regen_nc$bin)

regen_full = merge(regen_yc, regen_nc)

Idents(regen_full) = regen_full$celsr1_bin
yc_nc_bin_full = data.frame()
for (bin in c("High", "Medium", "Low")) {
  print(bin)
  bin_res = FindMarkers(regen_full, ident.1 = paste0("yc_", bin ), ident.1 = paste0("nc_", bin ))
  bin_res$gene = rownames(bin_res)
  bin_res$cluster = bin
  bin_res$bin = bin
  yc_nc_bin_full = rbind(yc_nc_bin_full, bin_res)
}
write.csv(yc_nc_bin_full, "yc_nc_bin_simple_deg_all.csv")

yc_nc_bin_full_sig = yc_nc_bin_full[which(yc_nc_bin_full$p_val_adj < 0.05),]
write.csv(yc_nc_bin_full_sig, "yc_nc_bin_simple_deg_sig.csv")