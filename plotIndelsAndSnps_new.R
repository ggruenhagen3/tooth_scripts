####################
# Fst File to Plot #
####################

# Load modules if installed,
# if not installed, use 'install.packages("package_name")'
library("Cairo")
library("stringr")
library("RColorBrewer")
library("ggplot2")
library("fdrtool")

# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
# path_to_snps <- "C:/Users/miles_000/Downloads/Research/Tooth/results/snps.fst"
# path_to_indels <- "C:/Users/miles_000/Downloads/Research/Tooth/results/indels.fst"
# path_to_snps <- "C:/Users/miles/Downloads/d_tooth/data/snp.fst"
# path_to_indels <- "C:/Users/miles/Downloads/d_tooth/data/indel.fst"
path_to_snps <- "~/research/tooth/data/snp.fst"
path_to_indels <- "~/research/tooth/data/indel.fst"

# Load the data
snps_data <- read.table(path_to_snps, sep="\t", header=TRUE)
indels_data <- read.table(path_to_indels, sep="\t", header=TRUE)
# Convert the negative fst scores to 0
snps_data[, c(5)][snps_data[, c(5)] < 0] <- 0
indels_data[, c(5)][indels_data[, c(5)] < 0] <- 0
# Make the indels be plotted below the SNPs by making all their Zfst score inverted
# Add a column that just has row number
snps_data$ID <- seq.int(nrow(snps_data))
indels_data$ID <- seq.int(nrow(indels_data))

# indels_data <- readRDS(file = "indels_data.rds")
# Only to be run first time, save the data, then load it
# ~20 minute runtime
for (i in 1:nrow(indels_data)) {
  indel_start = indels_data[i,2]
  for (j in i:nrow(snps_data)) {
    snps_start = snps_data[j,2]
    snps_id = snps_data[j,7]
    if (indel_start == snps_start && indels_data$CHROM[i] == snps_data$CHROM[j]) {
      indels_data[i,7] <- snps_id
      break
    }
  }
}
indels_data <- indels_data[indels_data$ID != 0,]
# saveRDS(indels_data, file = "indels_data.rds")
# Add a column for Zfst
snps_data$Zfst <- (( snps_data$WEIGHTED_FST - mean(snps_data$WEIGHTED_FST) ) / sd(snps_data$WEIGHTED_FST)) + 1
indels_data$Zfst <- (( indels_data$WEIGHTED_FST - mean(indels_data$WEIGHTED_FST) ) / sd(indels_data$WEIGHTED_FST)+1)

# Circle
# In the circle plot, make sure to change indels_data$Zfst to positive, but in the traditional plot, it should be negative
total = snps_data
colnames(total)[which(colnames(total) == "Zfst")] = "Snp"
total$Indel = indels_data$Zfst[match(snps_data$ID, indels_data$ID)]
total$LG <- sub("^NW.*", "Unplaced", total$CHROM)  # rename unplaced contigs
total$LG <- sub("^NC_027944.1", "Unplaced", total$LG) # rename unplaced contigs
total$LG <- sub("^NC_0", "", total$LG)
total$LG[which(! is.na(as.numeric(total$LG)) )] = as.numeric(total$LG[which(! is.na(as.numeric(total$LG)) )]) - 36779.1
total$LG[which(total$LG %in% c("21", "22"))] = as.character(as.numeric(total$LG[which(total$LG %in% c("21", "22"))]) + 1)
total$LG = factor(total$LG, levels = c(as.character(1:20), 22:23, "Unplaced"))
total$Snp[which(total$Snp < 0)] = 0
total$Indel[which(total$Indel < 0)] = 0
snp_cutoff = min(snps_data$Zfst[which(p.adjust(2*pnorm(-abs(snps_data$Zfst)), method="BH") < 0.05)])
indel_cutoff = min(indels_data$Zfst[which(p.adjust(2*pnorm(-abs(indels_data$Zfst)), method="BH") < 0.05)])
pal = rep(brewer.pal(n=8,name="Dark2"), 6)
circos.par("gap.after" = 2, cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.02, 0.02))
circos.initialize(total$LG, xlim = cbind(rep(0, 23), aggregate(ID ~ LG, total, length)[,2]))
circos.track(ylim = c(min(total$Snp), max(total$Snp)), track.height = 0.15, bg.border = "black", panel.fun = function(x, y) {
  start_id = min(total$ID[which(total$LG == CELL_META$sector.index)])
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.points(total$ID[which(total$LG == CELL_META$sector.index)] - start_id, total$Snp[which(total$LG == CELL_META$sector.index)], col = paste0(pal[CELL_META$sector.numeric.index], "90"), pch = 20)
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(2),
              CELL_META$sector.index, facing = "downward", niceFacing = TRUE,
              adj = c(0.5, 0.5), cex = 0.6)
  circos.segments(CELL_META$cell.xlim[1], snp_cutoff, CELL_META$cell.xlim[2], snp_cutoff, col = "gray60", lty = 2)
})
circos.track(ylim = c(min(total$Indel, na.rm = T), max(total$Indel, na.rm = T)), track.height = 0.15, bg.border = "black", panel.fun = function(x, y) {
  start_id = min(total$ID[which(total$LG == CELL_META$sector.index)])
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.points(total$ID[which(total$LG == CELL_META$sector.index)] - start_id, total$Indel[which(total$LG == CELL_META$sector.index)], col = paste0(pal[CELL_META$sector.numeric.index], "90"), pch = 20)
  circos.segments(CELL_META$cell.xlim[1], indel_cutoff, CELL_META$cell.xlim[2], indel_cutoff, col = "gray60", lty = 2)
})

# Traditional
# In the traditional plot, make sure to change indels_data$Zfst to negative, but in the circle plot, it should be psotive
total <- rbind(snps_data, indels_data)

total$LG <- sub("^NW.*", "2", total$CHROM)  # rename unplaced contigs
total$LG <- sub("^NC_027944.1", "2", total$LG) # rename unplaced contigs
total$LG <- sub("^NC_0", "", total$LG)
# total$LG = as.numeric(as.vector(total$LG))
# total$col = as.numeric(as.vector(lg))

# Plot the data
palette(brewer.pal(n=8,name="Dark2"))
image_name <- "UMD2a_Bicuspid_vs_Tricuspid_Zfst_Separated.png"
#svg("UMD2a_Bicuspid_vs_Tricuspid_Zfst_Separated.svg", width = 1350, height = 401)
Cairo(image_name, type="png", , width = 12, height = 4, units = 'in', res=300)
# palette("default")
plot(total$ID, total$Zfst, col=as.numeric(as.character(total$LG)), xlab="Linkage Group", ylab=expression("Z"[fst]), xaxt="n", yaxt="n", bty="n", main="", pch=20, cex=0.30)
abline(h = 0)
title(main = expression("UMD2a Bicuspid vs Tricuspid Z"[fst]))

text(2, 15, labels="SNPs", cex=0.8)
text(2, -12, labels="Indels", cex=0.8)

text(2, 12, labels="10kb bins", cex=0.8)
text(2, -15, labels="10kb bins", cex=0.8)

yaxisat <- c(min(total$Zfst), 0, max(total$Zfst))
yaxislabels <- c("high", "low", "high")
axis(2, at=yaxisat, labels=yaxislabels)

# FDR Line
# snp_cutoff   = min(snps_data$Zfst[which(p.adjust(fdrtool(snps_data$Zfst, plot=FALSE)$pval, method="bonferroni") < 0.05)])
# indel_cutoff = max(indels_data$Zfst[which(p.adjust(fdrtool(indels_data$Zfst, cutoff.method = "pct0", plot=FALSE)$pval, method="bonferroni") < 0.05)])
snp_cutoff = min(snps_data$Zfst[which(p.adjust(2*pnorm(-abs(snps_data$Zfst)), method="BH") < 0.05)])
indel_cutoff = max(indels_data$Zfst[which(p.adjust(2*pnorm(-abs(indels_data$Zfst)), method="BH") < 0.05)])
abline(h=snp_cutoff, lty="9999", col="gray28")
abline(h=indel_cutoff, lty="9999", col="gray28")

xaxisat <- c(0, 3864, 7122, 10841, 13887, 17501, 21472, 27958, 30354, 32449, 35681, 38921, 42326, 45524, 49300, 52739, 56205, 59777, 62721, 65312, 68282, 71749, 75949)
xaxislabels <- as.character(1:20)
xaxislabels <- c(xaxislabels, c("22", "23", "unplaced"))
axis(1, at=xaxisat,labels=xaxislabels, cex.axis=0.75)
#dev.copy(png,"/mnt/c/Users/miles_000/Downloads/Research/Tooth/scripts/fst/BHLHE40_Bicuspid_vs_Tricuspid_Fst.png",width=8, height=10, units="in",res=750)
dev.off()

print(length(which(p.adjust(2*pnorm(-abs(snps_data$Zfst)), method="BH") < 0.05)))
print(length(which(p.adjust(2*pnorm(-abs(indels_data$Zfst)), method="BH") < 0.05)))


# Pit v Castle

fst = fst_rc
fst$id = 1:nrow(fst)
# fst$Zfst <- (( fst$WEIR_AND_COCKERHAM_FST - mean(fst$WEIR_AND_COCKERHAM_FST) ) / sd(fst$WEIR_AND_COCKERHAM_FST)) + 1
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST < 0)] = 0
fst$Zfst <- (( fst$WEIGHTED_FST - mean(fst$WEIGHTED_FST) ) / sd(fst$WEIGHTED_FST)) + 1

fst$LG = lgConverter(fst$CHROM)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
fst$LG[which(fst$CHROM == "unplaced")] = "unplaced"

test = sample(1:nrow(fst), 5000)
my_breaks = which(! duplicated(fst$LG))

image_name <- "~/scratch/brain/fst/rc.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Zfst, color = LG)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$LG), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
dev.off()
system(paste0("rclone copy ", image_name, " dropbox:BioSci-Streelman/George/Brain/fst/"))

image_name <- "~/scratch/brain/fst/pc_hist.png"
png(image_name, type="cairo", width = 6, height = 4, units = 'in', res=300)
hist(fst_pc_unbin$WEIR_AND_COCKERHAM_FST, main = "Pit-Castle", xlab = "FST", breaks = 50)
dev.off()
system(paste0("rclone copy ", image_name, " dropbox:BioSci-Streelman/George/Brain/fst/"))

# Pit v Castle Tang

fst = read.table("C:/Users/miles/Downloads/fst_pc_tang_10kb.windowed.weir.fst", sep="\t", header=TRUE)
fst$id = 1:nrow(fst)
# fst$Zfst <- (( fst$WEIR_AND_COCKERHAM_FST - mean(fst$WEIR_AND_COCKERHAM_FST) ) / sd(fst$WEIR_AND_COCKERHAM_FST)) + 1
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST < 0)] = 0
fst$Zfst <- (( fst$WEIGHTED_FST - mean(fst$WEIGHTED_FST) ) / sd(fst$WEIGHTED_FST)) + 1

fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)

test = sample(1:nrow(fst), 5000)
my_breaks = which(! duplicated(fst$CHROM))

image_name <- "C:/Users/miles/Downloads/tang_pc.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Zfst, color = CHROM)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$CHROM), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Pit vs Castle in Tanganyka")
print(p)
dev.off()

#***************************************************************************************
# MC vs Pit ============================================================================
#***************************************************************************************
pcrc2030 = read.csv("C:/Users/miles/Downloads/brain/data/markers/pcrc_FST20_30_LG11_evolution_genes_031821.csv")
pat = read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
pcrc2030[,c("LG", "START", "END")] = pat[match(pcrc2030$mzebra, pat$V2), c("V3", "V4", "V5")]
pcrc2030$START = as.numeric(pcrc2030$START)
pcrc2030$END   = as.numeric(pcrc2030$END)

fst = read.table("C:/Users/miles/Downloads/JTS09_pit_10kb.fst", header = T)
image_name = "C:/Users/miles/Downloads/JTS09_pit_10kb_highlight.png"
fst$id = 1:nrow(fst)
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST < 0)] = 0
fst$Zfst <- (( fst$WEIGHTED_FST - mean(fst$WEIGHTED_FST) ) / sd(fst$WEIGHTED_FST)) + 1

fst$LG = lgConverter(fst$CHROM)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
fst$LG[which(fst$CHROM == "unplaced")] = "unplaced"
fst$p = 2*pnorm(-abs(fst$Zfst))

test = sample(1:nrow(fst), 5000)
my_breaks = which(! duplicated(fst$LG))

fst$inHigh = F
fst$inHighGene = ""
fst$inHighStart = 0
fst$inHighEnd = 0
for (i in 1:nrow(pcrc2030)) {
  this_bool = fst$LG == pcrc2030$LG[i] & fst$BIN_START >= pcrc2030$START[i] & fst$BIN_END <= pcrc2030$END[i]
  this_idx = which(! fst$inHigh & this_bool )
  fst$inHigh[this_idx] = T
  fst$inHighGene[this_idx] = pcrc2030$mzebra[i]
  fst$inHighStart[this_idx] = pcrc2030$START[i]
  fst$inHighEnd[this_idx] = pcrc2030$END[i]
}
fst = fst[order(fst$inHigh),]
fst$LG = factor(fst$LG, levels = c(unique(fst$LG), "special"))
fst$LG[which(fst$inHigh)] = "special"
my_cols = lighten(rep(paste0(brewer.pal(n=8,name="Dark2"), "70"), 6), 0.2)
my_cols = my_cols[1:length(levels(fst$LG))]
my_cols[length(my_cols)] = "red"

# test$inHigh2 = test$inHigh
# test$inHigh2[which(! (test$inHighGene %in% c("mrpl13", "adamts16", "LOC101472004") | test$BIN_START >= 20740001 ) )] = F
png(image_name, type="cairo", width = 5, height = 4, units = 'in', res=300)
# p = ggplot(fst, aes(id, Zfst, color = LG)) + geom_point(size = 1) + theme_classic() + scale_color_manual(values = my_cols) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = levels(fst$LG)[1:22], expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("")
p = ggplot(fst[which(fst$CHROM == "NC_036790"),], aes(id, Zfst, color = inHigh)) + geom_point(size = 1.5) + theme_classic() + scale_color_manual(values = my_cols[c(11, 23)]) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = levels(fst$LG)[1:22], expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("")
print(p)
dev.off()

fst_mc_v_pit$CHROM_START = paste0(fst_mc_v_pit$CHROM, fst_mc_v_pit$BIN_START)
fst_mc_v_castle$CHROM_START = paste0(fst_mc_v_castle$CHROM, fst_mc_v_castle$BIN_START)
fst_mc_merge = list(fst_mc_v_pit, fst_mc_v_castle) %>% reduce(inner_join, by = "CHROM_START", suffix = c('_mc_v_pit', '_mc_v_castle'))
fst_mc_merge$Zfst_dif = fst_mc_merge$Zfst_mc_v_pit - fst_mc_merge$Zfst_mc_v_castle
my_breaks = which(! duplicated(fst_mc_merge$LG_mc_v_pit))
ggplot(fst_mc_merge, aes(x = id_mc_v_pit, y = Zfst_dif, color = LG_mc_v_pit)) + geom_point() + scale_color_manual(values = c(rep(brewer.pal(8, 'Dark2'), 6)[1:22], 'red') ) + scale_x_continuous(breaks = my_breaks, labels = unique(fst_mc_merge$LG_mc_v_pit), expand=c(0,0)) + ylab("Difference in Zfst (MC v Pit - MC v Castle)") + xlab("")

fst_rvp$CHROM_START = paste0(fst_rvp$CHROM, fst_rvp$BIN_START)
fst_rvc$CHROM_START = paste0(fst_rvc$CHROM, fst_rvc$BIN_START)
fst_r_merge = list(fst_rvp[which(fst_rvp$CHROM != "unplaced"),], fst_rvc[which(fst_rvc$CHROM != "unplaced"),]) %>% reduce(inner_join, by = "CHROM_START", suffix = c('_rvp', '_rvc'))
# fst_r_merge = merge(fst_rvp, fst_rvc, by = "CHROM_START", suffixes = c('_rvp', '_rvc'))
fst_r_merge$Zfst_dif = fst_r_merge$Zfst_rvp - fst_r_merge$Zfst_rvc
my_breaks = which(! duplicated(fst_r_merge$LG_rvp))
ggplot(fst_r_merge, aes(x = id_rvp, y = Zfst_dif, color = LG_rvp)) + geom_point() + scale_color_manual(values = c(rep(brewer.pal(8, 'Dark2'), 6)[1:22], 'red') ) + scale_x_continuous(breaks = my_breaks, labels = unique(fst_r_merge$LG_rvp), expand=c(0,0)) + ylab("Difference in Zfst (Rock v Pit - Rock v Castle)") + xlab("")


#***********************************
fst_mc_v_pit = fst
fst_pvc = fst
fst_rvc = fst

fst_mc_v_pit = fst_mc_v_pit[,c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST", "id", "Zfst", "LG", "p")]
fst_pvc = fst_pvc[,c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST", "id", "Zfst", "LG", "p")]
fst_rvc = fst_rvc[,c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST", "id", "Zfst", "LG", "p")]

fst_mc_v_pit$CHROM_START = paste0(fst_mc_v_pit$CHROM, fst_mc_v_pit$BIN_START)
fst_pvc$CHROM_START = paste0(fst_pvc$CHROM, fst_pvc$BIN_START)
fst_rvc$CHROM_START = paste0(fst_rvc$CHROM, fst_rvc$BIN_START)

library('tidyverse')
fst_merge = list(fst_mc_v_pit, fst_pvc, fst_rvc) %>% reduce(inner_join, by = "CHROM_START", suffix = c('_mc_v_pit', '_pvc'))
colnames(fst_merge)[(ncol(fst_merge)-9):ncol(fst_merge)] = paste0(tail(colnames(fst_merge), 10), '_rvc')

library("psych")
fst_merge$hmp = harmonic.mean(x=t(fst_merge[,c("p_pvc", 'p_mc_v_pit', 'p_rvc')]))
fst_merge$neg_log_hmp = -log10(fst_merge$hmp)
# fst_merge$num = as.numeric(reshape2::colsplit(fst_merge$LG_mc_v_pit, "LG", c('1', '2'))[,2])
# fst_merge = fst_merge[order(fst_merge$num),]
# fst_merge$new_id = 1:nrow(fst_merge)
my_breaks = which(! duplicated(fst$LG))
ggplot(fst_merge, aes(x = id_pvc, y = neg_log_hmp, color = LG_mc_v_pit)) + geom_point() + scale_color_manual(values = c(rep(brewer.pal(8, 'Dark2'), 6)[1:22], 'red') ) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$LG), expand=c(0,0))

fst_merge$q_mc_v_pit = p.adjust(fst_merge$p_mc_v_pit, method = "BH")
fst_merge$q_pvc      = p.adjust(fst_merge$p_pvc,      method = "BH")
fst_merge$q_rvc      = p.adjust(fst_merge$p_rvc,      method = "BH")

all_sig_gene = c()
for (gene in unique(fst_merge$inHighGene)) {
  this_gene_df = fst_merge[which(fst_merge$inHighGene == gene),]
  is_mc_v_pit_sig = any(this_gene_df$q_mc_v_pit < 0.05)
  is_pvc_sig = any(this_gene_df$q_pvc < 0.05)
  is_pvc_sig = any(this_gene_df$q_rvc < 0.05)
  if (is_mc_v_pit_sig && is_pvc_sig && is_rvc_sig) { all_sig_gene = c(all_sig_gene, gene) }
}

#***************************************************************************************
# Zack PC RC Plot ======================================================================
#***************************************************************************************
pc_fst = read.table("C:/Users/miles/Downloads/pit_castle_10kb.windowed.weir.fst", header = T)
rc_fst = read.table("C:/Users/miles/Downloads/rock_castle_10kb.windowed.weir.fst", header = T)

pc_fst$id = 1:nrow(pc_fst)
pc_fst$WEIGHTED_FST[which(pc_fst$WEIGHTED_FST < 0)] = 0
pc_fst$Zfst <- (( pc_fst$WEIGHTED_FST - mean(pc_fst$WEIGHTED_FST) ) / sd(pc_fst$WEIGHTED_FST)) + 1
pc_fst$LG = lgConverter(pc_fst$CHROM)
pc_fst$LG_BIN_START = paste0(pc_fst$LG, "_", pc_fst$BIN_START)
pc_fst = pc_fst[which(pc_fst$LG != "na"),]

rc_fst$id = 1:nrow(rc_fst)
rc_fst$WEIGHTED_FST[which(rc_fst$WEIGHTED_FST < 0)] = 0
rc_fst$Zfst <- (( rc_fst$WEIGHTED_FST - mean(rc_fst$WEIGHTED_FST) ) / sd(rc_fst$WEIGHTED_FST)) + 1
rc_fst$LG = lgConverter(rc_fst$CHROM)
rc_fst$LG_BIN_START = paste0(rc_fst$LG, "_", rc_fst$BIN_START)
rc_fst = rc_fst[which(rc_fst$LG != "na"),]

pcrc_fst = inner_join(pc_fst, rc_fst, by = "LG_BIN_START", suffix = c("_PC", "_RC"))
# write.csv(pcrc_fst, "C:/Users/miles/Downloads/pcrc_fst.csv")
# write.table(pcrc_fst_closest[,c(1, 2, 3)], "C:/Users/miles/Downloads/pcrc_fst_closest.bed", sep = "\t", row.names = F, col.names = F, quote = F)
pcrc_fst_closest = read.table("C:/Users/miles/Downloads/pcrc_fst_closest.bed", sep = "\t")
pcrc_fst_closest$gene = colsplit(pcrc_fst_closest$V12, ";", c("a", "b"))[,1]
pcrc_fst$gene = colsplit(pcrc_fst_closest$gene, "gene_id ", c("a", "b"))[,2]
pcrc_fst$gene_dist = pcrc_fst_closest$V13

pcrc_fst$Zfst_tot = pcrc_fst$Zfst_PC + pcrc_fst$Zfst_RC
pcrc_fst$isPCRC2030 = pcrc_fst$gene %in% pcrc2030$mzebra
ggplot(pcrc_fst, aes(x = Zfst_PC, y = Zfst_RC, color = Zfst_tot)) + geom_point() + geom_text_repel(data = pcrc_fst[which(pcrc_fst$Zfst_tot > 12 & pcrc_fst$isPCRC2030),], aes( label = gene)) + scale_color_gradientn(colors = plasma(100))
pcrc_fst = pcrc_fst[order(pcrc_fst$isPCRC2030),]
ggplot(pcrc_fst, aes(x = Zfst_PC, y = Zfst_RC, color = isPCRC2030)) + geom_point() + geom_text_repel(data = pcrc_fst[which(pcrc_fst$Zfst_tot > 12 & pcrc_fst$isPCRC2030),], aes( label = gene))

mc_pit$comp = "pit"
mc_castle$comp = "castle"
mc_all = rbind(mc_pit[which(mc_pit$inHigh),], mc_castle[which(mc_castle$inHigh),])
# ggplot(mc_all, aes(x = comp, y = Zfst)) + geom_boxplot() + geom_point(position = position_jitter())+ xlab("") + ggtitle("ZFst of BINs within PCRC2030 Genes")
mc_all = rbind(mc_pit, mc_castle)
ggplot(mc_all, aes(x = comp, y = Zfst, color = inHigh)) + geom_boxplot() + xlab("") + guides(color = guide_legend(title = "PCRC2030 Gene BIN"))

lg11_peak_start = min(pcrc_fst$BIN_START_PC[which(pcrc_fst$LG_PC == "LG11" & pcrc_fst$Zfst_PC >= 5.5)])
lg11_peak_end   = max(pcrc_fst$BIN_START_PC[which(pcrc_fst$LG_PC == "LG11" & pcrc_fst$Zfst_PC >= 5.5)])
mc_all$cat = "all"
mc_all$cat[which(mc_all$inHigh)] = "PCRC2030 Gene Bin"
mc_all$cat[which(mc_all$LG == "LG11" & mc_all$BIN_START >= lg11_peak_start & mc_all$BIN_START <= lg11_peak_end)] = "LG11 Peak"
ggplot(mc_all, aes(x = comp, y = Zfst, color = cat)) + geom_boxplot() + xlab("") + guides(color = guide_legend(title = ""))

#=======================================================================================
# Plot PCA =============================================================================
#=======================================================================================
fst = read.table("~/scratch/msc/pace_backup/snp_no_window.fst", sep = "\t", header = T)
vcf = read.table("~/scratch/msc/pace_backup/whole_filter_rm_SNP_biallelic_08_31_20.vcf", sep = "\t")
fst_pos = paste0(fst[,1], "_", fst[,2])
vcf_pos = paste0(vcf[,1], "_", vcf[,2])
sample_names = c("1619", "1818", "1860", "1863", "1912", "1983", "2162", "2241", "2277", "2298", "2302", "2319", "2320", "2332", "403", "404",  "493", "494", "495")
fst[,sample_names] = vcf[match(fst_pos, vcf_pos), 10:28]
write.table(fst, "~/scratch/msc/pace_backup/biallelic_snp_no_window_w_vcf.fst", sep = "\t", row.names = F)

# Start Here
sample_df = data.frame(samples = sample_names,
                       groups = c("bi", "bi", "bi", "bi", "bi", "bi", "tri", "tri", "bi", "tri", "tri", "tri", "bi", "tri", "tri", "tri", "tri", "tri", "tri"),
                       species = c("MZ", "MG", "MP", "MP", "PC", "PC", "PN", "PN", "ML", "LF", "LT", "LF", "ML", "PN", "LF", "LF", "LF", "LF", "LF"))
# PCA Method from: https://comppopgenworkshop2019.readthedocs.io/en/latest/contents/03_pca/pca.html
fst_simple = fst[1:3]
for (i in sample_names) {
  allele1 = substr(fst[,i], 0, 1)
  allele2 = substr(fst[,i], 3, 3)
  
  isRef1 = allele1 == 0
  isHomo = allele1 == allele2
  isMissing1 = allele1 == "."
  
  score = rep(-1, nrow(fst_simple)) # no score should be -1, if any is at the end, something went wrong
  score[which(isMissing1)] = 9
  score[which(!isMissing1 & isHomo & isRef1)] = 2
  score[which(!isMissing1 & isHomo & !isRef1)] = 0
  score[which(!isMissing1 & !isHomo )] = 1
  
  fst_simple = cbind(fst_simple, score)
}
colnames(fst_simple) = colnames(fst)
res.pca <- prcomp(t(fst_simple[,sample_names]))

fst_simple_lg5 = fst_simple[which(fst_simple$CHROM == "NC_036784.1"),]
res.pca.lg5 = prcomp(t(fst_simple_lg5[,sample_names]))

fst_simple_bhlhe40 = fst_simple[which(fst_simple$POS > (7528132 - 25000) & fst_simple$POS < (7531739 + 25000)),]
res.pca.bhlhe40 = prcomp(t(fst_simple_bhlhe40[,sample_names]))

library(factoextra)
png("~/scratch/d_tooth/results/ts_ms_pca_all_biallelic_snp.png", width = 500, height = 500)
fviz_pca_ind(res.pca, col.ind = sample_df$groups, palette = c("#00AFBB", "#FC4E07"), addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Groups", repel = FALSE, max.overlaps = Inf)
dev.off()
png("~/scratch/d_tooth/results/ts_ms_pca_all_biallelic_snp_species.png", width = 500, height = 500)
fviz_pca_ind(res.pca, col.ind = sample_df$species, addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Groups", repel = FALSE, max.overlaps = Inf)
dev.off()

png("~/scratch/d_tooth/results/ts_ms_pca_bhlhe40_biallelic_snp.png", width = 500, height = 500)
fviz_pca_ind(res.pca.bhlhe40, col.ind = sample_df$groups, palette = c("#00AFBB", "#FC4E07"), addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Groups", repel = FALSE, max.overlaps = Inf)
dev.off()
png("~/scratch/d_tooth/results/ts_ms_pca_bhlhe40_biallelic_snp_species.png", width = 500, height = 500)
fviz_pca_ind(res.pca.bhlhe40, col.ind = sample_df$species, addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Groups", repel = FALSE, max.overlaps = Inf)
dev.off()

png("~/scratch/d_tooth/results/ts_ms_pca_lg5_biallelic_snp.png", width = 500, height = 500)
fviz_pca_ind(res.pca.lg5, col.ind = sample_df$groups, palette = c("#00AFBB", "#FC4E07"), addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Groups", repel = FALSE, max.overlaps = Inf)
dev.off()
png("~/scratch/d_tooth/results/ts_ms_pca_lg5_biallelic_snp_species.png", width = 500, height = 500)
fviz_pca_ind(res.pca.lg5, col.ind = sample_df$species, addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Groups", repel = FALSE, max.overlaps = Inf)
dev.off()
