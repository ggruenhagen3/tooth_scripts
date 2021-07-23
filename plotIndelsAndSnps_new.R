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

fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)

test = sample(1:nrow(fst), 5000)
my_breaks = which(! duplicated(fst$CHROM))

image_name <- "~/scratch/brain/fst/rc.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Zfst, color = CHROM)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$CHROM), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
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
