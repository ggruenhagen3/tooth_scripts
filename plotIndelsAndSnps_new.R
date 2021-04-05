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
path_to_snps <- "C:/Users/miles/Downloads/d_tooth/data/snp.fst"
path_to_indels <- "C:/Users/miles/Downloads/d_tooth/data/indel.fst"


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
indels_data$Zfst <- -(( indels_data$WEIGHTED_FST - mean(indels_data$WEIGHTED_FST) ) / sd(indels_data$WEIGHTED_FST)+1)

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
