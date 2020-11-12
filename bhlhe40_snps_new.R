ggplot(snp_fst_data, aes(x=POS, y=FST, color=FST, alpha=FST+0.01)) + geom_point() + xlab("Position on LG5") + ylab("Fst") + theme_classic() + scale_x_continuous(expand = c(0,0)) + scale_color_gradient(low = "#56B1F7", high = "#132B43") + theme(legend.position="none")


# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
library("Cairo")
library("TigR")
path_to_snp <- "C:/Users/miles/Downloads/rna/data/snps_no_window.fst"
start <- 5594151
stop <- 8133628

# Load the data
snp_fst_data <- read.table(path_to_snp, sep="\t", header=TRUE)
# Rename the columns
names(snp_fst_data) <- c("LG", "POS", "FST")
snp_fst_data <- snp_fst_data[snp_fst_data$POS > start,]
snp_fst_data <- snp_fst_data[snp_fst_data$POS < stop,]
snp_fst_data <- snp_fst_data[snp_fst_data$LG == "NC_036784.1",]

# Convert the negative fst scores to 0
snp_fst_data[, c(3)][snp_fst_data[, c(3)] < 0] <- 0

# Custom Coloring Scheme
palette(c("grey68", "grey50", "grey0"))
snp_fst_data$COL <- snp_fst_data$FST
snp_fst_data$COL[snp_fst_data$FST > .99] <- 3
snp_fst_data$COL[snp_fst_data$FST < .99] <- 2
snp_fst_data$COL[snp_fst_data$FST < .50] <- 1

snp_fst_data$COL <- colorRampPalette(c("#f0f00c", "#ff0000"))(10)[as.numeric(cut(snp_fst_data$FST,breaks = 10))]
snp_fst_data$alpha <- seq(29,99,7)[as.numeric(cut(snp_fst_data$FST,breaks = 10))]
#snp_fst_data$alpha <- sapply(1:nrow(snp_fst_data), function(x) substrRight(format(round(snp_fst_data[x,c("FST")], 2), nsmall = 2),2) )
#snp_fst_data$alpha[which(as.numeric(as.vector(snp_fst_data$alpha)) < 0.4)] = "40"
## snp_fst_data$alpha[which(snp_fst_data$alpha == "00")] = "01"
#snp_fst_data$alpha[which(format(round(snp_fst_data$FST, 2), nsmall = 2) == "0.00")] <- "00"
snp_fst_data$COL <- paste0(snp_fst_data$COL, snp_fst_data$alpha)
snp_fst_data$COL[which(snp_fst_data$COL == "NANA")] <- "white"
## snp_fst_data$COL[which(format(round(snp_fst_data$FST, 2), nsmall = 2) == "1.00")] <- "black"

# Plot the data
png_name <- "BHLHE40_Bicuspid_vs_Tricuspid_Fst_SNP_Zoom_out.png"
Cairo(png_name, type="png", , width = 13, height = 4, units = 'in', res=300)
par(fig=c(0,1,0,0.8))
plot(snp_fst_data$POS, snp_fst_data$FST, col=snp_fst_data$COL, xlab="Position on LG5", ylab=expression("Fst"), main="", pch=20, xaxs="i")
# xaxisat <- c(start, 6000000, 6500000, 7000000, 7500000, 8000000, stop)
# axis(1, at=xaxisat, labels=xaxisat)
#Plot the genes
par(fig=c(0,1,0.35,1), new=TRUE)
plot(c(start, stop), c(1,1) , col="white", main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n", xaxs="i")
title(main = expression("BHLHE40 Bicuspid vs Tricuspid Fst - SNPs"), line = -0.4)
new_coord <- list( c(5594151, 5792593, "CNTN4"), c(5922372, 6050031, "CNTN4"), c(6080665, 6084765, "106675611"), c(6128021, 6128652, "112434873"), c(6137631, 6220665, "CHL1"), c(6352281, 6406729, "CNTN6"), c(6558824, 6657918, "PDZRN3"), c(6676147, 6684020, "PPP4R2"), c(6686902, 6704710, "GXYLT2"), c(6702319, 6752270, "SHQ1"), c(6780372, 6797823, "RYBP"), c(6864703, 6866489, "PROK1"), c(6869876, 6872788, "GPR27"), c(6877764, 6886369, "EIF4E"), c(6894717, 7040306, "FOXP1"),  c(7127935, 7139474, "101474366"), c(7167098,  7170405, "101474077"), c(7200706, 7235571, "MITFA"), c(7247226, 7289133, "FRMD4B"), c(7293110, 7297473, "ARL6IP5"), c(7297541, 7306936, "TRNT1"), c(7307116, 7310492, "AVPR2"), c(7405790, 7419704, "LRRN1"), c(7432720, 7438156, "SETMAR"), c(7445848, 7526086, "ITPR1"), c(7528132, 7531739, "bhlhe40"), c( 7543994, 7547291, ""), c (7549522, 7551078, ""), c(7555195, 7558074, ""), c(7567775, 7570522, ""), c(7574807, 7576428, ""), c(7586405, 7589692, ""), c(7591989, 7620486, ""), c(7629619, 7635457, "NTN4"), c(7635595, 7638445, "CHCHD4"), c(7645665, 7650070, "SEC61A1"), c(7748815, 7753100, "OGFR"), c(7754716, 7758783, "OGFR"), c(7759408, 7779041, "DIDO1"), c(7779313, 7801172, "ZMYND8"), c(7801154, 7836689, "NCOA3"), c(7841798, 7843092, "ID1"), c(7884943, 7918096, "SULF2"), c(7920378, 7950460, "ARHGAP39"), c(7952476, 7963638, "PHF20"), c(7964194, 7999246, "NDRG3"), c(8000818, 8004267, "SLA2"), c(8005017, 8017305, "TRPC4AP"), c(8017356, 8031751, "MYH7B"), c(8048347, 8051809, "SNAI1"), c(8058293, 8062175, "GGT7"), c(8063344, 8069444, "CYLD"), c(8071003, 8082439, "TOP1"), c(8083263, 8096859, "ZHX3"), c(8098606, 8102254, "FAM83D"), c(8119598, 8121529, "MAFB"), c(8123574, 8125956, "PREX1"), c(8126167, 8130526, "ADNP"), c(8131101, 8133628, "BCL2L1"))

xstart = 0
ystart = 0
ystop = 1
for (pair in new_coord) {
  genestart <- pair[1]
  genestop  <- pair[2]
  name      <- pair[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(ystart, ystop, ystop, ystart)
  polygon(x, y, col="#45454580", border=NA)
  if (name == "PROK1" || name == "NTN4" ) {
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(0.8), labels = name, srt=55, cex=0.5)
  }
  else {
    if ( name == "PPP4R2" || name == "GPR27" || name == "CHCHD4" || name == "ARL6IP5" || name == "FRMD4B") {
      text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(1), labels = name, srt=55, cex=0.5)
    }
    else {
      if (name == "OGFR" || name == "DIDO1" || name == "ZMYND8" || name == "NCOA3" || name == "ID1" || name == "SULF2" || name == "ARHGAP39" || name == "PHF20" || name == "NDRG3" || name == "SLA2" || name == "TRPC4AP" || name == "MYH7B" || name == "SNAI1" || name == "GGT7" || name == "CYLD" || name == "TOP1" || name == "ZHX3" || name == "FAM83D" || name == "MAFB" || name == "PREX1" || name == "ADNP" || name == "BCL2L1") {
        # Do Nothing
      }
      else {
        if (name == "AVPR2") {
          text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(0.75), labels = name, srt=55, cex=0.5) 
        }
        else {
          if (name == "bhlhe40") {
            text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(0.75), labels = name, srt=55, cex=0.425, font=2) 
          }
          else {
            text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(0.8), labels = name, srt=55, cex=0.5) 
          }
        }
      }
    }
  }
}

# Area for the genes at the end
text(7921221 ,as.numeric(1.35), labels = "OGFR,OGFR,DIDO1,ZMYND8,NCOA3,ID1,SULF2,", cex=0.5)
text(7921221, as.numeric(1.25), labels = "ARHGAP39,PHF20,NDRG3,SLA2,TRPC4AP,", cex=0.5)
text(7921221, as.numeric(1.15), labels = "MYH7B,SNAI1,GGT7,CYLD,TOP1,ZHX3,", cex=0.5)
text(7921221, as.numeric(1.05), labels = "FAM83D,MAFB,PREX1,ADNP,BCL2L1", cex=0.5)
segments(7746815, ystop, 7726000, 1.3)
segments(8133628, ystop, 8112813, 1.3)

# lncRNA's
text(7582240, as.numeric(1.35), labels = "lncRNAs", cex=0.5)
segments(7543994, ystop, 7582240, 1.3)
segments(7620486, ystop, 7582240, 1.3)

dev.off()