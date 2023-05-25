library(Cairo)

# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
# path_to_fst <- "C:/Users/miles_000/Downloads/Research/Tooth/results/x/x_w.fst"
# path_to_var <- "C:/Users/miles_000/Downloads/Research/Tooth/results/x/x_no_w.fst"
path_to_fst = "C:/Users/miles/Downloads/d_tooth/data/x_w.fst"
path_to_var = "C:/Users/miles/Downloads/d_tooth/data/x_no_w.fst"

# Load the binned data
fst_data <- read.table(path_to_fst, sep="\t", header=TRUE)
# Convert the negative fst scores to 0
fst_data[, c(5)][fst_data[, c(5)] < 0] <- 0
# Add a column that just has row number
fst_data$ID <- seq.int(nrow(fst_data))
# Add a column that for LG numbers
fst_data$LG <- sub("^NW.*", "36802.1", fst_data$CHROM)
fst_data$LG <- sub("^NC_027944.1", "36802.1", fst_data$LG)
fst_data$LG <- sub("^NC_0", "", fst_data$LG)
# Add a column for Zfst
fst_data$Zfst <- ( fst_data$WEIGHTED_FST - mean(fst_data$WEIGHTED_FST) ) / sd(fst_data$WEIGHTED_FST)
fst_data <- fst_data[,c("ID","Zfst", "LG", "CHROM", "BIN_START","BIN_END", "WEIGHTED_FST")]
lg5 <- fst_data[fst_data$CHROM == "NC_036784.1",]

# Load the unbinned data
var_fst_data <- read.table(path_to_var, sep="\t", header=TRUE)
# Rename the columns
names(var_fst_data) <- c("LG", "POS", "FST")
var_fst_data <- var_fst_data[var_fst_data$LG == "NC_036784.1",]
# Convert the negative fst scores to 0
var_fst_data[, c(3)][var_fst_data[, c(3)] < 0] <- 0

# Custom Color Palette
palette(c("grey68", "grey50", "grey0", "mediumpurple1", "mediumpurple4", "grey0", "cadetblue2", "cadetblue4", "grey0", "maroon2", "maroon4", "grey0"))
lg5$COL <- lg5$Zfst
lg5$COL[lg5$Zfst > 2.5] <- 2
lg5$COL[lg5$Zfst < 2.5] <- 1
var_fst_data$COL <- var_fst_data$FST
var_fst_data$COL[var_fst_data$FST > .99] <- 3
var_fst_data$COL[var_fst_data$FST < .99] <- 2
var_fst_data$COL[var_fst_data$FST < .50] <- 1

# Find important/high Zfst points
high <- subset(lg5, lg5$Zfst > 5)
high$Zfst = format(round(high$Zfst, 2), nsmall = 2)

# Constants for plots
lg5ystop <- 15
lg5Range <- 17501 - 13888
lg5xstart <- 13888
lg5xstop <- 17501
bufferImportantPoints <- 25
xplotsystart <- 0.45
xplotcex <- 0.6
genesystart <- 1.05
genesystop <- 1.25
geneslabely <- 1.15
geneslabelcex <- 0.4
subplotstop <- 2.0

# Constants for Plot 3
plot3col <- rgb(199,21,133, max=255, alpha=175)
plot3xstart <- 33000001
plot3xstop <- 36170001
plot3xstartID <- subset(lg5, BIN_START==plot3xstart)[1,1]
plot3xstopID <- subset(lg5, BIN_START==plot3xstop)[1,1]
plot3 <- var_fst_data[var_fst_data$POS > plot3xstart,]
plot3 <- plot3[plot3$POS < plot3xstop,]
lg5$COL[lg5$BIN_START >= plot3xstart & lg5$BIN_END <= plot3xstop] <- 11
plot3[,4] <- plot3[, 4] + 9
# Genes for Plot 3
plot3genes <- list( c(33000284, 33004528, "DNASE1L1"), c(33020839, 33022458, "GZMA"), c(33051343, 33053133, "GZMA"), c(33060440, 33075984, "DNASE1L1"), c(33076920, 33086420, "SEC13"), c(33086940, 33104562, "CAND2"), c(33113841, 33116502, "ADORA3"), c(33135724, 33146117, "NCR2"), c(33170391, 33180872, "NCR2"), c(33200828, 33206154, "ADORA3"), c(33246299, 33316457, "ANGPT4"), c(33329161, 33373791, "STK4"), c(33381629, 33405900, "KCNS1"), c(33508179, 33536598, "MATN4"), c(33536326, 33560184, "RBPJL"), c(33916796, 33973294, "DAG1"), c(34369205, 34379167, "NICN1"), c(34642256, 34643870, "TRIM35"), c(34675168, 34836681, "GRM7"), c(34871935, 34873272, "TRIM39"), c(34893168, 35004251, "RNF123"), c(35005643, 35015323, "ALDH4A1"), c(35044233, 35048954, "PPP1R15B"), c(35051777, 35055989, "PIK3C2B"), c(35059764, 35089134, "KCNQ2"), c(35117209, 35127472, "DDX27"), c(35131780, 35166595, "PLTP"), c(35168038, 35177215, "KHK"), c(35211660, 35300827, "PCBP4"), c(35326453, 35331785, "RBM15B"), c(35331891, 35339781, "MANF"), c(35341377, 35346427, "C20orf111"), c(35385545, 35400256, "SDC4"), c(35428748, 35438611, "GDF5"), c(35507962, 35525146, "UQCC"), c(35549711, 35565594, "FAM83C"), c(35569932, 35591346, "EML6"), c(35620272, 35705767, "MMP24"), c(35726737, 35731720, "EIF6"), c(35740551, 35748347, "MAP1LC3A"), c(35750710, 35753783, "DYNLRB1"), c(35755299, 35772139, "SSUH2"), c(35773337, 35879838, "SEMA3E"), c(35880419, 35890888, "SMPD4"), c(35893896, 35914333, "VPRBP"), c(35920834, 36170178, "DOCK3"), c(36131652, 36139062, "BTN1A1"))
plot3hgncgenes_noclutter <- list( c(33000284, 33004528, "DNASE1L1"), c(33020839, 33022458, "GZMA"), c(33051343, 33053133, "GZMA"), c(33060440, 33075984, "DNASE1L1"), c(33076920, 33086420, "SEC13"), c(33086940, 33104562, "CAND2"), c(33113841, 33116502, "ADORA3"), c(33135724, 33146117, "NCR2"), c(33170391, 33180872, "NCR2"), c(33200828, 33206154, "ADORA3"), c(33246299, 33316457, "ANGPT4"), c(33329161, 33373791, "STK4"), c(33381629, 33405900, "KCNS1"), c(33508179, 33536598, "MATN4"), c(33536326, 33560184, "RBPJL"), c(33916796, 33973294, "DAG1"), c(34369205, 34379167, "NICN1"), c(34642256, 34643870, "TRIM35"), c(34675168, 34836681, "GRM7"), c(34871935, 34873272, "TRIM39"), c(34893168, 35004251, "RNF123"), c(35005643, 35015323, "ALDH4A1"), c(35044233, 35048954, "PPP1R15B"), c(35059764, 35089134, "KCNQ2"), c(35117209, 35127472, "DDX27"), c(35131780, 35166595, "PLTP"), c(35168038, 35177215, "KHK"), c(35211660, 35300827, "PCBP4"), c(35326453, 35331785, "RBM15B"), c(35331891, 35339781, "MANF"), c(35341377, 35346427, "C20orf111"), c(35385545, 35400256, "SDC4"), c(35428748, 35438611, "GDF5"), c(35507962, 35525146, "UQCC"), c(35549711, 35565594, "FAM83C"), c(35569932, 35591346, "EML6"), c(35620272, 35705767, "MMP24"), c(35726737, 35731720, "EIF6"), c(35740551, 35748347, "MAP1LC3A"), c(35750710, 35753783, "DYNLRB1"), c(35755299, 35772139, "SSUH2"), c(35773337, 35879838, "SEMA3E"), c(35880419, 35890888, "SMPD4"), c(35893896, 35914333, "VPRBP"), c(35920834, 36170178, "DOCK3"), c(36131652, 36139062, "BTN1A1"))


png_name <- "x_top3LG5_3.png"
Cairo(png_name, type="png", , width = 13, height = 8, units = 'in', res=300)
plot(plot3$POS, plot3$FST, col = plot3$COL, xaxs="i", xlim=c(33000000,36170000), ylim=c(0,subplotstop), xlab="Position on LG5 (Mbp)", ylab="Fst", xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
axis(1, at=c(34000000,35000000), label=c(34,35), cex=xplotcex)
axis(1, at=c(33000000,36170000), label=c(33, 36.17), col.axis="black", cex = xplotcex)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1), label=c(0,0.2,0.4,0.6,0.8,1), cex =xplotcex)
# genes
previousgenestop <- 0
labelbuffer <- 0.1
cummulativelabelbuffer <- 0
for (gene in plot3genes) {
  genestart <- gene[1]
  genestop  <- gene[2]
  name      <- gene[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(genesystart, genesystop, genesystop, genesystart)
  polygon(x, y, col=plot3col, border=NA)
  if (as.numeric(genestart) - previousgenestop < 4000) {
    cummulativelabelbuffer <- cummulativelabelbuffer + labelbuffer
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=55, cex=geneslabelcex)
    if (cummulativelabelbuffer >= labelbuffer*7) { cummulativelabelbuffer <- 0 }
  }
  else {
    cummulativelabelbuffer <- 0
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=55, cex=geneslabelcex) 
  }
  previousgenestop <- as.numeric(genestop)
}
dev.off()