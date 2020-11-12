library(Cairo)

# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
# path_to_fst <- "C:/Users/miles_000/Downloads/Research/Tooth/results/fst.vcf.windowed.weir.fst"
# path_to_var <- "C:/Users/miles_000/Downloads/Research/Tooth/results/no_window.fst"
path_to_fst = "C:/Users/miles/Downloads/d_tooth/data/whole_filter_rm_window_08_31_20.fst"
path_to_var = "C:/Users/miles/Downloads/d_tooth/data/whole_filter_rm_no_window_08_31_20.fst"

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
my_colors <- c("grey68", "grey50", "grey0", "gold2", "gold4", "grey0", "chartreuse2", "chartreuse4", "grey0", "cadetblue2", "cadetblue4", "grey0")
my_colors_2 <- c("grey68", "grey50", "grey0", "gold2", "gold4", "grey0", "chartreuse2", "green4", "grey0", "cadetblue2", "blue4", "grey0")
palette(my_colors)
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

# Constants for Plot 1
plot1col <- rgb(255, 255, 0, max=255, alpha=175)
plot1xstart <- 4140001
plot1xstop <- 7440001
plot1xstartID <- subset(lg5, BIN_START==plot1xstart)[1,1]
plot1xstopID <- subset(lg5, BIN_START==plot1xstop)[1,1]
lg5$COL[lg5$BIN_START >= plot1xstart & lg5$BIN_END <= plot1xstop] <- 5
plot1 <- var_fst_data[var_fst_data$POS > plot1xstart,]
plot1 <- plot1[plot1$POS < plot1xstop,]
plot1[,4] <- plot1[, 4] + 3
plot1$COL <- colorRampPalette(c(my_colors_2[4], my_colors_2[5]))(10)[as.numeric(cut(plot1$FST,breaks = 10))]
plot1$alpha <- seq(29,99,7)[as.numeric(cut(plot1$FST,breaks = 10))]
plot1$COL <- paste0(plot1$COL, plot1$alpha)
plot1$COL[which(plot1$COL == "NANA")] <- "white"
plot1$COL[which(plot1$FST == 1)] <- "black"
# Genes for Plot 1
plot1genes <- list( c(4159928, 4171786, "SLC45A3"), c(4175323, 4183757, "KLF15"), c(4193028, 4218295, "ELK4"), c(4218399, 4234891, "MFSD4"), c(4238528, 4309859, "CDK18"), c(4281108, 4283214, "LOC105941108"), c(4405299, 4407285, "ZBED4"), c(4417800, 4417871, "trnaq-cug"), c(4442733, 4473382, "LOC106675000"), c(4520496, 4597123, "KLHL17"), c(4618584, 4653072, "SLC41A1"), c(4640404, 4643164, "LOC112434816"), c(4654310, 4660810, "EFCC1"), c(4686659, 4705584, "ETNK2"), c(4716042, 4767253, "SOX13"), c(4807699, 4832038, "PBLD"), c(4832860, 4836914, "SNRPE"), c(4837514, 4848576, "CDC42"), c(4859650, 4885996, "WNT4"), c(5055328, 5058235, "CCDC23"), c(5058614, 5063063, "FBXO2"), c(5065750, 5117199, "ITGA5"), c(5126330, 5139991, "BIN2"), c(5145998, 5157055, "TXNRD3"), c(5172266, 5191290, "CHST13"), c(5197030, 5255185, "PHACTR3"), c(5257224, 5261272, "ACOT8"), c(5261507, 5265925, "SRSF6"), c(5267946, 5272821, "IFT52"), c(5274456, 5283440, "MYBL2"), c(5284687, 5296315, "RPN2"), c(5299543, 5302127, "GHRH"), c(5307024, 5320582, "CDK5RAP1"), c(5323503, 5327027, "LOC112434867"), c(5329156, 5332516, "LOC112434797"), c(5332951, 5334446, "LOC112434811"), c(5334837, 5345548, "LOC112434886"), c(5364137, 5364732, "LOC112434885"), c(5377530, 5378575, "CDK5RAP1"), c(5378643, 5380908, "BPI"), c(5433728, 5437330, "LOC112434812"), c(5451075, 5457038, "BPI"), c(5488563, 5589138, "EPB41L1"), c(5594151, 5792593, "CNTN4"), c(5922372, 6050031, "CNTN4"), c(6080665, 6084765, "LOC106675611"), c(6128021, 6128652, "LOC112434873"), c(6137631, 6220665, "CHL1"), c(6352281, 6406729, "CNTN6"), c(6558824, 6657918, "PDZRN3"), c(6676147, 6684020, "PPP4R2"), c(6686902, 6704710, "GXYLT2"), c(6702319, 6752270, "SHQ1"), c(6780372, 6797823, "RYBP"), c(6864703, 6866489, "PROK1"), c(6869876, 6872788, "GPR27"), c(6877764, 6886369, "EIF4E"), c(6894717, 7040306, "FOXP1"), c(7127935, 7139474, "LOC101474366"), c(7167098, 7170405, "LOC101474077"), c(7200706, 7235571, "MITF"), c(7247226, 7289133, "FRMD4B"), c(7293110, 7297473, "ARL6IP5"), c(7297541, 7306936, "TRNT1"), c(7307116, 7310492, "AVPR2"), c(7405790, 7419704, "LRRN1"), c(7432720, 7438156, "LOC101470671"))
plot1hgncgenes <- list( c(4159928, 4171786, "SLC45A3"), c(4175323, 4183757, "KLF15"), c(4193028, 4218295, "ELK4"), c(4218399, 4234891, "MFSD4"), c(4238528, 4309859, "CDK18"), c(4405299, 4407285, "ZBED4"), c(4520496, 4597123, "KLHL17"), c(4618584, 4653072, "SLC41A1"), c(4654310, 4660810, "EFCC1"), c(4686659, 4705584, "ETNK2"), c(4716042, 4767253, "SOX13"), c(4807699, 4832038, "PBLD"), c(4832860, 4836914, "SNRPE"), c(4837514, 4848576, "CDC42"), c(4859650, 4885996, "WNT4"), c(5055328, 5058235, "CCDC23"), c(5058614, 5063063, "FBXO2"), c(5065750, 5117199, "ITGA5"), c(5126330, 5139991, "BIN2"), c(5145998, 5157055, "TXNRD3"), c(5172266, 5191290, "CHST13"), c(5197030, 5255185, "PHACTR3"), c(5257224, 5261272, "ACOT8"), c(5261507, 5265925, "SRSF6"), c(5267946, 5272821, "IFT52"), c(5274456, 5283440, "MYBL2"), c(5284687, 5296315, "RPN2"), c(5299543, 5302127, "GHRH"), c(5307024, 5320582, "CDK5RAP1"), c(5377530, 5378575, "CDK5RAP1"), c(5378643, 5380908, "BPI"), c(5451075, 5457038, "BPI"), c(5488563, 5589138, "EPB41L1"), c(5594151, 5792593, "CNTN4"), c(5922372, 6050031, "CNTN4"), c(6137631, 6220665, "CHL1"), c(6352281, 6406729, "CNTN6"), c(6558824, 6657918, "PDZRN3"), c(6676147, 6684020, "PPP4R2"), c(6686902, 6704710, "GXYLT2"), c(6702319, 6752270, "SHQ1"), c(6780372, 6797823, "RYBP"), c(6864703, 6866489, "PROK1"), c(6869876, 6872788, "GPR27"), c(6877764, 6886369, "EIF4E"), c(6894717, 7040306, "FOXP1"), c(7200706, 7235571, "MITF"), c(7247226, 7289133, "FRMD4B"), c(7293110, 7297473, "ARL6IP5"), c(7297541, 7306936, "TRNT1"), c(7307116, 7310492, "AVPR2"), c(7405790, 7419704, "LRRN1"))


# Constants for Plot 2
plot2col <- rgb(0, 255, 0, max=255, alpha=175)
plot2xstart <- 20830001
plot2xstop <- 23330001
plot2xstartID <- subset(lg5, BIN_START==plot2xstart)[1,1]
plot2xstopID <- subset(lg5, BIN_START==plot2xstop)[1,1]
plot2 <- var_fst_data[var_fst_data$POS > plot2xstart,]
plot2 <- plot2[plot2$POS < plot2xstop,]
lg5$COL[lg5$BIN_START >= plot2xstart & lg5$BIN_END <= plot2xstop] <- 8
plot2[,4] <- plot2[, 4] + 6
plot2$COL <- colorRampPalette(c(my_colors_2[7], my_colors_2[8]))(10)[as.numeric(cut(plot2$FST,breaks = 10))]
plot2$alpha <- seq(29,99,7)[as.numeric(cut(plot2$FST,breaks = 10))]
plot2$COL <- paste0(plot2$COL, plot2$alpha)
plot2$COL[which(plot2$COL == "NANA")] <- "white"
plot2$COL[which(plot2$FST == 1)] <- "black"
# Genes for Plot 2
plot2genes <- list( c(20834882, 20850871, "MST1"), c(20851959, 20860716, "HMGA2"), c(20862667, 20867100, "C3orf18"), c(20868713, 20878609, "CALML5"), c(20874338, 20885219, "GNB1"), c(20886603, 20890367, "TARDBP"), c(20890796, 20896104, "APITD1-CORT"), c(20897518, 20898638, "RBP7"), c(20902260, 20909067, "H6PD"), c(20916392, 20929391, "EPHA2"), c(20931552, 20932439, "HSPB7"), c(20932855, 20935509, "RSG1"), c(20937240, 20941904, "ARHGEF19"), c(20942464, 20952607, "DDI2"), c(20960080, 20963036, "TNFRSF14"), c(20968913, 20971181, "TNFRSF14"), c(20979768, 20983303, "TNFRSF8"), c(20984263, 20989684, "ARL8B"), c(20990983, 21002192, "EDEM1"), c(21004514, 21008253, "GNAT1"), c(21011309, 21057398, "SEMA3F"), c(21110849, 21125353, "RBM6"), c(21126195, 21142379, "CAMKV"), c(21153520, 21166796, "MST1R"), c(21168897, 21201665, "UBA7"), c(21203302, 21205876, "FAM212A"), c(21207278, 21256866, "GNAI2"), c(21259492, 21289854, "SLC38A5"), c(21292979, 21303573, "RBM5"), c(21304427, 21306177, "BRK1"), c(21306980, 21308722, "GPX2"), c(21358351, 21386821, "GRM2"), c(21387925, 21395349, "RRP9"), c(21395501, 21401387, "LOC101473799"), c(21402097, 21406337, "GPR61"), c(21413051, 21415557, "CAV2"), c(21415847, 21499740, "SRGAP3"), c(21434319, 21435164, "LOC112434822"), c(21510146, 21520938, "PTGIS"), c(21521602, 21563132, "KCNB1"), c(21579834, 21588703, "STAU1"), c(21588867, 21603528, "LOC101477452"), c(21606334, 21625857, "CAMK1"), c(21626239, 21641350, "FAM3C"), c(21642142, 21677641, "C3orf67"), c(21755505, 21873152, "PTPRG"), c(21784439, 21785956, "LOC112434823"), c(21876080, 21942932, "CADPS"), c(21957430, 21967085, "SYNPR"), c(21969464, 21971537, "RHO"), c(21980328, 21985726, "LOC101482996"), c(21987096, 21989317, "RHO"), c(21996068, 21999492, "RHO"), c(22019318, 22025722, "SLC6A12"), c(22030503, 22052030, "SLC6A12"), c(22054364, 22067083, "FKBP5"), c(22067365, 22080218, "SRSF3"), c(22080420, 22098858, "LOC101485711"), c(22099451, 22104060, "TRUB2"), c(22104184, 22108329, "HDHD3"), c(22109427, 22119679, "TGFA"), c(22122493, 22147554, "SLC2A11"), c(22155087, 22168063, "PATZ1"), c(22169498, 22185057, "EIF4ENIF1"), c(22185313, 22195622, "SFI1"), c(22195620, 22205506, "TF"), c(22207494, 22211667, "TNFSF14"), c(22212879, 22214141, "TRIM35"), c(22215012, 22224663, "GMPPB"), c(22235076, 22242643, "PDE12"), c(22243203, 22252218, "ATP6AP1L"), c(22254250, 22256232, "RPS23"), c(22258147, 22269495, "ASB14"), c(22269597, 22306245, "IP6K1"), c(22307606, 22321531, "ccdc66"), c(22324196, 22351316, "FAM208A"), c(22356741, 22369609, "ARHGEF3"), c(22361993, 22364980, "LOC105940873"), c(22390748, 22423769, "IL17RD"), c(22430051, 22436841, "KLHL10"), c(22439258, 22454910, "APPL1"), c(22459858, 22463558, "ARF1"), c(22464287, 22469039, "ARF1"), c(22487270, 22494315, "NPTX1"), c(22513818, 22519177, "NPTX1"), c(22524523, 22542365, "DENND6A"), c(22550919, 22633366, "SLMAP"), c(22638861, 22702632, "FLNB"), c(22703807, 22706879, "TRIM35"), c(22709578, 22712114, "TRIM35"), c(22712897, 22715744, "TRIM35"), c(22718507, 22730927, "OSGIN2"), c(22757753, 22762951, "CCDC51"), c(22762995, 22766354, "TMA7"), c(22767930, 22781328, "DNASE1L1"), c(22781582, 22784288, "TRIM35"), c(22788211, 22791494, "TRIM35"), c(22809538, 22813106, "CBLN3"), c(22817362, 22822347, "C1QC"), c(22825851, 22830064, "COL10A1"), c(22835447, 22844312, "LOC105940861"), c(22856509, 22936343, "MTSS1"), c(23000434, 23116782, "RALY"), c(23131658, 23148552, "ASIP"), c(23226948, 23229218, "LOC112434877"), c(23238420, 23248456, "AHCY"), c(23254545, 23264567, "CHMP4B"), c(23267333, 23268583, "C20orf85"), c(23272013, 23275019, "LOC101472532"), c(23277144, 23283453, "LOC101482013"), c(23295242, 23300604, "PCK1"))
plot2hgncgenes <- list( c(20834882, 20850871, "MST1"), c(20851959, 20860716, "HMGA2"), c(20862667, 20867100, "C3orf18"), c(20868713, 20878609, "CALML5"), c(20874338, 20885219, "GNB1"), c(20886603, 20890367, "TARDBP"), c(20890796, 20896104, "APITD1-CORT"), c(20897518, 20898638, "RBP7"), c(20902260, 20909067, "H6PD"), c(20916392, 20929391, "EPHA2"), c(20931552, 20932439, "HSPB7"), c(20932855, 20935509, "RSG1"), c(20937240, 20941904, "ARHGEF19"), c(20942464, 20952607, "DDI2"), c(20960080, 20963036, "TNFRSF14"), c(20968913, 20971181, "TNFRSF14"), c(20979768, 20983303, "TNFRSF8"), c(20984263, 20989684, "ARL8B"), c(20990983, 21002192, "EDEM1"), c(21004514, 21008253, "GNAT1"), c(21011309, 21057398, "SEMA3F"), c(21110849, 21125353, "RBM6"), c(21126195, 21142379, "CAMKV"), c(21153520, 21166796, "MST1R"), c(21168897, 21201665, "UBA7"), c(21203302, 21205876, "FAM212A"), c(21207278, 21256866, "GNAI2"), c(21259492, 21289854, "SLC38A5"), c(21292979, 21303573, "RBM5"), c(21304427, 21306177, "BRK1"), c(21306980, 21308722, "GPX2"), c(21358351, 21386821, "GRM2"), c(21387925, 21395349, "RRP9"), c(21402097, 21406337, "GPR61"), c(21413051, 21415557, "CAV2"), c(21415847, 21499740, "SRGAP3"), c(21510146, 21520938, "PTGIS"), c(21521602, 21563132, "KCNB1"), c(21579834, 21588703, "STAU1"), c(21606334, 21625857, "CAMK1"), c(21626239, 21641350, "FAM3C"), c(21642142, 21677641, "C3orf67"), c(21755505, 21873152, "PTPRG"), c(21876080, 21942932, "CADPS"), c(21957430, 21967085, "SYNPR"), c(21969464, 21971537, "RHO"), c(21987096, 21989317, "RHO"), c(21996068, 21999492, "RHO"), c(22019318, 22025722, "SLC6A12"), c(22030503, 22052030, "SLC6A12"), c(22054364, 22067083, "FKBP5"), c(22067365, 22080218, "SRSF3"), c(22099451, 22104060, "TRUB2"), c(22104184, 22108329, "HDHD3"), c(22109427, 22119679, "TGFA"), c(22122493, 22147554, "SLC2A11"), c(22155087, 22168063, "PATZ1"), c(22169498, 22185057, "EIF4ENIF1"), c(22185313, 22195622, "SFI1"), c(22195620, 22205506, "TF"), c(22207494, 22211667, "TNFSF14"), c(22212879, 22214141, "TRIM35"), c(22215012, 22224663, "GMPPB"), c(22235076, 22242643, "PDE12"), c(22243203, 22252218, "ATP6AP1L"), c(22254250, 22256232, "RPS23"), c(22258147, 22269495, "ASB14"), c(22269597, 22306245, "IP6K1"), c(22324196, 22351316, "FAM208A"), c(22356741, 22369609, "ARHGEF3"), c(22390748, 22423769, "IL17RD"), c(22430051, 22436841, "KLHL10"), c(22439258, 22454910, "APPL1"), c(22459858, 22463558, "ARF1"), c(22464287, 22469039, "ARF1"), c(22487270, 22494315, "NPTX1"), c(22513818, 22519177, "NPTX1"), c(22524523, 22542365, "DENND6A"), c(22550919, 22633366, "SLMAP"), c(22638861, 22702632, "FLNB"), c(22703807, 22706879, "TRIM35"), c(22709578, 22712114, "TRIM35"), c(22712897, 22715744, "TRIM35"), c(22718507, 22730927, "OSGIN2"), c(22757753, 22762951, "CCDC51"), c(22762995, 22766354, "TMA7"), c(22767930, 22781328, "DNASE1L1"), c(22781582, 22784288, "TRIM35"), c(22788211, 22791494, "TRIM35"), c(22809538, 22813106, "CBLN3"), c(22817362, 22822347, "C1QC"), c(22825851, 22830064, "COL10A1"), c(22856509, 22936343, "MTSS1"), c(23000434, 23116782, "RALY"), c(23131658, 23148552, "ASIP"), c(23238420, 23248456, "AHCY"), c(23254545, 23264567, "CHMP4B"), c(23267333, 23268583, "C20orf85"), c(23295242, 23300604, "PCK1"))
plot2hgncgenes_noclutter <- list( c(20834882, 20850871, "MST1"), c(20851959, 20860716, "HMGA2"), c(20862667, 20867100, "C3orf18"), c(20868713, 20878609, "CALML5"), c(20874338, 20885219, "GNB1"), c(20886603, 20890367, "TARDBP"), c(20890796, 20896104, "APITD1-CORT"), c(20897518, 20898638, "RBP7"), c(20902260, 20909067, "H6PD"), c(20916392, 20929391, "EPHA2"), c(20931552, 20932439, "HSPB7"), c(20932855, 20935509, "RSG1"), c(20937240, 20941904, "ARHGEF19"), c(20942464, 20952607, "DDI2"), c(20960080, 20963036, "TNFRSF14"), c(20968913, 20971181, "TNFRSF14"), c(20979768, 20983303, "TNFRSF8"), c(20984263, 20989684, "ARL8B"), c(20990983, 21002192, "EDEM1"), c(21004514, 21008253, "GNAT1"), c(21011309, 21057398, "SEMA3F"), c(21110849, 21125353, "RBM6"), c(21126195, 21142379, "CAMKV"), c(21153520, 21166796, "MST1R"), c(21168897, 21201665, "UBA7"), c(21203302, 21205876, "FAM212A"), c(21207278, 21256866, "GNAI2"), c(21259492, 21289854, "SLC38A5"), c(21292979, 21303573, "RBM5"), c(21304427, 21306177, "BRK1"), c(21306980, 21308722, "GPX2"), c(21358351, 21386821, "GRM2"), c(21387925, 21395349, "RRP9"), c(21402097, 21406337, "GPR61"), c(21413051, 21415557, "CAV2"), c(21415847, 21499740, "SRGAP3"), c(21510146, 21520938, "PTGIS"), c(21521602, 21563132, "KCNB1"), c(21579834, 21588703, "STAU1"), c(21606334, 21625857, "CAMK1"), c(21626239, 21641350, "FAM3C"), c(21642142, 21677641, "C3orf67"), c(21755505, 21873152, "PTPRG"), c(21876080, 21942932, "CADPS"), c(21957430, 21967085, "SYNPR"), c(21987096, 21989317, "RHO"), c(21996068, 21999492, "RHO"), c(22019318, 22025722, "SLC6A12"), c(22099451, 22104060, "TRUB2"), c(22104184, 22108329, "HDHD3"), c(22109427, 22119679, "TGFA"), c(22122493, 22147554, "SLC2A11"), c(22155087, 22168063, "PATZ1"), c(22169498, 22185057, "EIF4ENIF1"), c(22185313, 22195622, "SFI1"), c(22195620, 22205506, "TF"), c(22207494, 22211667, "TNFSF14"), c(22212879, 22214141, "TRIM35"), c(22215012, 22224663, "GMPPB"), c(22235076, 22242643, "PDE12"), c(22243203, 22252218, "ATP6AP1L"), c(22254250, 22256232, "RPS23"), c(22258147, 22269495, "ASB14"), c(22269597, 22306245, "IP6K1"), c(22324196, 22351316, "FAM208A"), c(22390748, 22423769, "IL17RD"), c(22439258, 22454910, "APPL1"), c(22459858, 22463558, "ARF1"), c(22464287, 22469039, "ARF1"), c(22487270, 22494315, "NPTX1"), c(22513818, 22519177, "NPTX1"), c(22550919, 22633366, "SLMAP"), c(22638861, 22702632, "FLNB"), c(22703807, 22706879, "TRIM35"), c(22709578, 22712114, "TRIM35"), c(22712897, 22715744, "TRIM35"), c(22718507, 22730927, "OSGIN2"), c(22757753, 22762951, "CCDC51"), c(22781582, 22784288, "TRIM35"), c(22788211, 22791494, "TRIM35"), c(22809538, 22813106, "CBLN3"), c(22817362, 22822347, "C1QC"), c(22825851, 22830064, "COL10A1"), c(22856509, 22936343, "MTSS1"), c(23000434, 23116782, "RALY"), c(23131658, 23148552, "ASIP"), c(23238420, 23248456, "AHCY"), c(23254545, 23264567, "CHMP4B"), c(23267333, 23268583, "C20orf85"), c(23295242, 23300604, "PCK1"))

# Constants for Plot 3
plot3col <- rgb(0,0,255, max=255, alpha=175)
plot3xstart <- 29010001
plot3xstop <- 31860001
plot3xstartID <- subset(lg5, BIN_START==plot3xstart)[1,1]
plot3xstopID <- subset(lg5, BIN_START==plot3xstop)[1,1]
plot3 <- var_fst_data[var_fst_data$POS > plot3xstart,]
plot3 <- plot3[plot3$POS < plot3xstop,]
lg5$COL[lg5$BIN_START >= plot3xstart & lg5$BIN_END <= plot3xstop] <- 11
plot3[,4] <- plot3[, 4] + 9
plot3$COL <- colorRampPalette(c(my_colors_2[10], my_colors_2[11]))(10)[as.numeric(cut(plot3$FST,breaks = 10))]
plot3$alpha <- seq(29,99,7)[as.numeric(cut(plot3$FST,breaks = 10))]
plot3$COL <- paste0(plot3$COL, plot3$alpha)
plot3$COL[which(plot3$COL == "NANA")] <- "white"
plot3$COL[which(plot3$FST == 1)] <- "black"
# Genes for Plot 3
plot3genes <- list( c(29018627, 29031263, "TSPAN2"), c(29032240, 29037994, "SLC25A22"), c(29042036, 29043571, "TSHB"), c(29045224, 29051904, "SLC5A8"), c(29052226, 29053513, "LOC101470491"), c(29053615, 29062173, "SYCP1"), c(29062794, 29071727, "SLC16A1"), c(29184340, 29185334, "LOC112434878"), c(29214203, 29328256, "FAM19A1"), c(29359912, 29375773, "PPM1J"), c(29379069, 29395463, "RHOC"), c(29395536, 29413061, "MOV10"), c(29414653, 29419791, "CAPZA3"), c(29420373, 29430461, "NCOA5"), c(29430975, 29641597, "SLC12A5"), c(29531108, 29531533, "LOC112434829"), c(29652542, 29654603, "LOC112434828"), c(29663579, 29672209, "NAT8L"), c(29677012, 29683286, "NAT8L"), c(29684697, 29687054, "NAT8L"), c(29718584, 29725342, "CYP24A1"), c(29727067, 29727139, "trnat-ugu"), c(29728172, 29740385, "NFATC2"), c(29741561, 29744926, "MANBAL"), c(29745986, 29763944, "SRC"), c(29813777, 29830148, "APCDD1L"), c(29840983, 29865691, "VAPB"), c(29868540, 29883678, "RAB22A"), c(29883888, 29899549, "PPP4R1L"), c(29902421, 29903778, "BCDIN3D"), c(29904615, 29906699, "COX14"), c(29906784, 29932421, "CERS5"), c(29934473, 29938894, "GPD1"), c(29940371, 29958942, "SMARCD1"), c(29967228, 30148945, "ASIC1"), c(30174097, 30177067, "CELA1"), c(30186073, 30188748, "CELA1"), c(30198248, 30210267, "GLS2"), c(30218632, 30222930, "BLOC1S1"), c(30223430, 30236778, "ITCH"), c(30239122, 30240596, "PXMP4"), c(30262928, 30268093, "TGM5"), c(30273165, 30278130, "TGM5"), c(30278372, 30302770, "EYA2"), c(30330324, 30353577, "RTFDC1"), c(30335291, 30343819, "GCNT7"), c(30357862, 30370582, "TFAP2C"), c(30517986, 30520530, "LOC112434853"), c(30529982, 30533197, "NAT8L"), c(30535792, 30539154, "NAT8L"), c(30540495, 30547148, "NAT8L"), c(30552104, 30560761, "NAT8L"), c(30569634, 30571773, "LOC101474941"), c(30630306, 30631157, "LOC112434798"), c(30668260, 30856329, "CDH22"), c(30888951, 30890442, "LOC112434869"), c(30890450, 30891445, "LOC112434870"), c(30954415, 30968875, "ZNFX1"), c(30969535, 30976403, "TBC1D20"), c(30976960, 30984300, "RBCK1"), c(30984928, 31032677, "NECAB3"), c(31039994, 31076618, "CBFA2T2"), c(31085947, 31117388, "SNTA1"), c(31129990, 31235221, "TOX2"), c(31251324, 31268116, "ST7L"), c(31272612, 31281722, "WNT2B"), c(31287438, 31302590, "CTTNBP2NL"), c(31331517, 31430416, "KCND3"), c(31439195, 31452853, "FAM212B"), c(31487428, 31509064, "RAP1A"), c(31510178, 31515467, "TIMM17A"), c(31516891, 31522221, "LMOD1"), c(31524221, 31534161, "SHISA4"), c(31536278, 31552145, "IPO9"), c(31554567, 31592401, "GNAS"), c(31598733, 31602506, "TUBB1"), c(31604014, 31609561, "SLMO2"), c(31612394, 31619444, "AURKA"), c(31622370, 31627134, "SLC32A1"), c(31654783, 31679849, "ARHGAP40"), c(31704958, 31796462, "PPP1R16B"), c(31794119, 31795542, "LOC105940844"), c(31806412, 31817286, "FAM210B"), c(31839663, 31854111, "NDUFA4L2"))
plot3hgncgenes <- list( c(29018627, 29031263, "TSPAN2"), c(29032240, 29037994, "SLC25A22"), c(29042036, 29043571, "TSHB"), c(29045224, 29051904, "SLC5A8"), c(29053615, 29062173, "SYCP1"), c(29062794, 29071727, "SLC16A1"), c(29214203, 29328256, "FAM19A1"), c(29359912, 29375773, "PPM1J"), c(29379069, 29395463, "RHOC"), c(29395536, 29413061, "MOV10"), c(29414653, 29419791, "CAPZA3"), c(29420373, 29430461, "NCOA5"), c(29430975, 29641597, "SLC12A5"), c(29663579, 29672209, "NAT8L"), c(29677012, 29683286, "NAT8L"), c(29684697, 29687054, "NAT8L"), c(29718584, 29725342, "CYP24A1"), c(29728172, 29740385, "NFATC2"), c(29741561, 29744926, "MANBAL"), c(29745986, 29763944, "SRC"), c(29813777, 29830148, "APCDD1L"), c(29840983, 29865691, "VAPB"), c(29868540, 29883678, "RAB22A"), c(29883888, 29899549, "PPP4R1L"), c(29902421, 29903778, "BCDIN3D"), c(29904615, 29906699, "COX14"), c(29906784, 29932421, "CERS5"), c(29934473, 29938894, "GPD1"), c(29940371, 29958942, "SMARCD1"), c(29967228, 30148945, "ASIC1"), c(30174097, 30177067, "CELA1"), c(30186073, 30188748, "CELA1"), c(30198248, 30210267, "GLS2"), c(30218632, 30222930, "BLOC1S1"), c(30223430, 30236778, "ITCH"), c(30239122, 30240596, "PXMP4"), c(30262928, 30268093, "TGM5"), c(30273165, 30278130, "TGM5"), c(30278372, 30302770, "EYA2"), c(30330324, 30353577, "RTFDC1"), c(30335291, 30343819, "GCNT7"), c(30357862, 30370582, "TFAP2C"), c(30529982, 30533197, "NAT8L"), c(30535792, 30539154, "NAT8L"), c(30540495, 30547148, "NAT8L"), c(30552104, 30560761, "NAT8L"), c(30668260, 30856329, "CDH22"), c(30954415, 30968875, "ZNFX1"), c(30969535, 30976403, "TBC1D20"), c(30976960, 30984300, "RBCK1"), c(30984928, 31032677, "NECAB3"), c(31039994, 31076618, "CBFA2T2"), c(31085947, 31117388, "SNTA1"), c(31129990, 31235221, "TOX2"), c(31251324, 31268116, "ST7L"), c(31272612, 31281722, "WNT2B"), c(31287438, 31302590, "CTTNBP2NL"), c(31331517, 31430416, "KCND3"), c(31439195, 31452853, "FAM212B"), c(31487428, 31509064, "RAP1A"), c(31510178, 31515467, "TIMM17A"), c(31516891, 31522221, "LMOD1"), c(31524221, 31534161, "SHISA4"), c(31536278, 31552145, "IPO9"), c(31554567, 31592401, "GNAS"), c(31598733, 31602506, "TUBB1"), c(31604014, 31609561, "SLMO2"), c(31612394, 31619444, "AURKA"), c(31622370, 31627134, "SLC32A1"), c(31654783, 31679849, "ARHGAP40"), c(31704958, 31796462, "PPP1R16B"), c(31806412, 31817286, "FAM210B"), c(31839663, 31854111, "NDUFA4L2"))
plot3hgncgenes_noclutter <- list( c(29018627, 29031263, "TSPAN2"), c(29032240, 29037994, "SLC25A22"), c(29042036, 29043571, "TSHB"), c(29045224, 29051904, "SLC5A8"), c(29053615, 29062173, "SYCP1"), c(29062794, 29071727, "SLC16A1"), c(29214203, 29328256, "FAM19A1"), c(29359912, 29375773, "PPM1J"), c(29379069, 29395463, "RHOC"), c(29395536, 29413061, "MOV10"), c(29414653, 29419791, "CAPZA3"), c(29420373, 29430461, "NCOA5"), c(29430975, 29641597, "SLC12A5"), c(29663579, 29672209, "NAT8L"), c(29718584, 29725342, "CYP24A1"), c(29728172, 29740385, "NFATC2"), c(29741561, 29744926, "MANBAL"), c(29745986, 29763944, "SRC"), c(29813777, 29830148, "APCDD1L"), c(29840983, 29865691, "VAPB"), c(29868540, 29883678, "RAB22A"), c(29883888, 29899549, "PPP4R1L"), c(29902421, 29903778, "BCDIN3D"), c(29904615, 29906699, "COX14"), c(29906784, 29932421, "CERS5"), c(29934473, 29938894, "GPD1"), c(29940371, 29958942, "SMARCD1"), c(29967228, 30148945, "ASIC1"), c(30174097, 30177067, "CELA1"), c(30186073, 30188748, "CELA1"), c(30198248, 30210267, "GLS2"), c(30218632, 30222930, "BLOC1S1"), c(30223430, 30236778, "ITCH"), c(30239122, 30240596, "PXMP4"), c(30262928, 30268093, "TGM5"), c(30278372, 30302770, "EYA2"), c(30330324, 30353577, "RTFDC1"), c(30335291, 30343819, "GCNT7"), c(30357862, 30370582, "TFAP2C"), c(30529982, 30533197, "NAT8L"), c(30668260, 30856329, "CDH22"), c(30954415, 30968875, "ZNFX1"), c(30969535, 30976403, "TBC1D20"), c(30976960, 30984300, "RBCK1"), c(30984928, 31032677, "NECAB3"), c(31039994, 31076618, "CBFA2T2"), c(31085947, 31117388, "SNTA1"), c(31129990, 31235221, "TOX2"), c(31251324, 31268116, "ST7L"), c(31272612, 31281722, "WNT2B"), c(31287438, 31302590, "CTTNBP2NL"), c(31331517, 31430416, "KCND3"), c(31487428, 31509064, "RAP1A"), c(31510178, 31515467, "TIMM17A"), c(31516891, 31522221, "LMOD1"), c(31524221, 31534161, "SHISA4"), c(31536278, 31552145, "IPO9"), c(31554567, 31592401, "GNAS"), c(31598733, 31602506, "TUBB1"), c(31604014, 31609561, "SLMO2"), c(31612394, 31619444, "AURKA"), c(31622370, 31627134, "SLC32A1"), c(31654783, 31679849, "ARHGAP40"), c(31704958, 31796462, "PPP1R16B"), c(31806412, 31817286, "FAM210B"), c(31839663, 31854111, "NDUFA4L2"))

lg5$COL[lg5$Zfst > 5] <- 3

# Actual Plotting ============================================================

# Plot LG5
png_name <- "top3LG5.png"
Cairo(png_name, type="png", , width = 15, height = 7, units = 'in', res=300)
par(fig=c(0,1,0,0.7))
plot(lg5$ID, lg5$Zfst, col=lg5$COL, xlab="Position on LG5 (Mbp)", ylab=expression("Z"[fst]), xaxs="i",xaxt="n", yaxt="n", ylim=c(0, lg5ystop), main="", pch=20, bty="n" )
# Label the important points
text(high[,1:2], labels=paste0(format(round(high$BIN_START/1000000, 2)), " Mbp"), pos=4, cex= 0.75)
#title(main = expression("LG5 UMD2a Bicuspid vs Tricuspid Z"[fst]))
axis(2, at=c(0,3,6,9), label=c(0,3,6,9))
axis(1, at=c(13888,14088, 14288, 14488, 14688,14888,15088, 15288, 15488, 15688, 15886,16086,16286,16486,16686,16886,17086,17286, 17501), label=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36.17))
text(lg5xstart+100, 6, pos=3, labels=c("10 kb bins"), cex=0.8)

# Line segments to illustrate that the new plots are "blown up"
sep <- 1.25
# ## Segments for Plot 1
# segments(plot1xstartID, 6, lg5xstart + 0, lg5ystop-sep, lwd=2)
# segments(plot1xstopID, 6, lg5xstart + (lg5Range/3) - bufferImportantPoints, lg5ystop-sep, lwd=2)
# ## Segments for Plot 2
# segments(plot2xstartID, 6, lg5xstart + (lg5Range/3) + bufferImportantPoints, lg5ystop-sep, lwd=2)
# segments(plot2xstopID, 6, lg5xstart + 2*(lg5Range/3) - bufferImportantPoints, lg5ystop-sep, lwd=2)
# ## Segments for Plot 3
# segments(plot3xstartID, 6, lg5xstart + 2*(lg5Range/3) + bufferImportantPoints, lg5ystop-sep, lwd=2)
# segments(plot3xstopID, 6, lg5xstop, lg5ystop-sep, lwd=2)

# PLOT 1
par(fig=c(0,0.38,xplotsystart,1), new=TRUE)
plot(plot1$POS, plot1$FST, col = plot1$COL, xaxs="i", ylim=c(0,subplotstop), xlab="", ylab="Fst", xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
axis(2, at=c(0,0.5,1), label=c(0,"",1), cex =xplotcex)
axis(1, at=c(5000000,6000000,7000000), label=c(5,6,7), cex = xplotcex)
axis(1, at=c(4140163,7439752), label=c(4.1, 7.4), col.axis="black", cex = xplotcex)
axis(1, at=c(5390000,6190000), label=c("5.39","6.19"), cex = xplotcex, padj = -1)
# genes
previousgenestop <- 0
labelbuffer <- 0.1
cummulativelabelbuffer <- 0
for (gene in plot1hgncgenes) {
  genestart <- gene[1]
  genestop  <- gene[2]
  name      <- gene[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(genesystart, genesystop, genesystop, genesystart)
  polygon(x, y, col=plot1col, border=NA)
  if (as.numeric(genestart) - previousgenestop < 10000) {
    if (cummulativelabelbuffer == 0) { cummulativelabelbuffer<- cummulativelabelbuffer + labelbuffer }
    cummulativelabelbuffer <- cummulativelabelbuffer + labelbuffer
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=0, cex=geneslabelcex) 
  }
  else {
    cummulativelabelbuffer <- 0
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=55, cex=geneslabelcex) 
  }
  previousgenestop <- as.numeric(genestop)
}

# PLOT 2
par(fig=c(0.31,0.688,xplotsystart,1), new=TRUE)
plot(plot2$POS, plot2$FST, col = plot2$COL, xaxs="i", ylim=c(0,subplotstop), xlab="", ylab="", xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
axis(1, at=c(21000000,22000000,23000000), label=c(21,22,23), cex = xplotcex)
axis(1, at=c(20830003,23329917), label=c("",""), col.axis="black", cex = xplotcex)
axis(1, at=c(22080000), label=c("22.08"), cex = xplotcex, padj = -1)
# genes
previousgenestop <- 0
labelbuffer <- 0.1
cummulativelabelbuffer <- 0
for (gene in plot2hgncgenes_noclutter) {
  genestart <- gene[1]
  genestop  <- gene[2]
  name      <- gene[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(genesystart, genesystop, genesystop, genesystart)
  polygon(x, y, col=plot2col, border=NA)
  if (as.numeric(genestart) - previousgenestop < 10000) {
    if (cummulativelabelbuffer == 0) { cummulativelabelbuffer<- cummulativelabelbuffer + labelbuffer }
    cummulativelabelbuffer <- cummulativelabelbuffer + labelbuffer
    if (genestart < 20900000) { text(as.numeric(genestart)+45000 + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=0, cex=geneslabelcex)  }
    else { text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=0, cex=geneslabelcex) }
    
  }
  else {
    cummulativelabelbuffer <- 0
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=55, cex=geneslabelcex) 
  }
  previousgenestop <- as.numeric(genestop)
}

# PLOT 3
par(fig=c(0.618,1,xplotsystart,1), new=TRUE)
plot(plot3$POS, plot3$FST, col = plot3$COL, xaxs="i", ylim=c(0,subplotstop), xlab="", ylab="", xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
axis(1, at=c(30000000,31000000), label=c(30,31), cex=xplotcex)
axis(1, at=c(29010050,31859939), label=c(29, 31.9), col.axis="black", cex = xplotcex)
axis(1, at=c(30260000,30610000), label=c("30.26","30.61"), cex = xplotcex, padj = -1)
axis(1, at=30610000, label="30.61", cex = xplotcex, padj = -1)
# genes
previousgenestop <- 0
labelbuffer <- 0.1
cummulativelabelbuffer <- 0
for (gene in plot3hgncgenes_noclutter) {
  genestart <- gene[1]
  genestop  <- gene[2]
  name      <- gene[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(genesystart, genesystop, genesystop, genesystart)
  polygon(x, y, col=plot3col, border=NA)
  if (as.numeric(genestart) - previousgenestop < 10000) {
    if (cummulativelabelbuffer == 0) { cummulativelabelbuffer<- cummulativelabelbuffer + labelbuffer }
    cummulativelabelbuffer <- cummulativelabelbuffer + labelbuffer
    
    if (genestart < 29110030) { text(as.numeric(genestart)+2000 + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=0, cex=geneslabelcex) }
    else { text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=0, cex=geneslabelcex) }
  }
  else {
    cummulativelabelbuffer <- 0
    text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=55, cex=geneslabelcex) 
  }
  previousgenestop <- as.numeric(genestop)
}

dev.off()

##############################
# Top 3 for Every Chromosome #
##############################
# Helper Functions #
creatSubPlot = function(plot_data, plot_data_start, plot_data_end, plot_genes, plot_col, plot_num) {
  my_ylab=""
  if (plot_num==1) { my_ylab="Fst" }
  plot(plot_data$POS, plot_data$FST, col = plot_data$COL, xaxs="i", ylim=c(0,subplotstop), xlim = c(plot_data_start, plot_data_end), xlab="", ylab=my_ylab, xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
  if (plot_num==1) { axis(2, at=c(0,0.5,1), label=c(0,"",1), cex =xplotcex) }
  axis(1, at=seq(plot_data_start+1, plot_data_end-1, by=zoombuffer), label=format(round(seq(plot_data_start, plot_data_end, by=zoombuffer)/1e6, 2), nsmall = 2), cex = xplotcex)
  # genes
  previousx <- 0
  # One loop for plotting the boxes and one loop for plotting names.
  # If both were done at the same time, the boxes would overlap the names.
  for (row in 1:nrow(plot_genes)) {
    genestart <- plot_genes$start[row]
    genestop  <- plot_genes$stop[row]
    genelength = genestop - genestart
    name      <- plot_genes$gene_name[row]
    x <- as.numeric(c(genestart, genestart, genestop, genestop))
    y <- as.numeric(c(genesystart, genesystop, genesystop, genesystart))
    polygon(x, y, col=plot_col, border=NA)
  }
  for (row in 1:nrow(plot_genes)) {
    genestart <- plot_genes$start[row]
    genestop  <- plot_genes$stop[row]
    genelength = genestop - genestart
    name      <- plot_genes$gene_name[row]
    thisx = as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2
    if (name == "bhlhe40" || thisx > previousx + 20000 || thisx < previousx-20000 ) {
      thisx = as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2
      text(thisx,as.numeric(geneslabely), labels = name, srt=55, cex=geneslabelcex) 
      previousx <- thisx
    }
  }
}

library(Cairo)
# *** Specify the path to Fst File ***
path_to_fst = "C:/Users/miles/Downloads/d_tooth/data/whole_filter_rm_window_08_31_20.fst"
path_to_var = "C:/Users/miles/Downloads/d_tooth/data/whole_filter_rm_no_window_08_31_20.fst"
# path_to_fst = "C:/Users/miles/Downloads/d_tooth/data/no_1615_1621_1984_windowed.windowed.weir.fst"
# path_to_var = "C:/Users/miles/Downloads/d_tooth/data/no_1615_1621_1984_no_window.weir.fst"

# Load the gtf file
# gtf = read.table("C:/Users/miles/Downloads/d_tooth/data/genes.gff", sep="\t", header=F, stringsAsFactors = F)
# gtf = gtf[which(gtf[,3] == "gene"),]
gtf = read.table("C:/Users/miles/Downloads/brain/brain_scripts/full_ens_w_ncbi_gene_NC_sort.gtf", sep="\t", header=F, stringsAsFactors = F)
gtf = gtf[which(gtf[,3] == "gene" & gtf[,1] != "NC_027944.1"),]
gtf_gene_name <- c()
for (i in 1:nrow(gtf)) {
  start <- gregexpr(pattern ='gene_name', gtf$V9[i])[[1]]
  stop  <- gregexpr(pattern =';', substr(gtf$V9[i], start, nchar(gtf$V9[i])))[[1]][1]
  gene_name <- substr(gtf$V9[i], start+10, start+stop-2)
  if (start == -1) {
    gene_name <- substr(gtf$V9[i], start+10, start+stop)
  } else if ( grepl('gene_biotype lncRNA', gtf$V9[i], fixed = TRUE) ) {
    gene_name <- "lncRNA"
  } else if ( grepl('gene_biotype "miRNA"', gtf$V9[i], fixed = TRUE) ) {
    gene_name <- "miRNA"
  }
  gtf_gene_name <- c(gtf_gene_name, gene_name)
}
gtf$gene_name <- gtf_gene_name
colnames(gtf) <- c("LG", "source", "type", "start", "stop", "idk", "idk1", "idk2", "info", "gene_name")
gtf = gtf[which(! startsWith(gtf$gene_name, "ENSMZEG")),]
gene_gtf = gtf

# Load the binned data
fst_data <- read.table(path_to_fst, sep="\t", header=TRUE)
# Convert the negative fst scores to 0
fst_data[, c(5)][fst_data[, c(5)] < 0] <- 0
# Add a column that just has row number
fst_data$ID <- seq.int(nrow(fst_data))
# Add a column for Zfst
fst_data$Zfst <- ( fst_data$WEIGHTED_FST - mean(fst_data$WEIGHTED_FST) ) / sd(fst_data$WEIGHTED_FST)
fst_data <- fst_data[,c("ID","Zfst", "CHROM", "BIN_START","BIN_END", "WEIGHTED_FST")]

# Load the Unbinned data
var_fst_data <- read.table(path_to_var, sep="\t", header=TRUE)
# Rename the columns
names(var_fst_data) <- c("LG", "POS", "FST")
# Convert the negative fst scores to 0
var_fst_data[, c(3)][var_fst_data[, c(3)] < 0] <- 0

# Set a Color Palette
my_colors <- c("grey68", "grey50", "grey0", "gold2", "gold4", "grey0", "chartreuse2", "chartreuse4", "grey0", "cadetblue2", "cadetblue4", "grey0")
my_colors_2 <- c("grey68", "grey50", "grey0", "gold2", "gold4", "grey0", "chartreuse2", "green4", "grey0", "cadetblue2", "blue4", "grey0")
palette(my_colors)

# Constants for all plots
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
zoombuffer = 1e6

lgs = unique(gtf$LG)
lgs_name = c(1:20, 22, 23)
# lgs = c("NC_036784.1")
# lgs_name = c(5)
for (i in 1:length(lgs)) {
  # Load data for this LG
  print(i)
  lg = lgs[i]
  lg_name = lgs_name[i]
  binned = fst_data[fst_data$CHROM == lg,]
  unbinned <- var_fst_data[var_fst_data$LG == lg,]
  this_genes = gene_gtf[which(gene_gtf$LG == lg),]
  
  binned$COL <- binned$Zfst
  binned$COL[binned$Zfst > 2.5] <- 2
  binned$COL[binned$Zfst < 2.5] <- 1
  
  # Identify important high Zfst points
  # high <- subset(binned, binned$Zfst > 5)
  # high$Zfst = format(round(high$Zfst, 2), nsmall = 2)
  
  print("Finding Top 3 Regions")
  # Identify Top 3 Regions Using the Following Criteria
  # 3. # bins > 2.5 Zfst score in 10 10kb bins
  top_3 = c()
  window_size = 20
  for (j in 1:nrow(binned)) {
    if (j > window_size) {
      this_zfst = binned$Zfst[(j-window_size):j]
      top_3 = c(top_3, length(this_zfst[which(this_zfst > 2.5)]))
      # top_3 = c(top_3, sum(this_zfst))
    } else {
      top_3 = c(top_3, NA)
    }
  }
  top_3[1:window_size] = NA
  plot(binned$BIN_START, top_3, xaxt="n", xlab=paste("Position on LG", lg_name, "(Mbp)"), ylab=paste("Number of 10kb Bins out of", window_size, "w/ Zfst > 2.5"))
  axis(1, at=(seq(0,max(binned$BIN_START), by=1e6)), label=seq(0, max(binned$BIN_START)/1e6))
  # top_3 = binned$Zfst
  top_3_i = c()
  top_3_order_i = sort(top_3, decreasing = T, index.return = T)$ix
  for (j in 1:length(top_3_order_i)) {
    if (length(top_3_i) == 0 || length(top_3_i[which(top_3_order_i[j] < top_3_i-(zoombuffer/10000) | top_3_order_i[j] > top_3_i+(zoombuffer/10000))]) == length(top_3_i)) {
      top_3_i = c(top_3_i, top_3_order_i[j])
    }
  }
  head(top_3_i, 20)
  top_3_i = sort(top_3_i[1:3])
  top_3_i
  
  # Prepare the Data for the top 3 zoom-ins
  plot1_start = binned$BIN_START[top_3_i[1]] - zoombuffer
  plot1_end   = binned$BIN_END  [top_3_i[1]] + zoombuffer
  plot2_start = binned$BIN_START[top_3_i[2]] - zoombuffer
  plot2_end   = binned$BIN_END  [top_3_i[2]] + zoombuffer
  plot3_start = binned$BIN_START[top_3_i[3]] - zoombuffer
  plot3_end   = binned$BIN_END  [top_3_i[3]] + zoombuffer
  if (plot1_start < 0)                          { plot1_start = 0 }
  if (plot1_end > binned$BIN_END[nrow(binned)]) { plot1_end   = binned$BIN_END[nrow(binned)] }
  if (plot2_start < 0)                          { plot2_start = 0 }
  if (plot2_end > binned$BIN_END[nrow(binned)]) { plot2_end   = binned$BIN_END[nrow(binned)] }
  if (plot3_start < 0)                          { plot3_start = 0 }
  if (plot3_end > binned$BIN_END[nrow(binned)]) { plot3_end   = binned$BIN_END[nrow(binned)] }
  
  plot1 = unbinned[which(unbinned$POS > plot1_start & unbinned$POS < plot1_end),]
  plot2 = unbinned[which(unbinned$POS > plot2_start & unbinned$POS < plot2_end),]
  plot3 = unbinned[which(unbinned$POS > plot3_start & unbinned$POS < plot3_end),]
  
  plot1genes = gene_gtf[which(gene_gtf$LG == lg & gene_gtf$start >= plot1_start & gene_gtf$stop <= plot1_end),]
  plot2genes = gene_gtf[which(gene_gtf$LG == lg & gene_gtf$start >= plot2_start & gene_gtf$stop <= plot2_end),]
  plot3genes = gene_gtf[which(gene_gtf$LG == lg & gene_gtf$start >= plot3_start & gene_gtf$stop <= plot3_end),]
  
  binned$COL[which(binned$BIN_START >= plot1_start & binned$BIN_END <= plot1_end)] = 5
  binned$COL[which(binned$BIN_START >= plot2_start & binned$BIN_END <= plot2_end)] = 8
  binned$COL[which(binned$BIN_START >= plot3_start & binned$BIN_END <= plot3_end)] = 11
  
  plot1col <- rgb(255, 255, 0, max=255, alpha=175)
  plot1$COL <- colorRampPalette(c(my_colors_2[4], my_colors_2[5]))(10)[as.numeric(cut(plot1$FST,breaks = 10))]
  plot1$alpha <- seq(29,99,7)[as.numeric(cut(plot1$FST,breaks = 10))]
  plot1$COL <- paste0(plot1$COL, plot1$alpha)
  plot1$COL[which(plot1$COL == "NANA")] <- "white"
  plot1$COL[which(plot1$FST == 1)] <- "black"

  plot2col <- rgb(0, 255, 0, max=255, alpha=175)
  plot2$COL <- colorRampPalette(c(my_colors_2[7], my_colors_2[8]))(10)[as.numeric(cut(plot2$FST,breaks = 10))]
  plot2$alpha <- seq(29,99,7)[as.numeric(cut(plot2$FST,breaks = 10))]
  plot2$COL <- paste0(plot2$COL, plot2$alpha)
  plot2$COL[which(plot2$COL == "NANA")] <- "white"
  plot2$COL[which(plot2$FST == 1)] <- "black"
  
  plot3col <- rgb(0,0,255, max=255, alpha=175)
  plot3$COL <- colorRampPalette(c(my_colors_2[10], my_colors_2[11]))(10)[as.numeric(cut(plot3$FST,breaks = 10))]
  plot3$alpha <- seq(29,99,7)[as.numeric(cut(plot3$FST,breaks = 10))]
  plot3$COL <- paste0(plot3$COL, plot3$alpha)
  plot3$COL[which(plot3$COL == "NANA")] <- "white"
  plot3$COL[which(plot3$FST == 1)] <- "black"
  
  # Actual Plotting ============================================================
  # Binned Data
  png_name <- paste0("C:/Users/miles/Downloads/d_tooth/results/ts_ms/top3_", lg_name, ".png")
  Cairo(png_name, type="png", , width = 15, height = 7, units = 'in', res=300)
  par(fig=c(0,1,0,0.6))
  plot(binned$BIN_START, binned$Zfst, col=binned$COL, xlab=paste("Position on LG", lg_name, "(Mbp)"), ylab=expression("Z"[fst]), xaxs="i",xaxt="n", main="", pch=20, bty="n" )
  # Label the important points
  # text(x=high$BIN_START, y=high$Zfst, labels=paste0(format(round(high$BIN_START/1000000, 2)), " Mbp"), pos=4, cex= 0.75)
  #title(main = expression("LG5 UMD2a Bicuspid vs Tricuspid Z"[fst]))
  axis(1, at=(seq(0,max(binned$BIN_START), by=1e6)), label=seq(0, max(binned$BIN_START)/1e6))
  text(3e6, 6, pos=3, labels=c("10 kb bins"), cex=0.8)
  
  # Plot 1
  par(fig=c(0,0.36,xplotsystart,1), new=TRUE)
  creatSubPlot(plot_data = plot1, plot_data_start = plot1_start, plot_data_end = plot1_end, plot_genes = plot1genes, plot_col = plot1col, plot_num = 1)
  
  # Plot 2
  par(fig=c(0.31,0.688,xplotsystart,1), new=TRUE)
  creatSubPlot(plot_data = plot2, plot_data_start = plot2_start, plot_data_end = plot2_end, plot_genes = plot2genes, plot_col = plot2col, plot_num = 2)
  
  # Plot 3
  par(fig=c(0.64,1,xplotsystart,1), new=TRUE)
  creatSubPlot(plot_data = plot3, plot_data_start = plot3_start, plot_data_end = plot3_end, plot_genes = plot3genes, plot_col = plot3col, plot_num = 3)
  
  dev.off()
}