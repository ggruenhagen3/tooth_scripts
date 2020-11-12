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
my_colors = c("grey68", "grey50", "grey0", "mediumpurple1", "mediumpurple4", "grey0", "cadetblue2", "cadetblue4", "grey0", "maroon2", "maroon4", "grey0")
my_colors_2 = c("grey68", "grey50", "grey0", "mediumpurple1", "purple4", "grey0", "cadetblue2", "blue4", "grey0", "maroon1", "maroon4", "grey0")
palette(my_colors)
lg5$COL <- lg5$Zfst
lg5$COL[lg5$Zfst > 2.5] <- 2
lg5$COL[lg5$Zfst < 2.5] <- 1
var_fst_data$COL <- var_fst_data$FST
var_fst_data$COL[var_fst_data$FST > .99] <- 3
var_fst_data$COL[var_fst_data$FST < .99] <- 2
var_fst_data$COL[var_fst_data$FST < .50] <- 1


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
plot1col <- rgb(148, 69, 211, max=255, alpha=175)
plot1xstart <- 1
plot1xstop <- 3000001
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
# lol the CantFind gene is a protein coding gene and it's id is LOC112434809
plot1genes <- list( c(16846, 20191, "SCRT2"), c(30757, 33343, "TCF15"), c(44609, 59968, "CSNK2A1"), c(60181, 67435, "SLC2A10"), c(68107, 69721, "TP53RK"), c(70584, 79723, "TAPBPL"), c(81885, 98476, "SERINC3"), c(120523, 121935, "P2RY11"), c(130366, 144083, "STK35"), c(146097, 170886, "HELZ2"), c(319546, 324332, "FEZF2"), c(375886, 715617, "CADPS"), c(721736, 728136, "PHF7"), c(840383, 844854, "CantFind"), c(902798, 975585, "PLEKHA6"), c(991848, 999932, "APOBEC2"), c(1053014, 1225986, "NFASC"), c(1232575, 1235405, "ZFP41_HUMAN"), c(1236158, 1251098, "ZC3H11A"), c(1251792, 1294161, "SYPL2"), c(1305662, 1320213, "ATXN7L2"), c(1324400, 1328932, "CYB561D1"), c(1332427, 1349317, "ZBTB40"), c(1351057, 1378840, "USP48"), c(1388365, 1408067, "CAMK1G"), c(1423264, 1442012, "LAMB4"), c(1495541, 1515853, "LDLRAD2"), c(1554629, 1558345, "ELF3"), c(1613593, 1765287, "FHIT"), c(2029118, 2035792, "UBE2V2"), c(2049261, 2145746, "PLXNB1"), c(2149850, 2161566, "PTRH2"), c(2150052, 2153154, "PGBD4"), c(2163281, 2192606, "SETD5"), c(2200437, 2228961, "LHFPL4"), c(2234634, 2255017, "MTMR14"), c(2259236, 2463550, "CPNE9"), c(2466502, 2475280, "ZXDB"), c(2488301, 2499839, "KLF15"), c(2580508, 2609617, "ALDH1L1"), c(2588083, 2595849, "ZNF853"), c(2616813, 2619522, "FAM217B"), c(2621499, 2624628, "PDRG1"), c(2674969, 2696498, "EDN3"), c(2704560, 2723384, "ZNF831"), c(2728731, 2738977, "CSTF1"), c(2739768, 2757511, "CASS4"), c(2761465, 2788475, "FAM120A"), c(2790881, 2831947, "PHF2"), c(2901347, 2904564, "BARX1"), c(2922955, 2952074, "PTPDC1"), c(2962231, 2978168, "WNT7A"))
plot1hgncgenes <- list( c(16846, 20191, "SCRT2"), c(30757, 33343, "TCF15"), c(44609, 59968, "CSNK2A1"), c(60181, 67435, "SLC2A10"), c(68107, 69721, "TP53RK"), c(70584, 79723, "TAPBPL"), c(81885, 98476, "SERINC3"), c(120523, 121935, "P2RY11"), c(130366, 144083, "STK35"),c(319546, 324332, "FEZF2"), c(375886, 715617, "CADPS"), c(721736, 728136, "PHF7"), c(840383, 844854, "CantFind"), c(902798, 975585, "PLEKHA6"), c(991848, 999932, "APOBEC2"), c(1053014, 1225986, "NFASC"), c(1232575, 1235405, "ZFP41_HUMAN"), c(1236158, 1251098, "ZC3H11A"), c(1251792, 1294161, "SYPL2"), c(1305662, 1320213, "ATXN7L2"), c(1324400, 1328932, "CYB561D1"), c(1332427, 1349317, "ZBTB40"), c(1351057, 1378840, "USP48"), c(1388365, 1408067, "CAMK1G"), c(1423264, 1442012, "LAMB4"), c(1495541, 1515853, "LDLRAD2"), c(1554629, 1558345, "ELF3"), c(1613593, 1765287, "FHIT"), c(2029118, 2035792, "UBE2V2"), c(2049261, 2145746, "PLXNB1"), c(2149850, 2161566, "PTRH2"), c(2150052, 2153154, "PGBD4"), c(2163281, 2192606, "SETD5"), c(2200437, 2228961, "LHFPL4"), c(2234634, 2255017, "MTMR14"), c(2259236, 2463550, "CPNE9"), c(2466502, 2475280, "ZXDB"), c(2488301, 2499839, "KLF15"), c(2580508, 2609617, "ALDH1L1"), c(2616813, 2619522, "FAM217B"), c(2621499, 2624628, "PDRG1"), c(2674969, 2696498, "EDN3"), c(2704560, 2723384, "ZNF831"), c(2728731, 2738977, "CSTF1"), c(2739768, 2757511, "CASS4"), c(2761465, 2788475, "FAM120A"), c(2790881, 2831947, "PHF2"), c(2901347, 2904564, "BARX1"), c(2922955, 2952074, "PTPDC1"), c(2962231, 2978168, "WNT7A"))

# Constants for Plot 2
plot2col <- rgb(0,0,255, max=255, alpha=175)
plot2xstart <- 29010001
plot2xstop <- 31860001
plot2xstartID <- subset(lg5, BIN_START==plot2xstart)[1,1]
plot2xstopID <- subset(lg5, BIN_START==plot2xstop)[1,1]
plot2 <- var_fst_data[var_fst_data$POS > plot2xstart,]
plot2 <- plot2[plot2$POS < plot2xstop,]
lg5$COL[lg5$BIN_START >= plot2xstart & lg5$BIN_END <= plot2xstop] <- 7
plot2[,4] <- plot2[, 4] + 6
plot2$COL <- colorRampPalette(c(my_colors_2[7], my_colors_2[8]))(10)[as.numeric(cut(plot2$FST,breaks = 10))]
plot2$alpha <- seq(29,99,7)[as.numeric(cut(plot2$FST,breaks = 10))]
plot2$COL <- paste0(plot2$COL, plot2$alpha)
plot2$COL[which(plot2$COL == "NANA")] <- "white"
plot2$COL[which(plot2$FST == 1)] <- "black"
# Genes for Plot 2
plot2genes <- list( c(29018627, 29031263, "TSPAN2"), c(29032240, 29037994, "SLC25A22"), c(29042036, 29043571, "TSHB"), c(29045224, 29051904, "SLC5A8"), c(29052226, 29053513, "LOC101470491"), c(29053615, 29062173, "SYCP1"), c(29062794, 29071727, "SLC16A1"), c(29184340, 29185334, "LOC112434878"), c(29214203, 29328256, "FAM19A1"), c(29359912, 29375773, "PPM1J"), c(29379069, 29395463, "RHOC"), c(29395536, 29413061, "MOV10"), c(29414653, 29419791, "CAPZA3"), c(29420373, 29430461, "NCOA5"), c(29430975, 29641597, "SLC12A5"), c(29531108, 29531533, "LOC112434829"), c(29652542, 29654603, "LOC112434828"), c(29663579, 29672209, "NAT8L"), c(29677012, 29683286, "NAT8L"), c(29684697, 29687054, "NAT8L"), c(29718584, 29725342, "CYP24A1"), c(29727067, 29727139, "trnat-ugu"), c(29728172, 29740385, "NFATC2"), c(29741561, 29744926, "MANBAL"), c(29745986, 29763944, "SRC"), c(29813777, 29830148, "APCDD1L"), c(29840983, 29865691, "VAPB"), c(29868540, 29883678, "RAB22A"), c(29883888, 29899549, "PPP4R1L"), c(29902421, 29903778, "BCDIN3D"), c(29904615, 29906699, "COX14"), c(29906784, 29932421, "CERS5"), c(29934473, 29938894, "GPD1"), c(29940371, 29958942, "SMARCD1"), c(29967228, 30148945, "ASIC1"), c(30174097, 30177067, "CELA1"), c(30186073, 30188748, "CELA1"), c(30198248, 30210267, "GLS2"), c(30218632, 30222930, "BLOC1S1"), c(30223430, 30236778, "ITCH"), c(30239122, 30240596, "PXMP4"), c(30262928, 30268093, "TGM5"), c(30273165, 30278130, "TGM5"), c(30278372, 30302770, "EYA2"), c(30330324, 30353577, "RTFDC1"), c(30335291, 30343819, "GCNT7"), c(30357862, 30370582, "TFAP2C"), c(30517986, 30520530, "LOC112434853"), c(30529982, 30533197, "NAT8L"), c(30535792, 30539154, "NAT8L"), c(30540495, 30547148, "NAT8L"), c(30552104, 30560761, "NAT8L"), c(30569634, 30571773, "LOC101474941"), c(30630306, 30631157, "LOC112434798"), c(30668260, 30856329, "CDH22"), c(30888951, 30890442, "LOC112434869"), c(30890450, 30891445, "LOC112434870"), c(30954415, 30968875, "ZNFX1"), c(30969535, 30976403, "TBC1D20"), c(30976960, 30984300, "RBCK1"), c(30984928, 31032677, "NECAB3"), c(31039994, 31076618, "CBFA2T2"), c(31085947, 31117388, "SNTA1"), c(31129990, 31235221, "TOX2"), c(31251324, 31268116, "ST7L"), c(31272612, 31281722, "WNT2B"), c(31287438, 31302590, "CTTNBP2NL"), c(31331517, 31430416, "KCND3"), c(31439195, 31452853, "FAM212B"), c(31487428, 31509064, "RAP1A"), c(31510178, 31515467, "TIMM17A"), c(31516891, 31522221, "LMOD1"), c(31524221, 31534161, "SHISA4"), c(31536278, 31552145, "IPO9"), c(31554567, 31592401, "GNAS"), c(31598733, 31602506, "TUBB1"), c(31604014, 31609561, "SLMO2"), c(31612394, 31619444, "AURKA"), c(31622370, 31627134, "SLC32A1"), c(31654783, 31679849, "ARHGAP40"), c(31704958, 31796462, "PPP1R16B"), c(31794119, 31795542, "LOC105940844"), c(31806412, 31817286, "FAM210B"), c(31839663, 31854111, "NDUFA4L2"))
plot2hgncgenes <- list( c(29018627, 29031263, "TSPAN2"), c(29032240, 29037994, "SLC25A22"), c(29042036, 29043571, "TSHB"), c(29045224, 29051904, "SLC5A8"), c(29053615, 29062173, "SYCP1"), c(29062794, 29071727, "SLC16A1"), c(29214203, 29328256, "FAM19A1"), c(29359912, 29375773, "PPM1J"), c(29379069, 29395463, "RHOC"), c(29395536, 29413061, "MOV10"), c(29414653, 29419791, "CAPZA3"), c(29420373, 29430461, "NCOA5"), c(29430975, 29641597, "SLC12A5"), c(29663579, 29672209, "NAT8L"), c(29677012, 29683286, "NAT8L"), c(29684697, 29687054, "NAT8L"), c(29718584, 29725342, "CYP24A1"), c(29728172, 29740385, "NFATC2"), c(29741561, 29744926, "MANBAL"), c(29745986, 29763944, "SRC"), c(29813777, 29830148, "APCDD1L"), c(29840983, 29865691, "VAPB"), c(29868540, 29883678, "RAB22A"), c(29883888, 29899549, "PPP4R1L"), c(29902421, 29903778, "BCDIN3D"), c(29904615, 29906699, "COX14"), c(29906784, 29932421, "CERS5"), c(29934473, 29938894, "GPD1"), c(29940371, 29958942, "SMARCD1"), c(29967228, 30148945, "ASIC1"), c(30174097, 30177067, "CELA1"), c(30186073, 30188748, "CELA1"), c(30198248, 30210267, "GLS2"), c(30218632, 30222930, "BLOC1S1"), c(30223430, 30236778, "ITCH"), c(30239122, 30240596, "PXMP4"), c(30262928, 30268093, "TGM5"), c(30273165, 30278130, "TGM5"), c(30278372, 30302770, "EYA2"), c(30330324, 30353577, "RTFDC1"), c(30335291, 30343819, "GCNT7"), c(30357862, 30370582, "TFAP2C"), c(30529982, 30533197, "NAT8L"), c(30535792, 30539154, "NAT8L"), c(30540495, 30547148, "NAT8L"), c(30552104, 30560761, "NAT8L"), c(30668260, 30856329, "CDH22"), c(30954415, 30968875, "ZNFX1"), c(30969535, 30976403, "TBC1D20"), c(30976960, 30984300, "RBCK1"), c(30984928, 31032677, "NECAB3"), c(31039994, 31076618, "CBFA2T2"), c(31085947, 31117388, "SNTA1"), c(31129990, 31235221, "TOX2"), c(31251324, 31268116, "ST7L"), c(31272612, 31281722, "WNT2B"), c(31287438, 31302590, "CTTNBP2NL"), c(31331517, 31430416, "KCND3"), c(31439195, 31452853, "FAM212B"), c(31487428, 31509064, "RAP1A"), c(31510178, 31515467, "TIMM17A"), c(31516891, 31522221, "LMOD1"), c(31524221, 31534161, "SHISA4"), c(31536278, 31552145, "IPO9"), c(31554567, 31592401, "GNAS"), c(31598733, 31602506, "TUBB1"), c(31604014, 31609561, "SLMO2"), c(31612394, 31619444, "AURKA"), c(31622370, 31627134, "SLC32A1"), c(31654783, 31679849, "ARHGAP40"), c(31704958, 31796462, "PPP1R16B"), c(31806412, 31817286, "FAM210B"), c(31839663, 31854111, "NDUFA4L2"))
plot2hgncgenes_noclutter <- list( c(29018627, 29031263, "TSPAN2"), c(29032240, 29037994, "SLC25A22"), c(29042036, 29043571, "TSHB"), c(29045224, 29051904, "SLC5A8"), c(29053615, 29062173, "SYCP1"), c(29062794, 29071727, "SLC16A1"), c(29214203, 29328256, "FAM19A1"), c(29359912, 29375773, "PPM1J"), c(29379069, 29395463, "RHOC"), c(29395536, 29413061, "MOV10"), c(29414653, 29419791, "CAPZA3"), c(29420373, 29430461, "NCOA5"), c(29430975, 29641597, "SLC12A5"), c(29663579, 29672209, "NAT8L"), c(29718584, 29725342, "CYP24A1"), c(29728172, 29740385, "NFATC2"), c(29741561, 29744926, "MANBAL"), c(29745986, 29763944, "SRC"), c(29813777, 29830148, "APCDD1L"), c(29840983, 29865691, "VAPB"), c(29868540, 29883678, "RAB22A"), c(29883888, 29899549, "PPP4R1L"), c(29902421, 29903778, "BCDIN3D"), c(29904615, 29906699, "COX14"), c(29906784, 29932421, "CERS5"), c(29934473, 29938894, "GPD1"), c(29940371, 29958942, "SMARCD1"), c(29967228, 30148945, "ASIC1"), c(30174097, 30177067, "CELA1"), c(30186073, 30188748, "CELA1"), c(30198248, 30210267, "GLS2"), c(30218632, 30222930, "BLOC1S1"), c(30223430, 30236778, "ITCH"), c(30239122, 30240596, "PXMP4"), c(30262928, 30268093, "TGM5"), c(30278372, 30302770, "EYA2"), c(30330324, 30353577, "RTFDC1"), c(30335291, 30343819, "GCNT7"), c(30357862, 30370582, "TFAP2C"), c(30529982, 30533197, "NAT8L"), c(30668260, 30856329, "CDH22"), c(30954415, 30968875, "ZNFX1"), c(30969535, 30976403, "TBC1D20"), c(30976960, 30984300, "RBCK1"), c(30984928, 31032677, "NECAB3"), c(31039994, 31076618, "CBFA2T2"), c(31085947, 31117388, "SNTA1"), c(31129990, 31235221, "TOX2"), c(31251324, 31268116, "ST7L"), c(31272612, 31281722, "WNT2B"), c(31287438, 31302590, "CTTNBP2NL"), c(31331517, 31430416, "KCND3"), c(31487428, 31509064, "RAP1A"), c(31510178, 31515467, "TIMM17A"), c(31516891, 31522221, "LMOD1"), c(31524221, 31534161, "SHISA4"), c(31536278, 31552145, "IPO9"), c(31554567, 31592401, "GNAS"), c(31598733, 31602506, "TUBB1"), c(31604014, 31609561, "SLMO2"), c(31612394, 31619444, "AURKA"), c(31622370, 31627134, "SLC32A1"), c(31654783, 31679849, "ARHGAP40"), c(31704958, 31796462, "PPP1R16B"), c(31806412, 31817286, "FAM210B"), c(31839663, 31854111, "NDUFA4L2"))

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
plot3$COL <- colorRampPalette(c(my_colors_2[10], my_colors_2[11]))(10)[as.numeric(cut(plot3$FST,breaks = 10))]
plot3$alpha <- seq(29,99,7)[as.numeric(cut(plot3$FST,breaks = 10))]
plot3$COL <- paste0(plot3$COL, plot3$alpha)
plot3$COL[which(plot3$COL == "NANA")] <- "white"
plot3$COL[which(plot3$FST == 1)] <- "black"
# Genes for Plot 3
plot3genes <- list( c(33000284, 33004528, "DNASE1L1"), c(33020839, 33022458, "GZMA"), c(33051343, 33053133, "GZMA"), c(33060440, 33075984, "DNASE1L1"), c(33076920, 33086420, "SEC13"), c(33086940, 33104562, "CAND2"), c(33113841, 33116502, "ADORA3"), c(33135724, 33146117, "NCR2"), c(33170391, 33180872, "NCR2"), c(33200828, 33206154, "ADORA3"), c(33246299, 33316457, "ANGPT4"), c(33329161, 33373791, "STK4"), c(33381629, 33405900, "KCNS1"), c(33508179, 33536598, "MATN4"), c(33536326, 33560184, "RBPJL"), c(33916796, 33973294, "DAG1"), c(34369205, 34379167, "NICN1"), c(34642256, 34643870, "TRIM35"), c(34675168, 34836681, "GRM7"), c(34871935, 34873272, "TRIM39"), c(34893168, 35004251, "RNF123"), c(35005643, 35015323, "ALDH4A1"), c(35044233, 35048954, "PPP1R15B"), c(35051777, 35055989, "PIK3C2B"), c(35059764, 35089134, "KCNQ2"), c(35117209, 35127472, "DDX27"), c(35131780, 35166595, "PLTP"), c(35168038, 35177215, "KHK"), c(35211660, 35300827, "PCBP4"), c(35326453, 35331785, "RBM15B"), c(35331891, 35339781, "MANF"), c(35341377, 35346427, "C20orf111"), c(35385545, 35400256, "SDC4"), c(35428748, 35438611, "GDF5"), c(35507962, 35525146, "UQCC"), c(35549711, 35565594, "FAM83C"), c(35569932, 35591346, "EML6"), c(35620272, 35705767, "MMP24"), c(35726737, 35731720, "EIF6"), c(35740551, 35748347, "MAP1LC3A"), c(35750710, 35753783, "DYNLRB1"), c(35755299, 35772139, "SSUH2"), c(35773337, 35879838, "SEMA3E"), c(35880419, 35890888, "SMPD4"), c(35893896, 35914333, "VPRBP"), c(35920834, 36170178, "DOCK3"), c(36131652, 36139062, "BTN1A1"))
plot3hgncgenes_noclutter <- list( c(33000284, 33004528, "DNASE1L1"), c(33020839, 33022458, "GZMA"), c(33051343, 33053133, "GZMA"), c(33060440, 33075984, "DNASE1L1"), c(33076920, 33086420, "SEC13"), c(33086940, 33104562, "CAND2"), c(33113841, 33116502, "ADORA3"), c(33135724, 33146117, "NCR2"), c(33170391, 33180872, "NCR2"), c(33200828, 33206154, "ADORA3"), c(33246299, 33316457, "ANGPT4"), c(33329161, 33373791, "STK4"), c(33381629, 33405900, "KCNS1"), c(33508179, 33536598, "MATN4"), c(33536326, 33560184, "RBPJL"), c(33916796, 33973294, "DAG1"), c(34369205, 34379167, "NICN1"), c(34642256, 34643870, "TRIM35"), c(34675168, 34836681, "GRM7"), c(34871935, 34873272, "TRIM39"), c(34893168, 35004251, "RNF123"), c(35005643, 35015323, "ALDH4A1"), c(35044233, 35048954, "PPP1R15B"), c(35059764, 35089134, "KCNQ2"), c(35117209, 35127472, "DDX27"), c(35131780, 35166595, "PLTP"), c(35168038, 35177215, "KHK"), c(35211660, 35300827, "PCBP4"), c(35326453, 35331785, "RBM15B"), c(35331891, 35339781, "MANF"), c(35341377, 35346427, "C20orf111"), c(35385545, 35400256, "SDC4"), c(35428748, 35438611, "GDF5"), c(35507962, 35525146, "UQCC"), c(35549711, 35565594, "FAM83C"), c(35569932, 35591346, "EML6"), c(35620272, 35705767, "MMP24"), c(35726737, 35731720, "EIF6"), c(35740551, 35748347, "MAP1LC3A"), c(35750710, 35753783, "DYNLRB1"), c(35755299, 35772139, "SSUH2"), c(35773337, 35879838, "SEMA3E"), c(35880419, 35890888, "SMPD4"), c(35893896, 35914333, "VPRBP"), c(35920834, 36170178, "DOCK3"), c(36131652, 36139062, "BTN1A1"))


# Find important/high Zfst points
high <- subset(lg5, (lg5$Zfst > 5 & lg5$ID > plot1xstopID) | (lg5$Zfst > 5.9 & lg5$ID < plot1xstopID) )
high$Zfst = format(round(high$Zfst, 2), nsmall = 2)
high <- high[which(high$BIN_START != 35130001),]
lg5$COL[lg5$Zfst > 5 & lg5$ID > plot1xstopID ] <- 3
lg5$COL[lg5$Zfst > 5.9 & lg5$ID < plot1xstopID ] <- 3

# Actual Plotting ============================================================

# Plot LG5
png_name <- "x_top3LG5.png"
Cairo(png_name, type="png", , width = 15, height = 7, units = 'in', res=300)
par(fig=c(0,1,0,0.7))
plot(lg5$ID, lg5$Zfst, col=lg5$COL, xlab="Position on LG5 (Mbp)", ylab=expression("Z"[fst]), xaxs="i",xaxt="n", yaxt="n", ylim=c(0, lg5ystop), main="", pch=20, bty="n" )
# Label the important points
text(high[,1:2], labels=paste0(format(round(high$BIN_START/1000000, 2)), " Mbp"), pos=4, cex= 0.75)
#title(main = expression("LG5 UMD2a Bicuspid vs Tricuspid Z"[fst]))
axis(2, at=c(0,3,6,9), label=c(0,3,6,9))
axis(1, at=c(13888,14088, 14288, 14488, 14688,14888,15088, 15288, 15488, 15688, 15886,16086,16286,16486,16686,16886,17086,17286, 17501), label=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36.17))
text(lg5xstart+100, 10, pos=3, labels=c("10 kb bins"), cex=0.8)

# PLOT 1
par(fig=c(0,0.38,xplotsystart,1), new=TRUE)
plot(plot1$POS, plot1$FST, col = plot1$COL, xaxs="i", xlim=c(0,3000001), ylim=c(0,subplotstop), xlab="", ylab="Fst", xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
axis(2, at=c(0,0.5,1), label=c(0,"",1), cex =xplotcex)
axis(1, at=c(1000000, 2000000), label=c(1, 2), cex = xplotcex)
axis(1, at=c(0,3000001), label=c(0, 3), col.axis="black", cex = xplotcex)
axis(1, at=850000, label=".85", cex = xplotcex, padj = -1)
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
  if (as.numeric(genestart) - previousgenestop < 14000) {
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
axis(1, at=c(30000000,31000000), label=c(30,31), cex = xplotcex)
axis(1, at=c(29010050,31859939), label=c(29, 32), col.axis="black", cex = xplotcex)
axis(1, at=30250000, label="30.26", cex = xplotcex, padj = -1)
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
    
    if (genestart < 29110030) { text(as.numeric(genestart)+2000 + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(geneslabely + cummulativelabelbuffer), labels = name, srt=0, cex=geneslabelcex) }
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
plot(plot3$POS, plot3$FST, col = plot3$COL, xaxs="i", xlim=c(33000000, 36170000), ylim=c(0,subplotstop), xlab="", ylab="", xaxt = "n", yaxt = "n", bty="n", main="", cex = xplotcex)
axis(1, at=c(34000000,35000000), label=c(34,35), cex=xplotcex)
axis(1, at=c(33000000,36170000), label=c(33, 36.17), col.axis="black", cex = xplotcex)
axis(1, at=35100000, label=35.1, cex = xplotcex, padj = -1)
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
