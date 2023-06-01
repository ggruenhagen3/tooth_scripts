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

# Load the unbinned data
var_fst_data <- read.table(path_to_var, sep="\t", header=TRUE)
# Rename the columns
names(var_fst_data) <- c("LG", "POS", "FST")
# Convert the negative fst scores to 0
var_fst_data[, c(3)][var_fst_data[, c(3)] < 0] <- 0

output_folder = "C:/Users/miles/Downloads/d_tooth/results/ts_ms/talha/"
# gene_info = read.table("C:/Users/miles/Downloads/all_research/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)

buffer = 100e3
ish_genes_df = data.frame()
ish_genes = c("bmp2", "fst", "shh", "fgf3", "LOC101468418", "LOC101468567", "notch1", "LOC101467699", "hes1", "bhlhe40")
ish_genes_hgnc = read.csv("C:/Users/miles/Downloads/d_tooth/ts_ms_markers_to_check.txt", header = F, stringsAsFactors = F)[,1]
gene_info = read.table("C:/Users/miles/Downloads/all_research/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
ish_genes = gene_info$mzebra[match(ish_genes_hgnc, gene_info$human)]
ish_genes = ish_genes[which( ! is.na(ish_genes) )]
for (gene in ish_genes) {
  print(gene)
  gene_start = as.numeric(pat$V4[which(pat$V2 == gene)])
  gene_stop = as.numeric(pat$V5[which(pat$V2 == gene)])
  plot_start = gene_start - buffer
  plot_stop = gene_stop + buffer
  if (plot_start <= 0) { plot_start = 0 }
  # if (plot_start <= 0) { plot_start = 0 }
  lg = pat$V3[which(pat$V2 == gene)]
  lg_nc = as.numeric(substr(lg, 3, nchar(lg))) + 36779.1
  lg_nc = paste0("NC_0", as.character(lg_nc))
  if (lg == "LG22") { lg_nc = "NC_036800.1" }
  if (lg == "LG23") { lg_nc = "NC_036801.1" }
  if (substr(lg, 0,2) == "NW") { lg_nc = lg }
  
  data = var_fst_data[which( var_fst_data$LG == lg_nc & var_fst_data$POS < plot_stop & var_fst_data$POS > plot_start ),]
  data$COL = "gray50"
  data$COL[which(data$POS < gene_stop & data$POS > gene_start)] = "green4"
  data = data[order(data$COL),]
  
  data_b = fst_data[which( fst_data$CHROM == lg_nc & fst_data$BIN_START < plot_stop & fst_data$BIN_END > plot_start ),]
  data_b$COL = "gray50"
  data_b$COL[which(data_b$BIN_START < gene_stop & data_b$BIN_END > gene_start)] = "green4"
  # data_b = data_b[order(data_b$COL),]
  
  png_name = paste0(output_folder, gene, ".png")
  Cairo(png_name, type="png", width = 7, height = 4, res = 300, units = "in")
  plot(data$POS, data$FST, col = data$COL, ylim=c(0,1), xaxs="i", xlab = paste0("Position on ", lg), ylab = "FST", main = gene, pch = 20)
  dev.off()
  
  png_name = paste0(output_folder, gene, "_b.png")
  Cairo(png_name, type="png", width = 7, height = 4, res = 300, units = "in")
  plot(data_b$ID, data_b$Zfst, col = data_b$COL, xaxs="i", xlab = paste0("Position on ", lg), ylab = "Zfst", main = gene, pch = 16, xaxt = "n")
  text(data_b$ID[1] + 2, 1, pos=3, labels=c("10 kb bins"), cex=0.8)
  dev.off()
}

# High FST genes on LG5 and LG15
high_fst_genes = read.csv("C:/Users/miles/Downloads/d_tooth/results/ts_ms/var_fix_sig_bin_25kb_genes.txt", header = F, stringsAsFactors = F)[,1]
high_fst_genes_df = data.frame(genes = high_fst_genes, lg = pat$V3[match(high_fst_genes, pat$V2)], human = gene_info$human[match(high_fst_genes, gene_info$mzebra)])
high_fst_genes_5_15 = high_fst_genes_df[which(high_fst_genes_df$lg == "LG5" | high_fst_genes_df$lg == "LG15"),]
high_fst_genes_5_15 = high_fst_genes_5_15[order(high_fst_genes_5_15$lg, decreasing = T),]
write.csv(high_fst_genes_5_15, "C:/Users/miles/Downloads/d_tooth/results/ts_ms/high_fst_genes_lg5_lg15.csv", quote = F)

high_fst_genes_5_15 = read.csv("C:/Users/miles/Downloads/d_tooth/results/ts_ms/high_fst_genes_lg5_lg15.csv")
high_fst_genes_5 = high_fst_genes_5_15[which(high_fst_genes_5_15$lg == "LG5"),]
high_fst_genes_5$start = pat$V4[match(high_fst_genes_5$genes, pat$V2)]
high_fst_genes_5$stop = pat$V5[match(high_fst_genes_5$genes, pat$V2)]
bhlhe40_row = pat[which(pat$V2 == "bhlhe40"),]
high_fst_genes_5$dist_bhlhe40_start = as.numeric(as.vector(high_fst_genes_5$start)) - as.numeric(as.vector(bhlhe40_row$V4))
high_fst_genes_5$abs_dist_bhlhe40_start = as.numeric(as.vector(high_fst_genes_5$start)) - as.numeric(as.vector(bhlhe40_row$V4))
write.csv(high_fst_genes_5_15, "C:/Users/miles/Downloads/d_tooth/results/ts_ms/high_fst_genes_lg5_dist_from_bhlhe40.csv", quote = F)

ish_genes_df = gene_info[match(ish_genes[which(ish_genes %in% high_fst_genes)], gene_info$mzebra),]
ish_genes_df$LG = pat$V3[match(ish_genes_df$mzebra, pat$V2)]
write.csv(ish_genes_df, "C:/Users/miles/Downloads/d_tooth/results/ts_ms/talha_markers_in_high_fst_genes.csv")

# Plot some high FST genes that are close to bhlhe40
library("scales")
plot_genes = c("id1", "sema3f", "LOC101484715", "sox13", "sema3f", "bhlhe40")
col_pal = hue_pal()(length(plot_genes))
names(col_pal) = plot_genes
data = fst_data[which( fst_data$CHROM == "NC_036784.1" ),]
data$COL = "gray50"
for (i in 1:length(plot_genes)) {
  gene = plot_genes[i]
  gene_start = as.numeric(pat$V4[which(pat$V2 == gene)])
  gene_stop = as.numeric(pat$V5[which(pat$V2 == gene)])
  data$COL[which(data$BIN_START < gene_stop & data$BIN_END > gene_start)] = col_pal[i]
}

data2 = data[which(data$COL != "gray50"),]
data = data[which(data$COL == "gray50"),]
png_name = paste0("C:/Users/miles/Downloads/d_tooth/results/ts_ms/lg5_high_fst_genes_near_bhlhe40.png")
Cairo(png_name, type="png", width = 15, height = 5, res = 300, units = "in")
plot(data$ID, data$Zfst, col = data$COL, xaxs="i", xlab = paste0("Position on LG5"), ylab = "Zfst", main = "High FST Genes Near Bhlhe40", pch = 16, xaxt = "n")
label_pos = data$ID[which(data$BIN_END %% 1000000 == 0)]
axis(1, at=label_pos, label=1:length(label_pos))
points(data2$ID, data2$Zfst, col = data2$COL, pch = 20)
text(data$ID[1] + 100, 4, pos=3, labels=c("10 kb bins"), cex=0.8)
for (i in 1:length(col_pal)) {
  text(data2$ID[which(data2$COL == col_pal[i])], 3 + i*0.25, pos = 3, labels = c(names(col_pal[i])), cex = 0.8, col = col_pal[i])
}
dev.off()
