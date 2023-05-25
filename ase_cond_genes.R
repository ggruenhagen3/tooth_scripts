rna_path <- "C:/Users/miles/Downloads/brain/"
data <- read.csv(paste(rna_path, "/data/pnas.1810140115.sd02.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

gene_counts <- table(data[,1])
three <- names(gene_counts[which(gene_counts == 3)])
three <- data[which(data[,1] %in% three),]
baseline <- data.frame(genes <- c())
same_sign <- c()
for (i in 1:nrow(three)) {
  gene <- three[i,1]
  this_rows <- three[which(three[,1] == gene),]
  if (!gene %in% same_sign) {
    if (this_rows$rep_1_ase_ratio[1] < 0 && this_rows$rep_2_ase_ratio[1] < 0 &&
        this_rows$rep_1_ase_ratio[2] < 0 && this_rows$rep_2_ase_ratio[2] < 0 &&
        this_rows$rep_1_ase_ratio[3] < 0 && this_rows$rep_2_ase_ratio[3] < 0) {
      baseline <- rbind(baseline, this_rows)
      same_sign <- c(same_sign, gene)
    }
    if (this_rows$rep_1_ase_ratio[1] > 0 && this_rows$rep_2_ase_ratio[1] > 0 &&
        this_rows$rep_1_ase_ratio[2] > 0 && this_rows$rep_2_ase_ratio[2] > 0 &&
        this_rows$rep_1_ase_ratio[3] > 0 && this_rows$rep_2_ase_ratio[3] > 0) {
      baseline <- rbind(baseline, this_rows)
      same_sign <- c(same_sign, gene)
    }
  }
}

build_dep <- data.frame()
build <- data[which(data$condition == "building"),]
dig <- data[which(data$condition == "digging"),]
iso <- data[which(data$condition == "isolated"),]
for (i in 1:nrow(build)) {
  gene <- build[i,1]
  pos <- build$rep_1_ase_ratio > 0 && build$rep_2_ase_ratio > 0
  if (gene %in% iso[,1]) {
    if (pos && iso$rep_1_ase_ratio < 0 && iso$rep_2_ase_ratio < 0) {
      build_dep <- rbind(build_dep, build[i,])
    }
    if (! pos && iso$rep_1_ase_ratio > 0 && iso$rep_2_ase_ratio > 0) {
      build_dep <- rbind(build_dep, build[i,])
    }
  } else {
    build_dep <- rbind(build_dep, build[i,])
  }
}
bhve_dep <- build_dep
for (i in 1:nrow(dig)) {
  gene <- dig[i,1]
  pos <- dig$rep_1_ase_ratio > 0 && dig$rep_2_ase_ratio > 0
  if (gene %in% iso[,1]) {
    if (pos && iso$rep_1_ase_ratio < 0 && iso$rep_2_ase_ratio < 0) {
      bhve_dep <- rbind(bhve_dep, dig[i,])
    }
    if (! pos && iso$rep_1_ase_ratio > 0 && iso$rep_2_ase_ratio > 0) {
      bhve_dep <- rbind(bhve_dep, dig[i,])
    }
  } else {
    bhve_dep <- rbind(bhve_dep, dig[i,])
  }
}

df <- data.frame(data <- unique(baseline[,1]), bio <- rep("BASELINE", length(unique(baseline[,1]))))
write.table(df, file = "C:/Users/miles/Downloads/brain/data/markers/baseline.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
df <- data.frame(data <- unique(bhve_dep[,1]), bio <- rep("BHVE_DEP", length(unique(bhve_dep[,1]))))
write.table(df, file = "C:/Users/miles/Downloads/brain/data/markers/bhve_dep.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
df <- data.frame(data <- unique(build_dep[,1]), bio <- rep("BUILD_DEP", length(unique(build_dep[,1]))))
write.table(df, file = "C:/Users/miles/Downloads/brain/data/markers/build_dep.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
