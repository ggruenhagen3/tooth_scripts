# Species Tooth Shape
ts_df = data.frame( sample= c("403", "404", "493", "494", "495", "1615", "1619", "1621", "1818", "1860", "1863", "1912", "1983", "1984", "2162", "2241", "2277", "2298", "2302", "2319", "2320", "2332"),
                    label = c("403", "404", "493", "494", "495", "1615", "1619", "1621", "1818", "1860", "1863", "1912", "1983", "1984", "2162", "2241", "2277", "2298", "2302", "2319", "2320", "2332"),
                    round = c("1",   "1",   "2",   "1",   "2",   "2",    "2",    "1",    "2",    "1",    "1",    "2",    "1",    "1",    "2",    "2",    "1",    "2",    "2",    "2",    "2",    "1"),
                    ts=     c("tri", "tri", "tri", "tri", "tri", "bi",   "bi",    "bi",  "bi",   "bi",   "bi",   "bi",   "bi",   "bi",   "tri",  "tri",  "bi",   "tri",  "tri",  "tri",  "bi",   "tri"),
                    species=c("LF",  "LF",  "LF",  "LF",  "LF",  "MZ",   "MZ",   "MZ",   "MG",   "MP",   "MP",   "PC",   "PC",   "PC",   "PN",   "PN",   "ML",   "LF",   "LT",   "LF",   "ML",   "PN"))
palette(c("orange", "skyblue", "red", "blue", "green", "purple", "grey", "black", "pink"))


###################
# VCF2PHY - RAXML #
###################
library("ggplot2")
library("ape")
library("treeio")
library("ggtree")
library("Cairo")
library("dplyr")

# tree.nwk = read.newick("C:/Users/miles/Downloads/d_tooth/results/ts_ms/RAxML_bestTree.test_2")
tree.nwk = read.newick("C:/Users/miles/Downloads/d_tooth/results/ts_ms/RAxML_bestTree.vcf2phylip_08_31_20")
# ape_tree = read.tree("C:/Users/miles/Downloads/d_tooth/results/ts_ms/RAxML_bestTree.vcf2phylip_08_31_20")
tree2 <- full_join(tree.nwk, ts_df[which(ts_df$label %in% tree.nwk$tip.label),], by = "label")
print(ggtree(tree2, aes(col=species)) + geom_text(aes(x=branch, label=format(round(branch.length,2),nsmall=2)), vjust=-.8, color='gray48') + geom_tiplab() + geom_label(aes(x=branch, label=species), label.padding = unit(0.10, "lines")) + guides(color=FALSE))
# ggtree(tree.nwk) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

Cairo("C:/Users/miles/Downloads/d_tooth/results/ts_ms/vcf2phylip_08_31_20_phylo_species.png", type="png", , width = 14, height = 7, units = 'in', res=300)
print(ggtree(tree2, aes(col=species)) + geom_text(aes(x=branch, label=format(round(branch.length,2),nsmall=2)), vjust=-.8, color='gray48') + geom_tiplab() + geom_label(aes(x=branch, label=species), label.padding = unit(0.10, "lines")) + guides(color=FALSE))
dev.off()

Cairo("C:/Users/miles/Downloads/d_tooth/results/ts_ms/vcf2phylip_08_31_20_phylo.png", type="png", , width = 14, height = 7, units = 'in', res=300)
print(ggtree(tree2, aes(col=ts)) + geom_text(aes(x=branch, label=format(round(branch.length,2),nsmall=2)), vjust=-.8, color='gray48') + geom_tiplab() + geom_treescale() + geom_label(aes(x=branch, label=ts), label.padding = unit(0.10, "lines")) + guides(color=FALSE))
dev.off()

Cairo("C:/Users/miles/Downloads/d_tooth/results/ts_ms/vcf2phylip_08_31_20_phylo_test.png", type="png", , width = 14, height = 7, units = 'in', res=300)
print(ggtree(tree2, aes(col=round)) + geom_text(aes(x=branch, label=format(round(branch.length,2),nsmall=2)), vjust=-.8, color='gray48') + geom_tiplab() + geom_treescale() + geom_label(aes(x=branch, label=round), label.padding = unit(0.10, "lines")) + guides(color=FALSE))
dev.off()

#############
# SNPRelate #
#############
library("gdsfmt")
library("SNPRelate")
library("Cairo")
library("ggrepel")

#vcf to GDS
snpgdsVCF2GDS("C:/Users/miles/Downloads/filter_50.recode.vcf", "C:/Users/miles/Downloads/filter_50.recode.gds")
snpgdsSummary("C:/Users/miles/Downloads/filter_50.recode.gds")
genofile <- openfn.gds("C:/Users/miles/Downloads/filter_50.recode.gds")

#dendogram
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
dissMatrix  <-  snpgdsDiss(genofile , sample.id=NULL, snp.id=NULL, autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=10, verbose=TRUE)
snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)
ts_group = factor(ts_df$ts[match(snpHCluster$sample.id, ts_df$sample)])
species_group = factor(ts_df$species[match(snpHCluster$sample.id, ts_df$sample)])
cutTree <- snpgdsCutTree(snpHCluster, samp.group=species_group, z.threshold=15, outlier.n=5, n.perm = 5000,col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, verbose=TRUE)

Cairo("C:/Users/miles/Downloads/d_tooth/results/ts_ms/no_404_405_2320.png", type="png", , width = 14, height = 7, units = 'in', res=300)
# par(fig=c(0,0.7,0,1))
plot(cutTree$dendrogram, ylim=c(0,1.3), lwd=2, ps=6, cex=6, cex.lab=6, main="Tooth Samples")
# par(fig=c(0.7,1,0,1))
# legend("topleft", legend=levels(cutTree$samp.group), col=1:nlevels(cutTree$samp.group), pch=19, bty="n")
legend(x=0,y=1.3, legend=levels(cutTree$samp.group), col=1:nlevels(cutTree$samp.group), pch=19, bty="n", horiz=T, x.intersp = 0.4, text.width = 2.1)
dev.off()

#pca
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.gdsn(index.gdsn(genofile, "sample.id"))
pca <- snpgdsPCA(genofile, autosome.only = FALSE)
tab <- data.frame(sample.id = pca$sample.id,pop = factor(pop_code)[match(pca$sample.id, sample.id)],EV1 = pca$eigenvect[,1],EV2 = pca$eigenvect[,2],stringsAsFactors = FALSE)
par(fig=c(0,0.7,0,1))
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),xlab="eigenvector 2", ylab="eigenvector 1")
par(fig=c(0.7,1,0,1))
legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

tab$species = ts_df$species[match(tab$sample.id, ts_df$sample)]
tab$ts = ts_df$ts[match(tab$sample.id, ts_df$sample)]
ggplot(tab, aes(EV1, EV2, col=species)) + geom_point() + geom_text_repel(aes(label=sample.id),hjust=0, vjust=0) + scale_color_brewer(palette="Dark2")
