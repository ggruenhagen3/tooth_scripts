# Combine Epi and Mes
library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("monocle3")
library("EnhancedVolcano")

mes <- readRDS("C:/Users/miles/Downloads/rna/data/combined.rds")
epi <- readRDS("C:/Users/miles/Downloads/d_tooth/data/epi_full.rds")

mes$project <- "MES"
epi$project <- "EPI"

neo <- merge(mes, epi, merge.data = TRUE)
neo <- FindVariableFeatures(object = neo, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
neo <- ScaleData(object = neo, vars.to.regress = NULL)
neo <- RunPCA(neo, npcs = 30, verbose = FALSE) 
# All PC's significant
neo <- RunUMAP(neo, reduction = "pca", dims = 1:30)
neo <- FindNeighbors(neo, reduction = "umap", dims = 1:2)
neo <- FindClusters(neo, resolution = 0.25)
DimPlot(neo, reduction = "umap", label = TRUE)
DimPlot(neo, reduction = "umap", split.by = "project", label = TRUE)

num_epi_clusters <- as.numeric(tail(levels(epi@meta.data$seurat_clusters), n=1))
num_mes_clusters <- as.numeric(tail(levels(mes@meta.data$seurat_clusters), n=1))
for (epi_clust in 0:num_epi_clusters) {
  neo <- SetIdent(neo, cells=WhichCells(epi, idents = epi_clust), value=paste("epi", epi_clust, sep="_"))
}
for (mes_clust in 0:num_mes_clusters) {
  neo <- SetIdent(neo, cells=WhichCells(mes, idents = mes_clust), value=paste("mes", mes_clust, sep="_"))
}
neo$orig.cluster <- neo@active.ident
# saveRDS(neo, "C:/Users/miles/Downloads/d_tooth/data/epi_mes.rds")
# neo <- readRDS("C:/Users/miles/Downloads/d_tooth/data/epi_mes.rds")

# Cluster our tooth data with neo
tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
tj <- NormalizeData(tj, normalization.method = "LogNormalize", scale.factor = 10000)
# VlnPlot(neo_tj, features = "Wnt10a")
neo_backup <- neo
tj <- mzebraToMouseInPlace2(tj)

# neo <- mzebraToMouseInPlace(neo)
# new_tj_counts <- as.matrix(tj@assays$RNA@counts[rownames(neo)[which(rownames(neo) %in% rownames(tj))],])
# new_tj_data   <- tj@assays$RNA@data[rownames(neo)[which(rownames(neo) %in% rownames(tj))],]
# tj <- CreateSeuratObject(counts = new_tj_counts, project = "TJ")
# tj <- SetAssayData(object = tj, slot = 'data', new.data = new_tj_data)
# tj$seurat_clusters <- tj_backup$seurat_clusters

new_neo_counts <- as.matrix(neo@assays$RNA@counts[rownames(neo)[which(rownames(neo) %in% rownames(tj))],])
new_neo_data   <- neo@assays$RNA@data[rownames(neo)[which(rownames(neo) %in% rownames(tj))],]
neo <- CreateSeuratObject(counts = new_neo_counts, project = "NEO")
neo <- SetAssayData(object = neo, slot = 'data', new.data = new_neo_data)
neo$seurat_clusters <- neo_backup$seurat_clusters

tj$org <- "MZEBRA"
neo$org <- "MOUSE"
neo_tj <- merge(neo, tj, merge.data = TRUE)
neo_tj <- FindVariableFeatures(object = neo_tj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
neo_tj <- ScaleData(object = neo_tj, vars.to.regress = NULL)
neo_tj <- RunPCA(neo_tj, npcs = 30, verbose = FALSE) 
# All PC's significant
neo_tj <- RunUMAP(neo_tj, reduction = "pca", dims = 1:30)
neo_tj <- FindNeighbors(neo_tj, reduction = "umap", dims = 1:2)
neo_tj <- FindClusters(neo_tj, resolution = 0.25)
DimPlot(neo_tj, reduction = "umap", label = TRUE)
DimPlot(neo_tj, reduction = "umap", split.by = "org", label = TRUE)
# saveRDS(neo_tj, "C:/Users/miles/Downloads/d_tooth/data/neo_tj.rds")
# neo_tj <- readRDS("C:/Users/miles/Downloads/d_tooth/data/neo_tj.rds")

# DEG Across Org
Idents(neo_tj) <- neo_tj$org
mzebra.markers <- FindMarkers(neo_tj, ident.1 = "MZEBRA", ident.2 = "MOUSE", only.pos = TRUE)
mzebra.markers <- mzebra.markers[which(mzebra.markers$p_val_adj < .05),]
mzebra.markers$gene <- rownames(mzebra.markers)
mzebra.markers  <- mzebra.markers[,c(ncol(mzebra.markers), 1:(ncol(mzebra.markers)-1))]
write.table(mzebra.markers, file = paste(rna_path, "/results/mzebra_vs_mouse_mzebra_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

FeaturePlot(neo_tj, features = "ptchd4", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)

num_epi_clusters <- as.numeric(tail(levels(epi@meta.data$seurat_clusters), n=1))
num_mes_clusters <- as.numeric(tail(levels(mes@meta.data$seurat_clusters), n=1))
num_tj_clusters <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
Idents(tj) <- tj$seurat_clusters
for (epi_clust in 0:num_epi_clusters) {
  neo_tj <- SetIdent(neo_tj, cells=WhichCells(epi, idents = epi_clust), value=paste("epi", epi_clust, sep="_"))
}
for (mes_clust in 0:num_mes_clusters) {
  neo_tj <- SetIdent(neo_tj, cells=WhichCells(mes, idents = mes_clust), value=paste("mes", mes_clust, sep="_"))
}
for (tj_clust in 0:num_tj_clusters) {
  neo_tj <- SetIdent(neo_tj, cells=WhichCells(tj, idents = tj_clust), value=paste("tj", tj_clust, sep="_"))
}
neo_tj$orig.cluster <- neo_tj@active.ident
DimPlot(neo_tj, reduction = "umap", split.by = "org", label = TRUE)
######################
## Helper Functions ##
######################
library(biomaRt)
library(dplyr)
library(stringr)
# library(rentrez)

geneCap <- function(gene, gene_names) {
  # Gene the gene name in the right format
  gene_lower <- tolower(gene)
  gene_upper <- toupper(gene)
  gene_title <- str_to_title(gene)
  error <- FALSE
  if (gene_lower %in% gene_names) {
    gene <- gene_lower
  } else if (gene_upper %in% gene_names) {
    gene <- gene_upper
  } else if (gene_title %in% gene_names) {
    gene <- gene_title
  } else {
    error <- TRUE
  }
  
  return(c(gene, error))
}

validGenes <- function(genes, gene_names) {
  valid_genes <- c()
  for (gene in genes) {
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    if (! error) {
      valid_genes <- c(valid_genes, gene)
    }
  } # end gene for
  valid_genes <- unique(valid_genes)
  return(valid_genes)
} # end validGenes function


convertToHgncObj <- function(obj, organism) {
  # Convert a Seurat object to have all HGNC gene names
  print("Converting Genes Names...")
  genes <- rownames(obj)
  if (organism == "mouse") {
    dataset_name <- "mmusculus_gene_ensembl"
  } else if (organism == "mzebra") {
    dataset_name <- "mzebra_gene_ensembl"
  } else {
    stop("Organism not recognized. Options are mouse or mzebra.")
  }
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  org = useMart(biomart="ensembl", dataset=dataset_name)
  
  # DF to convert from org to HGNC
  all_hgnc <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = org, attributesL = c("external_gene_name"), martL = human, uniqueRows=T)
  
  # Initialize New Matricies
  print("Creating New Matrices...")
  new_counts_matrix <- as.matrix(obj@assays$RNA@counts)
  new_data_matrix   <- as.matrix(obj@assays$RNA@data)
  
  not_i <- 0
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  for (gene in genes) {
    if (gene %in% all_hgnc[,1]) {
      hgnc_gene <- all_hgnc[which(all_hgnc[,1] == gene),2]
      if (length(hgnc_gene) > 1) {
        multiple_hgnc <- multiple_hgnc + 1
        upper_hgnc_gene <- hgnc_gene[which(startsWith(tolower(hgnc_gene), gene))]
        if (length(upper_hgnc_gene) == 1) {
          hgnc_gene <- upper_hgnc_gene
        } else {
          bad_multiple_hgnc <- bad_multiple_hgnc + 1
          hgnc_gene <- hgnc_gene[1]
        } # end bad multiple
      } # end multiple
      
      rownames(new_counts_matrix)[which(rownames(new_counts_matrix) == gene)] <- hgnc_gene
      rownames(new_data_matrix)[which(rownames(new_data_matrix) == gene)] <- hgnc_gene
    } else {
      not_i <- not_i + 1
    } # end gene not an hgnc gene
  } # end gene for
  print(paste("Number of Genes not converted to HGNC:", not_i))
  print(paste("Number of Genes with multple HGNC:", multiple_hgnc))
  print(paste("Number of Genes with multple HGNC and non-ideal circumstance:", bad_multiple_hgnc))
  
  # Merge the duplicated rows
  print("Removing duplicated HGNC rows...")
  dup_genes <- rownames(new_counts_matrix)[which(duplicated(rownames(new_counts_matrix)))]
  dup_ind <- c()
  for (gene in dup_genes) {
    ind <- which(rownames(new_counts_matrix) == gene)
    
    new_counts_row <- new_counts_matrix[ind[1],]
    new_data_row   <- new_data_matrix[ind[1],]
    for (i in 2:length(ind)) {
      new_counts_row <- new_counts_row + new_counts_matrix[ind[i],]
      new_data_row   <- new_data_row   + new_data_matrix[ind[i],]
    }
    new_counts_matrix[ind[1],] <- new_counts_row
    new_data_matrix[ind[1],]   <- new_data_row

    # Delete the duplicated rows
    dup_ind <- c(dup_ind, ind[2:length(ind)])
  }
  # Delete all the duplicated rows at once
  new_counts_matrix <- new_counts_matrix[-dup_ind,]
  new_data_matrix   <- new_data_matrix[-dup_ind,]

  # Remove mzebra rows
  print("Removing old org rows...")
  ind <- which(rownames(new_counts_matrix) %in% all_hgnc[,2])
  new_counts_matrix <- new_counts_matrix[ind,]
  new_data_matrix   <- new_data_matrix[ind,]
  
  print("Creating New Seurat Object...")
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  obj_2$seurat_clusters <- obj$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj@meta.data)) {
    obj_2@meta.data[col] <- obj@meta.data[col]
  }
  
  return(obj_2)
}

convertToHgnc <- function(genes) {
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  zfin_genes    <- getLDS(attributes = c("zfin_id_symbol"), filters = "zfin_id_symbol", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  hgnc_genes    <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  all_hgnc <- rbind(ensembl_genes, setNames(zfin_genes, names(ensembl_genes)), setNames(hgnc_genes, names(ensembl_genes)))
  all_hgnc = all_hgnc[!duplicated(all_hgnc[,1]),]
  colnames(all_hgnc) <- c("mzebra", "hgnc")
  
  # all_hgnc <- unique(c(ensembl_genes[,2], zfin_genes[,2], hgnc_genes[,2]))
  return(all_hgnc)
}

hgncMzebra <- function(genes, gene_names) {
  pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  valid_genes <- validGenes(genes, gene_names)
  all_hgnc <- convertToHgnc(gene_names)
  ind <- match(genes,pat$V2)
  ind <- ind[! is.na(ind)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  
  df1 <- setNames(as.data.frame(found_names_hgnc), c("hgnc"))
  found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")
  # found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")$mzebra
  # ind_found_hgnc <- match(found_names_hgnc, all_hgnc$hgnc)
  # ind_found_hgnc <- ind_found_hgnc[! is.na(ind_found_hgnc)]
  # found_names_hgnc <- as.vector(all_hgnc[ind_found_hgnc,1])
  
  pseudo_hgnc <- toupper(genes)
  df1 <- setNames(as.data.frame(pseudo_hgnc), c("hgnc"))
  found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")
  # found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")$mzebra
  # ind_hgnc <- match(pseudo_hgnc, all_hgnc$hgnc)
  # ind_hgnc <- ind_hgnc[! is.na(ind_hgnc)]
  # found_mzebra <- as.vector(all_hgnc[ind_hgnc,1])
  
  found_mzebra <- found_mzebra[,2:1]
  found_names_hgnc <- found_names_hgnc[,2:1]
  good_df <- rbind(all_hgnc, setNames(found_names, names(all_hgnc)), setNames(found_mzebra, names(all_hgnc)), setNames(found_names_hgnc, names(all_hgnc)))
  # good <- c(valid_genes, found_names, found_names_hgnc, found_mzebra)
  good <- good[which(good != "")]
  good <- validGenes(good, gene_names)
  good <- unique(good)
  good <- sort(good)
  return(good)
}


# Helper Functions
mzebraToMouseInPlace <- function(obj) {
  genes <- rownames(obj@assays$RNA@counts)
  mouse  = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  all_mzebra <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = mouse, attributesL = c("external_gene_name"), martL = mzebra, uniqueRows=T)
  
  all_mzebra <- all_mzebra[!duplicated(all_mzebra[,2]),]
  all_mzebra <- all_mzebra[!duplicated(all_mzebra[,1]),]
  colnames(all_mzebra) <- c("mouse", "mzebra")
  
  new_counts_matrix <- as.matrix(obj@assays$RNA@counts[which(rownames(obj@assays$RNA@counts) %in% all_mzebra[,1]),])
  # new_counts_matrix <- new_counts_matrix[which(rownames(new_counts_matrix) %in% all_mzebra[,1]),]
  new_data_matrix   <- as.matrix(obj@assays$RNA@data[which(rownames(obj@assays$RNA@counts) %in% all_mzebra[,1]),])
  # new_data_matrix   <- new_data_matrix[which(rownames(new_data_matrix) %in% all_mzebra[,1]),]
  genes <- rownames(new_data_matrix)
  for (i in 1:nrow(all_mzebra)) {
    mouse_gene <- all_mzebra[i,1]
    mzebra_gene  <- all_mzebra[i,2]
    # print(mzebra_gene)
    # print(mouse_gene)
    # print(ind)
    ind <- which(genes == mouse_gene)
    rownames(new_counts_matrix)[ind] <- mzebra_gene
    rownames(new_data_matrix)[ind]   <- mzebra_gene
  }
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  obj_2$seurat_clusters <- obj$seurat_clusters
  
  return(obj_2)
}

mzebraToMouseInPlace2 <- function(obj) {
  genes <- rownames(obj@assays$RNA@counts)
  mz_to_mo_file <- read.table("C:/Users/miles/Downloads/all_research/mzebra_mouse.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
  colnames(mz_to_mo_file) <- c("Gene_stable_ID", "Gene_name", "Mouse_gene_stable_ID", "Mouse_gene_name")
  
  new_counts_matrix <- as.matrix(obj@assays$RNA@counts)
  new_data_matrix   <- as.matrix(obj@assays$RNA@data)
  mz_to_mo_df <- data.frame()
  i <- 0
  not_i <- 0
  for (gene in genes) {
    mo_rows <- mz_to_mo_file[which(mz_to_mo_file$Gene_stable_ID == gene | mz_to_mo_file$Gene_name == gene),]

    if (nrow(mo_rows) > 1) {

      mo_gene <- mo_rows$Mouse_gene_name[which(startsWith(tolower(mo_rows$Mouse_gene_name), tolower(substr(gene,1,2))))]
      if (! identical(mo_gene, character(0))) {
        mo_gene <- mo_gene[1]
      } else {
        mo_gene <- mo_rows[1,4]
      }
    } else if (nrow(mo_rows) == 1) {
      mo_gene <- mo_rows[1,4]
    } else {
      # cat(paste(gene, "not found", "\n"))
      not_i <- not_i + 1
      mo_gene <- ""
    }
    
    if( (! identical(mo_gene, character(0))) && mo_gene != "" ) {
      i <- i + 1
      mz_to_mo_df <- rbind(mz_to_mo_df, t(c(gene, mo_gene)))
      rownames(new_counts_matrix)[which(rownames(new_counts_matrix) == gene)] <- mo_gene
      rownames(new_data_matrix)[which(rownames(new_data_matrix) == gene)] <- mo_gene
    }
    
  }
  
  # Merge the duplicated rows
  n_occur <- data.frame(table(rownames(new_counts_matrix)))
  dup_genes <- n_occur$Var1[which(n_occur$Freq > 1)]
  for (gene in dup_genes) { 
    ind <- which(rownames(new_counts_matrix) == gene)
    # this_rows <- new_counts_matrix[ind,]
    
    new_row_counts <- new_counts_matrix[ind[1],]
    new_row_data   <- new_data_matrix[ind[1],]
    for (i in 2:length(ind)) {
      new_row_counts <- new_row_counts + new_counts_matrix[ind[i],]
      new_row_data   <- new_row_data   + new_data_matrix[ind[i],]
    }
    
    # Delete the duplicated rows
    new_counts_matrix <- new_counts_matrix[-ind[2:length(ind)],]
    new_data_matrix   <- new_data_matrix[-ind[2:length(ind)],]
  }
  
  # Remove mzebra rows
  ind <- which(rownames(new_counts_matrix) %in% mz_to_mo_df$V2)
  new_counts_matrix <- new_counts_matrix[ind,]
  new_data_matrix   <- new_data_matrix[ind,]
  
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  obj_2$seurat_clusters <- obj$seurat_clusters
  
  return(obj_2)
}