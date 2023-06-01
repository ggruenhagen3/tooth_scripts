#*******************************************************************************
# Load Libraries ===============================================================
#*******************************************************************************
message('Loading Libraries')
wdstr = substr(getwd(), 1, 12)
switch(wdstr,
       "C:/Users/mil" = { main_path = "C:/Users/miles/Downloads/";        },
       "/home/george" = { main_path = "~/research/"                       },
       "/storage/scr" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/hom" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/cod" = { main_path = "/storage/scratch1/6/ggruenhagen3/" })
switch(wdstr,
       "C:/Users/mil" = { tooth_path = "d_tooth/" },
       "/home/george" = { tooth_path = "tooth/"   },
       "/storage/scr" = { tooth_path = "d_tooth/" },
       "/storage/hom" = { tooth_path = "d_tooth/" },
       "/storage/cod" = { tooth_path = "d_tooth/" })
switch(wdstr,
       "C:/Users/mil" = { gene_info_path = "all_research/" },
       "/home/george" = { gene_info_path = "all_research/" },
       "/storage/scr" = { gene_info_path = "m_zebra_ref/"  },
       "/storage/hom" = { gene_info_path = "m_zebra_ref/"  },
       "/storage/cod" = { gene_info_path = "m_zebra_ref/"  })
full_path = paste0(main_path, tooth_path)
gene_info_path = paste0(main_path, gene_info_path)
brain_dir = paste0(main_path, "brain/")
data_dir  = paste0(full_path, "data/")
out_dir   = paste0(full_path, "results/")
source(paste0(brain_dir, "/brain_scripts/all_f.R"))
source(paste0(full_path, "/tooth_scripts/plk_glmmseq.R"))
setwd(out_dir)

library(glmmSeq)
library(Seurat)
library(lme4)
library(lmerTest)
library(parallel)
library(edgeR)
library(DESeq2)
library(glmGamPoi)
library(scran)
library(BiocParallel)

#*******************************************************************************
# Load Objects =================================================================
#*******************************************************************************
message('Loading Objects')
gene_info = read.csv(paste0(gene_info_path, "gene_info_3.csv"))
plk = readRDS(paste0(data_dir, "plkall_053023.rds"))
plk_subject = readRDS(paste0(data_dir, "plkall_subject_053023.rds"))

# Main =====
message('Starting glmmseq analysis')
obj = subset(plk_subject, cells = colnames(plk_subject)[which(plk_subject$exp == "plk3")])
obj$pair = obj$subject
n_pairs = length(unique(obj$pair))
for (this_clust in sort(unique(obj$seurat_clusters))) {
  this_cells = colnames(obj)[which(obj$seurat_clusters == this_clust)]
  if (length(unique(obj$pair[this_cells])) < n_pairs) {
    message(paste0("Not all pairs present in cluster ", this_clust))
  } else {
    message(paste0("Performing glmmSeq on cluster ", this_clust))
    res = fastGlmm(obj, this_cells, num_cores = 20, out_path = paste0("~/scratch/d_tooth/results/plk_glmmseq_plk3_clusters50/cluster_", this_clust, ".csv")) 
  }
}