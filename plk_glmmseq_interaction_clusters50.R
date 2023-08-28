#*******************************************************************************
# Load Libraries ===============================================================
#*******************************************************************************
.libPaths(c("/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/conda_envs/glmmseq/lib/R/library", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
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
library(SeuratObject)
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

#*******************************************************************************
# Load Objects =================================================================
#*******************************************************************************
obj = plk_subject
obj$pair = obj$subject
n_pairs = length(unique(obj$pair))
big_res = data.frame()
glmm_out_dir = "~/scratch/d_tooth/results/plk_glmmseq_interaction_clusters50/"

for (this_clust in sort(unique(obj$seurat_clusters))) {
  this_cells = colnames(obj)[which(obj$seurat_clusters == this_clust)]
  pair_count = table(obj$pair[this_cells])
  if (length(pair_count) < n_pairs || any(pair_count < 3)) {
    message(paste0("There are not 3 cells present from all pairs in cluster ", this_clust))
  } else {
    message(paste0("Performing glmmSeq on cluster ", this_clust))
    my_formula = ~ cond*exp + (1|subject)
    res = fastGlmm(obj, this_cells, my_formula = my_formula, my_formula2 = ~ subject, calc_pct = F, return_means=T, do_contrasts=T, num_cores = 24, out_path = paste0(glmm_out_dir, "cluster_", this_clust, ".csv")) 
    if (!is.null(res)) { 
      res$gene = rownames(res)
      res$cluster = this_clust
      big_res = rbind(res, big_res)
    } else { message("fastGlmm returned NULL.") }
  }
  message("===============================")
}
big_res$bh = p.adjust(big_res$P_exp, method = "BH")
big_res$hgnc = gene_info$seurat_name[match(big_res$gene, gene_info$seurat_name)]
write.csv(big_res, paste0(glmm_out_dir, "all.csv"))
deg_sig = big_res[which(big_res$bh < 0.05 & big_res$up_pct > 0.1),]
write.csv(big_res, paste0(glmm_out_dir, "sig.csv"))

