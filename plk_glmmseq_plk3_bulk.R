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
# source(paste0(brain_dir, "/brain_scripts/all_f.R"))
source(paste0(full_path, "/tooth_scripts/plk_glmmseq.R"))
setwd(out_dir)

.libPaths(c("/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/conda_envs/glmmseq/lib/R/library", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
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
message('Starting glmmseq analysis')
obj = subset(plk_subject, cells = colnames(plk_subject)[which(plk_subject$exp == "plk3")])
obj$pair = obj$subject
big_res = data.frame()
glmm_out_dir = "~/scratch/d_tooth/results/"

this_cells = colnames(obj)
message(paste0("Performing glmmSeq on bulk"))
my_formula = ~ cond + (1|subject)
res = fastGlmm(obj, this_cells, num_cores = 24, my_formula = my_formula, out_path = paste0(glmm_out_dir, "plk_glmmseq_plk3_bulk.csv")) 
