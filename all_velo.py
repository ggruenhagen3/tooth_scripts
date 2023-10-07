"""
RNA Velocity
"""
# cond activate velo3
import scvelo as scv
import scanpy
import loompy
import pandas as pd
from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join

# Things to run the first time (commented out to save time)
# mypath = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plk_velo/"
# loomfiles = [mypath+f for f in listdir(mypath) if isfile(join(mypath, f))]
# all_velo_file = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/all_velo.loom"
# loompy.combine(loomfiles, all_velo_file)

# Read in the velo files for all the data
all_velo_file = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/all_velo.loom"
adata = scv.read(all_velo_file) # I get errors with cache = True
adata.obs['file'] = adata.obs.index.str.split(':').str[0]
adata.obs['barcode'] = adata.obs.index.str.split(':').str[1].str.split('x').str[0]
file_to_sample_dict = dict(zip(['JTS15_cnt12_bcl_cr7', 'JTS15_clp12_bcl_cr7', 'JTS17_Cont03_1_2_bcl', 'JTS17_Plk03_1_2_bcl', 'JTS17_Cont01_1_2_bcl', 'JTS17_Plk01_1_2_bcl', 'JTS16_c12_bcl', 'JTS16_c34_bcl', 'JTS16_p12_bcl', 'JTS16_p34_bcl'], ['plk7_c12', 'plk7_p12', 'plk3_c12', 'plk3_p12', 'plk1_c12', 'plk1_p12', 'plk60_c12', 'plk60_c34', 'plk60_p12', 'plk60_p34']))
adata.obs['sample'] = adata.obs['file'].replace(file_to_sample_dict, regex=True)
adata.obs['my_barcode'] = adata.obs['sample'] + '_' + adata.obs['barcode']

# R code:
# plk.loom <- as.loom(plk, filename = "~/scratch/d_tooth/data/plkall_053023.loom", verbose = FALSE)
# plk.loom$close_all()

# Add seurat metadata in
seurat_data = scv.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plkall_053023.loom")
seurat_data.obs['my_barcode'] = seurat_data.obs.index.str.split('-').str[0]
adata = adata[adata.obs['my_barcode'].isin(seurat_data.obs['my_barcode']),]
adata.obs = pd.merge(adata.obs, seurat_data.obs, on='my_barcode', how='inner')

# Process Data and Perform RNA Velocity
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata, n_jobs = 20)  # had to downgrade numpy. "Downgrading from numpy-1.24.3 to numpy-1.23.5 fixes the problem."
scv.tl.umap(adata)
adata.write_h5ad("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/all_velo2.h5ad")

# had to change a line of code see: https://github.com/theislab/scvelo/issues/811
adata = scanpy.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/all_velo2.h5ad")
scv.pl.velocity_embedding_stream(adata, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/all_velo_stream2.svg")
# scv.pl.velocity_embedding_stream(adata, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/all_velo_stream.svg")

# Have to run first on all before you can plot them separately
plk60 = adata[adata.obs['exp'] == "plk60",]
plk1  = adata[adata.obs['exp'] == "plk1",]
plk3  = adata[adata.obs['exp'] == "plk3",]
plk7  = adata[adata.obs['exp'] == "plk7",]
scv.pl.velocity_embedding_stream(plk60, color='cond', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_velo_stream2_by_cond.svg")
scv.pl.velocity_embedding_stream(plk1,  color='cond', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_velo_stream2_by_cond.svg")
scv.pl.velocity_embedding_stream(plk3,  color='cond', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_velo_stream2_by_cond.svg")
scv.pl.velocity_embedding_stream(plk7,  color='cond', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_velo_stream2_by_cond.svg")

plk60_con = adata[(adata.obs['exp'] == "plk60") & (adata.obs['cond'] == 'con'),]
plk60_plk = adata[(adata.obs['exp'] == "plk60") & (adata.obs['cond'] == 'plk'),]
plk1_con  = adata[(adata.obs['exp'] == "plk1")  & (adata.obs['cond'] == 'con'),]
plk1_plk  = adata[(adata.obs['exp'] == "plk1")  & (adata.obs['cond'] == 'plk'),]
plk3_con  = adata[(adata.obs['exp'] == "plk3")  & (adata.obs['cond'] == 'con'),]
plk3_plk  = adata[(adata.obs['exp'] == "plk3")  & (adata.obs['cond'] == 'plk'),]
plk7_con  = adata[(adata.obs['exp'] == "plk7")  & (adata.obs['cond'] == 'con'),]
plk7_plk  = adata[(adata.obs['exp'] == "plk7")  & (adata.obs['cond'] == 'plk'),]
scv.pl.velocity_embedding_stream(plk60_con, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_con_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk60_plk, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_plk_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk1_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_con_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk1_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_plk_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk3_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_con_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk3_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_plk_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk7_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_con_velo_stream2.svg")
scv.pl.velocity_embedding_stream(plk7_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_plk_velo_stream2.svg")

scv.pl.velocity_embedding_grid(plk60_con, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_con_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk60_plk, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_plk_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk1_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_con_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk1_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_plk_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk3_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_con_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk3_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_plk_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk7_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_con_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")
scv.pl.velocity_embedding_grid(plk7_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_plk_velo_grid.svg", arrow_size = 2, arrow_length = 2, legend_loc = 'on data', legend_fontweight = "light")

scv.pl.velocity_embedding_grid(plk60_con, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_con_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk60_plk, color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk60_plk_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk1_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_con_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk1_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk1_plk_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk3_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_con_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk3_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk3_plk_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk7_con,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_con_velo_grid2.svg")
scv.pl.velocity_embedding_grid(plk7_plk,  color='seurat_clusters', save = "/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/plk7_plk_velo_grid2.svg")

# Old velocity code
# merged = scv.utils.merge(adata, seurat_data)
# scv.pp.filter_and_normalize(adata)
# scv.pp.moments(adata)
# scv.tl.velocity(adata, mode='stochastic')