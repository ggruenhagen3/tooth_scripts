import scanpy
import numpy as np
import pandas as pd
from numba import njit
import numba
from numba import prange
from sklearn.preprocessing import LabelBinarizer
from scipy import sparse
import seaborn as sns
import matplotlib.pyplot as plt

adata = scanpy.read_h5ad("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plk_div_081823.h5ad")
umap_coord = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plk_div_081823_umap.csv", index_col=0)
plk_subjects = ['plk7_1', 'plk7_2', 'plk7_NA', 'plk60_2', 'plk60_1', 'plk60_NA', 'plk60_3', 'plk60_4', 'plk1_2', 'plk1_1', 'plk1_NA', 'plk3_2', 'plk3_NA', 'plk3_1']
adata.obs['isPlk'] = adata.obs['subject'].isin(plk_subjects)
adata.obs['plk_col'] = adata.obs['X2_0.2_30_res.0.6'].astype("str")
adata.obs.loc[adata.obs['plk_col'] == "NA", "plk_col"] = "unassigned"
adata.obs['div_col'] = adata.obs['X1_0.2_40_res.2'].astype("str")
adata.obs.loc[adata.obs['div_col'] == "NA", "div_col"] = "unassigned"
# adata.obs['plk_col'] = "unassigned"
# adata.obs.loc[adata.obs['isPlk'], "plk_col"] = adata.obs.loc[adata.obs['isPlk'], "seurat_clusters"].astype("string")
# adata.obs['div_col'] = "unassigned"
# adata.obs.loc[~adata.obs['isPlk'], "div_col"] = adata.obs.loc[~adata.obs['isPlk'], "seurat_clusters"].astype("string")

plk_names = adata.obs.loc[  adata.obs['isPlk']].index
div_names = adata.obs.loc[~ adata.obs['isPlk']].index
dist = euclidean_distance2(umap_coord.loc[plk_names].to_numpy(), umap_coord.loc[div_names].to_numpy())
dist_sparse = sparse.csr_matrix(dist)
knn = dist_sparse
mean_dist = my_mapper(mz_col = 'plk_col', mm_col = 'div_col')
to_p = 1/mean_dist
to_p.loc[to_p > 0.03] = 0.03
sns.clustermap(to_p, annot=False, cmap = "viridis", figsize = [mean_dist.shape[1]/4.5, mean_dist.shape[0]/4.5], dendrogram_ratio=(.05), cbar_pos = None, zscore = 0)
plt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/div_to_plk_inverse_z.png")
mean_dist.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/div_plk_mean_dist.csv")

@njit(parallel=True)
def euclidean_distance2(coords1, coords2):
    c1_length, c2_length = len(coords1), len(coords2)
    out = np.empty(shape=(c1_length, c2_length), dtype=np.float64)
    for lat_ix in prange(c1_length):
        for lon_ix in prange(c2_length):
            out[lat_ix, lon_ix] = ( (coords1[lat_ix, 0] - coords2[lon_ix, 0]) ** 2 + (coords1[lat_ix, 1] - coords2[lon_ix, 1]) ** 2 ) ** 0.5
    return out


def my_mapper(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
    meta = adata.obs
    all_mz_cluster = meta[mz_col].unique()
    all_mm_cluster = meta[mm_col].unique()
    all_mz_cluster.sort()
    all_mm_cluster.sort()
    lb = LabelBinarizer(sparse_output=True)
    knn_mz  = lb.fit_transform(adata.obs.loc[adata.obs['isPlk'], mz_col]).T.dot(knn)
    knn_sum = pd.DataFrame(lb.fit_transform(adata.obs.loc[~adata.obs['isPlk'], mm_col]).T.dot(knn_mz.T).todense())
    mz_counts = meta[mz_col].value_counts().sort_index()
    mm_counts = meta[mm_col].value_counts().sort_index()
    mz_mm_counts = pd.DataFrame(np.array(mz_counts) * np.array(mm_counts).reshape(len(mm_counts), 1))
    knn_mean = knn_sum / mz_mm_counts
    knn_mean.index = all_mm_cluster
    knn_mean.columns = all_mz_cluster
    knn_mean = knn_mean.drop("unassigned", axis=0)
    knn_mean = knn_mean.drop("unassigned", axis=1)
    return(knn_mean)


# def square_to_condensed(i, j, n):
#     assert i != j, "no diagonal elements in condensed matrix"
#     if i < j:
#         i, j = j, i
#     return n*j - j*(j+1)//2 + i - 1 - j
#
# c = list(itertools.product(np.where(sm.samap.adata.obs['mz_good_names'] == "8.1_Glut")[0], np.where(sm.samap.adata.obs['mm_ClusterName'] == "TEGLU6")[0]))
# d, e = zip(*c)
# this_idx = square_to_condensed(d, e, sm.samap.adata.obsm['X_umap'].shape[0])
# with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
#     tmp = mp_pool.starmap(square_to_condensed, zip(d, e, sm.samap.adata.obsm['X_umap'].shape[0] * len(d)))
#
# this_idx = square_to_condensed(np.where(sm.samap.adata.obs['mz_good_names'] == "8.1_Glut"), np.where(sm.samap.adata.obs['mm_ClusterName'] == "TEGLU6"), sm.samap.adata.obsm['X_umap'].shape[0])
