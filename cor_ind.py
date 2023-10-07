import pandas
import numpy as np
import scanpy
import random
import argparse
import h5py
import multiprocessing
from pathlib import Path
from itertools import repeat

global data_mat
global gene_labels
global cluster_labels
global do_abs
global cluster_set

def parseArgs():
    parser = argparse.ArgumentParser(description='Shuffle cluster labels and see if a gene has significantly greater node strength in the real vs perm.')
    parser.add_argument("dataset", metavar="dataset", type=str, help="Dataset (ct, cj, mi, mie, mim, hm)")
    parser.add_argument('perm_num', metavar='perm_num', type = int, help='The current permutation number. This is used for the seed.')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='The number of permutations to complete.')
    parser.add_argument("-o", "--output_folder", help="Output Folder", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/py_ns/",
                        const="/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/py_ns/")
    parser.add_argument("-n", "--no_perm", help="Do no permutations?", action="store_true")
    parser.add_argument("-a", "--do_abs", help="Take the absolute value of the correlations?", action="store_true")
    parser.add_argument("-r", "--cor_only", help="Only find correlations in the data?", action="store_true")
    args = parser.parse_args()
    return args.perm_num, args.num_perm, args.output_folder, args.no_perm, args.dataset, args.do_abs, args.cor_only

def corOnlyAndWrite(this_idx, output_path):
    """
    Given idexes of cells, create a matrix and find correlations only
    :param this_idx: Indexes of columns
    :param output_path: Output path of h5 correlation matrix file
    :return success: Function completed? True/False
    """
    cor = pandas.DataFrame(data=sparse_corrcoef(data_mat[:, this_idx].todense()), index=gene_labels, columns=gene_labels)
    if do_abs:
        print("Taking absolute value of correlations")
        cor = cor.abs()
    else:
        print("NOT taking absolute value of correlations. Using raw values.")
    h5f = h5py.File(output_path, 'w')
    h5f.create_dataset('name', data=cor)
    h5f.close()
    return True

def myShuffle(this_list):
    """
    Shuffle and return the list. (I made this function bc I don't like the how regular shuffle returns the list).
    :param this_list: input list (unshuffled)
    :return this_list: shuffled list
    """
    random.shuffle(this_list)
    return(this_list)

def corAndNodeStrength(this_idx):
    """
    Given idexes of cells, create a matrix, find correlations, and find node strengths
    :param i: i
    :param cluster: Cluster
    :return ns: Node Strength
    """
    # Find Correlations
    cor = pandas.DataFrame(data = sparse_corrcoef(data_mat[:, this_idx].todense()), index = gene_labels, columns = gene_labels)
    if do_abs:
        print("Taking absolute value of correlations")
        cor = cor.abs()
    else:
        print("NOT taking absolute value of correlations. Using raw values.")
    # Find Node Strength
    ns = cor.sum(axis=1)
    return ns

def sparse_corrcoef(A, B=None):
    """
    Find correlations in sparse matrix
    """
    if B is not None:
        A = sparse.vstack((A, B), format='csr')
    A = A.astype(np.float64)
    n = A.shape[1]
    # Compute the covariance matrix
    rowsum = A.sum(1)
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)
    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    coeffs = C / np.sqrt(np.outer(d, d))
    return coeffs

def permuteLabels(num_perm):
    """
    Permute cluster labels and find the indexes of the matrix belong to each.
    :param num_perm: Number of permutations
    :return mat_idx: dictionary of lists of column indexes that belong to BHVE and CTRL for each permutation
    """
    label_dict = {} # key is permutation #, value is a smaller dictionary
    for i in range(0, num_perm):
        perm_label = myShuffle(cluster_labels)
        small_dict = {}  # key is the cluster and value is indexes of the cluster
        for cluster in cluster_set:
            # Idx of labels equal to this cluster
            small_dict[cluster] = np.hstack(np.argwhere(perm_label == cluster))
        label_dict[i] = small_dict
    return label_dict

def permuteLabels2(num_perm):
    """
    Permute cluster labels and find the indexes of the matrix belong to each.
    :param num_perm: Number of permutations
    :return mat_idx: dictionary of lists of column indexes that belong to BHVE and CTRL for each permutation
    """
    label_dict = {} # key is permutation #, value is a smaller dictionary
    for cluster in cluster_set:
        perm_label = myShuffle(cluster_labels)
        small_dict = {}  # key is the cluster and value is indexes of the cluster
        for i in range(0, num_perm):
            # Idx of labels equal to this cluster
            small_dict[i] = np.hstack(np.argwhere(perm_label == cluster))
        label_dict[cluster] = small_dict
    return label_dict

def main():
    # Start the timer
    start_time = time.perf_counter()

    # Read Inputs
    global data_mat
    global gene_labels
    global cluster_labels
    global do_abs
    global cluster_set
    perm_num, num_perm, output_folder, no_perm, dataset, do_abs, cor_only = parseArgs()

# Find gene-gene correlations in all individuals
sct = scanpy.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plkall_sct_040523.h5ad")
sct = sct[sct.obs['subject_num'] != "NA",]
gene_labels = sct.var_names
cluster_labels = sct.obs['seurat_clusters']

# Version 1?? Weird merge problem right here
sct.obs['subject_cond'] = sct.obs['subject'] + '_' + sct.obs['cond']
df_list = []
for ind in np.sort(sct.obs['subject_cond'].unique()):
    print(ind)
    sct_ind = sct[sct.obs['subject_cond'] == ind,]
    sct_ind = sct[sct.obs['my_sub'] == ind,]
    this_mat = sct_ind.X.T
    cor = pandas.DataFrame(data=sparse_corrcoef(this_mat.todense()), index=sct.var_names, columns=sct.var_names)
    # le_zero_mask = np.array(cor[cor.le()] > 0)[0]
    # rows = cor.nonzero()[0][le_zero_mask]
    # cols = cor.nonzero()[1][le_zero_mask]
    # cor[rows, cols] = NaN
    df = cor.where(np.triu(np.ones(cor.shape)).astype(np.bool))
    df = df.stack().reset_index()
    df.columns = ['gene1', 'gene2', 'cor']
    df = df[df['cor'] > 0]
    df['id'] = df['gene1'] + '_' + df['gene2']
    df.index = df['id']
    df = df.drop(['gene1', 'gene2', 'id'], axis = 1)
    df_list.append(df['cor'])

# Version 2?? Weird merge problem right here
sct.obs['my_sub'] = sct.obs['subject'] + '_' + sct.obs['cond']
df_list = []
# dict_list = []
for ind in np.sort(sct.obs['my_sub'].unique()):
    print(ind)
    sct_ind = sct[sct.obs['my_sub'] == ind,]
    this_mat = sct_ind.X.T
    cor = pandas.DataFrame(data=sparse_corrcoef(this_mat.todense()), index=sct.var_names, columns=sct.var_names)
    # le_zero_mask = np.array(cor[cor.le()] > 0)[0]
    # rows = cor.nonzero()[0][le_zero_mask]
    # cols = cor.nonzero()[1][le_zero_mask]
    # cor[rows, cols] = NaN
    df = cor.where(np.triu(np.ones(cor.shape)).astype(np.bool))
    df = df.stack().reset_index()
    df.columns = ['gene1', 'gene2', 'cor']
    df = df[df['cor'] > 0]
    df['id'] = df['gene1'] + '_' + df['gene2']
    df.index = df['id']
    df = df.drop(['gene1', 'gene2', 'id'], axis = 1)
    df_list.append(df['cor'])

df2 = pandas.concat(df_list, axis=1, join = "inner")
df2.columns = np.sort(sct.obs['my_sub'].unique())
df2.to_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/all_ind_cor.csv')

df_list = []
plk60_inds = ['plk60_1_con', 'plk60_1_plk', 'plk60_2_con', 'plk60_2_plk', 'plk60_3_con', 'plk60_3_plk', 'plk60_4_con', 'plk60_4_plk']
for ind in plk60_inds:
    print(ind)
    sct_ind = sct[(sct.obs['subject_cond'] == ind) & (sct.obs['seurat_clusters'] == 19),]
    this_mat = sct_ind.X.T
    this_mat_bin = sct_ind.raw.X.T
    nonzero_idx = np.nonzero(this_mat_bin)
    this_mat_bin[nonzero_idx[0], nonzero_idx[1]] = 1
    this_mat_bin_counts = this_mat_bin.sum(axis=1)
    # this_mat_bin_counts[np.where(this_mat_bin_counts > 5)]
    # sct_ind.var_names[np.where(this_mat_bin_counts > 5)[0]]
    this_mat = this_mat[np.where(this_mat_bin_counts > 5)[0],]
    cor = pandas.DataFrame(data=sparse_corrcoef(this_mat.todense()), index=sct_ind.var_names[np.where(this_mat_bin_counts > 5)[0]], columns=sct_ind.var_names[np.where(this_mat_bin_counts > 5)[0]])
    df = cor.where(np.triu(np.ones(cor.shape)).astype(np.bool))
    df = df.stack().reset_index()
    df.columns = ['gene1', 'gene2', 'cor']
    df['id'] = df['gene1'] + '_' + df['gene2']
    df.index = df['id']
    df = df.drop(['gene1', 'gene2', 'id'], axis=1)
    df_list.append(df['cor'])

df2 = pandas.concat(df_list, axis=1, join="inner")
df2.columns = np.sort(plk60_inds)
df2.to_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/plk60_in_19_ind_cor.csv')

    cluster_set = list(set(cluster_labels))

    # Change folder name based on input
    base_name = dataset + "_perm_" + str(perm_num)
    if no_perm:
        base_name = dataset + "_real"
    if do_abs:
        base_name = base_name + "_abs"

    # Set random seed so all the permutations are different
    print("Seed = " + str(perm_num))
    random.seed(perm_num)

    # If necessary, find the correlations only (and don't do node strength stuff) and do no permutations
    if cor_only:
        print("Finding Correlations Only")
        base_name = base_name + "_cor"
        cor_output = output_folder + "/" + base_name + ".h5"
        cor_success = corOnlyAndWrite(np.array(range(0, data_mat.shape[1])), cor_output)
        print("Done")
        return

    # Permute Cluster Labels
    if no_perm:
        print("Not permuting Data")
        mat_idx = {}
        small_dict = {}  # key is the cluster and value is indexes of the cluster
        for cluster in cluster_set:
            # Idx of labels equal to this cluster
            small_dict[cluster] = np.hstack(np.argwhere(cluster_labels == cluster))
        mat_idx[0] = small_dict
    else:
        print("Permuting Data " + str(num_perm) + " times.")
        mat_idx = permuteLabels(num_perm)
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")

    # Set up dataframes that store results
    all_cluster_df = {}
    for cluster in cluster_set:
        all_cluster_df[cluster] = pandas.DataFrame(index = gene_labels, columns = list(range(1,num_perm+1)))

    # Find Correlations and Node Strengths
    print("Finding Correlations")
    for i in range(0, num_perm):
        print("Perm: " + str(i))
        perm_start = time.perf_counter()
        this_idx_list = mat_idx[i].values()
        with multiprocessing.Pool(len(cluster_set)) as pool:
            pool_ns = pool.map(corAndNodeStrength, this_idx_list)
            print(len(pool_ns))
            print(len(pool_ns[0]))
            # pool_ns = pool.starmap(corAndNodeStrength, zip(repeat(i), cluster_set))
            for j in range(0, len(cluster_set)):
                all_cluster_df[cluster_set[j]][i+1] = pool_ns[j]
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - perm_start:0.4f} seconds")
    # print("Finding Correlations")
    # for cluster in cluster_set:
    #     perm_start = time.perf_counter()
    #     this_idx_list = mat_idx[cluster]  # TODO
    #     with multiprocessing.Pool(24) as pool:
    #         pool_ns = pool.map(corAndNodeStrength, this_idx_list)
    #     print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - perm_start:0.4f} seconds")
    print(f"All Done. Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")

    # Write results for each cluster to a file
    for cluster in cluster_set:
        Path(output_folder + "/" + base_name + "/").mkdir(parents=True, exist_ok=True)
        all_cluster_df[cluster].to_csv(output_folder + "/" + base_name + "/" + cluster + "_" + str(perm_num) + ".csv")

if __name__ == '__main__':
    main()