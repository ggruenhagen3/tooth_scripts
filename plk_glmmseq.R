#' Fast GlmmSeq Analysis
#'
#' @author George Gruehnagen (Georgia Institute of Technology) and Zachary Johnson (Emory)
#'
#' @param obj Seurat object
#' @param cells cells to run the analysis on
#' @param subset.genes.by genes to subset results by
#' @param sequenced.genes the gene_list paramater from dndscv, which is a list of genes to restrict the analysis (use for targeted sequencing studies)
#' @param ... other parameters passed to dncdscv, all defaults from dndscv are used except max_muts_per_gene_per_sample is set to Infinity
#'
#' @return coselens returns a list containing six objects: 1) "substitutions", a summary table of gene-level conditional selection for single-nucleotide substitutions (including missense, nonsense, and essential splice site mutations); 2) "indels", same for small indels; (3) "missense_sub", same for missense substitutions only; (4) "truncating_sub", same for truncating substitutions only (nonsense and essential splice site mutations); (5) "overall_mut", a summary table with the combined analysis of single-nucleotide substitutions and indels; and (6) "dndscv", a list of objects with the complete output of (non-conditional) selection analyses separately run on each group of samples, as provided by the dndscv package. The first table should be sufficient for most users. If interested in indels, please note that the indel analysis uses a different null model that makes the test for conditional selection notably less sensitive than in the case of substitutions. Such lower sensitivity also extends to the "overall_mut" table. The dataframes (1-5) contain the following:
#' @return - gene_name: name of the gene that was tested for conditional selection
#' @return - num.driver.group1: estimate of the number of drivers per sample per gene in group 1
#' @return - num.driver.group2: estimate of the number of drivers per sample per gene in group 2
#' @return - Delta.Nd: absolute difference in the average number of driver mutations per sample (group 1 minus group 2)
#' @return - classification: classification of conditional selection. The most frequent classes are strict dependence (drivers only in group 1), facilitation (drivers more frequent in group 1), independence, inhibition (drivers less frequent in group 1), and strict inhibition (drivers absent from group 1). If negative selection is present, other possibilities are strict dependence with sign change (drivers positively selected in group 1 but negatively selected in group 2), strict inhibition with sign change (drivers positively selected in group 2 but negatively selected in group 1), aggravation (purifying selection against mutations becomes stronger in group 1), and relaxation (purifying selection against mutations becomes weaker in group 1).
#' @return - dependency: dependency index, measuring the association between the grouping variable (group 1 or 2) and the average number of drivers observed in a gene. It serves as a quantitative measure of the qualitative effect described in "classification". In the most common cases, a value of 1 indicates strict dependence or inhibition (drivers only observed in one group) and a value of 0 (or NA) indicates independence.
#' @return - pval: p-value for conditional selection
#' @return - qval: q-value for conditional selection using Benjamini-Hochberg correction of false discovery rate.
#'
#' @return The "dndscv" list contains two objects. Please, read the documentation of the dndscvpackage for further information about th
#' @return - dndscv_group1: output of dndscv for group 1
#' @return - dndscv_group2: output of dndscv for group 2
#'
#' @export
fastGlmm = function(obj, cells, my_formula, num_cores = 24, out_path = NULL) {
  this_counts = obj@assays$RNA@counts[,cells]
  this_meta   = data.frame(obj@meta.data[cells,])
  if (!"pair" %in% colnames(this_meta)) { message("Required metadata column 'pair' not found in metadata. Exiting."); return(NULL); }
  
  gene_not_present_in_pairs = lapply(unique(this_meta$pair), function(x) rowSums(this_counts[,which(this_meta$pair == x)]) == 0)
  gene_not_present_in_pairs = Reduce(`+`, gene_not_present_in_pairs)
  genes_present_in_all_pairs = names(gene_not_present_in_pairs[which(gene_not_present_in_pairs == 0)])
  this_counts = this_counts[genes_present_in_all_pairs,]
  
  this_max_cluster_size = max(table(this_meta$seurat_clusters)) # TODO make clusters as a variable
  multicoreParam <- MulticoreParam(workers = num_cores)
  dds = DESeqDataSetFromMatrix(countData = this_counts, colData = this_meta, design = ~ cond)
  size_factors = calculateSumFactors(this_counts, BPPARAM = multicoreParam, max.cluster.size = this_max_cluster_size, clusters = NULL, ref.clust = NULL, positive = TRUE, scaling = NULL,  min.mean = NULL, subset.row = NULL)
  sizeFactors(dds) = size_factors
  dds = estimateDispersions(dds, fitType = "glmGamPoi", useCR = TRUE, maxit = 100, weightThreshold = 0.01, quiet = FALSE, modelMatrix = NULL, minmu = 1e-06)
  disp = as.matrix(mcols(dds))
  disp = disp[,11]
  names(disp) = genes_present_in_all_pairs
  
  results <- glmmSeq(my_formula, id = "subject",
                     countdata = this_counts,
                     metadata = this_meta,
                     dispersion = disp,
                     sizeFactors = size_factors,
                     removeSingles=FALSE,
                     progress=TRUE,
                     cores = num_cores)
  results_df = data.frame(summary(results))
  if(!is.null(out_path)) { write.csv(results_df, out_path) }
  return(results_df)
}