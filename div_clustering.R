obj_path = "C:/Users/miles/Downloads/div_filtered_081423.rds"
obj = readRDS(obj_path)
library("Seurat")
min.dist = 0.2
n.neighbors = 30
res = 2
raw.assay = "RNA"
obj = SCTransform(obj, assay = raw.assay, vars.to.regress = "sample", verbose = F)
obj@active.assay = "SCT"
obj = RunPCA(obj, dim = 50, verbose = F)
obj = RunUMAP(obj, dims=1:50, min.dist=min.dist, spread=1, n.neighbors=n.neighbors, n.epochs=1000, metric = "euclidean")
this.graph.name = paste0("X", paste(2, min.dist, n.neighbors, sep = "_"))
obj = FindNeighbors(obj, reduction="umap", k.param=n.neighbors, dims=1:2, n.trees=500, prune.SNN=0, graph.name = this.graph.name)
obj = FindClusters(obj, resolution=res, algorithm = 2, graph.name = this.graph.name)

pdf("~/scratch/d_tooth/results/r6_plk_test.pdf", width = 7, height = 7)
DimPlot(plk2, label = T)
dev.off()

merged = ProjectUMAP(query=div, reference=plk2, query.reduction="pca", reference.reduction="pca", reduction.model = "umap")

anchors <- FindTransferAnchors(reference = plk2, query = div, normalization.method = "SCT", reference.reduction = "pca", dims = 1:50)
merged = MapQuery(anchorset = anchors, query = div, reference = plk2, refdata = list(celltype.l1 = "seurat_clusters", celltype.l2 = "seurat_clusters", predicted_ADT = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")
