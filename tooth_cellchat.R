# Read Input ===================================================================
# this.run = 1; do.down = T; is.real = F; num.perms = 100;
# this.run = 1; do.down = F; is.real = T; num.perms = 1; ind = 0; cond.use = "all";
args = commandArgs(trailingOnly=TRUE)
this.run  = as.numeric(args[1])
do.down   = as.logical(args[2])
is.real   = as.logical(args[3])
cond.use  = as.character(argrs[4])
num.perms = as.numeric(args[5])
if (length(args) == 6) { ind = as.numeric(args[6]) } else { ind = 0 }
set.seed(this.run)
message(paste0("Initializng run with parameters: this.run=", this.run, ", do.down=", do.down, ", is.real=", is.real, ", cond.use=", cond.use, ", num.perms=", num.perms, ", ind=", ind, "."))

# Load Libraries ===============================================================
suppressMessages(library('CellChat',  quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('patchwork', quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('stringr',   quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('parallel',  quietly = T, warn.conflicts = F, verbose = F))
options(stringsAsFactors = FALSE)
source("~/scratch/bcs/bcs_scripts/bcs_f.R")

# Load Data ====================================================================
genePopFnc = function(x) {
  if (genePops$level[x] != "goi") {
    combined$this = switch(genePops$level[x], "primary" = combined$seuratclusters15, "secondary" = combined$seuratclusters53, "goi" = combined)
    this.cells = colnames(combined)[which(combined$this == genePops$cluster[x])]
  } else { this.cells = colnames(combined) }
  this.cells = this.cells[which(combined@assays$RNA@counts[genePops$mzebra[x], this.cells] > 0)]
  return(subset(combined, cells = this.cells))
}

message("Loading the object...")
# setwd("~/scratch/d_tooth/cellchat/")
gene_info = read.table("~/scratch/m_zebra_ref/gene_info_2.txt", header = T, stringsAsFactors = F)
combined = readRDS("~/scratch/d_tooth/data/plk_120922.rds")
combined$label = paste0("cluster.", combined$seurat_clusters)
message("Done.")

# Condition ====================================================================
if (cond.use == "all") {
  message("Using cells from both conditions")
} else if (cond.use == "plk") {
  message("Using cells from pluck")
  combined = subset(combined, cells = colnames(combined)[which(combined$cond == "plk")])
} else if (cond.use == "con") {
  message("Using cells from control")
  combined = subset(combined, cells = colnames(combined)[which(combined$cond == "con")])
}

# Downsample ===================================================================
if (do.down) {
  message("Doing Downsampling...")
  down_sample_fnc = function() { unlist(lapply(unique(combined$label), function(x) sample(colnames(combined)[which(combined$label == x)], 49) )) }
  down_cell_list  = mclapply(1:num.perms, function(x) down_sample_fnc(), mc.cores = 24)
  message("Done.")
} else { message("Not Downsampling.") }

# Permutations =================================================================
getLabels = function(x) {
  this.labels = c()
  if (is.real) {
    if (do.down) { this.labels = combined$label[down_cell_list[[x]]]         } else { this.labels = combined$label         }
  } else {
    if (do.down) { this.labels = sample(combined$label[down_cell_list[[x]]]) } else { this.labels = sample(combined$label) }
  }
  return (this.labels)
}

if (is.real) { message("Gather Cell Labels...") } else { message("Permuting Data...") }
label_list = mclapply(1:num.perms, function(x) getLabels(x), mc.cores = 24)
message("Done.")

# Human Object =================================================================
message("Creating a Human Object...")
gene_info_3 = read.csv("~/scratch/m_zebra_ref/gene_info_3.csv")
data.input = as.matrix(combined@assays$RNA@data[gene_info_3$seurat_name[which(!is.na(gene_info_3$one_to_one_human))],])
rownames(data.input) = gene_info_3$one_to_one_human[which(!is.na(gene_info_3$one_to_one_human))]
message("Done.")

rm(combined) # delete original Seurat object to save memory

# Cell Chat ====================================================================
CellChatWeights = function(x) {
  if (do.down) { this.cells = down_cell_list[[x]] } else { this.cells = colnames(data.input) }
  this.meta = data.frame(label = label_list[[x]], row.names = this.cells)
  
  cellchat = createCellChat(object = data.input[,this.cells], meta = this.meta, group.by = "label")
  cellchat = addMeta(cellchat, meta = this.meta)
  cellchat = setIdent(cellchat, ident.use = "label")
  
  cellchat@DB = CellChatDB.human
  cellchat = subsetData(cellchat)
  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)
  
  cellchat = computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1, population.size = F)
  df.net_lig_recept = subsetCommunication(cellchat) 
  cellchat = aggregateNet(cellchat)
  cellchat = filterCommunication(cellchat, min.cells = 10)
  net_weight = data.frame(cellchat@net$weight)
  net_weight_vect = unlist(net_weight)
  name_rep = rep(rownames(net_weight), ncol(net_weight))
  names(net_weight_vect) = paste0(name_rep, ".", sort(name_rep))
  return(net_weight_vect)
}

message("Running cellchat (this while take awhile)...")
# if (do.down) { num.parallel.jobs = 10 } else { num.parallel.jobs = 6 }
# if (do.down) { num.parallel.jobs = 10 } else { num.parallel.jobs = 1 }
if (do.down) { 
  num.parallel.jobs = 10
  if (is.real) { num.parallel.jobs = 5 } 
} else { num.parallel.jobs = 2 }
message(paste0("Using ", num.parallel.jobs, " cores."))
# onerun = suppressMessages(CellChatWeights(1))
sink(file="~/scratch/brain/cellchat_sink.txt")
# run_outs = mclapply(1:num.perms, function(x) suppressMessages(CellChatWeights(x)), mc.cores = num.parallel.jobs)
run_outs = Mclapply(1:num.perms, function(x) suppressMessages(CellChatWeights(x)), mc.cores = num.parallel.jobs)
sink()

n.success = length(run_outs)
if (n.success != num.perms) { message(paste0("Not all runs were successful (", (num.perms - n.success), "/", num.perms, ")")) }
out = as.data.frame(do.call('cbind', run_outs))
colnames(out) = paste0("run", 1:n.success)
out[, c("clust1", "clust2")] = reshape2::colsplit(names(run_outs[[1]]), "\\.", c("1", "2"))
out = out[, c(n.success+1, n.success+2, 1:n.success)]
message("Done.")

# Save Output ==================================================================
message("Writing Output...")
todays.date = stringr::str_split(Sys.Date(), pattern = "-")[[1]]
todays.date = paste0(todays.date[2], todays.date[3], substr(todays.date[1], 3, 4))
out.str = paste0("~/scratch/d_tooth/results/cellchat/output1/cellchat_", ifelse(do.down, "downsampled_", "full_"), ifelse(is.real, "real_", "perm_"), num.perms, "nruns_run", this.run, ".csv")
if (ind > 0) { out.str = paste0("~/scratch/d_tooth/results/cellchat/output1/cellchat_", ifelse(do.down, "downsampled_", "full_"), ifelse(is.real, "real_", "perm_"), num.perms, "nruns_run", this.run, "_ind", ind, ".csv") }
write.csv(out, out.str)
message("Done.")
message("All Done.")
