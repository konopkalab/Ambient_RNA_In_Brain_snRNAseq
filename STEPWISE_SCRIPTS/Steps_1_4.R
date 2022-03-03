rm(list = ls())
library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)

#####
## STEP 01: LOAD RAW COUNT MATRIX
#####

# Load raw count matrix.
rawMat = readRDS('DATA/p56_rawmat.RDS')

# Create seurat object
rawMatSeur = CreateSeuratObject(counts = rawMat, min.cells = 10)
rawMatSeur$orig.ident = gsub('.*_YL', 'YL', colnames(rawMatSeur))


####
## STEP 02: CLUSTER EXCESS NUMBER OF CELL BARCODES
####

# Keep 20k cells per sample (~10k cells were targeted per sample)
samps = names(table(rawMatSeur$orig.ident))
keepL = list()
for(i in 1:length(samps)){

	subSeur = subset(rawMatSeur, subset = orig.ident == samps[i])
	keepL[[i]] = subSeur[[]][order(subSeur$nCount_RNA, decreasing = T)[1:20000], ] %>% rownames
}

keepV = unlist(keepL)
filtMatSeur = subset(rawMatSeur, cells = keepV)

# Dimensionality reduction and plotting without any batch correction
options(future.globals.maxSize = 4000 * 1024^2)
filtMatSeur = SCTransform(filtMatSeur, ncells = 1000)
filtMatSeur = RunPCA(filtMatSeur, verbose = FALSE)
filtMatSeur = RunUMAP(filtMatSeur, dims = 1:30)

pdf("P56_NoBatchCorrect.pdf", width = 12, height = 10)
DimPlot(filtMatSeur, group.by = 'orig.ident', label = T)
dev.off()

# Harmony batch correction 
library(harmony)
filtMatSeur_corr = RunHarmony(object = filtMatSeur, group.by.vars = 'orig.ident', assay.use="SCT")

# After batch correction
filtMatSeur_corr = RunUMAP(filtMatSeur_corr, dims = 1:30, reduction = 'harmony')
filtMatSeur_corr = FindNeighbors(filtMatSeur_corr, dims = 1:30, verbose = FALSE, reduction = 'harmony')
filtMatSeur_corr = FindClusters(filtMatSeur_corr, verbose = FALSE, resolution = 0.5)

stackedbarplot(filtMatSeur_corr[[]], 'seurat_clusters', 'orig.ident', 'P56_AfterBatchCorrect_Stack', wd = 15)

pdf("P56_AfterBatchCorrect_UMAP.pdf")
DimPlot(filtMatSeur_corr, group.by = 'orig.ident', label = T, pt.size = 0.1) + NoLegend()
dev.off()

pdf("P56_CLUSTERS.pdf")
DimPlot(filtMatSeur_corr, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

# Save (optional)
saveRDS(filtMatSeur_corr, 'DATA/p56_seurat_clustered.RDS')

filtMatSeur_corr = readRDS('DATA/p56_seurat_clustered.RDS')

####
## STEP 03: FIND CLUSTERS THAT MOSTLY CONTAIN CELL BARCODES
####

# Annotate cell barcodes by whether they belong to top 10k or bottom 10k by UMI count per sample
tmpmeta = filtMatSeur_corr[[]]
tmpmeta2 = tmpmeta %>% add_column(barc = rownames(tmpmeta)) %>% group_by(orig.ident) %>% dplyr::arrange(desc(nFeature_RNA)) %>% slice_head(n=10000) %>% as.data.frame
tmpmeta$cb_type = ifelse(rownames(tmpmeta) %in% tmpmeta2$barc, 'Top', 'Bottom')

# Per cluster, find the ratio of cell barcodes that belong to bottom 10k
clDF = tmpmeta$seurat_clusters %>% table %>% names %>% data.frame(clusters = .)
clDF$size = tmpmeta[,'seurat_clusters'] %>% table
clDF$nfSize = tmpmeta[tmpmeta$cb_type == 'Bottom', 'seurat_clusters'] %>% table
clDF$ratio = clDF$nfSize / clDF$size
clDF$clRatio = clDF$size / ncol(filtMatSeur_corr)

# Find clusters that have >75% of cells from the bottom 10k cell barcodes
ambCl = clDF[clDF$ratio > 0.75, 'clusters']
tmpmeta$clus_amb = as.character(tmpmeta$seurat_clusters)
tmpmeta$clus_amb = ifelse(tmpmeta$clus_amb %in% ambCl, 'Ambient_Cluster', 'Others')

# Plot ambient RNA clusters in UMAP
filtMatSeur_corr$clus_amb = tmpmeta$clus_amb
q = cbind(filtMatSeur_corr@reductions$umap@cell.embeddings, cb_type = filtMatSeur_corr$clus_amb) %>% as.data.frame
q$UMAP_1 = as.numeric(q$UMAP_1)
q$UMAP_2 = as.numeric(q$UMAP_2)

pdf("P56_Ambient_UMAP.pdf",)
ggscatter(q, x = 'UMAP_1', y = 'UMAP_2', alpha = 0.01, color = 'cb_type', palette = c('red', 'grey'))
dev.off()

####
## STEP 04_01: FIND INTRONIC READ RATIO
####

# Load count matrix computed only with intronic reads
hmatFinal = readRDS('DATA/p56_rawmat_intronic.RDS')

# Match the barcodes and calculate intronic read ratio
hmatFinal = hmatFinal[, colnames(filtMatSeur_corr)]
allMat = filtMatSeur_corr@assays$RNA@counts
filtMatSeur_corr$intronRat = colSums(hmatFinal) / colSums(allMat)

# Plot intronic read ratio with increasing UMI
meta = filtMatSeur_corr[[]]
meta$countInt = '>2000'

levsL = list()
increm = 100
for(i in seq(100,1900,increm)){
	meta[meta$nCount_RNA > i & meta$nCount_RNA < i+increm, 'countInt'] = paste0(i, '-', i+increm)
	levsL[[i]] = paste0(i, '-', i+increm)
}

meta$countInt = factor(meta$countInt, levels = c(unlist(levsL), '>2000'))

pdf("IntronicRatio_Intervals.pdf", width = 8)
ggviolin(meta, x = 'countInt', y = 'intronRat', color = 'grey', fill = 'grey', ylim = c(0,1)) +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(45) +
ylab('Intronic_Reads / All_Reads') + xlab('UMI Count') + NoLegend()
dev.off()

# Plot histogram of intronic read ratio
tmp = filtMatSeur_corr[[]]
metaD = tmp
metaD = metaD[metaD$clus_amb != 'Others',]
metaD$log10UMI = log10(metaD$nCount_RNA)

pdf('Histogram_IntronicRatio.pdf', width = 9)
gghistogram(metaD, x = "intronRat", add = "mean", rug = TRUE, fill = "clus_amb", palette = c("#00AFBB", "#E7B800"), bins = 100) +
 theme(text=element_text(size=20, face = 'bold')) +
 ylab('Count') + xlab('Intronic Read Ratio')
dev.off()

saveRDS(filtMatSeur_corr, 'DATA/p56_seurat_clustered.RDS')

# Assign high intronic and low intronic ambient RNAs based on the histogram
cutoff = 0.68
tmpmeta = filtMatSeur_corr[[]]
tmpmeta$clus_amb2 = 'Others'
tmpmeta[tmpmeta$intronRat < cutoff & tmpmeta$clus_amb == 'Ambient_Cluster', 'clus_amb2'] = 'Ambient_Low_Intronic'
tmpmeta[tmpmeta$intronRat >= cutoff & tmpmeta$clus_amb == 'Ambient_Cluster', 'clus_amb2'] = 'Ambient_High_Intronic'
filtMatSeur_corr$clus_amb = tmpmeta$clus_amb2


####
## STEP 04_02: FIND MARKERS FOR AMBIENT TYPES
####

# PSEUDOBULK DGE
library(scran)
library(scater)
library(SingleCellExperiment)

# Convert cells in ambient cluster types SCE object
DefaultAssay(filtMatSeur_corr) = 'RNA'
filtMatSeur_corr = NormalizeData(filtMatSeur_corr)
filtMatSeur_corr_sub = subset(filtMatSeur_corr, subset = clus_amb != 'Others')
seurmergSCE = as.SingleCellExperiment(filtMatSeur_corr_sub)

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups = colData(seurmergSCE)[, c("clus_amb", "orig.ident")]
pseudoSCE = sumCountsAcrossCells(seurmergSCE, groups)

# Create experimental design
pseudoSCE$orig.ident = factor(pseudoSCE$orig.ident)
pseudoSCE$clus_amb = factor(pseudoSCE$clus_amb)
pseudoSCE$TOTEST = factor(pseudoSCE$clus_amb)
pseudoSCE$tmpid = 'All'

# DGE analysis
resBulkL = pseudoBulkDGE(pseudoSCE, 
   label=pseudoSCE$tmpid,
   condition=pseudoSCE$TOTEST,
   design=~TOTEST,
   coef="TOTESTAmbient_Low_Intronic", method = 'edgeR')

resBulk = resBulkL[[1]]
resBulk = resBulk[!is.na(resBulk[,1]),]

resBulkSign_LowInt = resBulk[resBulk$logFC > 0.3 & resBulk$FDR < 0.05,]
resBulkSign_HighInt = resBulk[resBulk$logFC < -0.3 & resBulk$FDR < 0.05,]

resBulkSign_LowInt = resBulkSign_LowInt[order(resBulkSign_LowInt$FDR, decreasing = T),]
resBulkSign_HighInt = resBulkSign_HighInt[order(resBulkSign_HighInt$FDR, decreasing = T),]

resBulkSign_LowInt$Gene = rownames(resBulkSign_LowInt)
resBulkSign_HighInt$Gene = rownames(resBulkSign_HighInt)

rio::export(resBulkSign_LowInt, 'DATA/ambMark_PSEUDOBULK_LOWINTRON.xlsx')
rio::export(resBulkSign_HighInt, 'DATA/ambMark_PSEUDOBULK_HIGHINTRON.xlsx')


