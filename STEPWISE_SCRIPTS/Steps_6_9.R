rm(list = ls())
library(dplyr)
library(plyr)
library(Seurat)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(tidyr)
library(tidyverse)
library(harmony)
library(ggrepel)
library(GeneOverlap)

####
## STEP 06_01: CREATE SEURAT OBJECT AFTER CELLBENDER
####

# Samples
sp = c('YL01_S1', 'YL08_S8', 'YL10_S10', 'YL12_S12')

# Set your working directory
setwd('path/to/workdir')

hmatL = list()
for(i in 1:length(sp)){

	hmat = Read10X_h5(paste0('CELLBENDER/', sp[i], '/CellBender_out_filtered.h5'))
	colnames(hmat) = gsub('-1', paste0('_',sp[i]), colnames(hmat))
	
	hmatL[[i]] = hmat

}

# Merge matrices using genes detected in all samples
cmnGns = Reduce(intersect, lapply(hmatL, function(x){rownames(x)}))
cmnGns = cmnGns[!grepl('^NA\\.', cmnGns)]
cmnGns = cmnGns[cmnGns != '']
hmatL2 = lapply(hmatL, function(x){x[cmnGns,]})
hmatFinal = Reduce(cbind, hmatL2)


# Create seurat object
afterCBSeur = CreateSeuratObject(counts = hmatFinal,
				min.cells = 10,
				min.features = 0)
afterCBSeur$orig.ident = gsub('.*_YL', 'YL', colnames(afterCBSeur))

####
## STEP 06_02: CLUSTER AND ANNOTATE MAJOR CELL TYPES
####

afterCBSeur = subset(afterCBSeur, subset = nFeature_RNA > 0)

# Plot without any batch correction
options(future.globals.maxSize = 4000 * 1024^2)
afterCBSeur = SCTransform(afterCBSeur, ncells = 1000)
afterCBSeur = RunPCA(afterCBSeur, verbose = FALSE)

# Harmony batch correction
afterCBSeur_corr = RunHarmony(object = afterCBSeur, group.by.vars = 'orig.ident', assay.use="SCT")

# Re-compute the UMAP using batch corrected PCAs
afterCBSeur_corr = RunUMAP(afterCBSeur_corr, dims = 1:30, reduction = 'harmony')
afterCBSeur_corr = FindNeighbors(afterCBSeur_corr, dims = 1:30, verbose = FALSE, reduction = 'harmony')
DefaultAssay(afterCBSeur_corr) = 'SCT'
afterCBSeur_corr = FindClusters(afterCBSeur_corr, verbose = FALSE, resolution = 0.1)

pdf("P56_CB_HARMONY_SAMPLES.pdf")
DimPlot(afterCBSeur_corr, group.by = 'orig.ident', label = T, pt.size = 0.1) + NoLegend()
dev.off()

pdf("P56_CB_HARMONY_CLUSTERS.pdf")
DimPlot(afterCBSeur_corr, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()


# Check marker genes
DefaultAssay(afterCBSeur_corr) = 'RNA'
afterCBSeur_corr = NormalizeData(afterCBSeur_corr)

pdf("P56_CB_VLNPLOT.pdf", width = 15)
VlnPlot(afterCBSeur_corr, features = c('Mog', 'Mbp', 'Pcdh15', 'Ptprz1', 'Slc1a3', 'Slc1a2', 'Apbb1ip', 'Slc17a7', 'Gad1', 'Flt1', 'Cobll1'), pt.size = 0) + NoLegend()
dev.off()

pdf("P56_CB_FEATURE_COUNT_PER_CLUSTER.pdf", width = 15)
ggboxplot(afterCBSeur_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA') + NoLegend()
dev.off()

# Broadly annotate based on marker genes
mapnames = setNames(c('Excitatory', 'Excitatory', 'Astrocyte', 'Excitatory', 'Inhibitory', 'Microglia',
			'Oligodendrocyte', 'Inhibitory', 'Excitatory', 'OPC', 'Astrocyte',
			'Excitatory', 'Excitatory', 'Endothelia', 'Potential_Doublet', 'Unknown'),
		      c(0:15))

afterCBSeur_corr[["broadannot"]] = mapnames[afterCBSeur_corr[["seurat_clusters"]][,1]]
afterCBSeur_corr$broadannot = factor(afterCBSeur_corr$broadannot)

# Remove unknown, endothelia and potential doublet
afterCBSeur_corr = subset(afterCBSeur_corr, subset = broadannot %in% c('Unknown', 'Potential_Doublet', 'Endothelia'), invert = T)
afterCBSeur_corr = RunUMAP(afterCBSeur_corr, dims = 1:30, reduction = 'harmony')

pdf('P56_CB_ANNOTATED.pdf')
DimPlot(afterCBSeur_corr, group.by = 'broadannot', label = T, raster = T) + NoLegend()
dev.off()

# Optional: re-cluster and re-annotate

# Save after annotation
saveRDS(afterCBSeur_corr, 'DATA/afterCBSeur_annotated_CellBender.RDS')


####
## STEP 07: SUBCLUSTER EACH MAJOR CELL TYPE (ONLY ASTROCYTES ARE SHOWN HERE)
####

# Keep only given cell type
p56CBAst = subset(afterCBSeur_corr, subset = broadannot == 'Astrocyte' & nFeature_RNA > 0)
p56Ast = subset(afterCBSeur, cells = colnames(p56CBAst))

# Dimensionality reduction
options(future.globals.maxSize = 4000 * 1024^2)
p56Ast = SCTransform(p56Ast, ncells = 1000)
p56Ast = RunPCA(p56Ast, verbose = FALSE)

# Harmony batch correction 
p56Ast_corr = RunHarmony(object = p56Ast, group.by.vars = 'orig.ident', assay.use="SCT")

# Compute UMAP after batch correction
p56Ast_corr = RunUMAP(p56Ast_corr, dims = 1:15, reduction = 'harmony')
p56Ast_corr = FindNeighbors(p56Ast_corr, dims = 1:15, verbose = FALSE, reduction = 'harmony')
p56Ast_corr = FindClusters(p56Ast_corr, verbose = FALSE, resolution = 0.5)

pdf("P56_HARMONY_SAMPLES_CB_AST.pdf")
DimPlot(p56Ast_corr, group.by = 'orig.ident', pt.size =1, label = T)
dev.off()

pdf("P56_HARMONY_CLUSTERS_CB_AST.pdf")
DimPlot(p56Ast_corr, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

# Find cluster markers
DefaultAssay(p56Ast_corr) = 'RNA'
p56Ast_corr = NormalizeData(p56Ast_corr)
p56Ast_corrM = FindAllMarkers(p56Ast_corr, only.pos = T, logfc.threshold = 0.25)
p56Ast_corrM = p56Ast_corrM[p56Ast_corrM$p_val_adj < 0.05,]
p56Ast_corrM = p56Ast_corrM %>% group_by(cluster) %>% slice_head(n=100) %>% as.data.frame

p56Ast_corrL = p56Ast_corrM[, c('cluster', 'gene')]
colnames(p56Ast_corrL) = c('CellType', 'Gene')
p56Ast_corrL = split(p56Ast_corrL, p56Ast_corrL$CellType)
p56Ast_corrL = lapply(p56Ast_corrL, function(x){x$Gene})


####
## STEP 08: RUN ENRICHMENT BETWEEN SUBCLUSTER MARKERS AND AMBIENT RNA MARKERS
####

# Load ambient markers
resBulkSign_LowInt = rio::import('DATA/ambMark_PSEUDOBULK_LOWINTRON.xlsx')
resBulkSign_HighInt = rio::import('DATA/ambMark_PSEUDOBULK_HIGHINTRON.xlsx')

# Take top 500 most abundant markers
resBulkSign_LowIntM = resBulkSign_LowInt[order(resBulkSign_LowInt$logCPM, decreasing = T)[1:500], 'Gene']
resBulkSign_HighIntM = resBulkSign_HighInt[order(resBulkSign_HighInt$logCPM, decreasing = T)[1:500], 'Gene']
ambList = union(resBulkSign_LowIntM, resBulkSign_HighIntM)

# Run enrichment by fisher's exact test
resgom = newGOM(list(ambList), p56Ast_corrL, genome.size=13000)
oddsr = getMatrix(resgom, name="odds.ratio")
pvals = getMatrix(resgom, name="pval")
jac = getMatrix(resgom, name="Jaccard")

# Reshape and combine to plot
oddsrM = melt(oddsr)
pvalsM = melt(pvals)
pvalsM$log10FDR = -log10(p.adjust(as.numeric(pvalsM$value)))

toplot = cbind(oddsrM, log10FDR = pvalsM$log10FDR)
toplot$Var1 = rownames(toplot)
toplot$type = 'Ambient_Markers'
toplot$FDRplot = formatC(p.adjust(as.numeric(pvalsM$value)), format = "e", digits = 2)

pdf('P56_AMBIENT_ENRICH_ASTROCYTES.pdf')
ggscatter(toplot, x = 'log10FDR', y = 'value', color = 'black', size = 6) +
 theme_classic() +
 xlab(expression('-log'[10]*'(FDR)')) + ylab('Odds Ratio') +
 geom_label_repel(data = toplot[toplot$log10FDR > 30,], aes(label = Var1), nudge_y = 0.02, nudge_x = -0.05, fontface = 'bold', size = 6) +
 theme(text=element_text(size=20, face = 'bold'))
dev.off()


####
## STEP 09: REMOVE HIGHLY ENRICHED SUBCLUSTERS
####

# Remove and save
p56Ast_corr = subset(p56Ast_corr, subset = seurat_clusters == 6, invert = T)
saveRDS(p56Ast_corr, 'DATA/p56_Ast_cellbendered_cleaned.RDS')





