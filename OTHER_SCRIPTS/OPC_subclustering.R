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

# Note: this script is only for NSD1 dataset.
# But other datasets are processed with the same pipeline.

####
## CLUSTER NSD1
####

# Load CellBender + Clean OPCs.
seurOPC = readRDS('NSD1_OPC_Seurat.RDS')

# Plot without any batch correction
options(future.globals.maxSize = 4000 * 1024^2)

seurOPC = SCTransform(seurOPC, ncells = 1000)
seurOPC <- RunPCA(seurOPC, verbose = FALSE)
seurOPC <- RunUMAP(seurOPC, dims = 1:15)

# Harmony batch correction
library(harmony)
seurOPC_corr <- RunHarmony(object = seurOPC, group.by.vars = 'orig.ident', assay.use="SCT")

# re-compute the UMAP using batch corrected PCAs
seurOPC_corr <- RunUMAP(seurOPC_corr, dims = 1:5, reduction = 'harmony')
seurOPC_corr <- FindNeighbors(seurOPC_corr, dims = 1:5, verbose = FALSE, reduction = 'harmony')
seurOPC_corr <- FindClusters(seurOPC_corr, verbose = FALSE, resolution = 0.5)

pdf("NSD1_HARMONY_SAMPLES_OPC.pdf")
DimPlot(seurOPC_corr, group.by = 'orig.ident', pt.size =2, label = T)
dev.off()

pdf("NSD1_HARMONY_CLUSTERS_OPC.pdf")
DimPlot(seurOPC_corr, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

pdf("NSD1_CLUSTER_DEPTH_OPC.pdf")
ggboxplot(seurOPC_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA', ylim = c(0,3000), outlier.shape = NA) +
rotate_x_text(90)
dev.off()

# Find cluster markers
DefaultAssay(seurOPC_corr) = 'RNA'
seurOPC_corr = NormalizeData(seurOPC_corr)
nsd1M = FindAllMarkers(seurOPC_corr, only.pos = T, logfc.threshold = 0.25)

# Annotate
seurOPC_corr$newannot = ifelse(seurOPC_corr$seurat_clusters == 2, 'COP', 'OPC')

pdf("NSD1_HARMONY_ANNOTATION_OPC_FINAL.pdf")
DimPlot(seurOPC_corr, group.by = 'newannot', label = T, raster = T, label.size = 8, cols = c('lightgreen','lightblue')) + NoLegend()
dev.off()

# Find and save final markers
nsd1M = FindMarkers(seurOPC_corr, only.pos = F, logfc.threshold = 0, ident.1 = 'COP', min.pct = 0, group.by = 'newannot')
nsd1M$FDR = p.adjust(nsd1M$p_val, method = 'fdr')
nsd1M$Gene = rownames(nsd1M)

rio::export(nsd1M, 'NSD1_COP_Markers.xlsx')


