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
library(celda)
library(SoupX)

#####
## PREPARE THE DATA
#####

# Read NeuN- SDs
hodgeSeur = readRDS('hodgeSeur.rds')
bakkenSeur = readRDS('bakkenSeur.RDS')
sadSeur = readRDS('sadSeur.RDS')

# Read datasets
velmPub = readRDS('velmSeur.RDS')
velmNOCB = readRDS('seurat_raw.RDS')
velmCB = readRDS('seurat_raw_CellBender.RDS')

# Add sample ids and match to published barcode ids
srrMeta = read.table('SRA_Table_10X.txt', sep = '\t', header = T)
srrMeta$Library.Name = gsub('_Nova', '', srrMeta$Library.Name)

meta = velmNOCB[[]]
meta$Sample = srrMeta[match(meta$orig.ident, srrMeta$Run), 'Library.Name']
velmNOCB$Sample = meta$Sample
velmNOCB$authBarc = paste0(gsub('_.*', '', colnames(velmNOCB)), '-1_', velmNOCB$Sample)

meta = velmCB[[]]
meta$Sample = srrMeta[match(meta$orig.ident, srrMeta$Run), 'Library.Name']
velmCB$Sample = meta$Sample
velmCB$authBarc = paste0(gsub('_.*', '', colnames(velmCB)), '-1_', velmCB$Sample)

# Keep same barcodes in all datasets
cmnBarcs = Reduce(intersect, list(velmCB$authBarc, velmNOCB$authBarc, colnames(velmPub)))
seurmergNOCB_Filt = subset(velmNOCB, subset = authBarc %in% cmnBarcs)
seurmergCB_Filt = subset(velmCB, subset = authBarc %in% cmnBarcs)
seurpubFilt = subset(velmPub, cells = cmnBarcs)

# Carry annotations
pubmeta = seurpubFilt[[]]
seurmergNOCB_Filt$cluster = pubmeta[match(seurmergNOCB_Filt$authBarc, rownames(pubmeta)), 'cluster']
seurmergCB_Filt$cluster = pubmeta[match(seurmergCB_Filt$authBarc, rownames(pubmeta)), 'cluster']


####
## RUN DECONTX AND SOUPX
####

# Run decontx
seurmergNOCB_SCE <- as.SingleCellExperiment(seurmergNOCB_Filt)
decontSCE <- decontX(seurmergNOCB_SCE)

saveRDS(decontSCE, 'decontX_SCE.RDS')

# Create the same seurat object on decontaminated matrix
decontSeur = CreateSeuratObject(decontXcounts(decontSCE), meta.data = seurmergNOCB_Filt[[]])

# Create soupX object
toc = seurmergNOCB_Filt@assays$RNA@counts # filtered
tod = velmNOCB@assays$RNA@counts # raw
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(seurmergNOCB_Filt$cluster, rownames(seurmergNOCB_Filt[[]])))

# Run soupX
sc = autoEstCont(sc)
out = adjustCounts(sc)

# Create the same seurat object on soupX decontaminated matrix
soupSeur = CreateSeuratObject(out, meta.data = seurmergNOCB_Filt[[]])
saveRDS(soupSeur, 'soupSeur.RDS')

####
## MERGE GLIA
####

# Give the same labels in all datasets
decontSeur$cluster = decontSeur$newannot
decontSeur$cluster = gsub('Oligodendrocytes', 'Oli', decontSeur$cluster)
decontSeur$cluster = gsub('OPC', 'OPC', decontSeur$cluster)
decontSeur$cluster = gsub('Microglia', 'Mic', decontSeur$cluster)
decontSeur$cluster = gsub('AST-FB|AST-PP', 'Ast', decontSeur$cluster)

soupSeur$cluster = gsub('Oligodendrocytes', 'Oli', soupSeur$cluster)
soupSeur$cluster = gsub('OPC', 'OPC', soupSeur$cluster)
soupSeur$cluster = gsub('Microglia', 'Mic', soupSeur$cluster)
soupSeur$cluster = gsub('AST-FB|AST-PP', 'Ast', soupSeur$cluster)

seurmergNOCB_Filt$cluster = gsub('Oligodendrocytes', 'Oli', seurmergNOCB_Filt$cluster)
seurmergNOCB_Filt$cluster = gsub('OPC', 'OPC', seurmergNOCB_Filt$cluster)
seurmergNOCB_Filt$cluster = gsub('Microglia', 'Mic', seurmergNOCB_Filt$cluster)
seurmergNOCB_Filt$cluster = gsub('AST-FB|AST-PP', 'Ast', seurmergNOCB_Filt$cluster)

seurmergCB_Filt$cluster = gsub('Oligodendrocytes', 'Oli', seurmergCB_Filt$cluster)
seurmergCB_Filt$cluster = gsub('OPC', 'OPC', seurmergCB_Filt$cluster)
seurmergCB_Filt$cluster = gsub('Microglia', 'Mic', seurmergCB_Filt$cluster)
seurmergCB_Filt$cluster = gsub('AST-FB|AST-PP', 'Ast', seurmergCB_Filt$cluster)


hodgeSeur$cluster = hodgeSeur$cluster
bakkenSeur$cluster = bakkenSeur$cluster_label
sadSeur$cluster = sadSeur$newannot

hodgeSeur$cluster = gsub('Oligo_L1-6_OPALIN_Oligo', 'Oli', hodgeSeur$cluster)
hodgeSeur$cluster = gsub('OPC_L1-6_PDGFRA_OPC', 'OPC', hodgeSeur$cluster)
hodgeSeur$cluster = gsub('Micro_L1-3_TYROBP_Micro', 'Mic', hodgeSeur$cluster)
hodgeSeur$cluster = gsub('Astro_L1-2_FGFR3_GFAP|Astro_L1-6_FGFR3_SLC14A1', 'Ast', hodgeSeur$cluster)

bakkenSeur$cluster = gsub('Oligo L1-6 OPALIN', 'Oli', bakkenSeur$cluster)
bakkenSeur$cluster = gsub('Micro L1-6 TYROBP', 'Mic', bakkenSeur$cluster)
bakkenSeur$cluster = gsub('OPC L1-6 PDGFRA', 'OPC', bakkenSeur$cluster)
bakkenSeur$cluster = gsub('Astro L1-2 FGFR3 GFAP|Astro L1-6 FGFR3 SLC14A1', 'Ast', bakkenSeur$cluster)

sadSeur$cluster = gsub('Oligo', 'Oli', sadSeur$cluster)
sadSeur$cluster = gsub('Micro', 'Mic', sadSeur$cluster)
sadSeur$cluster = gsub('OPC', 'OPC', sadSeur$cluster)
sadSeur$cluster = gsub('Astro_A|Astro_B', 'Ast', sadSeur$cluster)


# Subset
gliaNoRem = subset(seurmergNOCB_Filt, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))
gliaCB = subset(seurmergCB_Filt, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))
gliaDECONTX = subset(decontSeur, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))
gliaSOUPX = subset(soupSeur, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))

hodgeSeur$cluster = hodgeSeur$cluster
bakkenSeur$cluster = bakkenSeur$cluster_label
sadSeur$cluster = sadSeur$newannot

hodgeSeur$cluster = gsub('Oligo_L1-6_OPALIN_Oligo', 'Oli', hodgeSeur$cluster)
hodgeSeur$cluster = gsub('OPC_L1-6_PDGFRA_OPC', 'OPC', hodgeSeur$cluster)
hodgeSeur$cluster = gsub('Micro_L1-3_TYROBP_Micro', 'Mic', hodgeSeur$cluster)
hodgeSeur$cluster = gsub('Astro_L1-2_FGFR3_GFAP|Astro_L1-6_FGFR3_SLC14A1', 'Ast', hodgeSeur$cluster)

bakkenSeur$cluster = gsub('Oligo L1-6 OPALIN', 'Oli', bakkenSeur$cluster)
bakkenSeur$cluster = gsub('Micro L1-6 TYROBP', 'Mic', bakkenSeur$cluster)
bakkenSeur$cluster = gsub('OPC L1-6 PDGFRA', 'OPC', bakkenSeur$cluster)
bakkenSeur$cluster = gsub('Astro L1-2 FGFR3 GFAP|Astro L1-6 FGFR3 SLC14A1', 'Ast', bakkenSeur$cluster)

sadSeur$cluster = gsub('Oligo', 'Oli', sadSeur$cluster)
sadSeur$cluster = gsub('Micro', 'Mic', sadSeur$cluster)
sadSeur$cluster = gsub('OPC', 'OPC', sadSeur$cluster)
sadSeur$cluster = gsub('Astro_A|Astro_B', 'Ast', sadSeur$cluster)

gliaHODGE = subset(hodgeSeur, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))
gliaBAKKEN = subset(bakkenSeur, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))
gliaSAD = subset(sadSeur, subset = cluster %in% c('Ast', 'Oli', 'OPC', 'Mic'))


# Merge all
gliaNoRem$type = 'No_Removal'
gliaCB$type = 'CB'
gliaDECONTX$type = 'DecontX'
gliaSOUPX$type = 'SoupX'
gliaHODGE$type = 'NeuN-_SD1'
gliaBAKKEN$type = 'NeuN-_SD2'
gliaSAD$type = 'NeuN-_SD3'

gliaNoRem@meta.data = gliaNoRem[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]
gliaCB@meta.data = gliaCB[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]
gliaDECONTX@meta.data = gliaDECONTX[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]
gliaSOUPX@meta.data = gliaSOUPX[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]
gliaHODGE@meta.data = gliaHODGE[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]
gliaBAKKEN@meta.data = gliaBAKKEN[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]
gliaSAD@meta.data = gliaSAD[[c('nCount_RNA', 'nFeature_RNA', 'cluster', 'type')]]

gliaMERG = Reduce(merge, list(gliaNoRem, gliaCB, gliaDECONTX, gliaSOUPX, gliaHODGE, gliaBAKKEN, gliaSAD))

gliaMERG$type = factor(gliaMERG$type, levels = c('NeuN-_SD1', 'NeuN-_SD2', 'NeuN-_SD3', 'CB', 'DecontX', 'SoupX', 'No_Removal'))

saveRDS(gliaMERG, 'gliaMERGED_RemovalComparison.RDS')
#gliaMERG = readRDS('gliaMERGED_RemovalComparison.RDS')


####
## AMBIENT MARKER ENRICHMENT
####

velmAmb = rio::import('~/project/pr4/VELMESHEV_2019/SEURAT/ambMark_PSEUDOBULK.xlsx')
lakeAmb = rio::import('~/project/pr4/LAKE_2018/SEURAT/ambMark_PSEUDOBULK.xlsx')

sortedL = union( velmAmb$Gene[1:500], lakeAmb$Gene[1:500] )

gliaMERG[["percent_ambient"]] <- PercentageFeatureSet(object = gliaMERG, features = intersect(sortedL, rownames(gliaMERG)))
gliaMERG$type = gsub('CB', 'CellBender', gliaMERG$type)
gliaMERG$type = factor(gliaMERG$type, levels = c('NeuN-_SD1', 'NeuN-_SD2', 'NeuN-_SD3', 'CellBender', 'DecontX', 'SoupX', 'No_Removal'))

meta = gliaMERG[[]]
tmpMeta = meta %>% dplyr::filter(type == 'CellBender' & !(is.na(percent_ambient)) ) %>%
	group_by(cluster) %>% dplyr::summarize(mean = median(percent_ambient)) %>% 
	as.data.frame
meta = merge(tmpMeta, meta, by = 'cluster')


pdf('NSD1_Ambient_Removal_Comparison_Glia_Faceted.pdf', width = 12, height = 5)
ggboxplot(meta, x = 'type', y = 'percent_ambient', color = 'black', outlier.shape = NA, ylim = c(0,30)) +
facet_wrap(~cluster, nrow = 1) +
theme_classic() +
theme(text = element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
geom_hline(aes(yintercept = mean), color = 'red', linetype = 'dashed') +
NoLegend() +
xlab('') + ylab('Ambient RNA Percentage')
dev.off()

# Plot some ambient RNAs across datasets and ambient removal analyses
DefaultAssay(gliaMERG) = 'RNA'
gliaMERG = NormalizeData(gliaMERG)
subOli = subset(gliaMERG, subset = cluster == 'Oli')

pdf('NSD1_AMBIENT_COMPARISON_VlnPlot.pdf', width = 8)
VlnPlot(subOli, features = c('MOG', 'MOBP', 'MBP', 'SYT1', 'CSMD1', 'KCNIP4'), pt.size = 0, group.by = 'type')
dev.off()


####
## EXTRACT FOR FINAL CORRELATION
####

sd1_gliaMERG = readRDS('SD1/gliaMERGED_RemovalComparison.RDS')
sd2_gliaMERG = readRDS('SD2/gliaMERGED_RemovalComparison.RDS')
nsd1_gliaMERG = readRDS('NSD1/gliaMERGED_RemovalComparison.RDS')
nsd2_gliaMERG = readRDS('NSD2/gliaMERGED_RemovalComparison.RDS')


sd1_gliaMERG[["percent_ambient"]] <- PercentageFeatureSet(object = sd1_gliaMERG, features = intersect(sortedL, rownames(sd1_gliaMERG)))
sd2_gliaMERG[["percent_ambient"]] <- PercentageFeatureSet(object = sd2_gliaMERG, features = intersect(sortedL, rownames(sd2_gliaMERG)))
nsd1_gliaMERG[["percent_ambient"]] <- PercentageFeatureSet(object = nsd1_gliaMERG, features = intersect(sortedL, rownames(nsd1_gliaMERG)))
nsd2_gliaMERG[["percent_ambient"]] <- PercentageFeatureSet(object = nsd2_gliaMERG, features = intersect(sortedL, rownames(nsd2_gliaMERG)))


sd1_gliaMERG$type = gsub('CB', 'CellBender', sd1_gliaMERG$type)
sd2_gliaMERG$type = gsub('CB', 'CellBender', sd2_gliaMERG$type)
nsd1_gliaMERG$type = gsub('CB', 'CellBender', nsd1_gliaMERG$type)
nsd2_gliaMERG$type = gsub('CB', 'CellBender', nsd2_gliaMERG$type)

sd1_meta = sd1_gliaMERG[[]]
sd2_meta = sd2_gliaMERG[[]]
nsd1_meta = nsd1_gliaMERG[[]]
nsd2_meta = nsd2_gliaMERG[[]]

# Combine the medians of ambient RNA percentage per dataset and type of analysis
ctypes = names(table(sd1_gliaMERG$cluster))
dfL = list()
for(i in 1:length(ctypes)){

	neun1 = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'NeuN-_SD1', 'percent_ambient'] %>% na.omit %>% median
	neun2 = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'NeuN-_SD2', 'percent_ambient'] %>% na.omit %>% median
	neun3 = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'NeuN-_SD3', 'percent_ambient'] %>% na.omit %>% median

	sd1_nr = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'No_Removal', 'percent_ambient'] %>% na.omit %>% median
	sd1_cb = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'CellBender', 'percent_ambient'] %>% na.omit %>% median
	sd1_dc = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'DecontX', 'percent_ambient'] %>% na.omit %>% median
	sd1_sp = sd1_meta[sd1_meta$cluster == ctypes[i] & sd1_meta$type == 'SoupX', 'percent_ambient'] %>% na.omit %>% median

	sd2_nr = sd2_meta[sd2_meta$cluster == ctypes[i] & sd2_meta$type == 'No_Removal', 'percent_ambient'] %>% na.omit %>% median
	sd2_cb = sd2_meta[sd2_meta$cluster == ctypes[i] & sd2_meta$type == 'CellBender', 'percent_ambient'] %>% na.omit %>% median
	sd2_dc = sd2_meta[sd2_meta$cluster == ctypes[i] & sd2_meta$type == 'DecontX', 'percent_ambient'] %>% na.omit %>% median
	sd2_sp = sd2_meta[sd2_meta$cluster == ctypes[i] & sd2_meta$type == 'SoupX', 'percent_ambient'] %>% na.omit %>% median

	nsd1_nr = nsd1_meta[nsd1_meta$cluster == ctypes[i] & nsd1_meta$type == 'No_Removal', 'percent_ambient'] %>% na.omit %>% median
	nsd1_cb = nsd1_meta[nsd1_meta$cluster == ctypes[i] & nsd1_meta$type == 'CellBender', 'percent_ambient'] %>% na.omit %>% median
	nsd1_dc = nsd1_meta[nsd1_meta$cluster == ctypes[i] & nsd1_meta$type == 'DecontX', 'percent_ambient'] %>% na.omit %>% median
	nsd1_sp = nsd1_meta[nsd1_meta$cluster == ctypes[i] & nsd1_meta$type == 'SoupX', 'percent_ambient'] %>% na.omit %>% median

	nsd2_nr = nsd2_meta[nsd2_meta$cluster == ctypes[i] & nsd2_meta$type == 'No_Removal', 'percent_ambient'] %>% na.omit %>% median
	nsd2_cb = nsd2_meta[nsd2_meta$cluster == ctypes[i] & nsd2_meta$type == 'CellBender', 'percent_ambient'] %>% na.omit %>% median
	nsd2_dc = nsd2_meta[nsd2_meta$cluster == ctypes[i] & nsd2_meta$type == 'DecontX', 'percent_ambient'] %>% na.omit %>% median
	nsd2_sp = nsd2_meta[nsd2_meta$cluster == ctypes[i] & nsd2_meta$type == 'SoupX', 'percent_ambient'] %>% na.omit %>% median

	vals = c(sd1_nr, sd1_cb, sd1_dc, sd1_sp, sd2_nr, sd2_cb, sd2_dc, sd2_sp, nsd1_nr, nsd1_cb, nsd1_dc, nsd1_sp, nsd2_nr, nsd2_cb, nsd2_dc, nsd2_sp, neun1, neun2, neun3)
	dataset = c(rep('SD1', 4), rep('SD2', 4), rep('NSD1', 4), rep('NSD2', 4), 'NeuN-_SD1', 'NeuN-_SD2', 'NeuN-_SD3')
	type = c(rep(c('Not_Removed', 'CellBender', 'DecontX', 'SoupX'), 4), 'NeuN-_SD1', 'NeuN-_SD2', 'NeuN-_SD3')
	dfL[[i]] = data.frame(vals = vals, ctypes = ctypes[i], dataset = dataset, type = type)

}

dfAll = do.call(rbind, dfL)
dfAll$type = gsub('NeuN-.*', 'NeuN-_SD', dfAll$type)
dfAll$type = gsub('Not_Removed', 'No_Removal', dfAll$type)

dfAll$type = factor(dfAll$type, levels = c('NeuN-_SD', 'CellBender', 'DecontX', 'SoupX', 'No_Removal'))
comps = list(c('NeuN-_SD', 'CellBender'), c('NeuN-_SD', 'DecontX'), c('NeuN-_SD', 'SoupX'), c('NeuN-_SD', 'No_Removal'))

pdf('Ambient_Tool_Comparison_Summary.pdf')
ggviolin(dfAll, x = 'type', y = 'vals', color = 'black', outlier.shape = NA, add = 'dotplot', ylim = c(0,20)) +
theme_classic() +
theme(text = element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
NoLegend() +
xlab('') + ylab('Ambient RNA Percentage') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', method.args = list(alternative = 'less'), label.y = c(14:18), size = 5)
dev.off()


