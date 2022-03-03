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

#####
## READ DATASETS
#####

# Load seurat object that contains frontal cortex and occipital cortex
seurpub = readRDS('sd1_seurat_fcx_occ.RDS')

# Keep only oligo lineage
sd1OL = subset(seurpub, subset = cluster %in% c('OPC', 'Oli') & nFeature_RNA > 0)
sd1OL$id = gsub('.*_', '', colnames(sd1OL))
sd1OL$id = paste0(sd1OL$orig.ident, '_', sd1OL$id)

####
## CLUSTER
####

# Plot without any batch correction
options(future.globals.maxSize = 4000 * 1024^2)
sd1OL = SCTransform(sd1OL, ncells = 1000)
sd1OL <- RunPCA(sd1OL, verbose = FALSE)
sd1OL <- RunUMAP(sd1OL, dims = 1:15)

# Harmony batch correction
library(harmony)
sd1OL_corr <- RunHarmony(object = sd1OL, group.by.vars = 'orig.ident', assay.use="SCT")

# re-compute the UMAP using batch corrected PCAs
sd1OL_corr <- RunUMAP(sd1OL_corr, dims = 1:15, reduction = 'harmony')
sd1OL_corr <- FindNeighbors(sd1OL_corr, dims = 1:15, verbose = FALSE, reduction = 'harmony')
sd1OL_corr <- FindClusters(sd1OL_corr, verbose = FALSE, resolution = 0.5)

pdf("SD1_HARMONY_SAMPLES_OPC_OL.pdf")
DimPlot(sd1OL_corr, group.by = 'orig.ident', pt.size =1, label = T)
dev.off()

pdf("SD1_HARMONY_CLUSTERS_OPC_OL.pdf")
DimPlot(sd1OL_corr, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

pdf("SD1_CLUSTER_DEPTH_OPC_OL.pdf")
ggboxplot(sd1OL_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA', ylim = c(0,1600), outlier.shape = NA) +
rotate_x_text(90)
dev.off()

saveRDS(sd1OL_corr, 'sd1Seur_opc_ol.RDS')

# Find cluster markers
DefaultAssay(sd1OL_corr) = 'RNA'
sd1OL_corr = NormalizeData(sd1OL_corr)
sd1OL_corrM = FindAllMarkers(sd1OL_corr, only.pos = T, logfc.threshold = 0.25)

####
## DESTINY -- OPC_OL
####

library(destiny)
sd1OL_corr = readRDS('sd1Seur_opc_ol.RDS')

pubVisOli = subset(sd1OL_corr, subset = region == 'occ')
pubVisOli$orig.ident = factor(pubVisOli$orig.ident)
pubVisOli = NormalizeData(pubVisOli)
marks = FindMarkers(pubVisOli, ident.1 = 'OPC', group.by = 'cluster')

## Run Destiny
ct <- pubVisOli@assays$RNA@data[rownames(marks),]
dm <- DiffusionMap(t(as.matrix(ct)), n_pcs = 100, k = 100)

saveRDS(dm, 'Destiny_dm_OPC_OL_PUB.RDS')
dm = readRDS('Destiny_dm_OPC_OL_PUB.RDS')

pubdf <- data.frame(DC1 = -eigenvectors(dm)[,1],
	          DC2 = eigenvectors(dm)[,2],
	          color = pubVisOli$cluster[colnames(ct)])

minDC1 = min(pubdf$DC1)
maxDC1 = max(pubdf$DC1)
minDC2 = min(pubdf$DC2)
maxDC2 = max(pubdf$DC2)

# Plot custom
pubVisOli = NormalizeData(pubVisOli)
pubVisOli = ScaleData(pubVisOli)

# Plot genes
expvec = pubVisOli@assays$RNA@scale.data[c('GPR17', 'FYN', 'BCAS1', 'TFEB', 'ENPP6'),]
toplot = cbind(pubdf[,1:2], t(expvec))
toplot = melt(toplot, id.vars = c('DC1', 'DC2'))
toplot[toplot$value > 2, 'value'] = 2
toplot = toplot[order(toplot$value),]

pdf('SD1_PUB_OL_TRAJECTORY_FEATUREPLOT_AMBIENTs.pdf', width = 18, height = 10)
ggscatter(toplot, x = 'DC1', y = 'DC2', color = 'value') +
scale_color_gradient2(midpoint = 0, low = 'blue', high = 'red', mid = 'white') +
theme(text=element_text(size=20, face = 'bold'), legend.position = 'right') +
facet_wrap(~variable) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2))
dev.off()

# Take 200 nuclei from the middle
midpoint = sum(max(pubdf$DC1), min(pubdf$DC1))/2

incdf = pubdf[pubdf$DC1 > midpoint,]
decdf = pubdf[pubdf$DC1 < midpoint,]
toOL = incdf[order(incdf$DC1)[1:100],]
toOPC = decdf[order(decdf$DC1, decreasing = T)[1:100],]

imOli = c(rownames(toOL), rownames(toOPC))

# Plot clusters with borders
toplot = cbind(pubdf[,1:2], pubVisOli$cluster)
toplot = melt(toplot, id.vars = c('DC1', 'DC2'))

pdf('SD1_PUB_OL_TRAJECTORY_FEATUREPLOT_CLUSTERs.pdf')
ggscatter(toplot, x = 'DC1', y = 'DC2', color = 'value', palette = c('red', 'blue')) +
theme(text=element_text(size=20, face = 'bold')) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2)) +
geom_vline(xintercept = c(min(toOPC$DC1),max(toOL$DC1)), linetype = 'dashed', color = 'darkgreen')
dev.off()

# Plot with annotation
Idents(pubVisOli) = pubVisOli$cluster
opcNuc = WhichCells(pubVisOli, idents = 'OPC')
olNuc = WhichCells(pubVisOli, idents = 'Oli')

pubVisOli$dmannot = ifelse(colnames(pubVisOli) %in% imOli, 'Middle',
			ifelse(colnames(pubVisOli) %in% opcNuc, 'OPC', 'Oli'))

pubVisOli$log10UMI = log10(pubVisOli$nCount_RNA)
pubVisOli$dmannot = factor(pubVisOli$dmannot, levels = c('OPC', 'Middle', 'Oli'))
pdf("SD1_TRAJECTORY_UMI_COMPARISON.pdf")
ggviolin(pubVisOli[[]], x = 'dmannot', y = 'log10UMI', color = 'white', fill = 'grey', outlier.shape = NA) +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
ylab(expression('log'[10]*'(UMI)')) + xlab('') + NoLegend()
dev.off()

# Find markers of three zones
DefaultAssay(pubVisOli) = 'RNA'
pubVisOli = NormalizeData(pubVisOli)

Idents(pubVisOli) = pubVisOli$dmannot
opcOlMarks = FindAllMarkers(pubVisOli, only.pos = T, group.by = 'dmannot')
opcOlMarks = split(opcOlMarks, opcOlMarks$cluster) %>% lapply(., function(x){x$gene})
names(opcOlMarks) = paste0('OPC-OL_', names(opcOlMarks))

# Ratio of COPs in transition zone
sd1COP = readRDS('sd1OL_CellBender_OPC_CLEANED_CLUSTERED.RDS')
sd1COP = subset(sd1COP, subset = newannot == 'COP')

ratTransition = sum(pubVisOli[[]][pubVisOli$dmannot == 'Middle', 'id'] %in% colnames(sd1COP)) / sum(pubVisOli$dmannot %in% c('Middle', 'OPC'))
ratOPC = sum(pubVisOli[[]][pubVisOli$dmannot == 'OPC', 'id'] %in% colnames(sd1COP)) /  sum(pubVisOli$dmannot %in% c('Middle', 'OPC'))

tots = sum(pubVisOli$dmannot %in% c('Middle', 'OPC'))
vals = c(ratOPC, ratTransition)
vars = c('COP_In_OPC', 'COP_In_Transition')

toplot = data.frame(vals = vals, vars = vars, tots = tots)
toplot$counts = toplot[,1] * toplot[,3]

pdf('COP_DEPLETION_IN_TRANSIENT.pdf', width = 10, height = 3)
ggbarplot(toplot, x = 'vars', y = 'vals', fill = 'vars', color = 'white') +
theme(text=element_text(size=20, face = 'bold')) +
xlab('') + ylab('Ratio of COPs') + NoLegend() +
coord_flip()
dev.off()

# Chi-square test
prop.test(toplot[,4], toplot[,3])$p.value #0.02584862

####
## DESTINY -- OPC-ASTROCYTES
####

# Filter and normalize
sd1AST = subset(seurpub, subset = cluster %in% c('OPC', 'Ast') & nFeature_RNA > 0)
sd1AST = NormalizeData(sd1AST)
pubVisAst = subset(sd1AST, subset = region == 'occ')
marks = FindMarkers(pubVisAst, ident.1 = 'OPC', group.by = 'cluster')

## Run Destiny
ct <- pubVisAst@assays$RNA@data[rownames(marks),]
dm <- DiffusionMap(t(as.matrix(ct)), n_pcs = 100, k = 100)
saveRDS(dm, 'Destiny_dm_OPC_AST_PUB.RDS')

pubdf <- data.frame(DC1 = eigenvectors(dm)[,1],
	          DC2 = eigenvectors(dm)[,2],
	          color = pubVisAst$cluster[colnames(ct)])

minDC1 = min(pubdf$DC1)
maxDC1 = max(pubdf$DC1)
minDC2 = min(pubdf$DC2)
maxDC2 = max(pubdf$DC2)

# PLOT CUSTOM
DefaultAssay(pubVisAst) = 'RNA'
pubVisAst = NormalizeData(pubVisAst)
pubVisAst = ScaleData(pubVisAst)

# Plot genes
expvec = pubVisAst@assays$RNA@scale.data[c('SYT1', 'SNAP25', 'GRIN2A', 'PCDH15', 'SLC1A3'),]
toplot = cbind(pubdf[,1:2], t(expvec))
toplot = melt(toplot, id.vars = c('DC1', 'DC2'))
toplot[toplot$value > 2, 'value'] = 2
toplot = toplot[sample(1:nrow(toplot), nrow(toplot)),]

pdf('SD1_PUB_AST_TRAJECTORY_FEATUREPLOT_AMBIENTs.pdf', width = 16, height = 10)
ggscatter(toplot, x = 'DC1', y = 'DC2', color = 'value') +
scale_color_gradient2(midpoint = 0, low = 'blue', high = 'red', mid = 'white') +
theme(text=element_text(size=20, face = 'bold')) +
facet_wrap(~variable) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2))
dev.off()


# Take 200 nuclei from the middle
midpoint = sum(max(pubdf$DC1), min(pubdf$DC1))/2

incdf = pubdf[pubdf$DC1 > midpoint,]
decdf = pubdf[pubdf$DC1 < midpoint,]
toOL = incdf[order(incdf$DC1)[1:100],]
toOPC = decdf[order(decdf$DC1, decreasing = T)[1:100],]

imOli = c(rownames(toOL), rownames(toOPC))

# Plot clusters with borders
toplot = cbind(pubdf[,1:2], pubVisAst$cluster)
toplot = melt(toplot, id.vars = c('DC1', 'DC2'))

pdf('SD1_PUB_AST_TRAJECTORY_FEATUREPLOT_CLUSTERs.pdf')
ggscatter(toplot, x = 'DC1', y = 'DC2', color = 'value', palette = c('red', 'blue')) +
theme(text=element_text(size=20, face = 'bold')) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2)) +
geom_vline(xintercept = c(min(toOPC$DC1),max(toOL$DC1)), linetype = 'dashed', color = 'darkgreen')
dev.off()

# Find markers of these nuclei
Idents(pubVisAst) = pubVisAst$cluster
opcNuc = WhichCells(pubVisAst, idents = 'OPC')
olNuc = WhichCells(pubVisAst, idents = 'Ast')

pubVisAst$dmannot = ifelse(colnames(pubVisAst) %in% imOli, 'Middle',
			ifelse(colnames(pubVisAst) %in% opcNuc, 'OPC', 'Ast'))

DefaultAssay(pubVisAst) = 'RNA'
pubVisAst = NormalizeData(pubVisAst)
Idents(pubVisAst) = pubVisAst$dmannot

opcAstMarks = FindAllMarkers(pubVisAst, only.pos = T, group.by = 'dmannot')
opcAstMarks = split(opcAstMarks, opcAstMarks$cluster) %>% lapply(., function(x){x$gene})
names(opcAstMarks) = paste0('OPC-AST_', names(opcAstMarks))




####
## DESTINY -- ASTROCYTES-MICROGLIA
####

# Filter and normalize
sd1AST = subset(seurpub, subset = cluster %in% c('Mic', 'Ast') & nFeature_RNA > 0)
sd1AST = NormalizeData(sd1AST)
pubVisAst = subset(sd1AST, subset = region == 'occ')
marks = FindMarkers(pubVisAst, ident.1 = 'Ast', group.by = 'cluster')

## Run Destiny
ct <- pubVisAst@assays$RNA@data[rownames(marks),]
dm <- DiffusionMap(t(as.matrix(ct)), n_pcs = 100, k = 100)
saveRDS(dm, 'Destiny_dm_MIC_AST_PUB.RDS')


pubdf <- data.frame(DC1 = eigenvectors(dm)[,1],
	          DC2 = eigenvectors(dm)[,2],
	          color = pubVisAst$cluster[colnames(ct)])

minDC1 = min(pubdf$DC1)
maxDC1 = max(pubdf$DC1)
minDC2 = min(pubdf$DC2)
maxDC2 = max(pubdf$DC2)

# PLOT CUSTOM
DefaultAssay(pubVisAst) = 'RNA'
pubVisAst = NormalizeData(pubVisAst)
pubVisAst = ScaleData(pubVisAst)


# Take 200 nuclei from the middle
midpoint = sum(max(pubdf$DC1), min(pubdf$DC1))/2

incdf = pubdf[pubdf$DC1 > midpoint,]
decdf = pubdf[pubdf$DC1 < midpoint,]
toOL = incdf[order(incdf$DC1)[1:100],]
toMIC = decdf[order(decdf$DC1, decreasing = T)[1:100],]

imOli = c(rownames(toOL), rownames(toMIC))

# Plot clusters with borders
toplot = cbind(pubdf[,1:2], pubVisAst$cluster)
toplot = melt(toplot, id.vars = c('DC1', 'DC2'))

pdf('SD1_PUB_ASTMIC_TRAJECTORY_FEATUREPLOT_CLUSTERs.pdf')
ggscatter(toplot, x = 'DC1', y = 'DC2', color = 'value', palette = c('red', 'blue')) +
theme(text=element_text(size=20, face = 'bold')) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2)) +
geom_vline(xintercept = c(min(toMIC$DC1),max(toOL$DC1)), linetype = 'dashed', color = 'darkgreen')
dev.off()

# Find markers of these nuclei
Idents(pubVisAst) = pubVisAst$cluster
micNuc = WhichCells(pubVisAst, idents = 'Mic')
olNuc = WhichCells(pubVisAst, idents = 'Ast')

pubVisAst$dmannot = ifelse(colnames(pubVisAst) %in% imOli, 'Middle',
			ifelse(colnames(pubVisAst) %in% micNuc, 'Mic', 'Ast'))

DefaultAssay(pubVisAst) = 'RNA'
pubVisAst = NormalizeData(pubVisAst)
Idents(pubVisAst) = pubVisAst$dmannot

micAstMarks = FindAllMarkers(pubVisAst, only.pos = T, group.by = 'dmannot')
micAstMarks = split(micAstMarks, micAstMarks$cluster) %>% lapply(., function(x){x$gene})
names(micAstMarks) = paste0('MIC-AST_', names(micAstMarks))


#####
## ENRICHMENT
#####

# Run enrichment with ambient markers
nsd1Amb = rio::import('ExtraNuclear_ambMark.xlsx') # From Table S2
sd1Amb = rio::import('Nuclear_ambMark.xlsx') # From Table S2

ambList = union(sd1Amb$Gene[1:500], nsd1Amb$Gene[1:500])
ambList = list(ambList) 

# Lake et al immature oli markers
sd1Marks = rio::import('Lake_Cluster_Markers.xlsx')
ambList[[2]] = sd1Marks[sd1Marks[,2] == 'Immature_Oli', 1]
names(ambList) = c('Ambient_Markers', 'Immature_Oli')

library(GeneOverlap)
resgom = newGOM(ambList, c(opcOlMarks, opcAstMarks, micAstMarks), genome.size=13000)
oddsr = getMatrix(resgom, name="odds.ratio")
pvals = getMatrix(resgom, name="pval")
jac = getMatrix(resgom, name="Jaccard")

oddsrM = melt(oddsr)
pvalsM = melt(pvals)
pvalsM$FDR = p.adjust(as.numeric(pvalsM$value))

toplot = cbind(oddsrM, FDR = pvalsM$FDR)
toplot$value = ifelse(toplot$value > 1000, 1000, toplot$value)
toplot$is_top = ifelse(toplot$value > 3 & toplot$FDR > 3, 'Top', 'Others')
toplot$FDR = formatC(toplot$FDR, format = "e", digits = 2)

toplot$trajectory = gsub('_.*', '', toplot$Var2)
toplot$log10FDR = -log10(p.adjust(as.numeric(pvalsM$value)))
toplot$Var2 = factor(toplot$Var2, levels = rev(c('OPC-OL_OPC', 'OPC-OL_Middle', 'OPC-OL_Oli', 'OPC-AST_OPC', 'OPC-AST_Middle', 'OPC-AST_Ast', 'MIC-AST_Ast', 'MIC-AST_Middle', 'MIC-AST_Mic')))

pdf("ENRICH_TRAJECTORY_GROUPS_SD1.pdf", width = 14, height = 6)
ggplot(toplot, aes(Var2, Var1, fill = log10FDR))+
 geom_tile(color = "white") +
 rotate_x_text(90) +
 theme_classic() +
 scale_fill_gradient2(midpoint = 1, high = 'red', low = 'white') +
 geom_text(aes(label = FDR), fontface = 'bold') +
 ylab('') + xlab('') +
 theme(text=element_text(size=20, face = 'bold')) +
 labs(fill = expression('-log'[10]*'(FDR)')) +
 facet_wrap(~trajectory, scales = 'free') +
 coord_flip() + rotate_x_text(45)
dev.off()

####
## CELLBENDER CLEAN ON TRAJECTORY
####

## OPC-OLI ##
sd1CB_cleanOPC = readRDS('sd1OL_CB_OPC_CLEANED.RDS')
sd1CB_cleanOLI = readRDS('sd1OL_CB_OLI_CLEANED.RDS')

sd1CB_clean = merge(sd1CB_cleanOPC, sd1CB_cleanOLI)
sd1CB_clean$region = gsub('[0-9]', '', sd1CB_clean$orig.ident)

# Plot Diffusion Map highlighting removed cell barcodes
dm_opc_ol = readRDS('Destiny_dm_OPC_OL_PUB.RDS')
eigenDF = eigenvectors(dm_opc_ol) %>% as.data.frame
rownames(eigenDF) = gsub('Oli_', '', rownames(eigenDF)) %>% gsub('OPC_', '', .)
eigenDF$Filtered = 'Not_Filtered'
eigenDF[!(rownames(eigenDF) %in% colnames(sd1CB_clean)), 'Filtered'] = 'Filtered'

pubdf <- data.frame(DC1 = -eigenDF[,1],
	          DC2 = eigenDF[,2],
	          color = eigenDF$Filtered)

minDC1 = min(pubdf$DC1)
maxDC1 = max(pubdf$DC1)
minDC2 = min(pubdf$DC2)
maxDC2 = max(pubdf$DC2)

pdf('SD1_TRAJECTORY_CELLBENDER_REMOVED_OL_OPC.pdf')
ggscatter(pubdf, x = 'DC1', y = 'DC2', color = 'color', palette = c('lightblue', 'black', 'red')) +
theme(text=element_text(size=20, face = 'bold')) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2))
dev.off()

## OPC-AST ##
sd1CB_cleanOPC = readRDS('sd1OL_CB_OPC_CLEANED.RDS')
sd1CB_cleanAST = readRDS('sd1OL_CB_AST_CLEANED.RDS')

sd1CB_clean = sd1CB_cleanOPC
sd1CB_clean = merge(sd1CB_cleanOPC, sd1CB_cleanAST)
sd1CB_clean$region = gsub('[0-9]', '', sd1CB_clean$orig.ident)

# Plot Diffusion Map highlighting removed cell barcodes
dm_opc_ast = readRDS('Destiny_dm_OPC_AST_PUB.RDS')
eigenDF = eigenvectors(dm_opc_ast) %>% as.data.frame
rownames(eigenDF) = gsub('Ast_', '', rownames(eigenDF)) %>% gsub('OPC_', '', .)
eigenDF$Filtered = ifelse(rownames(eigenDF) %in% colnames(sd1CB_clean), 'Not_Filtered', 'Filtered')

pubdf <- data.frame(DC1 = eigenDF[,1],
	          DC2 = eigenDF[,2],
	          color = eigenDF$Filtered)

minDC1 = min(pubdf$DC1)
maxDC1 = max(pubdf$DC1)
minDC2 = min(pubdf$DC2)
maxDC2 = max(pubdf$DC2)

pdf('SD1_TRAJECTORY_CELLBENDER_REMOVED_AST_OPC.pdf')
ggscatter(pubdf, x = 'DC1', y = 'DC2', color = 'color', palette = c('black', 'lightblue')) +
theme(text=element_text(size=20, face = 'bold')) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2))
dev.off()


## MIC-AST ##
sd1CB_cleanMIC = readRDS('sd1OL_CB_MIC_CLEANED.RDS')
sd1CB_cleanAST = readRDS('sd1OL_CB_AST_CLEANED.RDS')

sd1CB_clean = merge(sd1CB_cleanMIC, sd1CB_cleanAST)
sd1CB_clean$region = gsub('[0-9]', '', sd1CB_clean$orig.ident)

# Plot Diffusion Map highlighting removed cell barcodes
dm_mic_ast = readRDS('Destiny_dm_MIC_AST_PUB.RDS')
eigenDF = eigenvectors(dm_mic_ast) %>% as.data.frame
rownames(eigenDF) = gsub('Ast_', '', rownames(eigenDF)) %>% gsub('Mic_', '', .)
eigenDF$Filtered = ifelse(rownames(eigenDF) %in% colnames(sd1CB_clean), 'Not_Filtered', 'Filtered')

pubdf <- data.frame(DC1 = eigenDF[,1],
	          DC2 = eigenDF[,2],
	          color = eigenDF$Filtered)

minDC1 = min(pubdf$DC1)
maxDC1 = max(pubdf$DC1)
minDC2 = min(pubdf$DC2)
maxDC2 = max(pubdf$DC2)

pdf('SD1_TRAJECTORY_CELLBENDER_REMOVED_AST_MIC.pdf')
ggscatter(pubdf, x = 'DC1', y = 'DC2', color = 'color', palette = c('black', 'lightblue')) +
theme(text=element_text(size=20, face = 'bold')) +
xlim (c(minDC1, maxDC1)) +
ylim(c(minDC2, maxDC2))
dev.off()




