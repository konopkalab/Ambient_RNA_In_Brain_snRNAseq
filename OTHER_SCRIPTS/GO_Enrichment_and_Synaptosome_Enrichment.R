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
source("utility_functions.R")

####
## LOAD AMBIENT MARKERS
####

nsd1Amb = rio::import('ExtraNuclear_ambMark.xlsx') # From Table S2
sd1Amb = rio::import('Nuclear_ambMark.xlsx') # From Table S2

####
## GO ENRICHMENTS
####

# Load SD1 dataset for background in gene ontology enrichment
seurmerg_corr = readRDS('SD1_Seurat.RDS')

# Run Enrichment
ambExtraN = nsd1Amb$Gene
ambN = sd1Amb$Gene
uni = rownames(seurmerg_corr)

print(length(ambExtraN))
print(length(ambN))
print(length(uni))

resExtraN = GOenrich(ambExtraN, uni)
resN = GOenrich(ambN, uni)

resExtraN$type = 'ExtraNuclear_Ambient'
resN$type = 'Nuclear_Ambient'

rnadf = rbind(resExtraN, resN)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)

rio::export(rnadf, 'GO_Enrichment.xlsx')

# PLOT
library(showtext)

rnadf = rio::import('GO_Enrichment.xlsx')

# Plot 5 non-redundant categories
rnadf = rnadf[rnadf$GO == 'CC',]
pltgns = rnadf[rnadf$ID %in% c('GO:0044391', 'GO:0098798', 'GO:0098793', 'GO:0045211', 'GO:0034702', 'GO:1902495'),]

pltgns$log10FDR = -log10(pltgns$p.adjust)
pltgns$Description = factor(pltgns$Description, levels = unique(pltgns[order(pltgns$log10FDR), 'Description']))

pdf('GO_ENRICH_AMBIENTS.pdf', width = 25)
print( ggscatter(pltgns, x = 'log10FDR', y = 'Description', size = 12,
	color = 'red', xlim = c(1,max(pltgns$log10FDR))) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
theme_bw() +
facet_wrap(~type, nrow = 1, scales = "free") +
theme(text=element_text(size=30, face = 'bold')) + ylab('') +
xlab(expression('-log'[10]*'(FDR)')))
dev.off()

####
## HAFNER ET AL. OVERLAP
####

nsd1Amb = rio::import('ExtraNuclear_ambMark.xlsx') # From Table S2
sd1Amb = rio::import('Nuclear_ambMark.xlsx') # From Table S2

ambExtraN = nsd1Amb[order(nsd1Amb$logCPM, decreasing = T)[1:500], 'Gene']
ambN = sd1Amb[order(sd1Amb$logCPM, decreasing = T)[1:500], 'Gene']

hafnEnr = rio::import('HAFNER_ENR_HUMAN_SYMBOLS.xlsx')
hafnDep = rio::import('HAFNER_DEP_HUMAN_SYMBOLS.xlsx')

# RUN ENRICHMENT
library(GeneOverlap)
hafnerL = list(hafnEnr[,1], hafnDep[,1])
names(hafnerL) = c('vGLUT_Enriched', 'vGLUT_Depleted')

ambL = list(ambExtraN, ambN)
names(ambL) = c('ExtraNuclear_Ambient', 'Nuclear_Ambient')

resgom = newGOM(hafnerL, ambL, genome.size = nrow(seurmerg_corr))
oddsr = getMatrix(resgom, name="odds.ratio")
pvals = getMatrix(resgom, name="pval")
jac = getMatrix(resgom, name="Jaccard")

oddsrM = melt(oddsr)
pvalsM = melt(pvals)
pvalsM$FDR = p.adjust(as.numeric(pvalsM$value))

toplot = cbind(oddsrM, FDR = pvalsM$FDR)
toplot$value = ifelse(toplot$value > 1000, 1000, toplot$value)
toplot$FDR = formatC(toplot$FDR, format = "e", digits = 2)
toplot$log10FDR = -log10(p.adjust(as.numeric(pvalsM$value)))

pdf("HAFNER_ENRICHMENT.pdf", width = 10)
ggplot(toplot, aes(Var2, Var1, fill = log10FDR))+
 geom_tile(color = "white") +
 rotate_x_text(90) +
 theme_classic() +
 scale_fill_gradient2(midpoint = 3, high = 'red', low = 'white') +
 geom_text(aes(label = FDR), fontface = 'bold', size = 6) +
 ylab('') + xlab('') +
 theme(text=element_text(size=20, face = 'bold')) +
 labs(fill = expression('-log'[10]*'(FDR)'))
dev.off()

####
## TOP EXPRESSED NEURONAL GENES OVERLAP
####

nsd1Seur = readRDS('NSD1_Seurat.RDS')
sd1Seur = readRDS('SD1_Seurat.RDS')

# Keep only neurons
nsd1Seur = subset(nsd1Seur, subset = cluster %in% c('IN-PV', 'IN-SST', 'IN-SV2C', 'IN-VIP', 'L2/3', 'L4', 'L5/6', 'L5/6-CC'))
sd1Seur = subset(sd1Seur, subset = cluster %in% c('Ex1', 'Ex2', 'Ex3a', 'Ex3b', 'Ex3d', 'Ex3e', 'Ex4', 'Ex5b', 'Ex5a', 'Ex6a', 'Ex6b', 'Ex8', 'In1a', 'In1b', 'In1c', 'In2', 'In3', 'In4a', 'In4b', 'In6a', 'In6b', 'In7', 'In8'))

# Find top expressed genes in neurons
Idents(nsd1Seur) = 'NSD1_Neurons'
Idents(sd1Seur) = 'SD1_Neurons'
DefaultAssay(nsd1Seur) = 'RNA'
DefaultAssay(sd1Seur) = 'RNA'
nsd1Seur = NormalizeData(nsd1Seur)
sd1Seur = NormalizeData(sd1Seur)

velmTop = rowMeans(nsd1Seur@assays$RNA) %>% sort %>% tail(500)
lakeTop = rowMeans(sd1Seur@assays$RNA) %>% sort %>% tail(500)

topExp = intersect(names(velmTop), names(lakeTop))

# Plot VENN
library(VennDiagram)
library(ggvenn)
set1 <- ambN
set2 <- ambExtraN
set3 <- topExp

setlist = list(set1, set2, set3)
cols = c('blue', 'orange', 'darkgreen')
nms = c('Nuclear_Ambient', 'ExtraNuclear_Ambient', 'Top_Expressed')
names(setlist) = nms

pdf('Ambient_TopExp_Overlap.pdf')
ggvenn(setlist, fill_color = cols, stroke_size = 0, set_name_size = 5,text_size = 13, digits = 0, stroke_alpha = 0.5, show_percentage = F)
dev.off()











