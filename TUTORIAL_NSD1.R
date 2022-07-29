rm(list = ls())
library(dplyr)
library(Seurat)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(DropletQC)
source('FUNCTIONS.R')

####
## STEP 00: DOWNLOAD DATASETS
####

# Raw count matrix (NSD1): https://cloud.biohpc.swmed.edu/index.php/s/WaNwndJmoG7TYwJ
# Glial cells after running CellBender (NSD1): https://cloud.biohpc.swmed.edu/index.php/s/krMDmEowmKcp6Hj

####
## STEP 01: LOAD RAW COUNT MATRIX AND CALCULATE INTRONIC READ RATIO
####

# We will use dropletQC to calculate ratio of reads mapping to introns.
# dropletQC requires the cell barcode to be in a separate tag in BAM file.

# If the preprocessing was done in a way that contains cell barcode...
# within the sequence id of the bam file,
# please refer to the script 'CB_HEADER_TO_TAG_IN_BAM.sh'

# Use dropletQC to calculate nuclear fraction (i.e intronic read ratio)

# Load the raw seurat object (without cell calling)
sampObj = readRDS('seurat_raw.RDS')

# Subset a sample
#sampObj = subset(sampObj, subset = orig.ident %in% c('SRR9262946'))

# Set the paths for the nuclear fraction function
gtfPath = '~/workdir/reference_genomes/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/genes/genes.gtf'
bamPath = '~/workShared/ForEmre/'

# Loop over samples and calculate nuclear fraction per cell barcode
samps = unique(sampObj$orig.ident)
meta = sampObj[[]]
nucFrac = list()
for(i in 1:length(samps)){

	# Extract cell barcodes
	barcs = gsub('_.*', '', rownames(meta[meta$orig.ident == samps[i],]) )

	# Loop over samples and calculate nuclear fraction per cell barcode
	nucFrac[[i]] = nuclear_fraction_annotation(
		annotation_path = gtfPath,
		bam = paste0(bamPath, samps[i], '_With_CB.bam'),
		barcodes = barcs,
		cores = 23, verbose = T)
}

nucFrac2 = lapply(1:length(nucFrac),
		function(x){rownames(nucFrac[[x]]) = paste0(rownames(nucFrac[[x]]), '_', samps[x]); nucFrac[[x]]}) %>%
		do.call(rbind, .)

# Add intronic read ratio to the metadata of the seurat object
sampObj$intronRat = nucFrac2[colnames(sampObj), 1]

# Plot intronic read ratio as a function of UMI
# Notice the cell barcodes with decent UMI number but low intronic read ratio
meta = sampObj@meta.data
meta$log10UMI = log10(meta$nCount_RNA)

pdf('Log10UMI_by_intronicReadRatio.pdf')
ggscatter(meta, x = 'log10UMI', y = 'intronRat') +
ylab('Intronic Read Ratio')
dev.off()

####
## STEPS 02 and 03: FIND AMBIENT CLUSTERS
####

# We will now find cell barcode clusters that are largely called as empty-droplets.
# These will be called ambient clusters.
# This will help us evaluate percentage of any group of annotated cells in the ambient clusters
# as well as find ambient RNA markers / most abundant genes in ambient clusters.

# Keep cells with at least 1 UMI. This shouldn't filter anything for most datasets.
sampObj = subset(sampObj, subset = nCount_RNA > 0)

# Find non-empty cell barcodes and ambient clusters
sampObjAmb = ambClusterFind(sampObj, batchCorrect = T)

# Run UMAP to visualize ambient clusters
sampObjAmb = RunUMAP(sampObjAmb, dims = 1:30)

# Visualize ambient clusters
pdf('Excess_AllClusters_UMAP.pdf')
DimPlot(sampObjAmb, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
	theme(text=element_text(size=20, face = 'bold')) + NoLegend() + ggtitle('')
dev.off()

pdf('Excess_AmbientClusters_UMAP.pdf', width = 12)
DimPlot(sampObjAmb, group.by = 'is_ambient_cluster', label = T, raster = T, label.size = 7) +
	theme(text=element_text(size=20, face = 'bold')) + NoLegend() + ggtitle('') +
	facet_wrap(~is_ambient_cluster)
dev.off()


####
## STEP 04: FIND AMBIENT CLUSTER MARKERS / MOST ABUNDANT GENES
####

# Now that we have ambient clusters, we can now find what genes mark / are most abundant
# in these clusters. We will set sorted = F since the data is not nuclei sorted. This
# will split ambient clusters into two by their intronic read ratio.

# We will find the marker genes for both low-intronic ambient clusters and high-intronic ambient clusters.

# This function returns a list with two data frames that has ambient RNA markers
# One data frame for low, one for high intronic cell barcodes.
# We are setting the logfc low since ambient clusters have many dropouts, reducing the fold change.

abundanceL = ambMarkFind(sampObjAmb, sorted = F, type = 'Markers', logfc = 0.1)
highMarks = abundanceL[["High_Intronic_Ambient_Markers"]]
lowMarks = abundanceL[["Low_Intronic_Ambient_Markers"]]

# Keep the significant markers
highMarks = highMarks[highMarks$p_val_adj < 0.05,] %>% rownames
lowMarks = lowMarks[lowMarks$p_val_adj < 0.05,] %>% rownames

ambGenes = union(highMarks, lowMarks)

# Save ambient genes
saveRDS(ambGenes, 'ambGenes.RDS')

# Find percentage of the ambient RNA marker genes per cell barcode
sampObjAmb[["percent_ambient_highIntronic"]] = PercentageFeatureSet(object = sampObjAmb, features = highMarks)
pdf('High_Intronic_Ambient_Percentage_Across_Annotations.pdf')
ggboxplot(sampObjAmb[[]], x = 'seurat_clusters', y = 'percent_ambient_highIntronic', outlier.shape = NA, ylim = c(0,30)) +
xlab('') + ylab('Nuclear Ambient Marker Percentage') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90)
dev.off()

sampObjAmb[["percent_ambient_LowIntronic"]] = PercentageFeatureSet(object = sampObjAmb, features = lowMarks)
pdf('Low_Intronic_Ambient_Percentage_Across_Annotations.pdf')
ggboxplot(sampObjAmb[[]], x = 'seurat_clusters', y = 'percent_ambient_LowIntronic', outlier.shape = NA, ylim = c(0,30)) +
xlab('') + ylab('Non-nuclear Ambient Marker Percentage') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90)
dev.off()

####
## STEP 05: RUN CELLBENDER
####

# Run CellBender with the default parameters. Please refer to: https://github.com/broadinstitute/CellBender

####
## STEP 06: ANNOTATE MAJOR CELL TYPES AFTER CELLBENDER
####

# Since this dataset is already annotated, we are only going to...
# create seurat object after CellBender.

# For general clustering and annotation of cell types, please refer to: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

# Combine matrices after running CellBender separately per sample 
sp = unique(sampObj$orig.ident)
hmatL = list()
for(i in 1:length(sp)){
	hmat = Read10X_h5(paste0('CELLBENDER/', sp[i], '/CellBender_out_filtered.h5'))
	colnames(hmat) = gsub('-1', paste0('_', sp[i]), colnames(hmat))
	hmatL[[i]] = hmat

}

cmnGns = Reduce(intersect, lapply(hmatL, function(x){rownames(x)}))
hmatL2 = lapply(hmatL, function(x){x[cmnGns,]})
hmatFinal = Reduce(cbind, hmatL2)

# Create and save the seurat object
seurmerg = CreateSeuratObject(counts = hmatFinal)
seurmerg$orig.ident = gsub('.*_', '', colnames(seurmerg))

saveRDS(seurmerg, 'CELLBENDER/seurat_combined_cellbender.RDS')


####
## STEP 07: FIND SUBCLUSTERS ENRICHED IN AMBIENT RNAs
####

# I recommend running this part after ambient RNA removal with CellBender.
# Primary use of subCLEAN is to find small clusters of cells that still carry...
#...high amounts of ambient RNA contamination after ambient RNA removal.

# The subCLEAN function will give subclusters and QC plots.
# The user can then determine which subclusters to remove.

# Read the seurat object that contains glial cell types.
afterCB = readRDS('NSD1_CellBender_Glia.RDS')

# Run the subcluster clean function
afterCB_AmbientMarked = subCLEAN(object = afterCB, group.by = 'cluster',
		key = 'Oligodendrocytes', ambientMarkers = ambGenes, batchCorrect = T)

# The function returns two values; seurat object with subclustering info...
# and a table of subcluster enrichment with ambient RNA markers


