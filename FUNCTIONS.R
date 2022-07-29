# Ambient Marker Finder
ambClusterFind = function(rawObj, finalCBs = NULL, ambCutoffCoef = 1.5, batchCorrect = T, timesExcess = 2, res = 0.5, npcs = 30, ncores = 8){

	if(batchCorrect == T & length(unique(rawObj$orig.ident)) == 1){
		stop('Only one batch detected in orig.ident. Please set batchCorrect to False or add the batch information to orig.ident')
	}

	# If there are cells without any UMI, remove them
	rawObj = subset(rawObj, subset = nCount_RNA > 0)

	## Keep excess number of cells

	# Find non-empty droplets if final list of cell barcodes is not specified

	if(is.null(finalCBs)){

		suppressMessages(require(DropletUtils))
		suppressMessages(require(parallel))

		print('Finding non-empty droplets with DropletUtils')

		samples = unique(rawObj$orig.ident)
		cellsL = mclapply(1:length(samples), mc.cores = ncores, function(x){
			print(samples[x])
			rawSub = subset(rawObj, subset = orig.ident == samples[x])	
			out = emptyDrops(rawSub@assays$RNA@counts)
			out = out[!(is.na(out$FDR)), ]
			out[out$FDR < 0.05,] %>% rownames
		})
		finalCBs = unlist(cellsL)
	}


	rawObj$is_empty = ifelse(colnames(rawObj) %in% finalCBs, 'Non_empty', 'Empty')
	
	if(prop.table(table(rawObj$is_empty))['Empty'] < 0.5){
		print('More than half of the cell barcodes were not filtered')
		print('We will keep all cell barcodes to find the ambient clusters')

		softFiltObj = rawObj
	} else{

		# Keep n times more than filtered
		keepL = list()
		rawMeta = rawObj@meta.data
		for(i in 1:length(samples)){
			rawMetaSub = rawMeta[rawMeta$orig.ident == samples[i],]
			rawMetaSubFinalCount = sum(rawMetaSub$is_empty == 'Non_empty')
			keepMax = max(rawMetaSubFinalCount * timesExcess, nrow(rawMetaSub))
			keepL[[i]] = rawMetaSub[order(rawMetaSub$nCount_RNA, decreasing = T)[1:keepMax], ] %>% rownames
		}

		# Soft filter
		softFiltObj = subset(rawObj, cells = unlist(keepL))
	}


	print('Running Clustering With Excess Barcodes')
	softFiltObj = NormalizeData(softFiltObj)
	softFiltObj = FindVariableFeatures(softFiltObj, selection.method = "vst", nfeatures = 2000)
	softFiltObj = ScaleData(softFiltObj, features = rownames(softFiltObj))
	softFiltObj = RunPCA(softFiltObj, verbose = F)

	# Batch correction
	if(batchCorrect == T){

	print('Running batch correction')
	suppressMessages(require(harmony))

	softFiltObj = RunHarmony(object = softFiltObj, group.by.vars = 'orig.ident', assay.use="RNA")
	softFiltObj = FindNeighbors(softFiltObj, dims = 1:npcs, verbose = FALSE, reduction = 'harmony')
	softFiltObj = FindClusters(softFiltObj, verbose = FALSE, resolution = res)
	} else{
	softFiltObj = FindNeighbors(softFiltObj, dims = 1:npcs, verbose = FALSE, reduction = 'pca')
	softFiltObj = FindClusters(softFiltObj, verbose = FALSE, resolution = res)
	}

	## Find Ambient Clusters
	filtObjMeta = softFiltObj[[]]
	clDF = filtObjMeta$seurat_clusters %>% table %>% names %>% data.frame(clusters = .)
	clDF$size = filtObjMeta[,'seurat_clusters'] %>% table
	clDF$emptySize = filtObjMeta[filtObjMeta$is_empty == 'Empty', 'seurat_clusters'] %>% table
	clDF$ratio = clDF$emptySize / clDF$size

	# Calculate cutoff to call a cluster among ambient clusters
	# Note that the largest cluster will automatically be considered the ambient cluster
	cutoff = prop.table(table(rawObj$is_empty))['Empty'] * ambCutoffCoef
	clDF_Ambient = clDF[clDF$ratio > cutoff | clDF$size == max(clDF$size),]

	filtObjMeta$is_ambient_cluster = 'Non_Ambient'
	filtObjMeta[filtObjMeta$seurat_clusters %in% clDF_Ambient$clusters, 'is_ambient_cluster'] = 'Ambient'
	softFiltObj@meta.data = filtObjMeta

	return(softFiltObj)
}

ambMarkFind = function(seuratObject, intronicCutoff = c(0.5, 0.7), sorted = F, type, logfc = 0.25){

	## Find Ambient Cluster Markers
	if(sorted == F){

		ambientObj = subset(seuratObject, subset = is_ambient_cluster == 'Ambient')
		ambientMeta = ambientObj@meta.data

		# Find high intronic and low intronic ambient cell barcodes
		ambientMeta$ambient_intronic_type = 'Unassigned'
		ambientMeta[ambientMeta$intronRat > intronicCutoff[2], 'ambient_intronic_type'] = 'High_Intronic'
		ambientMeta[ambientMeta$intronRat < intronicCutoff[1], 'ambient_intronic_type'] = 'Low_Intronic'

		# Plot histogram
		pdf('Ambient_Clusters_Histogram.pdf', width = 12)
		print( gghistogram(ambientMeta, x = "intronRat", add = "mean", rug = TRUE, fill = "ambient_intronic_type", palette = c("#00AFBB", "#E7B800", "gray")) +
		 theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
		 ylab('Count') + xlab('Intronic_Read_Ratio') )
		dev.off()

		if(type == 'Abundance'){

			# Find most abundant genes of top and bottom ambient cell barcodes
			highCB = rownames(ambientMeta[ambientMeta$ambient_intronic_type == 'High_Intronic',])
			highMat = ambientObj@assays$RNA@counts[, highCB] %>% rowSums %>% prop.table * 100

			lowCB = rownames(ambientMeta[ambientMeta$ambient_intronic_type == 'Low_Intronic',])
			lowMat = ambientObj@assays$RNA@data[, lowCB] %>% rowSums %>% prop.table * 100

			toReturn = list(highMat, lowMat)
			names(toReturn) = c('High_Intronic_Ambient_GeneAbundance', 'Low_Intronic_Ambient_GeneAbundance')
			return(toReturn)
		}

		if(type == 'Markers'){

			# Assign meta data
			ambientObj@meta.data = ambientMeta

			# Find marker genes
			highMarks = FindMarkers(ambientObj, ident.1 = 'High_Intronic', ident.2 = 'Unassigned', group.by = 'ambient_intronic_type', logfc.threshold = logfc, only.pos = T)
			lowMarks = FindMarkers(ambientObj, ident.1 = 'Low_Intronic', ident.2 = 'Unassigned', group.by = 'ambient_intronic_type', logfc.threshold = logfc, only.pos = T)

			toReturn = list(highMarks, lowMarks)
			names(toReturn) = c('High_Intronic_Ambient_Markers', 'Low_Intronic_Ambient_Markers')
			return(toReturn)
		}


	} else {

		if(type == 'Abundance'){
			ambientObj = subset(seuratObject, subset = is_ambient_cluster == 'Ambient')
			ambMat = ambientObj@assays$RNA@data %>% rowSums %>% prop.table * 100
			return(ambMat)
		}

		if(type == 'Markers'){

			# Find marker genes
			cbMarkers = FindMarkers(seuratObject, ident.1 = 'Ambient',
					ident.2 = 'Non_Ambient', group.by = 'is_ambient_cluster', logfc.threshold = logfc)
			ambEnriched = cbMarkers[cbMarkers[,2] > 0,]
			ambDepleted = cbMarkers[cbMarkers[,2] < 0,]

			toReturn = list(ambEnriched, ambDepleted)
			names(toReturn) = c('Ambient_Enriched_Markers', 'Ambient_Depleted')
			return(toReturn)
		}
	}

}



# Subcluster cleaning function
subCLEAN = function(object, group.by, key, ambientMarkers, batchCorrect = T, res = 0.5, npcs = 10){

	if(batchCorrect == T & length(unique(object$orig.ident)) == 1){
		stop('Only one batch detected in orig.ident. Please set batchCorrect to False or add the batch information to orig.ident')
	}

	# If there are cells without any UMI, remove them
	object = subset(object, subset = nCount_RNA > 0)

	# Subset cells
	Idents(object) = object[[group.by]]
	subcells = WhichCells(object, idents = key)
	subObj = subset(object, cells = subcells)

	# Set default assay to RNA for percentage features
	DefaultAssay(subObj) = 'RNA'
	subObj = NormalizeData(subObj)

	# Calculate ambient marker fraction
	subObj[["percent_ambient"]] = PercentageFeatureSet(object = subObj, features = intersect(ambientMarkers, rownames(subObj)))

	## Subcluster
	print('Running subclustering')
	
	subObj = SCTransform(subObj, ncells = 1000, verbose = F)
	subObj = RunPCA(subObj, verbose = FALSE)
	subObj = RunUMAP(subObj, dims = 1:npcs)

	## Batch correction
	if(batchCorrect == T){

	print('Running batch correction')
	require(harmony)

	subObj = RunHarmony(object = subObj, group.by.vars = 'orig.ident', assay.use="SCT")
	subObj = RunUMAP(subObj, dims = 1:npcs, reduction = 'harmony')
	subObj = FindNeighbors(subObj, dims = 1:npcs, verbose = FALSE, reduction = 'harmony')
	subObj = FindClusters(subObj, verbose = FALSE, resolution = res)
	} else{
	subObj = FindNeighbors(subObj, dims = 1:npcs, verbose = FALSE, reduction = 'pca')
	subObj = FindClusters(subObj, verbose = FALSE, resolution = res)
	}

	## Find Subcluster Markers
	print('Finding cluster markers')
	DefaultAssay(subObj) = 'RNA'
	subObj = NormalizeData(subObj)
	subObjMarkers = FindAllMarkers(subObj)
	subObjMarkers = subObjMarkers[subObjMarkers$p_val_adj < 0.05,]
	subObjMarkers_List = split(subObjMarkers, subObjMarkers$cluster)
	subObjMarkers_List = lapply(subObjMarkers_List, function(x){x$gene})

	## Run Enrichment with Ambient Markers
	print('Enrichment with ambient RNA markers')
	require(GeneOverlap)

	# Find total number of expressed genes for background
	bcg = sum(rowMeans(subObj@assays$RNA@data) > 0.05)

	resgom = newGOM(list(ambientMarkers), subObjMarkers_List, genome.size = bcg)
	oddsr = getMatrix(resgom, name="odds.ratio")
	pvals = getMatrix(resgom, name="pval")

	## Plot Enrichment with Ambient Markers
	pvals_log10FDR = -log10(p.adjust(pvals, method = 'fdr'))
	toplot = data.frame(oddsRatio = oddsr, pval = pvals, log10FDR = pvals_log10FDR)
	toplot$clusters = rownames(toplot)
	toplot$type = 'Ambient_Markers'
	toplot$FDRplot = formatC(p.adjust(as.numeric(toplot$pval)), format = "e", digits = 0)
	toplot$oddsRatio = round(toplot$oddsRatio, digits = 2)

	enrich_plot = ggplot(toplot, aes(clusters, type, fill = log10FDR))+
			 geom_tile(color = "white") +
			 rotate_x_text(90) +
			 theme_classic() +
			 scale_fill_gradient2(midpoint = 2, high = 'red', low = 'white') +
			 geom_text(aes(label = oddsRatio), fontface = 'bold') +
			 xlab('Clusters') + ylab('') +
			 theme(text=element_text(size=20, face = 'bold')) +
			 labs(fill = expression('-log'[10]*'(FDR)')) +
			 coord_flip()
			
	## Plot Ambient Marker Percentage
	ambPerc_plot = ggboxplot(subObj[[]], x = 'seurat_clusters', y = 'percent_ambient', color = 'black', outlier.shape = NA)+
			 rotate_x_text(90) +
			 theme_classic() +
			 xlab('Clusters') + ylab('Percent Ambient Markers') +
			 theme(text=element_text(size=20, face = 'bold'))
			
	## Sub clustering UMAP
	umap_plot = DimPlot(subObj, group.by = 'seurat_clusters', label = T,
			raster = T, label.size = 7) +
			theme(text=element_text(size=20, face = 'bold')) + NoLegend() +
			ggtitle('')

	pdf(paste0('Ambient_QC_Plots_Cluster_', key, '.pdf'), width = 20)
	print(umap_plot + enrich_plot + ambPerc_plot)
	dev.off()

	return(list(subObj, toplot))
}



