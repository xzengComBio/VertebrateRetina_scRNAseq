suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(clustree))
suppressMessages(library(Canek))

# set the working directory
source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Xu2020Development/')
set.seed(123)

# load the raw data
sce.list <- readRDS('results/rds/Xu2020Development_raw.rds')

hpf48 <- sce.list[c(5,6)]
hpf48 <- merge(hpf48$`48hpf1`,hpf48$`48hpf2`)

hpf48 <- hpf48 %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(hpf48) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     hpf48 <- FindClusters(hpf48,resolution = i, verbose= FALSE)
}

hpf48_clstree <- clustree(hpf48@meta.data, prefix = 'RNA_snn_res.') + NoLegend()
cls1 <- myPlotUMAP_clusters(hpf48, groups = 'RNA_snn_res.0.1')
cls2 <- myPlotUMAP_clusters(hpf48, groups = 'RNA_snn_res.0.2')

myplot <- hpf48_clstree | (cls1/cls2)
ggsave('results/analyze_by_stages/hpf48_clusters_umap.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 20,
    	 		height = 18)

feature_genes <- c('fabp11a','npm1a', 'her9',  # cluster1 RPC
				  'fabp7a', 'her4.1', 'her4.2', # cluster2 RPC
				   'dla', 'atoh7', 'neurod4', # cluster3 RPC
				   'rbpms2b', 'stmn2b', # RGC
				   'tfap2a', 'bhlhe22', #AC
				   'vsx1', 'otx2b', # BC
				   'rem1','onecut1', #HC
				   'nr2e3', 'otx5', # PR
				   'rlbp1a', 'efhd1') # MG

DefaultAssay(hpf48) <- 'SCT'
feature_genes_plot <- FeaturePlot(hpf48, features = feature_genes, ncol = 4, reduction = 'umap')
ggsave('results/analyze_by_stages/hpf48_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 30,
    	 		height = 28)

hpf48_filtered <- subset(hpf48, RNA_snn_res.0.1 != '7')

DefaultAssay(hpf48_filtered) <- 'RNA'
hpf48_filtered <- hpf48_filtered %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE,assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(hpf48_filtered) <- 'RNA'
resolutions <- c(0.2,0.4,0.6, 1)
for (i in resolutions){
     hpf48_filtered <- FindClusters(hpf48_filtered,resolution = i, verbose= FALSE)
}

DefaultAssay(hpf48_filtered) <- 'SCT'
feature_genes_plot <- FeaturePlot(hpf48_filtered, features = feature_genes, ncol = 4,reduction = 'umap')
ggsave('results/analyze_by_stages/hpf48_filtered_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 30,
    	 		height = 28)

hpf48_clstree <- clustree(hpf48_filtered@meta.data, prefix = 'RNA_snn_res.') + NoLegend()
cls1 <- myPlotUMAP_clusters(hpf48_filtered, groups = 'RNA_snn_res.0.2')
cls2 <- myPlotUMAP_clusters(hpf48_filtered, groups = 'RNA_snn_res.0.4')
myplot <- hpf48_clstree | (cls1/cls2)
ggsave('results/analyze_by_stages/hpf48_filtered_clusters_umap.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 20,
    	 		height = 18)

hpf48_filtered$cellType <- 'others'
hpf48_filtered$cellType[hpf48_filtered$RNA_snn_res.0.2 %in% c('1','4')] <- 'BC'
hpf48_filtered$cellType[hpf48_filtered$RNA_snn_res.0.2 == '2'] <- 'AC'
hpf48_filtered$cellType[hpf48_filtered$RNA_snn_res.0.2 == '3'] <- 'RGC'
hpf48_filtered$cellType[hpf48_filtered$RNA_snn_res.0.2 %in% c('5','6')] <- 'PR'
hpf48_filtered$cellType[hpf48_filtered$RNA_snn_res.0.2 == '7'] <- 'HC'

hpf48_filtered$cellType[hpf48_filtered$RNA_snn_res.0.4  %in% c('1','6', '11', '12','13')] <- 'RPC'

count_matrix <- hpf48_filtered@assays$RNA@counts
rlbp1_plus <- count_matrix['rlbp1a',]
rlbp1_selected_cells <- names(rlbp1_plus[rlbp1_plus > 0])
cluster11 <- colnames(subset(hpf48_filtered, RNA_snn_res.0.4 == '11'))
MG_select <- intersect(rlbp1_selected_cells, cluster11)
hpf48_filtered$cellType[MG_select] <- 'MG'


hpf48_filtered$module <- 'others'
hpf48_filtered$module[hpf48_filtered$RNA_snn_res.0.4 == '13'] <- '1'
hpf48_filtered$module[hpf48_filtered$RNA_snn_res.0.4 == '1'] <- '2'
hpf48_filtered$module[hpf48_filtered$RNA_snn_res.0.4 %in% c('6','12')] <- '3'

hpf48_filtered_selected <- subset(hpf48_filtered, cellType != 'others')



hpf48_filtered_selected <- hpf48_filtered_selected %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE,assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

hpf48_filtered_selected <- RunUMAP(hpf48_filtered_selected,dims = 1:20,reduction = 'canek', n.neighbors = 70, min.dist = 1,verbose = FALSE)

saveRDS(hpf48_filtered_selected, 'results/rds/hpf48_filtered.rds')

#--------------------------------------------------72hpf----------------------------------------------------
hpf72 <- sce.list[c(7,8)]
hpf72 <- merge(hpf72$`72hpf1`,hpf72$`72hpf2`)

hpf72 <- hpf72 %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(hpf72) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     hpf72 <- FindClusters(hpf72,resolution = i, verbose= FALSE)
}

feature_genes <- c(	   'hmgb2b', 'stmn1a',
				   'rbpms2b', 'stmn2b', # RGC
				   'tfap2a', 'bhlhe22', #AC
				   'vsx1', 'otx2b', # BC
				   'rem1','onecut1', #HC
				   'nr2e3', 'otx5', # PR
				   'rlbp1a', 'efhd1') # MG

hpf72_filtered <- subset(hpf72, RNA_snn_res.0.4 != '18')
hpf72_filtered$cellType <- 'BC'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 %in% c('0', '15', '10', '11')] <- 'AC'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 == '1'] <- 'RGC'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 %in% c('4', '16')] <- 'RPC'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 %in% c('6','12')] <- 'PR Precur.'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 == '7'] <- 'HC'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 == '9'] <- 'Late RPC'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 == '13'] <- 'Cone'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 == '14'] <- 'MG'
hpf72_filtered$cellType[hpf72_filtered$RNA_snn_res.0.4 == '17'] <- 'Rod'


saveRDS(hpf72_filtered,'results/rds/hpf72_filtered.rds')


