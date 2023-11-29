suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(clustree))

# set the working directory
source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Xu2020Development/')
set.seed(123)

# load the raw data
sce.list <- readRDS('results/rds/Xu2020Development_raw.rds')

#--------------------------------------------------24hpf----------------------------------------------------

hpf24 <- sce.list$`24hpf`

# conventional analysis pipeline
hpf24 <- hpf24 %>%
		 SCTransform(variable.features.n = 1000, vars.to.regress = c('S.Score', 'G2M.Score', 'percent.mt',verbose = FALSE)) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

# clustering
resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
     hpf24 <- FindClusters(hpf24,resolution = i, verbose= FALSE)
}



# check the relationship between clusters in different resolutions
hpf24_clstree <- clustree(hpf24@meta.data, prefix = 'SCT_snn_res.') + NoLegend()

cls1 <- myPlotUMAP_clusters(hpf24, groups = 'SCT_snn_res.0.05')
cls2 <- myPlotUMAP_clusters(hpf24, groups = 'SCT_snn_res.0.1')

myplot <- hpf24_clstree | (cls1/cls2)
ggsave('results/analyze_by_stages/hpf24_clusters_umap.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 20,
    	 		height = 18)
#check the marker genes in the paper
feature_genes <- c('fabp11a','npm1a', 'her9', 'fabp7a', 'her4.1', 'her4.2', 'dla', 'atoh7', 'neurod4')
feature_genes_plot <- FeaturePlot(hpf24, features = feature_genes, ncol = 3)

ggsave('results/analyze_by_stages/hpf24_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 22,
    	 		height = 20)

# filter out low quality cells
hpf24_filtered <- subset(hpf24,SCT_snn_res.0.05 == '0')

# reanalysis
DefaultAssay(hpf24_filtered) <- 'RNA'

hpf24_filtered <- hpf24_filtered %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
     hpf24_filtered <- FindClusters(hpf24_filtered,resolution = i, verbose= FALSE)
}

hpf24_filtered_clstree <- clustree(hpf24@meta.data, prefix = 'SCT_snn_res.') + NoLegend()
cls1 <- myPlotUMAP_clusters(hpf24_filtered, groups = 'SCT_snn_res.0.2')
cls2 <- myPlotUMAP_clusters(hpf24_filtered, groups = 'SCT_snn_res.0.4')
myplot <- hpf24_filtered_clstree | (cls1/cls2)

ggsave('results/analyze_by_stages/hpf24_filtered_clustree.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 22,
    	 		height = 20)

feature_genes_plot <- FeaturePlot(hpf24_filtered, features = feature_genes, ncol = 3)
ggsave('results/analyze_by_stages/hpf24_filtered_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 22,
    	 		height = 20)

hpf24_filtered$module <- '1'
hpf24_filtered$module[hpf24_filtered$SCT_snn_res.0.4 %in% c('0','1','4')] <- '2'
hpf24_filtered$module[hpf24_filtered$SCT_snn_res.0.4 == '5'] <- '3'

Idents(hpf24_filtered) <- hpf24_filtered$module
saveRDS(hpf24_filtered, 'results/rds/hpf24_filtered.rds')


#--------------------------------------------------36hpf----------------------------------------------------
hpf36 <- sce.list$`36hpf`

hpf36 <- hpf36 %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     hpf36 <- FindClusters(hpf36,resolution = i, verbose= FALSE)
}

hpf36_clstree <- clustree(hpf36@meta.data, prefix = 'SCT_snn_res.') + NoLegend()

cls1 <- myPlotUMAP_clusters(hpf36, groups = 'SCT_snn_res.0.1')
cls2 <- myPlotUMAP_clusters(hpf36, groups = 'SCT_snn_res.0.2')

myplot <- hpf36_clstree | (cls1/cls2)
ggsave('results/analyze_by_stages/hpf36_clusters_umap.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 20,
    	 		height = 18)

feature_genes <- c('fabp11a','npm1a', 'her9', 'fabp7a', 'her4.1', 'her4.2', 'dla', 'atoh7', 'neurod4', 'elavl3', 'stmn1b', 'hmgb3a')
feature_genes_plot <- FeaturePlot(hpf36, features = feature_genes, ncol = 3)
ggsave('results/analyze_by_stages/hpf36_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 22,
    	 		height = 24)
hpf36_filtered <- subset(hpf36, SCT_snn_res.0.2 %in% c('0','1','2','4','5','7'))

DefaultAssay(hpf36_filtered) <- 'RNA'
hpf36_filtered <- hpf36_filtered %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     hpf36_filtered <- FindClusters(hpf36_filtered,resolution = i, verbose= FALSE)
}
feature_genes_plot <- FeaturePlot(hpf36_filtered, features = feature_genes, ncol = 3)
ggsave('results/analyze_by_stages/hpf36_filtered_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 22,
    	 		height = 24)

hpf36_clstree <- clustree(hpf36_filtered@meta.data, prefix = 'SCT_snn_res.') + NoLegend()

cls1 <- myPlotUMAP_clusters(hpf36_filtered, groups = 'SCT_snn_res.0.1')
cls2 <- myPlotUMAP_clusters(hpf36_filtered, groups = 'SCT_snn_res.0.2')

myplot <- hpf36_clstree | (cls1/cls2)
ggsave('results/analyze_by_stages/hpf36_filtered_clusters_umap.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 20,
    	 		height = 18)

# Rerun hpf_filtered
hpf36_filtered <- subset(hpf36_filtered, SCT_snn_res.0.2 %in% c('0','1','2','4','5'))

DefaultAssay(hpf36_filtered) <- 'RNA'
hpf36_filtered <- hpf36_filtered %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     hpf36_filtered <- FindClusters(hpf36_filtered,resolution = i, verbose= FALSE)
}


feature_genes_plot <- FeaturePlot(hpf36_filtered, features = feature_genes, ncol = 3)
ggsave('results/analyze_by_stages/hpf36_filtered_feature_genes_plot.pdf',
    	 		feature_genes_plot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 22,
    	 		height = 24)

hpf36_clstree <- clustree(hpf36_filtered@meta.data, prefix = 'SCT_snn_res.') + NoLegend()

cls1 <- myPlotUMAP_clusters(hpf36_filtered, groups = 'SCT_snn_res.0.1')
cls2 <- myPlotUMAP_clusters(hpf36_filtered, groups = 'SCT_snn_res.0.2')

myplot <- hpf36_clstree | (cls1/cls2)
ggsave('results/analyze_by_stages/hpf36_filtered_clusters_umap.pdf',
    	 		myplot,
    	 		dpi = 200,
    	 		units = 'cm',
    	 		width = 20,
    	 		height = 18)

hpf36_filtered$cellType <- 'RPC'
hpf36_filtered$cellType[hpf36_filtered$SCT_snn_res.0.2 == '3'] <- 'RGC'

hpf36_filtered$module <- '2'
hpf36_filtered$module[hpf36_filtered$SCT_snn_res.0.2 == '0'] <- '3'
hpf36_filtered$module[hpf36_filtered$SCT_snn_res.0.2 == '3'] <- 'RGC'
hpf36_filtered$module[hpf36_filtered$SCT_snn_res.0.2 == '4'] <- '1'

saveRDS(hpf36_filtered, 'results/rds/hpf36_filtered.rds')

