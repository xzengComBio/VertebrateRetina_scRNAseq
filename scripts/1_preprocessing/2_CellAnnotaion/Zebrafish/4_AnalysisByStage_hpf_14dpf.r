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

dpf14 <- sce.list[c(5,6)]
dpf14 <- merge(dpf14$`14dpf1`,dpf14$`14dpf2`)

dpf14 <- dpf14 %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(dpf14) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     dpf14 <- FindClusters(dpf14,resolution = i, verbose= FALSE)
}

feature_genes <- c(	   'hmgb2b', 'stmn1a',
				   'rbpms2b', 'stmn2b', # RGC
				   'tfap2a', 'bhlhe22', #AC
				   'vsx1', 'otx2b', # BC
				   'rem1','onecut1', #HC
				   'nr2e3', 'otx5', # PR
				   'rlbp1a', 'efhd1') # MG

DefaultAssay(dpf14) <- 'SCT'
feature_genes_plot <- FeaturePlot(dpf14, features = feature_genes, ncol = 4, reduction = 'umap')

dpf14_filtered <- subset(dpf14,RNA_snn_res.0.6 != '11' & RNA_snn_res.0.6 != '18')

dpf14_filtered$cellType <- 'RPC'
dpf14_filtered$cellType[dpf14_filtered$RNA_snn_res.0.2  %in% c('2','13')] <- 'BC'
dpf14_filtered$cellType[dpf14_filtered$RNA_snn_res.0.2  %in% c('3','7')] <- 'PR'
dpf14_filtered$cellType[dpf14_filtered$RNA_snn_res.0.2  == '4'] <- 'RGC'
dpf14_filtered$cellType[dpf14_filtered$RNA_snn_res.0.2  == '6'] <- 'AC'
dpf14_filtered$cellType[dpf14_filtered$RNA_snn_res.0.2  == '9'] <- 'HC'

count_matrix <- dpf14_filtered@assays$RNA@counts
rlbp1_plus <- count_matrix['rlbp1a',]
rlbp1_selected_cells <- names(rlbp1_plus[rlbp1_plus > 0])
cluster11 <- colnames(subset(dpf14_filtered, RNA_snn_res.0.6 == '13'))
MG_select <- intersect(rlbp1_selected_cells, cluster11)
dpf14_filtered$cellType[MG_select] <- 'MG'


pr <- subset(dpf14_filtered, cellType == 'PR')

pr <- pr %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(pr) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     pr <- FindClusters(pr,resolution = i, verbose= FALSE)
}

dpf14_filtered$cellType[names(pr$RNA_snn_res.0.2[pr$RNA_snn_res.0.2 == '1'])] <- 'Cone'
dpf14_filtered$cellType[names(pr$RNA_snn_res.0.2[pr$RNA_snn_res.0.2 == '4'])] <- 'Rod'

dpf14_filtered$cellType[dpf14_filtered$cellType == 'PR'] <- 'PR Precur.'


dpf14_filtered <- RunUMAP(dpf14_filtered,dims = 1:20,reduction = 'canek', n.neighbors = 70, min.dist = 1,verbose = FALSE)
saveRDS(dpf14_filtered, 'results/rds/dpf14_filtered.rds')

