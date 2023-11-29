suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(clustree))
suppressMessages(library(Canek))

# set the working directory
source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Xu2020Development/')
set.seed(123)


hpf24 <- readRDS('results/rds/hpf24_filtered.rds')
hpf36 <- readRDS('results/rds/hpf36_filtered.rds')
hpf48 <- readRDS('results/rds/hpf48_filtered.rds')
hpf72 <- readRDS('results/rds/hpf72_filtered.rds')
dpf14 <- readRDS('results/rds/dpf14_filtered.rds')

hpf24$cellType <-  "RPC"

zebrafish_merge_list <- list(hpf24 = hpf24,
						hpf36 = hpf36,
						hpf48 = hpf48,
						hpf72 = hpf72,
						dpf14 = dpf14)

zebrafish_merge <- Reduce(merge, zebrafish_merge_list)

DefaultAssay(zebrafish_merge) <- 'RNA'
zebrafish_merge <- zebrafish_merge %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

zebrafish_merge$cellType[is.na(zebrafish_merge$cellType)] <- 'RPC'

 zebrafish_merge <- RunUMAP(zebrafish_merge, dims = 1:20,reduction = 'canek', verbose = FALSE, min.dist = 0.6, n.neighbors = 70)

DefaultAssay(zebrafish_merge) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     zebrafish_merge <- FindClusters(zebrafish_merge,resolution = i, verbose= FALSE)
}

zebrafish_merge$newCellType <- 'RPC'
zebrafish_merge$newCellType[zebrafish_merge$RNA_snn_res.0.1 == '1'] <- 'BC'
zebrafish_merge$newCellType[(zebrafish_merge$RNA_snn_res.0.1 == '3') & (zebrafish_merge$cellType %in% c('Cone', 'Rod', 'PR'))] <- 'PR'
zebrafish_merge$newCellType[zebrafish_merge$RNA_snn_res.0.1 %in% c('4', '7')] <- 'AC'
zebrafish_merge$newCellType[zebrafish_merge$RNA_snn_res.0.1 == '5'] <- 'RGC'
zebrafish_merge$newCellType[zebrafish_merge$RNA_snn_res.0.1 == '6'] <- 'HC'
zebrafish_merge$newCellType[(zebrafish_merge$RNA_snn_res.0.1 == '8') & (zebrafish_merge$cellType == 'MG') ] <- 'MG'


pr <- subset(zebrafish_merge, newCellType == 'PR')

DefaultAssay(pr) <- 'RNA'
pr <- pr %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

pr <-  RunUMAP(pr, dims = 1:20,reduction = 'canek', verbose = FALSE,min.dist = 0.6)
DefaultAssay(pr) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     pr <- FindClusters(pr,resolution = i, verbose= FALSE)
}

pr$pr_type <- 'PR Precurs.'
pr$pr_type[pr$RNA_snn_res.0.2 == '6'] <- 'Rod'
pr$pr_type[pr$RNA_snn_res.0.2 == '4'] <- 'Cone'

saveRDS(pr, 'results/rds/all_pr.rds')

zebrafish_merge$newCellType[names(pr$pr_type[pr$pr_type == 'PR'])] <- 'PR Precurs.'
zebrafish_merge$newCellType[names(pr$pr_type[pr$pr_type == 'Rod'])] <- 'Rod'
zebrafish_merge$newCellType[names(pr$pr_type[pr$pr_type == 'Cone'])] <- 'Cone'

saveRDS(zebrafish_merge, 'results/rds/merge_all.rds')

#---- ## generate metadata for RNA velocity

zebrafish_merge <- readRDS('results/rds/merge_all.rds')
filtered_cells <- read.table('../../velocyto/results/zebrafish/merge/filtered_cells.txt', sep = '\t')
zebrafish_merge_filtered <- subset(zebrafish_merge, cells = filtered_cells$V1)

zebrafish_merge_filtered <- zebrafish_merge_filtered %>% 
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)



output_file <- c('../../velocyto/Xu2020Development/merge/')
write.csv(Embeddings(zebrafish_merge_filtered, reduction = 'umap'), 
		file = paste0(output_file,'cell_embeddings.csv'),
		quote = FALSE)

write.csv(Cells(zebrafish_merge_filtered), 
		file = paste0(output_file,'cellID_obs.csv'),
		row.names = FALSE,
		quote = FALSE)

write.csv(zebrafish_merge_filtered@meta.data[,c('RNA_snn_res.0.2','newCellType', 'module','Phase', 'S.Score', 'G2M.Score')],
		file = paste0(output_file,'metadata.csv'),
		quote = FALSE)

saveRDS(zebrafish_merge_filtered,paste0(output_file, 'filtered_merge.rds'))


