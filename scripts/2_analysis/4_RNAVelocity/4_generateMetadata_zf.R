suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(clustree))
suppressMessages(library(Canek))

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

