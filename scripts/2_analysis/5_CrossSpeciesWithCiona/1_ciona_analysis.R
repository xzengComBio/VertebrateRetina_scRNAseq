suppressMessages(library(Seurat))
suppressMessages(library(harmony))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(clustree))
suppressMessages(library(readxl))
suppressMessages(library(do))
suppressMessages(library(igraph))
suppressMessages(library(future))
suppressMessages(library(SeuratDisk))
suppressMessages(library(harmony))
suppressMessages(library(Canek))

source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/')
set.seed(123)

options(future.globals.maxSize = 12000 * 1024^2)

#-----------------------------------------

sample_paths <- paste0(list.dirs('raw_data', recursive = FALSE), '/filtered_feature_bc_matrix/')
sample_paths <- sample_paths[c(19,20,21,22,23,24)]

names(sample_paths) <- c('LTB2.1', 'LTB2.1', 'LTB2.2',paste0('lv','.',c(3,4,1)))
ciona_counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE, gene.column = 1)

lv_ner_anno <- read.table('Cao2019_supp/myCNS.LTB2.lv.MSTRG.txt', sep = '\t',header = TRUE)
colnames(lv_ner_anno) <- c('barcodes','cellType')
lv_ner_anno$newBarcodes <- str_split(lv_ner_anno$barcodes, pattern = '-', simplify = TRUE)[,1]

select_barcodes <- intersect(lv_ner_anno$newBarcodes, colnames(ciona_counts))
ner_cells <- ciona_counts[,select_barcodes]
colnames(ner_cells) <- paste0('ciona_', colnames(ner_cells))
saveRDS(ner_cells, 'result/larve/larva_ner/larva_ner_forIntegration.rds')

ner_cells <- CreateSeuratObject(ner_cells)
ner_cells$sample <- str_split(colnames(ner_cells), pattern = '_', simplify = TRUE)[,2]

rownames(lv_ner_anno) <- paste0('ciona_',lv_ner_anno$newBarcodes)
lv_ner_anno <- lv_ner_anno[,-1]
ner_cells <- AddMetaData(ner_cells, lv_ner_anno)

ner_cells <- SCTransform(ner_cells, method = "glmGamPoi")
ner_cells <- RunPCA(ner_cells,verbose = FALSE)

p1 <- DimPlot(ner_cells, reduction = 'pca', group.by = 'sample', cols = plotColor) +
			theme_bw() +
			theme(panel.grid.major=element_blank(),
				  panel.grid.minor=element_blank(),
				  axis.title = element_text(face = "bold",size = rel(1))) + 
			ggtitle('Before integration')

ner_cells <- RunCanek(ner_cells, "sample")
ner_cells<- ScaleData(ner_cells)
ner_cells <- RunPCA(ner_cells)

p2 <-  DimPlot(ner_cells, reduction = 'pca', group.by = 'sample', cols = plotColor) +
			theme_bw() +
			theme(panel.grid.major=element_blank(),
				  panel.grid.minor=element_blank(),
				  axis.title = element_text(face = "bold",size = rel(1))) + 
			ggtitle('After integration')

ner_cells <- RunUMAP(ner_cells, dims = 1:20, verbose = FALSE)
saveRDS(ner_cells,'result/larve/larva_ner/larva_ner_processed.rds')


