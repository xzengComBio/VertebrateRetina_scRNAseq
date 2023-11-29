suppressMessages(library(Seurat))
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
suppressMessages(library(Canek))

source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Yamagata2021eLife/')
set.seed(123)

options(future.globals.maxSize = 12000 * 1024^2)

#--------------------------------------Analysis of E18 --------------------------------
chicken <- readRDS('results/rds/Yamagata2021eLife_filtered_processed.rds')
E18 <- subset(chicken, stage=='E18')

DefaultAssay(E18) <- 'RNA'
E18 <- SCTransform(E18, method = "glmGamPoi")
E18 <- RunPCA(E18)

# Integration
E18 <- RunCanek(E18, "orig.ident")
E18<- ScaleData(E18)
E18 <- RunPCA(E18)

E18 <- FindNeighbors(E18, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.2,0.4, 0.6)
for (i in resolutions){
     E18 <- FindClusters(E18,resolution = i, verbose= FALSE)
}
E18 <- RunUMAP(E18, dims = 1:20)


E18$newCellType <- 'Others'

E18_flitered <- subset(E18, cellType != 'Others')

saveRDS(E18_flitered,'results/rds/E18_filtered.rds')
#--------------------------------------Analysis of E16 --------------------------------
#chicken <- readRDS('results/rds/Yamagata2021eLife_filtered_processed.rds')
E16 <- subset(chicken, stage=='E16')

DefaultAssay(E16) <- 'RNA'
E16 <- subset(E16, cellType != 'Others')

E16 <- SCTransform(E16, method = "glmGamPoi")
E16 <- RunPCA(E16)

# Integration
E16 <- RunCanek(E16, "orig.ident")
E16<- ScaleData(E16)
E16 <- RunPCA(E16)


E16 <- FindNeighbors(E16, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.02,0.05,0.2)
for (i in resolutions){
     E16 <- FindClusters(E16,resolution = i, verbose= FALSE)
}
E16 <- RunUMAP(E16, dims = 1:20)

saveRDS(E16,'results/rds/E16_filtered.rds')
#--------------------------------------Analysis of E12 --------------------------------
#chicken <- readRDS('results/rds/Yamagata2021eLife_filtered_processed.rds')
E12 <- subset(chicken, stage=='E12')

DefaultAssay(E12) <- 'RNA'

E12 <- E12 %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

E12 <- FindNeighbors(E12, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.02,0.05,0.2)
for (i in resolutions){
     E12 <- FindClusters(E12,resolution = i, verbose= FALSE)
}


E12$cellType <- 'AC'
E12$cellType[E12$SCT_snn_res.0.2 == '0'] <- 'RGC'
E12$cellType[E12$SCT_snn_res.0.2 == '3'] <- 'OL'
E12$cellType[E12$SCT_snn_res.0.2 == '7'] <- 'MG'
E12$cellType[E12$SCT_snn_res.0.2 == '11'] <- 'HC'
E12$cellType[E12$SCT_snn_res.0.2 == '12'] <- 'BC'
E12$cellType[E12$SCT_snn_res.0.2 == '15'] <- 'PR'

saveRDS(E12,'results/rds/E12_filtered.rds')

#--------------------------------------Merge all annotaion --------------------------------

chicken$cellType[names(E12$cellType[E12$cellType == 'AC'])] <- 'AC'
chicken$cellType[names(E12$cellType[E12$cellType == 'RGC'])] <- 'RGC'
chicken$cellType[names(E12$cellType[E12$cellType == 'OL'])] <- 'OL'
chicken$cellType[names(E12$cellType[E12$cellType == 'MG'])] <- 'MG'
chicken$cellType[names(E12$cellType[E12$cellType == 'HC'])] <- 'HC'
chicken$cellType[names(E12$cellType[E12$cellType == 'BC'])] <- 'BC'
chicken$cellType[names(E12$cellType[E12$cellType == 'PR'])] <- 'PR Precurs.'

cellTypes_keep <- c('AC', 'BC', 'Cone', 'HC', 'MG', 'PR Precurs.', 'RGC', 'Rod')
sce_merge_final_filtered <- subset(chicken, cellType %in% cellTypes_keep)

DefaultAssay(sce_merge_final_filtered) <- 'RNA'
sce_merge_final_filtered <- sce_merge_final_filtered %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

sce_merge_final_filtered <- RunUMAP(sce_merge_final_filtered,dims = 1:20,reduction = 'canek', verbose = FALSE,n.neighbors=50,min.dist = 0.4)

saveRDS(sce_merge_final_filtered,'results/rds/Yamagata2021eLife_filtered_processed_final.rds')


#-----------------------------generate loom files for SCENIC --------------------------------------

chicken_retina_seurat <- readRDS('results/rds/Yamagata2021eLife_filtered_processed.rds')
chicken_counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE, gene.column = 1)
sce <- CreateSeuratObject(chicken_counts,names.field = 1, min.cells = 3, min.features = 200, project = 'chick retina')

sce <- subset(sce, cells = colnames(chicken_retina_seurat))

E18_names <- names(chicken_retina_seurat$stage[chicken_retina_seurat$stage == 'E18'])
E18 <- subset(sce, cells = E18_names)
as.loom(E18, filename ='other_species/Yamagata2021eLife/results/rds/E18.loom',overwrite = TRUE)


E12_names <- names(chicken_retina_seurat$stage[chicken_retina_seurat$stage == 'E12'])
E12 <- subset(sce, cells = E12_names)
E12_downsample <- E12[,sample(colnames(E12), size = 10000, replace = FALSE)]
as.loom(E12_downsample, filename ='other_species/Yamagata2021eLife/results/rds/E12.loom',overwrite = TRUE)

E16_names <- names(chicken_retina_seurat$stage[chicken_retina_seurat$stage == 'E16'])
E16 <- subset(sce, cells = E16_names)
as.loom(E16, filename ='other_species/Yamagata2021eLife/results/rds/E16.loom',overwrite = TRUE)
