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
sample_paths <- paste0(list.dirs('count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
names(sample_paths) <- c('E18A', 'E18B', 'E18C', 'E18D', 
						'dRGC1','dRGC2','E12A', 'E12B',
						'E12C', 'E12D', 'vRGC1', 'vRGC2')

counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE)
sce <- CreateSeuratObject(counts,names.field = 1, min.cells = 3, min.features = 200, project = 'chick retina')

mito_genes <- c("ND1", "ENSGALG00000043768", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")

sce[['percent.mt']] <- PercentageFeatureSet(sce, features = mito_genes)
sce <- subset(sce, subset = nFeature_RNA >= 600 & percent.mt <= 10)

E12_count_matrix <- read.csv('supp/GSE159107_E12chick_count.matrix.csv')
E18_count_matrix <- read.csv('supp/GSE159107_E18chick_count.matrix.csv')
E16_count_matrix <- read.csv('supp/GSE159107_E16chick_count.matrix.csv')

cells_use <- c(colnames(E12_count_matrix)[-1], colnames(E16_count_matrix)[-1], colnames(E18_count_matrix)[-1])


cells_use <- str_replace_all(cells_use, c('\\.1' = '',
								  'ChickenE12S1' = 'E12A',
								  'ChickenE12S2' = 'E12B',
								  'ChickenE12S3' = 'E12C',
								  'ChickenE12S4' = 'E12D',
								  'ChickendRGC1' = 'dRGC1',
								  'ChickendRGC2' = 'dRGC2',
								  'ChickenvRGC1' = 'vRGC1',
								  'ChickenvRGC2' = 'vRGC2',
								  'Chicken1' = 'E18'))

cells_use <- intersect(colnames(sce), cells_use)
sce_use <- subset(sce, cells = cells_use)

sce_use <- SplitObject(sce_use, split.by = 'orig.ident')

saveRDS(sce_use,'results/rds/Yamagata2021eLife_filtered.rds')

### Run 2_Integration.sh

sce_merge <- readRDS('results/rds/Yamagata2021eLife_filtered_integrated.rds')


sce_merge <- FindNeighbors(sce_merge, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.2,0.4, 0.6, 0.8,1)
for (i in resolutions){
     sce_merge <- FindClusters(sce_merge,resolution = i, verbose= FALSE)
}

sce_merge <- RunUMAP(sce_merge, dims = 1:20)

sce_merge$stage <- 'E18'
sce_merge$stage[sce_merge$orig.ident == "E12A"] <- "E12"
sce_merge$stage[sce_merge$orig.ident == "E12B"] <- "E12"
sce_merge$stage[sce_merge$orig.ident == "E12C"] <- "E12"
sce_merge$stage[sce_merge$orig.ident == "E12D"] <- "E12"

sce_merge$stage[sce_merge$orig.ident == "dRGC1"] <- "E16"
sce_merge$stage[sce_merge$orig.ident == "dRGC2"] <- "E16"
sce_merge$stage[sce_merge$orig.ident == "vRGC1"] <- "E16"
sce_merge$stage[sce_merge$orig.ident == "vRGC2"] <- "E16"


cell_anno <- read.csv('supp/SCP1159/metadata/Chick_retina_atlas_meta.csv')
cell_anno <- cell_anno[-1,]
cell_anno$new_name <-  str_replace_all(cell_anno$NAME, c('-1' = '',
                                                                'ChickenE12S1' = 'E12A',
                                                                'ChickenE12S2' = 'E12B',
                                                                'ChickenE12S3' = 'E12C',
                                                                'ChickenE12S4' = 'E12D',
                                                                'ChickendRGC1' = 'dRGC1',
                                                                'ChickendRGC2' = 'dRGC2',
                                                                'ChickenvRGC1' = 'vRGC1',
                                                                'ChickenvRGC2' = 'vRGC2',
                                                                'Chicken1' = 'E18'))
rownames(cell_anno) <- cell_anno$new_name
sce_merge$cellType <- 'Others'
sce_merge$anno <- 'Others'


cells_overlap <- intersect(colnames(sce_merge), rownames(cell_anno))
sce_merge$anno[cells_overlap] <- cell_anno[cells_overlap,]$Cluster

chick_AC <- read.csv('supp/SCP1159/cluster/Chick_AC_clusterfile.txt', sep = '\t')
chick_AC$new_name <- str_replace_all(chick_AC$NAME, c('-1'= ''))
sce_merge$cellType[colnames(sce_merge) %in% chick_AC$new_name] <- 'AC'

chick_BP <- read.csv('supp/SCP1159/cluster/Chick_BP_clusterfile.txt', sep = '\t')
chick_BP$new_name <- str_replace_all(chick_BP$NAME, c('-1'= ''))
sce_merge$cellType[colnames(sce_merge) %in% chick_BP$new_name] <- 'BC'

chick_HC <- read.csv('supp/SCP1159/cluster/Chick_HC_clusterfile.txt', sep = '\t')
chick_HC$new_name <- str_replace_all(chick_HC$NAME, c('-1'= ''))
sce_merge$cellType[colnames(sce_merge) %in% chick_HC$new_name] <- 'HC'

chick_MG <- read.csv('supp/SCP1159/cluster/Chick_MG_clusterfile.txt', sep = '\t')
chick_MG$new_name <- str_replace_all(chick_MG$NAME, c('-1'= ''))
sce_merge$cellType[colnames(sce_merge) %in% chick_MG$new_name] <- 'MG'

chick_PR <- read.csv('supp/SCP1159/cluster/Chick_PR_clusterfile.txt', sep = '\t')
chick_PR$new_name <- str_replace_all(chick_PR$NAME, c('-1'= ''))
sce_merge$cellType[colnames(sce_merge) %in% chick_PR$new_name] <- 'PR'

chick_RGC <- read.csv('supp/SCP1159/cluster/Chick_RGC_clusterfile.txt', sep = '\t')
chick_RGC$new_name <- str_replace_all(chick_RGC$NAME, c('-1'= '', 'Chicken'=''))
sce_merge$cellType[colnames(sce_merge) %in% chick_RGC$new_name] <- 'RGC'

sce_pr <- subset(sce_merge, cellType == 'PR')

sce_pr$cellType <- 'Cone'
sce_pr$cellType[sce_pr$anno %in% c('DevRods', 'Rods')] <- 'Rod'


sce_merge$cellType[names(sce_pr$cellType[sce_pr$cellType == 'Cone'])] <- 'Cone'
sce_merge$cellType[names(sce_pr$cellType[sce_pr$cellType == 'Rod'])] <- 'Rod'


saveRDS(sce_merge,'results/rds/Yamagata2021eLife_filtered_processed.rds')
