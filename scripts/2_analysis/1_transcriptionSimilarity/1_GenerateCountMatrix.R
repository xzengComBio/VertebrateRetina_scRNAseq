suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(do))
suppressMessages(library(future))
suppressMessages(library(SeuratDisk))


source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/')
set.seed(123)

options(future.globals.maxSize = 12000 * 1024^2)

#-----------------------------------------mouse sample info-------------------------------------------#

mouse_retina <- paste0(list.dirs('other_species/Clark2019neuron/count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
mouse_retina <- mouse_retina[-11]
names(mouse_retina) <- c('E11', 'E12', 'E14R1', 'E14R2', 'E16', 'E18R2', 'E18R3',
						 'P0', 'P14', 'P2R2', 'P5', 'P8R1', 'P8R2')

mouse_retina_counts <- Read10X(data.dir = mouse_retina, strip.suffix = TRUE, gene.column = 1)
mouse_retina_seurat <- readRDS('other_species/Clark2019neuron/results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')

mouse_retina_counts <- mouse_retina_counts[,colnames(mouse_retina_seurat)]
colnames(mouse_retina_counts) <- paste0('mouse-', colnames(mouse_retina_counts))
saveRDS(mouse_retina_counts, 'other_species/Clark2019neuron/results/rds/Clark2019neuron_forIntegration.rds')

#-----------------------------------------chicken sample info-------------------------------------------#
sample_paths <- paste0(list.dirs('other_species/Yamagata2021eLife/count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
names(sample_paths) <- c('E18A', 'E18B', 'E18C', 'E18D', 
						'dRGC1','dRGC2','E12A', 'E12B',
						'E12C', 'E12D', 'vRGC1', 'vRGC2')

chicken_counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE, gene.column = 1)
chicken_retina_seurat <- readRDS('other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_filtered_processed.rds')
chicken_counts <- chicken_counts[,colnames(chicken_retina_seurat)]
colnames(chicken_counts) <- paste0('chicken-', colnames(chicken_counts))
saveRDS(chicken_counts, 'other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_forIntegration.rds')

#-----------------------------------------zebrafish sample info-------------------------------------------#

sample_paths <- paste0(list.dirs('other_species/Xu2020Development/count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
names(sample_paths) <- c('14dpf1', '14dpf2', '24hpf', '36hpf', '48hpf1',
						 '48hpf2', '72hpf1', '72hpf2')

zebrafish_retina_counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE, gene.column = 1)
zebrafish_retina_seurat <- readRDS('other_species/Xu2020Development/results/rds/Xu2020Development_raw_processed.rds')

zebrafish_retina_counts <- zebrafish_retina_counts[,colnames(zebrafish_retina_seurat)]
colnames(zebrafish_retina_counts) <- paste0('zebrafish-', colnames(zebrafish_retina_counts))
saveRDS(zebrafish_retina_counts, 'other_species/Xu2020Development/results/rds/Xu2020Development_forIntegration.rds')

