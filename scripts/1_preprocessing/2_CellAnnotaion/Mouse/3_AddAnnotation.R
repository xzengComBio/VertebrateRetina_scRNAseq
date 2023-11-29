
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(clustree))
suppressMessages(library(readxl))
suppressMessages(library(do))
suppressMessages(library(igraph))
suppressMessages(library(SeuratDisk))
suppressMessages(library(homologene))

source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Clark2019neuron/')
set.seed(123)

options(future.globals.maxSize = 15000 * 1024^2)

#----------------------------------------------------------------------------------

sce_merge <- readRDS('results/rds/Clark2019neuron_raw_integrated.rds')

sce_merge <- FindNeighbors(sce_merge, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.2,0.4, 0.6, 0.8,1)
for (i in resolutions){
     sce_merge <- FindClusters(sce_merge,resolution = i, verbose= FALSE)
}

# read the annotation file
anno <- read.csv('supp/GSE118614_barcodes.tsv', sep = '\t')

anno$barcode <- str_replace_all(anno$barcode, c('\\.' = '_',
												'-1' = '', 
												'E14_rep1'='E14R1',
												'E14_rep2'='E14R2', 
												'E18_rep2' = 'E18R2',
												'E18_rep3'='E18R3', 
												'P2_rep2'='P2R2', 
												'P2_rep3' = 'P2R2',
												'E12_rep1'='E12',
												'P8_rep1'='P8R1',
												'P8_rep2'='P8R2'))

anno <- anno[anno$sample != 'P2_rep2',]
rownames(anno) <- anno$barcode

# keep the cells with annotation
annotated_cells <- intersect(anno$barcode, colnames(sce_merge))
sce_merge <- subset(sce_merge, cells = annotated_cells)

# add annotation
sce_merge$cellType <- anno[colnames(sce_merge),]$umap2_CellType

# extract retinal cells for further analysis
KeepcellTypes <- c('Amacrine Cells', 'Bipolar Cells', 'Cones', 'Early RPCs', 'Horizontal Cells', 'Late RPCs',
			     'Muller Glia', 'Neurogenic Cells', 'Photoreceptor Precursors', 'Retinal Ganglion Cells',
			     'Rods')
sce_merge <- subset(sce_merge, cellType %in% KeepcellTypes)

# calculate the cell cycle 
DefaultAssay(sce_merge) <- 'RNA'
mcc.s <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090)$`10090`
mcc.g2m <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090)$`10090`

sce_merge <- CellCycleScoring(sce_merge, s.features = mcc.s, g2m.features = mcc.g2m)

# run umap 
sce_merge <- RunUMAP(sce_merge, dims = 1:20)

saveRDS(sce_merge,'results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')


#----------------- generate loom file for SCENIC (by stage) -------
sce_finalKeep.list <- SplitObject(sce_merge, split.by = 'stage')

stage_name <- names(sce_finalKeep.list)
for (x in 1:length(sce_finalKeep.list)){
	temp <- sce_finalKeep.list[[x]]
	DefaultAssay(temp) <- 'RNA'
	as.loom(temp, filename = paste0('results/rds/',stage_name[x],'.loom'),overwrite = TRUE)
	}


#----------------- generate metadata for RNA velocity ---------------

output_file <- c('../../velocyto/Clark2019neuron/merge/')
write.csv(Embeddings(sce_merge, reduction = 'umap'), 
		file = paste0(output_file,'cell_embeddings.csv'),
		quote = FALSE)

write.csv(Cells(sce_merge), 
		file = paste0(output_file,'cellID_obs.csv'),
		row.names = FALSE,
		quote = FALSE)

write.csv(sce_merge@meta.data[,c('cellType','stage','Phase', 'S.Score', 'G2M.Score')],
		file = paste0(output_file,'metadata.csv'),
		quote = FALSE)

