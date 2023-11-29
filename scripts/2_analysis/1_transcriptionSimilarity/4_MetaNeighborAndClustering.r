source("~/Github/MetaNeighbor/2017-08-28-runMN-US.R")
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(homologene)
library(ComplexHeatmap)
library(biomaRt)
suppressMessages(library(harmony))
setwd('~/Desktop/Research/ciona/singleCell/')

source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')


mouse_pseudoCell <- readRDS('Integration/result/pseudoCell/mouse_all_pseudoCelll10_2.rds')
chicken_pseudoCell<-readRDS("Integration/result/pseudoCell/chicken_all_pseudoCell10_2.rds")
zebrafish_pseudoCell <- readRDS("Integration/result/pseudoCell/zebrafish_all_pseudoCell10_5.rds")

#---- ## read the metadata

mouse_metadata <-  as.data.frame(str_split(colnames(mouse_pseudoCell), pattern = '_', simplify = TRUE))[,c(1,2)]
mouse_metadata$Sample_ID <- colnames(mouse_pseudoCell)
colnames(mouse_metadata) <- c('Study_ID', 'Celltype', 'Sample_ID')
mouse_metadata <- mouse_metadata[,c(3,1,2)]

chicken_metadata <-  as.data.frame(str_split(colnames(chicken_pseudoCell), pattern = '_', simplify = TRUE))[,c(1,2)]
chicken_metadata$Sample_ID <- colnames(chicken_pseudoCell)
colnames(chicken_metadata) <- c('Study_ID', 'Celltype', 'Sample_ID')
chicken_metadata <- chicken_metadata[,c(3,1,2)]

zebrafish_metadata <-  as.data.frame(str_split(colnames(zebrafish_pseudoCell), pattern = '_', simplify = TRUE))[,c(1,2)]
zebrafish_metadata$Sample_ID <- colnames(zebrafish_pseudoCell)
colnames(zebrafish_metadata) <- c('Study_ID', 'Celltype', 'Sample_ID')
zebrafish_metadata <- zebrafish_metadata[,c(3,1,2)]

#------ Build the orth by homo
mouse_ensembl <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', mirror = 'asia', host = 'https://dec2021.archive.ensembl.org')
zebrafish_ensembl <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl', mirror = 'asia',host = 'https://dec2021.archive.ensembl.org')
chicken_ensembl <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'ggallus_gene_ensembl', mirror = 'asia',host = 'https://dec2021.archive.ensembl.org')

z2m <- getLDS(attributes = 'ensembl_gene_id', 
				 filters = 'ensembl_gene_id', 
				 mart = zebrafish_ensembl,
				 values = rownames(zebrafish_pseudoCell), 
				 attributesL ='ensembl_gene_id', 
				 martL = mouse_ensembl, 
				 uniqueRows = TRUE)
colnames(z2m) <- c('zebrafish_ensembl', 'mouse_ensembl')

z2m_unique <- z2m[(!duplicated(z2m$zebrafish_ensembl) & !duplicated(z2m$mouse_ensembl)),]


c2m <- getLDS(attributes = 'ensembl_gene_id', 
				 filters = 'ensembl_gene_id', 
				 mart = chicken_ensembl,
				 values = rownames(chicken_pseudoCell), 
				 attributesL ='ensembl_gene_id', 
				 martL = mouse_ensembl, 
				 uniqueRows = TRUE)

colnames(c2m) <- c('chicken_ensembl', 'mouse_ensembl')
c2m_unique <- c2m[(!duplicated(c2m$chicken_ensembl) & !duplicated(c2m$mouse_ensembl)),]

orth_df <- merge(c2m_unique,z2m_unique, by='mouse_ensembl')
orth_df <- orth_df[orth_df$mouse_ensembl %in% rownames(mouse_pseudoCell),]
orth_df$orth <- paste0('orth',1:length(orth_df$chicken))

write.table(orth_df,'Integration/result/ensembl_orth.txt', sep = '\t', col.names = TRUE, row.names = FALSE)


#---- ## integrate the pseudo-count-matrix

mouse_pseudoCell_orth <- mouse_pseudoCell[orth_df$mouse_ensembl,]
rownames(mouse_pseudoCell_orth) <- orth_df$orth

chicken_pseudoCell_orth <- chicken_pseudoCell[orth_df$chicken_ensembl,]
rownames(chicken_pseudoCell_orth) <- orth_df$orth

zebrafish_pseudoCell_orth <- zebrafish_pseudoCell[orth_df$zebrafish_ensembl,]
rownames(zebrafish_pseudoCell_orth) <- orth_df$orth

all_merge_matrix <- cbind(mouse_pseudoCell_orth,chicken_pseudoCell_orth)
all_merge_matrix <- cbind(all_merge_matrix,zebrafish_pseudoCell_orth)

metadata_list <- list( mouse = mouse_metadata,
					  chicken = chicken_metadata,
					  zebrafish = zebrafish_metadata)

metadata <- Reduce(function(x,y) rbind(x,y), metadata_list, accumulate = FALSE)

all_merge_matrix <- all_merge_matrix[,metadata$Sample_ID]
metadata$Celltype <- paste0(metadata$Study_ID, ' | ', metadata$Celltype)
merge_celltypes <- unique(metadata$Celltype)

#---- ## Run Metaneighbor
merge_var.genes <- get_variable_genes(all_merge_matrix,metadata)
# run MetaNeighbor
merge_celltype.default <- run_MetaNeighbor_US(merge_var.genes, all_merge_matrix, merge_celltypes, metadata)
pdf('Integration/result/MetaNeighbor/hvg_withoutCiona_default_homo.pdf', width = 12,height = 10)
ComplexHeatmap::Heatmap(merge_celltype.default,show_column_names = FALSE,show_row_dend = FALSE,name="AUC Score")
dev.off()

merge_var_df <- data.frame(hvg = merge_var.genes) 
write.table(merge_var_df, 'Integration/result/hvg.txt', sep = '\t', quote = FALSE,row.names = FALSE)

top_hits=get_top_hits(merge_celltype.default,metadata,threshold=0.9,filename="Integration/result/MetaNeighbor/hvg_withoutCiona_default_homo.NV_SRS_1.txt")


#---- ## Cluster pesudo-cell
merge_sce <- CreateSeuratObject(all_merge_matrix,names.field = 1, project = 'Integration')
merge_sce <- NormalizeData(merge_sce, normalization.method = "LogNormalize", scale.factor = 10000)

#merge_sce <- FindVariableFeatures(merge_sce, nfeatures = 500)
merge_sce <- ScaleData(merge_sce, features = rownames(merge_sce))
merge_sce <- RunPCA(merge_sce,features = merge_var.genes,verbose=FALSE)
merge_sce <- RunHarmony(merge_sce, group.by.vars = 'orig.ident')
merge_sce <- RunUMAP(merge_sce, reduction = 'harmony', dims = 1:20)

merge_sce <- FindNeighbors(merge_sce, reduction = 'harmony', dims = 1:20, verbose = FALSE)
resolutions <- c(0.001,0.005,0.01,0.05)
for (i in resolutions){
     merge_sce <- FindClusters(merge_sce,resolution = i, verbose= FALSE)
}

merge_sce$cellType <- metadata$Celltype
merge_sce$cluster <- str_split(metadata$Celltype, pattern = ' | ', simplify=TRUE)[,3]

p1 <- DimPlot(merge_sce, reduction = 'umap', group.by = 'cluster', cols = plotColor, label = TRUE) + 
			theme_bw() +
			theme(panel.grid.major=element_blank(),
				  panel.grid.minor=element_blank(),
				  axis.title = element_text(face = "bold",size = rel(1)),
				  plot.title = element_blank(),
				  legend.margin=margin(t = 0, unit='cm'))  + 
			NoLegend()

ggsave('Integration/result/integration_11_withoutCiona_800_umap.pdf', p1, dpi=200, units = 'cm', width = 15, height = 15)


col_use <- plotColor[c(1,4,5,3,11,6)]
p1 <- DimPlot(merge_sce, reduction = 'umap', group.by = 'RNA_snn_res.0.001', cols = col_use,label =FALSE,raster = TRUE) + 
			theme_bw() +
			theme(panel.grid.major=element_blank(),
				  panel.grid.minor=element_blank(),
				  plot.title = element_blank(),
				  legend.margin=margin(t = 0, unit='cm'),
				  axis.ticks.x = element_blank(),
				  axis.ticks.y = element_blank(),
				  axis.text.x  = element_blank(),
				  axis.text.y  = element_blank(),
				  panel.border = element_blank(),
				  axis.title= element_blank())

ggsave('Integration/result/integration_11_withoutCiona_392_clusters.pdf', p1, dpi=200, units = 'cm', width = 15, height = 12)



cellType_by_cluster <- table(merge_sce$cellType, merge_sce$RNA_snn_res.0.001)
cellType_by_cluster_df <- as.data.frame(cellType_by_cluster)
cellType_by_cluster_df <- as.data.frame(pivot_wider(cellType_by_cluster_df, names_from = Var2, values_from = Freq))

rownames(cellType_by_cluster_df) <- cellType_by_cluster_df$Var1
cellType_by_cluster_df <- cellType_by_cluster_df[,-1]

cellType_by_cluster_df <- apply(cellType_by_cluster_df, MARGIN = 1, FUN = function(x){round(x/sum(x), digits = 3)})

colnames(cellType_by_cluster_df) <- str_replace(colnames(cellType_by_cluster_df), 'chicken \\| ', 'C-')
colnames(cellType_by_cluster_df) <- str_replace(colnames(cellType_by_cluster_df), 'mouse \\| ', 'M-')
colnames(cellType_by_cluster_df) <- str_replace(colnames(cellType_by_cluster_df), 'zebrafish \\| ', 'Z-')


cellType_by_cluster_df <- cellType_by_cluster_df[,c('M-ERPC', 'M-LRPC', 'M-Neurog', 'M-MG', 'C-MG', 'Z-RPC', 'Z-MG',
										  'M-PRPrecurs', 'M-Cone', 'M-Rod', 'C-PRPrecurs','C-Rod', 'C-Cone', 'Z-PRPrecurs', 'Z-Cone','Z-Rod',
										  'M-AC', 'C-AC', 'Z-AC',
										  'M-RGC', 'C-RGC', 'Z-RGC',
										  'M-BC', 'C-BC', 'Z-BC',
										  'M-HC', 'C-HC', 'Z-HC')]

pdf('Integration/result/cellTypebyClusters_407.pdf', width =7,height = 3)

cys_heatmap <- Heatmap(cellType_by_cluster_df, 
                       cluster_columns = FALSE, 
                       cluster_rows = FALSE, 
                       column_names_side = 'top', 
                       row_names_side = 'left',
                       column_names_rot = 90,
                       column_names_centered = FALSE,
                       border_gp = gpar(col = "black", lwd = 2),
                       column_names_gp = gpar( fontface = "bold"),
                       row_names_gp = gpar(fontface = "bold"),
                       rect_gp = gpar(col= "white", lwd =0.1),
                       row_names_centered = FALSE,
                       heatmap_legend_param = list(title = 'Fraction of Cells within Cluster', direction = 'horizontal',title_position = 'lefttop'))

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()

