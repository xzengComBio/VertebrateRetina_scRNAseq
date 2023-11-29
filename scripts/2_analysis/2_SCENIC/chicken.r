library(SCopeLoomR)
library(SCENIC)
library(Seurat)

# For some of the plots:
library(KernSmooth)
library(clustree)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(stringr)
library(readxl)
library(circlize)
library(xlsx)
library(Canek)
setwd('~/Desktop/Research/ciona/SCENIC/results')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')
library(org.Gg.eg.db)
library(clusterProfiler)

library(igraph)
library(ggraph)
library(tidygraph)
# The process of regulon analysis
# 1. Read the final loom file of regulon and the seurat files(for annotation)
# 2. Estimate the threshold for each regulon. 
# 3. Filter out regulons based on 3 criterions:
#	- Regulons with mean of regulon acitivy across all cells lower than 0.01 
#	- Regulons with standard deviation lower than 0.01
#	- Regulons with TF express in lower than 10(default) cells
# 4. Change the regulon name (replace _(+) with ( #genes))
# 	- for chicken and zebrafish, replace the ENSEMBL name with gene symbol
# 5. Binarized the regulon activity
# 6. Plot the distribution of Binarized regulon activity 
# 7. Heatmap plot of regulon X CellType
# 8. Filter out regulons with mean of regulon activty acroos all cell type lower that 0.06 (default)
# 9. Heatmap plot of filtered regulon X CellType
# 10. Calculate the RSS (regulon specific score) for each cell type.

#-------------------------------------function-----------------------------------------------

myPreprocessingRegulonActivity <- function(regulon_info, seurat_object, stage){

	# 2. Filter out regulons based on 3 criterions
	print('Filtering out low quality regulons....')
	regulon_auc <- getAUC(regulon_info$regulonAUC)
	regulon_auc_stat <- calculateAUCInfo(regulon_auc)
	#regulon_auc_stat <- regulon_auc_stat[(regulon_auc_stat$Mean>0.01) & (regulon_auc_stat$SD>0.01), ]
	regulon_auc_stat <- filterRegulonsByLowExpressedTF(regulon_info,regulon_auc_stat, threshold=10)
	regulon_auc_filtered <- regulon_auc[rownames(regulon_auc_stat),]

	# 3. Change the regulon name (replace _(+) with ( #genes))
	tf_names <- str_split(rownames(regulon_auc_filtered), pattern = '_', simplify = TRUE)[,1]
	tf_info <- getTargetGeneInfo(tf_names, regulon_info)
	tf_info$gene_symbol <- chicken_gene_SymbolENSEMBL[tf_info$tf_name,]$Gene_symbol
	tf_info$new_names <- paste0(tf_info$gene_symbol, ' (', as.character(tf_info$target_gene_num), 'g)')
	rownames(tf_info) <- tf_info$tf_name

	# save the files
	regulon_info$seurat <- seurat_object
	regulon_info$filtered_regulons_auc <- regulon_auc_filtered
	regulon_info$tf_info <- tf_info

	print('saving regulon info....')
	saveRDS(regulon_info, paste0('R_downstream/chicken/',stage,'/regulon_info.rds'))

}

myCalRSS <- function(regulon_auc_filtered, metadata, stage,groupby='cellType'){

	#9. Calculate the RSS (regulon specific score) for each cell type.
	rss <- calcRSS(AUC=regulon_auc_filtered, cellAnnotation=metadata[, groupby])

	rss_list <- list()
	for (x in colnames(rss)){
		rss_list[[x]] <- names(sort(rss[,x],decreasing = TRUE)[1:10])
	}
	rss_df <- as.data.frame(rss_list)

	write.table(rss_df,paste0('R_downstream/chicken/',stage,'/rss.txt'), sep = '\t', quote = FALSE,row.names = FALSE)
	return(rss_df)

}

getTFNetworkDF <- function(rss, regulon_info, sample = 200,is_mouse = TRUE){

	if (is_mouse == TRUE){
	regulons <- regulon_info$regulons[str_replace_all(rss,' \\([0-9]+g\\)', '_(+)')]
	regulon_tf_names <- str_split(names(regulons), pattern = '_', simplify =  TRUE)[,1]
	target_num <- lengths(regulons)
	from_column <- rep(regulon_tf_names,target_num)

	regulons_df <- data.frame(from = from_column, to = unlist(regulons))
	
	}else{

		rownames(regulon_info$tf_info) <-regulon_info$tf_info$new_name
		regulons <- regulon_info$regulons[paste0(regulon_info$tf_info[rss,]$tf_name, '_(+)')]
		regulon_tf_names <- str_split(rss, pattern = ' ', simplify = TRUE)[,1]
		target_num <- lengths(regulons)
		from_column <- rep(regulon_tf_names,target_num)
		to_column <- unlist(regulons)
		replace_column <- regulon_tf_names
		names(replace_column) <- str_split(names(regulons), pattern='_', simplify=TRUE)[,1]
		to_column <- str_replace_all(unlist(regulons),replace_column)

		regulons_df <- data.frame(from = from_column, to = to_column)
	}

	
	keep_edges <- which(regulons_df$to %in% unique(regulons_df$from))
	regulons_df <- regulons_df[unique(c(sample(1:nrow(regulons_df), sample),keep_edges)),]

	regulons_df$line_width <- 1
	regulons_df$color <- 'grey'
	regulons_df$alpha <-  0.3
	if (sum(regulons_df$to %in% unique(regulons_df$from)) > 0){
		regulons_df[which(regulons_df$to %in% unique(regulons_df$from)),]$line_width <- 2
		regulons_df[which(regulons_df$to %in% unique(regulons_df$from)),]$color <- 'black'
		regulons_df[which(regulons_df$to %in% unique(regulons_df$from)),]$alpha <- 1
	}

	return(regulons_df)

}

createNetworkNode <- function(regulon_df){

	regulon_node <- data.frame(label = unique(c(regulon_df$from,regulon_df$to)),family = 'Others', node_color = 'lightgrey',node_side_color='lightgrey')
	rownames(regulon_node) <- regulon_node$label

	regulon_node[unique(regulon_df$from),]$family <- 'TF'
	regulon_node[unique(regulon_df$from),]$node_color <- '#9370DB'
	regulon_node[unique(regulon_df$from),]$node_side_color <- '#9370DB'
	return(regulon_node)
}


createAndPlotGraph <- function(edge_df,node_df){

	net_regulon <- graph_from_data_frame(d=edge_df,
																			 vertices=node_df,
																			 directed=TRUE)

	net_regulon <- as_tbl_graph(net_regulon) %>% 
  							 mutate(deg = centrality_degree(mode='out'),
          							group=group_infomap())

  p1 <- ggraph(net_regulon,layout = 'centrality',cent=deg) + 
    geom_node_point(aes(size = deg),fill = node_df$node_color,shape = 21,color = node_df$node_side_color,show.legend=FALSE)  + 
    geom_edge_link(aes(edge_width = 0.1 *line_width,edge_color=color, edge_alpha = alpha),arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm'),show.legend=FALSE) + 
    geom_node_text(aes(filter = family == 'TF',label=name),size=3) +
    theme_graph() + 
    scale_edge_width(range=c(0.1,0.3))

  return(p1)
}

createAndPlotGraph_star <- function(edge_df,node_df, no_red = FALSE){


	center <- names(sort(table(edge_df$from),decreasing = TRUE)[1])
	node_df <- node_df[c(center,sample(rownames(node_df),nrow(node_df))),]
	node_df <- node_df[!duplicated(node_df),]
	net_regulon <- graph_from_data_frame(d=edge_df,
										 vertices=node_df,
										 directed=TRUE)

	net_regulon <- as_tbl_graph(net_regulon) %>% 
  							 mutate(deg = centrality_degree(mode='out'),
          							group=group_infomap())
  	if ( no_red == FALSE){
  		p1 <- ggraph(net_regulon,layout = 'star') + 
    		  geom_node_point(aes(size = deg),fill = node_df$node_color,shape = 21,color = node_df$node_side_color,show.legend=FALSE)  + 
              geom_edge_link(aes(edge_width = 0.1 *line_width,edge_color=color, edge_alpha = alpha),arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm'),show.legend=FALSE) + 
              geom_node_text(aes(filter = family == 'TF',label=name),size=3) +
              theme_graph() + 
              scale_edge_width(range=c(0.1,0.3))
  	}
  	else{
  		 p1 <- ggraph(net_regulon,layout = 'star') + 
    			geom_node_point(aes(size = deg),fill = node_df$node_color,shape = 21,color = node_df$node_side_color,show.legend=FALSE)  + 
   				geom_edge_link(aes(edge_width = 0.1 *line_width, edge_alpha =alpha), edge_color='lightgrey',
   								   arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm'),show.legend=FALSE) + 
   				geom_node_text(aes(filter = family == 'TF',label=name),size=3) +
   				theme_graph() + 
    			scale_edge_width(range=c(0.1,0.3))
  	}

  return(p1)
}


#-------------------------------------------------------------------------------------------------
chicken_gene_SymbolENSEMBL <- read.table('../../TF/zebrafishAndChicken/gtf/chicken/geneSymbolENSEMBL.txt', sep = '\t')
colnames(chicken_gene_SymbolENSEMBL) <- c('ENSEMBL', 'Gene_symbol')
chicken_gene_SymbolENSEMBL[chicken_gene_SymbolENSEMBL$Gene_symbol == '',]$Gene_symbol <- chicken_gene_SymbolENSEMBL[chicken_gene_SymbolENSEMBL$Gene_symbol == '',]$ENSEMBL
rownames(chicken_gene_SymbolENSEMBL) <- chicken_gene_SymbolENSEMBL$ENSEMBL

chicken <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_filtered_processed_final.rds')

#--------------------------------------------------E18-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E18 <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/E18_filtered.rds')
replace_str <- c('Amacrine Cells' = 'AC', 'Bipolar Cells' = 'BC', 'Cones' = 'Cone', 'Horizontal Cells' = 'HC',
				'Muller Glia' = 'MG', 'Rods' = 'Rod')

E18$cellType <- str_replace_all(E18$CellTypes, replace_str)

anno_path <- '/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Yamagata2021eLife/supp/SCP1159/cluster/'
AC_anno <- read.table(paste0(anno_path,'Chick_AC_clusterfile.txt'), sep = '\t', header = TRUE)
BC_anno <- read.table(paste0(anno_path,'Chick_BP_clusterfile.txt'), sep = '\t', header = TRUE)
HC_anno <- read.table(paste0(anno_path,'Chick_HC_clusterfile.txt'), sep = '\t', header = TRUE)
PR_anno <- read.table(paste0(anno_path,'Chick_PR_clusterfile.txt'), sep = '\t', header = TRUE)

E18$newCellType <- E18$cellType

E18_AC <-  subset(E18, cellType == 'AC')

# GABA AC marker : SLC6A1(ENSGALG00000004921), GAD1, GAD2
# Glycine AC : SLC6A9

AC_anno$NAME <- str_split(AC_anno$NAME,pattern = '-', simplify = TRUE)[,1]

E18_AC$subType <- 'GABA AC'
E18_AC$subType[AC_anno[AC_anno$Cluster %in% c('48','20','31','58','37','2','15','1','29','22','47','10','21','34','42','6','13'),]$NAME] <- 'Glycine AC'

E18$newCellType[colnames(subset(E18_AC, subType == 'GABA AC'))] <- 'GABA AC'
E18$newCellType[colnames(subset(E18_AC, subType == 'Glycine AC'))] <- 'Glycine AC'

# ON BC marker : ISL1
# OFF BC marker :GRIK1(ENSGALG00000015835)
# ON&OFF BC 
# Rod BC : PRKCA
BC_anno$NAME <- str_split(BC_anno$NAME, pattern = '-', simplify = TRUE)[,1]
E18$newCellType[BC_anno[BC_anno$Cluster %in% c('3','9','22','4','14','18','13','21','1'),]$NAME] <- 'ON BC'
E18$newCellType[BC_anno[BC_anno$Cluster %in% c('2','5','8','6','7','12','11','15','20', '16','17'),]$NAME] <- 'OFF BC'
E18$newCellType[BC_anno[BC_anno$Cluster == '10',]$NAME] <- 'ON&OFF BC'
E18$newCellType[BC_anno[BC_anno$Cluster == '19',]$NAME] <- 'Rod BC'

# TYPE I HC : LHX1
# TYPE II HC : ISL1 
HC_anno$NAME <- str_split(HC_anno$NAME, pattern = '-', simplify = TRUE)[,1]
E18$newCellType[HC_anno[HC_anno$Cluster %in% c('1', '3'),]$NAME] <- 'T1 HC'
E18$newCellType[HC_anno[HC_anno$Cluster %in% c('2', '4','5'),]$NAME] <- 'T2 HC'


#PR
PR_anno$NAME <- str_split(PR_anno$NAME, pattern = '-', simplify = TRUE)[,1]
E18$newCellType[PR_anno[PR_anno$Cluster %in% c('DevRods', 'Rods'),]$NAME] <- 'Rod'
E18$newCellType[PR_anno[PR_anno$Cluster %in% c('DBCones1', 'DBCones2', 'DBCones3', 'DevDB13Cones', 'DevDB2Cones'),]$NAME] <- 'Double Cone'
E18$newCellType[PR_anno[PR_anno$Cluster == 'VioletCones',]$NAME] <- 'Violet Cone'
E18$newCellType[PR_anno[PR_anno$Cluster == 'RedCones',]$NAME] <- 'Red Cone'
E18$newCellType[PR_anno[PR_anno$Cluster == 'GreenCones',]$NAME] <- 'Green Cone'
E18$newCellType[PR_anno[PR_anno$Cluster == 'BlueCones',]$NAME] <- 'Blue Cone'
E18$newCellType[PR_anno[PR_anno$Cluster == 'DevNonDBCones',]$NAME] <- 'Dev. Cone'

E18_final_loom_path <- 'FromShirokane/Final/chickenE18_auc_final.loom'
E18_regulon_info <- getRegulonInfoFromLoom(E18_final_loom_path)
new_AUC_threshold <- AUCell_exploreThresholds(E18_regulon_info$regulonAUC,plotHist = FALSE, assignCells = TRUE)
E18_regulon_info$newRegulonAucThresholds <- getThresholdSelected(new_AUC_threshold)

myPreprocessingRegulonActivity(E18_regulon_info, E18, 'E18')
#--------------------------------------------------E16-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E16 <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/E16_filtered.rds')


E16_final_loom_path <- 'FromShirokane/Final/chickenE16_auc_final.loom'
E16_regulon_info <- getRegulonInfoFromLoom(E16_final_loom_path)

new_AUC_threshold <- AUCell_exploreThresholds(E16_regulon_info$regulonAUC[-24],plotHist = FALSE, assignCells = TRUE)
new_AUC_threshold <- getThresholdSelected(new_AUC_threshold)
add_regulons <- as.numeric(names(E16_regulon_info$regulonAucThresholds[24]))
names(add_regulons) <- 'ENSGALG00000003324_(+)'
new_AUC_threshold <- c(new_AUC_threshold[1:23], add_regulons,new_AUC_threshold[24:length(new_AUC_threshold)])
E16_regulon_info$newRegulonAucThresholds <- new_AUC_threshold

myPreprocessingRegulonActivity(E16_regulon_info, E16, 'E16')

#--------------------------------------------------E12-----------------------------------------------
E12 <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/E12_filtered.rds')

chicken_metadata <- chicken_metadata[chicken_metadata$stage == 'E12',]
E12 <- subset(E12, cells = rownames(chicken_metadata))

chicken_metadata <- chicken_metadata[colnames(E12),]
E12$cellType <- chicken_metadata$cellType

E12_final_loom_path <- 'FromShirokane/Final/chickenE12_auc_final.loom'
E12_regulon_info <- getRegulonInfoFromLoom(E12_final_loom_path)
new_AUC_threshold <- AUCell_exploreThresholds(E12_regulon_info$regulonAUC,plotHist = FALSE, assignCells = TRUE)
E12_regulon_info$newRegulonAucThresholds <- getThresholdSelected(new_AUC_threshold)

E12 <- subset(E12, cells = colnames(E12_regulon_info$exprMat))

myPreprocessingRegulonActivity(E12_regulon_info, E12, 'E12')



#--------------------------------------------------------------------------------------------------

stage_names <- c('E12', 'E16', 'E18')
seurat_list <- list(E12 = E12, E16 = E16, E18 = E18)


for (x in 1:3){


	regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[x],'/regulon_info.rds'))
	myBinarizedandPlotRegulonActivity(regulon_info, stage_names[x])

	metadata <- regulon_info$seurat@meta.data
	regulon_auc_filtered <- regulon_info$filtered_regulons_auc
	regulon_auc_filtered <- regulon_auc_filtered[,rownames(metadata)]

	regulonCellTypeActivity <- calculateCellGroupsActivity(regulon_auc_filtered
														   ,metadata,
														   'cellType')
	rownames(regulonCellTypeActivity) <- regulon_info$tf_info$new_names
	rownames(regulon_auc_filtered) <- regulon_info$tf_info$gene_symbol
	myPlotCellTypeHeatMap(regulonCellTypeActivity, stage_names[x])
	rss <- myCalRSS(regulon_auc_filtered,metadata, stage_names[x])
	print('-----------------------------------------------------')

}



#--------------------------------------------------------------------------------------------------
stage_names <- c('E12', 'E16', 'E18')

E12_regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[1],'/regulon_info.rds'))
E16_regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[2],'/regulon_info.rds'))
E18_regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[3],'/regulon_info.rds'))

regulon_num_by_stage <- data.frame(stage = stage_names, 
								   regulon = c(nrow(E12_regulon_info$tf_info),
								   			   nrow(E16_regulon_info$tf_info),
								   			   nrow(E18_regulon_info$tf_info)))

write.xlsx(regulon_num_by_stage,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'chicken_summary',append=TRUE,row.names = FALSE)

E12_excel <- E12_regulon_info$tf_info[,-c(2,5)]
E12_excel <- E12_excel[,c(1,3,2)]
write.xlsx(E12_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'chicken_E12',append=TRUE, row.names = FALSE)

E16_excel <- E16_regulon_info$tf_info[,-c(2,5)]
E16_excel <- E16_excel[,c(1,3,2)]
write.xlsx(E16_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'chicken_E16',append=TRUE, row.names = FALSE)

E18_excel <- E18_regulon_info$tf_info[,-c(2,5)]
E18_excel <- E18_excel[,c(1,3,2)]
write.xlsx(E18_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'chicken_E18',append=TRUE, row.names = FALSE)

#---------------------------------------------E12----------------------------------------------

E12_regulon_info$seurat <- E12_regulon_info$seurat %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(E12_regulon_info$seurat) <- 'RNA'

resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
	E12_regulon_info$seurat <- FindClusters(E12_regulon_info$seurat,resolution = i, verbose= FALSE)
}

E12_regulon_info$seurat$newCellType <-  E12_regulon_info$seurat$cellType


# GABA AC marker : SLC6A1(ENSGALG00000004921), GAD1, GAD2
# Glycine AC : SLC6A9

Glycine_AC <- colnames(subset(E12_regulon_info$seurat, cellType == 'AC' & RNA_snn_res.0.4 == 3))
E12_regulon_info$seurat$newCellType[Glycine_AC] <- 'Glycine AC'
E12_regulon_info$seurat$newCellType[E12_regulon_info$seurat$newCellType == 'AC'] <- 'GABA AC'

p1 <- myPlotUMAP_clusters(E12_regulon_info$seura,groups = 'newCellType')
ggsave('R_downstream/chicken/E12/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)


# ON BC marker : ISL1
# OFF BC marker :GRIK1(ENSGALG00000015835)
# ON&OFF BC 
# Rod BC : PRKCA

E12_metadata <- E12_regulon_info$seurat@meta.data
E12_regulon_auc <- E12_regulon_info$filtered_regulons_auc
rownames(E12_regulon_auc) <- E12_regulon_info$tf_info$new_names
E12_regulon_auc <- E12_regulon_auc[,rownames(E12_metadata)]

E12_metadata$cross_species <- E12_metadata$cellType
#E12_metadata[E12_metadata$cross_species %in% c('PR Precur.'),]$cross_species <- 'PR Precurs.'
#E18_metadata[E18_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
E12_cellactivity <- calculateCellGroupsActivity(E12_regulon_auc, E12_metadata, 'cross_species',scale = TRUE)
write.table(E12_cellactivity, 'R_downstream/chicken/E12/E12_cellactivity.txt', sep = '\t',quote = FALSE)

E12_rss <- myCalRSS(E12_regulon_auc,E12_metadata, 'E12',groupby='cellType')

for (x in colnames(E12_rss)){
    if (!file.exists(paste0('R_downstream/chicken/E12/rss_info/',x))){
        dir.create((paste0('R_downstream/chicken/E12/rss_info/',x)))
    }
    for (regulons in E12_rss[,x]){
        regulon_names <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(E12_regulon_info$seurat,
                                                             E12_regulon_auc,
                                                             regulon_names,
                                                             regulon_names,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/chicken/E12/rss_info/',x,'/',regulon_names,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}


rss_for_plot <- c(E12_rss$AC[c(1:5,8)], # RPC X 0
				  E12_rss$RGC[c(1,4,5,7)], #RPC X1
				  E12_rss$HC[c(1,2)],
				  E12_rss$BC[c(6,8)],
				  E12_rss$PR.Precurs.[c(3)],
				  E12_rss$MG[c(1,2,4,6,7,10)]) 
rss_for_plot <- unique(rss_for_plot)


E12_metadata$barcode <- rownames(E12_metadata)
E12_metadata_sample <- E12_metadata %>% group_by(cellType) %>% slice_sample(n=100)
E12_filtered_auc <-  E12_regulon_auc[rss_for_plot,E12_metadata_sample$barcode]
E12_filtered_auc_scaled <- t(scale(t(E12_filtered_auc), center = T, scale=T))

cols <- c(plotColor[5], plotColor[3], plotColor[6], plotColor[11], plotColor[4], plotColor[10])
names(cols) <- c('AC', 'RGC', 'HC','BC', 'PR Precurs.', 'MG')

top_anno <- HeatmapAnnotation(  cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(E12_metadata_sample$cellType, levels =c('AC', 'RGC', 'HC','BC', 'PR Precurs.', 'MG'))

pdf(file = paste0('R_downstream/chicken/E12/RegulonHeatmap.pdf'),width =3,height = 5)
cys_heatmap <- 	Heatmap(E12_filtered_auc_scaled, 
                        name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_title = NULL,
                        top_annotation = top_anno,
                        column_split = column_splits,
                        use_raster = TRUE, cluster_rows = FALSE,
                        border_gp = gpar(col = "black", lwd = 0.3),
                        heatmap_legend_param = list(title = 'Scaled Regulon Activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()
saveRDS(E12_regulon_info,'R_downstream/chicken/E12/regulon_info.rds')

E12_ac_edge <- getTFNetworkDF(E12_rss$AC[c(1:5,8)],E12_regulon_info, sample = 500, is_mouse = FALSE)
E12_ac_node <- createNetworkNode(E12_ac_edge)
pdf(file = paste0('R_downstream/chicken/E12/ac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E12_ac_edge,E12_ac_node)
dev.off()

E12_rgc_edge <- getTFNetworkDF(E12_rss$RGC[c(1,4,5,7)],E12_regulon_info, sample = 1000, is_mouse = FALSE)
E12_rgc_node <- createNetworkNode(E12_rgc_edge)
pdf(file = paste0('R_downstream/chicken/E12/rgc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E12_rgc_edge,E12_rgc_node)
dev.off()




# for (x in colnames(E12_rss)){

# 	E12_tf_info <- E12_regulon_info$tf_info[E12_regulon_info$tf_info$new_names %in% E12_rss[,x],]
#     E12_tf_target_list <- E12_regulon_info$regulons[paste0(E12_tf_info$tf_name, '_(+)')]

#     names(E12_tf_target_list) <- E12_tf_info$gene_symbol

#     E12_go <- compareCluster(E12_tf_target_list,
#                              fun=enrichGO,
#                              OrgDb = 'org.Gg.eg.db',
#                              keyType = 'ENSEMBL',
#                              ont = 'BP', 
#                              pAdjustMethod = 'BH',
#                              qvalueCutoff = 0.05)

#     saveRDS(E12_go,paste0('R_downstream/chicken/E12/',x,'_go.rds'))

#     pdf(file = paste0('R_downstream/chicken/E12/',x,'_go.pdf'),
#         width = 12,
#         height = 15)

#     dotplot(E12_go)
#     dev.off()
# }


#---------------------------------------------E16----------------------------------------------

E16_metadata <- E16_regulon_info$seurat@meta.data
E16_regulon_auc <- E16_regulon_info$filtered_regulons_auc
rownames(E16_regulon_auc) <- E16_regulon_info$tf_info$new_names
E16_regulon_auc <- E16_regulon_auc[,rownames(E16_metadata)]

#myBinarizedandPlotRegulonActivity(E16_regulon_info, 'E16')

p1 <- myPlotUMAP(E16_regulon_info$seurat)
ggsave('R_downstream/chicken/E16/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)

top10_regulons <- names(sort(rowSums(E16_regulon_auc),decreasing = TRUE)[1:10])
regulons_for_plot <- c(top10_regulons, 'ENSGALG00000027907 (36g)','EOMES (8g)',
						'SOX6 (22g)', 'SOX9 (54g)', 'TFAP2D (1445g)')

for (regulons in regulons_for_plot){
    regulon_names <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
    p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(E16_regulon_info$seurat,
                                                             E16_regulon_auc,
                                                             regulon_names,
                                                             regulon_names,
                                                             regulons,
                                                             'SCT')

    ggsave(paste0('R_downstream/chicken/E16/rss_info/',regulon_names,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
} 


regulons_for_plot <- unique(regulons_for_plot)
regulons_for_plot_df <- data.frame(RGC = regulons_for_plot)
write.table(regulons_for_plot_df, 'R_downstream/chicken/E16/rss.txt', col.names = TRUE,sep = '\t', quote = FALSE, row.names = FALSE)


E16_filtered_auc <-  E16_regulon_auc[regulons_for_plot,]
E16_filtered_auc_scaled <- t(scale(t(E16_filtered_auc), center = T, scale=T))

pdf(file = paste0('R_downstream/chicken/E16/RegulonHeatmap.pdf'),width =5,height = 5)
cys_heatmap <- 	Heatmap(E16_filtered_auc_scaled, 
                        name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        column_title = NULL,
                        heatmap_legend_param = list(title = 'Scaled Regulon Activity', direction = 'horizontal',title_position = 'lefttop'), 
                        column_gap = unit('0.1', 'mm'))

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()
saveRDS(E16_regulon_info,'R_downstream/chicken/E16/regulon_info.rds')


#---------------------------------------------E18----------------------------------------------

E18_metadata <- E18_regulon_info$seurat@meta.data
E18_regulon_auc <- E18_regulon_info$filtered_regulons_auc
rownames(E18_regulon_auc) <- E18_regulon_info$tf_info$new_names
E18_regulon_auc <- E18_regulon_auc[,rownames(E18_metadata)]


E18_metadata$cross_species <- E18_metadata$cellType
#E18_metadata[E18_metadata$cross_species %in% c('PR Precur.'),]$cross_species <- 'PR Precurs.'
#E18_metadata[E18_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
E18_cellactivity <- calculateCellGroupsActivity(E18_regulon_auc, E18_metadata, 'cross_species',scale = TRUE)
write.table(E18_cellactivity, 'R_downstream/chicken/E18/E18_cellactivity.txt', sep = '\t',quote = FALSE)



p1 <- myPlotUMAP_clusters(E18_regulon_info$seurat,groups = 'newCellType') + guides(color = guide_legend(override.aes = list(size=4), ncol=2))
ggsave('R_downstream/chicken/E18/umap.pdf', p1, dpi = 200, units = 'cm', width = 17, height = 10)

E18_rss <- myCalRSS(E18_regulon_auc,E18_metadata, 'E18',groupby='newCellType')

for (x in colnames(E18_rss)){
    if (!file.exists(paste0('R_downstream/chicken/E18/rss_info/',x))){
        dir.create((paste0('R_downstream/chicken/E18/rss_info/',x)))
    }
    for (regulons in E18_rss[,x]){
        regulon_names <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(E18_regulon_info$seurat,
                                                             E18_regulon_auc,
                                                             regulon_names,
                                                             regulon_names,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/chicken/E18/rss_info/',x,'/',regulon_names,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}


rss_for_plot <- c(
				  E18_rss$GABA.AC[c(1:4,7)],
				  E18_rss$Glycine.AC[c(1:4,9)],
				  E18_rss$Rod.BC[c(1)],
				  E18_rss$ON.BC[c(5,7)],
				  E18_rss$OFF.BC[c(1,5)],
				  E18_rss$T1.HC[c(1,2)],
				  E18_rss$T2.HC[c(1,2)],
				  E18_rss$Blue.Cone[c(1,2,5)],
				  E18_rss$Dev..Cone[c(2,3)],
				  E18_rss$Double.Cone[c(1,3,4)],
				  E18_rss$Green.Cone[c(1,3,4,5)],
				  E18_rss$Violet.Cone[c(1,3,4)],
				  E18_rss$Red.Cone[c(1,2,3)],
				  E18_rss$Rod[c(1,4,6)],
				  E18_rss$MG[c(3,4,5,6,9)]
	)

rss_for_plot <- unique(rss_for_plot)

E18_metadata$barcode <- rownames(E18_metadata)
#E18_metadata_sample <- E18_metadata %>% group_by(newCellType) %>% slice_sample(n=100)
E18_metadata_sample <- E18_metadata %>% group_by(cellType) %>% slice_sample(n=100)
E18_filtered_auc <-  E18_regulon_auc[rss_for_plot,E18_metadata_sample$barcode]
E18_filtered_auc_scaled <- t(scale(t(E18_filtered_auc), center = T, scale=T))

#cols <- c(plotColor[5],'#8a65d6', plotColor[11],'#80800d','#999900','#7a7a07', plotColor[6],'#87e087',
#	     '#d97373','Red','Green','Blue','Violet',plotColor[7],  plotColor[9],plotColor[10])

cols <- plotColor[c(5,11,6,7,9,10)]

#names(cols) <- c('GABA AC','Glycine AC','Rod BC','OFF BC','ON BC','ON&OFF BC','T1 HC','T2 HC', 
#	'Dev. Cone', 'Red Cone', 'Green Cone', 'Blue Cone', 'Violet Cone', 'Double Cone', 'Rod','MG')

names(cols) <- c('AC', 'BC', 'HC', 'Cone', 'Rod', 'MG')
top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

#column_splits <- factor(E18_metadata_sample$newCellType, levels = c('GABA AC','Glycine AC','Rod BC','OFF BC','ON BC','ON&OFF BC','T1 HC','T2 HC',
#	   'Dev. Cone', 'Red Cone', 'Green Cone', 'Blue Cone', 'Violet Cone', 'Double Cone', 'Rod','MG'))
column_splits <- factor(E18_metadata_sample$cellType, levels = c('AC', 'BC', 'HC', 'Cone', 'Rod','MG'))


row_loc <- c(2,3,7:11, 13,14,18,20:23)
row_label <- rownames(E18_filtered_auc_scaled)[row_loc]

row_anno <- rowAnnotation(
												foo = anno_mark(at = row_loc,
												 labels = row_label,
												 labels_gp = gpar(fontsize = 8)))

pdf(file = paste0('R_downstream/chicken/E18/RegulonHeatmap.pdf'),width =3,height = 3)

cys_heatmap <- 	Heatmap(E18_filtered_auc_scaled, 
                        name="Regulon rgctivity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_row_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_title = NULL,
                        right_annotation = row_anno,
                        top_annotation = top_anno,
                        column_split = column_splits,
                        use_raster = TRUE,
                        border_gp = gpar(col = "black", lwd = 0.3), cluster_rows = FALSE,
                        heatmap_legend_param = list(title = 'Scaled Regulon activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()


cone_regulons <- c(E18_rss$Blue.Cone[c(1,2,5)],
				  E18_rss$Dev..Cone[c(2,3)],
				  E18_rss$Double.Cone[c(1,3,4)],
				  E18_rss$Green.Cone[c(1,3,4,5)],
				  E18_rss$Violet.Cone[c(1,3,4)],
				  E18_rss$Red.Cone[c(1,2,3)])
cone_regulons <- unique(cone_regulons)
E18_cone_edge <- getTFNetworkDF(cone_regulons,E18_regulon_info, sample = 500, is_mouse = FALSE)
E18_cone_node <- createNetworkNode(E18_cone_edge)
pdf(file = paste0('R_downstream/chicken/E18/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E18_cone_edge,E18_cone_node, no_red = FALSE)
dev.off()


ac_regulons <- unique(c(E18_rss$GABA.AC[c(1:4,7)],E18_rss$Glycine.AC[c(1:4,9)]))

E18_ac_edge <- getTFNetworkDF(ac_regulons,E18_regulon_info, sample = 300, is_mouse = FALSE)
E18_ac_node <- createNetworkNode(E18_ac_edge)
pdf(file = paste0('R_downstream/chicken/E18/ac_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E18_ac_edge,E18_ac_node,no_red = TRUE)
dev.off()

bc_regulons <- c(E18_rss$Rod.BC[c(1)],
				  			 E18_rss$ON.BC[c(5,7)],
				  			 E18_rss$OFF.BC[c(1,5)])

bc_regulons <- unique(bc_regulons)

E18_bc_edge <- getTFNetworkDF(bc_regulons,E18_regulon_info, sample = 100, is_mouse = FALSE)
E18_bc_node <- createNetworkNode(E18_bc_edge)
pdf(file = paste0('R_downstream/chicken/E18/bc_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E18_bc_edge,E18_bc_node,no_red = TRUE)
dev.off()

E18_rgc_edge <- getTFNetworkDF(E18_rss$MG[c(3,4,5,6,9)],E18_regulon_info, sample = 1000, is_mouse = FALSE)
E18_rgc_node <- createNetworkNode(E18_rgc_edge)
pdf(file = paste0('R_downstream/chicken/E18/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_rgc_edge,E18_rgc_node)
dev.off()



E18_rgc_edge <- getTFNetworkDF(E18_rss$Rod[c(1,4,6)],E18_regulon_info, sample = 100, is_mouse = FALSE)
E18_rgc_node <- createNetworkNode(E18_rgc_edge)
pdf(file = paste0('R_downstream/chicken/E18/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E18_rgc_edge,E18_rgc_node, no_red = TRUE)
dev.off()


hc_regulons <- c(E18_rss$T1.HC[c(1,2)], E18_rss$T2.HC[c(1,2)])
hc_regulons <- unique(hc_regulons)

E18_hc_edge <- getTFNetworkDF(hc_regulons,E18_regulon_info, sample = 100, is_mouse = FALSE)
E18_hc_node <- createNetworkNode(E18_hc_edge)
pdf(file = paste0('R_downstream/chicken/E18/hc_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E18_hc_edge,E18_hc_node, no_red = TRUE)
dev.off()



saveRDS(E18_regulon_info,'R_downstream/chicken/E18/regulon_info.rds')

#---- network
E18_edge <- getTFNetworkDF(E18_rss$MG[c(3,4,5,6,9)],E18_regulon_info,sample=1000, is_mouse = FALSE)
E18_node <- createNetworkNode(E18_edge)
p1 <- createAndPlotGraph(E18_edge,E18_node)




