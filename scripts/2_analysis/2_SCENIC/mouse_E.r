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
library(clusterProfiler)
suppressMessages(library(Canek))

library(igraph)
library(ggraph)
library(tidygraph)


setwd('~/Desktop/Research/ciona/SCENIC/results')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')
set.seed(123)
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
myPreprocessingRegulonActivity <- function(final_loom_path, seurat_object, stage){


	regulon_info <- getRegulonInfoFromLoom(final_loom_path)
	print('Calculating AUC threshold for each regulon....')
	new_AUC_threshold <- AUCell_exploreThresholds(regulon_info$regulonAUC,plotHist = FALSE, assignCells = TRUE)
	regulon_info$newRegulonAucThresholds <- getThresholdSelected(new_AUC_threshold)

	# 2. Filter out regulons based on 3 criterions
	print('Filtering out low quality regulons....')
	regulon_auc <- getAUC(regulon_info$regulonAUC)
	regulon_auc_stat <- calculateAUCInfo(regulon_auc)
	regulon_auc_stat <- regulon_auc_stat[(regulon_auc_stat$Mean>0.01) & (regulon_auc_stat$SD>0.01), ]
	regulon_auc_stat <- filterRegulonsByLowExpressedTF(regulon_info,regulon_auc_stat, threshold=10)
	regulon_auc_filtered <- regulon_auc[rownames(regulon_auc_stat),]

	# 3. Change the regulon name (replace _(+) with ( #genes))
	tf_names <- str_split(rownames(regulon_auc_filtered), pattern = '_', simplify = TRUE)[,1]
	tf_info <- getTargetGeneInfo(tf_names, regulon_info)
	tf_info$new_names <- paste0(tf_info$tf_name, ' (', as.character(tf_info$target_gene_num), 'g)')
	rownames(tf_info) <- tf_info$tf_name

	# save the files
	regulon_info$seurat <- seurat_object
	regulon_info$filtered_regulons_auc <- regulon_auc_filtered
	regulon_info$tf_info <- tf_info

	print('saving regulon info....')
	saveRDS(regulon_info, paste0('R_downstream/mouse/',stage,'/regulon_info.rds'))

}


myCalRSS <- function(regulon_auc_filtered, metadata, stage, groupby){

	#9. Calculate the RSS (regulon specific score) for each cell type.
	rss <- calcRSS(AUC=regulon_auc_filtered, cellAnnotation=metadata[, groupby])

	rss_list <- list()
	for (x in colnames(rss)){
		rss_list[[x]] <- names(sort(rss[,x],decreasing = TRUE)[1:20])
	}
	rss_df <- as.data.frame(rss_list)

	write.table(rss_df,paste0('R_downstream/mouse/',stage,'/rss.txt'), sep = '\t', quote = FALSE,row.names = FALSE)
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

	regulons_df[which(regulons_df$to %in% unique(regulons_df$from)),]$line_width <- 2
	regulons_df[which(regulons_df$to %in% unique(regulons_df$from)),]$color <- 'black'
	regulons_df[which(regulons_df$to %in% unique(regulons_df$from)),]$alpha <- 1


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

#--------------------------------------------------regulon preprocessing-----------------------------------------------

#-------------E11------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E11 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/E11_filtered_lowcells.rds')

DefaultAssay(E11) <- 'Canek'
E11 <- RunUMAP(E11, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Early_RPC_1' = 'Early RPC 1', 'Early_RPC_2' = 'Early RPC 2', 'Early_RPC_3' = 'Early RPC 3',
					'Early_RPC_4' = 'Early RPC 4', 'Early_RPC_5' = 'Early RPC 5', 'Early_RPC_6' = 'Early RPC 6')
E11$cellType <- str_replace_all(E11$newCellType, replace_str)

#-------------E12------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E12 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/E12_filtered_lowcells.rds')

E12 <- subset(E12, newCellType %in% c('Early_RPC_1', 'Early_RPC_2', 'Early_RPC_3'))
DefaultAssay(E12) <- 'Canek'
E12 <- RunUMAP(E12, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Early_RPC_1' = 'Early RPC 1', 'Early_RPC_2' = 'Early RPC 2', 'Early_RPC_3' = 'Early RPC 3')
E12$cellType <- str_replace_all(E12$newCellType, replace_str)

#-------------E14
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E14 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/E14_filtered_lowcells.rds')

E14 <- subset(E14, cellType != 'Rods')
DefaultAssay(E14) <- 'Canek'
E14 <- RunUMAP(E14, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Neurogenic Cells' = 'Neurog.', 'Retinal Ganglion Cells' = 'RGC', 'Early RPCs' = 'Early RPC', 'Amacrine Cells' = 'AC',
				 'Photoreceptor Precursors' = 'PR Precur.', 'Late RPCs' = 'Late RPC', 'Horizontal Cells' = 'HC', 'Cones' = 'Cone')
E14$cellType <- str_replace_all(E14$cellType, replace_str)

#-------------E16
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E16 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/E16_filtered_lowcells.rds')
E16 <- subset(E16, cellType %in% c('Early RPCs', 'Horizontal Cells','Late RPCs', 'Neurogenic Cells', 'Photoreceptor Precursors', 'Retinal Ganglion Cells'))
DefaultAssay(E16) <- 'Canek'
E16 <- RunUMAP(E16, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Neurogenic Cells' = 'Neurog.', 'Retinal Ganglion Cells' = 'RGC', 'Early RPCs' = 'Early RPC',
				 'Photoreceptor Precursors' = 'PR Precur.', 'Late RPCs' = 'Late RPC', 'Horizontal Cells' = 'HC')
E16$cellType <- str_replace_all(E16$cellType, replace_str)

#-------------E18
# 1. Read the final loom file of regulon and the seurat files(for annotation)
E18 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/E18_filtered_lowcells.rds')

DefaultAssay(E18) <- 'Canek'
E18 <- RunUMAP(E18, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Cones' = 'Cone', 'Neurogenic Cells' = 'Neurog.', 'Rods' = 'Rod', 'Retinal Ganglion Cells' = 'RGC',
				 'Amacrine Cells' = 'AC', 'Photoreceptor Precursors' = 'PR Precur.', 'Late RPCs' = 'Late RPC')
E18$cellType <- str_replace_all(E18$cellType, replace_str)

#-------------P0
# 1. Read the final loom file of regulon and the seurat files(for annotation)
P0 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/P0_filtered_lowcells.rds')
DefaultAssay(P0) <- 'Canek'
P0 <- RunUMAP(P0, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Cones' = 'Cone', 'Neurogenic Cells' = 'Neurog.', 'Rods' = 'Rod',
				 'Amacrine Cells' = 'AC', 'Photoreceptor Precursors' = 'PR Precur.', 'Late RPCs' = 'Late RPC')
P0$cellType <- str_replace_all(P0$cellType, replace_str)

#-------------P2
# 1. Read the final loom file of regulon and the seurat files(for annotation)
P2 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/P2_filtered_lowcells.rds')
DefaultAssay(P2) <- 'Canek'
P2 <- RunUMAP(P2, dims = 1:20, reduction = 'pca')
replace_str <- c( 'Cones' = 'Cone', 'Rods' = 'Rod', 'Neurogenic Cells' = 'Neurog.',
				 'Amacrine Cells' = 'AC', 'Photoreceptor Precursors' = 'PR Precur.', 'Late RPCs' = 'Late RPC')
P2$cellType <- str_replace_all(P2$cellType, replace_str)

#-------------P5
# 1. Read the final loom file of regulon and the seurat files(for annotation)
P5 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/P5_filtered_lowcells.rds')
DefaultAssay(P5) <- 'Canek'
P5 <- RunUMAP(P5, dims = 1:20, reduction = 'pca')
replace_str <- c('Bipolar Cells' = 'BC', 'Cones' = 'Cone', 'Muller Glia' = 'MG', 'Rods' = 'Rod',
				 'Amacrine Cells' = 'AC', 'Photoreceptor Precursors' = 'PR Precur.', 'Late RPCs' = 'Late RPC')
P5$cellType <- str_replace_all(P5$cellType, replace_str)

#-------------P8
# 1. Read the final loom file of regulon and the seurat files(for annotation)
P8 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/P8_filtered_lowcells.rds')

P8<- subset(P8,cellType %in% c('Bipolar Cells', 'Rods','Cones', 'Muller Glia'))
DefaultAssay(P8) <- 'Canek'
P8 <- RunUMAP(P8, dims = 1:20, reduction = 'pca')
replace_str <- c('Bipolar Cells' = 'BC', 'Cones' = 'Cone', 'Muller Glia' = 'MG', 'Rods' = 'Rod')
P8$cellType <- str_replace_all(P8$cellType, replace_str)

#-------------P14
# 1. Read the final loom file of regulon and the seurat files(for annotation)
P14 <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/P14_filtered_lowcells.rds')

P14<- subset(P14,newCellType %in% c('Rod_PR_1', 'Muller_glia_1','Biopolar_cell_1','Cone_PR_1'))
P14 <- subset(P14,cellType != 'Early RPCs')
DefaultAssay(P14) <- 'Canek'
P14 <- RunUMAP(P14, dims = 1:20, reduction = 'pca')

replace_str <- c('Bipolar Cells' = 'BC', 'Cones' = 'Cone', 'Muller Glia' = 'MG', 'Rods' = 'Rod')
P14$cellType <- str_replace_all(P14$cellType, replace_str)

#----------------------------------------------------------------------------------------------------

stage_names <- c('E11', 'E12', 'E14', 'E16', 'E18', 'P0', 'P2', 'P5', 'P8', 'P14')

seurat_list <- list(E11 = E11, E12 = E12, E14 = E14, E16 = E16,
					E18 = E18, P0 = P0, P2 = P2, P5 = P5, P8 = P8, P14 = P14)

for (x in 1:10){

	print('-----------------------------------------------------')
	print(paste0('Processing ', stage_names[x],'...'))
	loom_file = paste0('FromShirokane/Final/mouse', stage_names[x], '_auc_final.loom')
	myPreprocessingRegulonActivity(loom_file,seurat_list[[x]], stage_names[x])
	print('-----------------------------------------------------')

}


# stastical analysis the number of regulons.
E11_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[1],'/regulon_info.rds'))
E12_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[2],'/regulon_info.rds'))
E14_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[3],'/regulon_info.rds'))
E16_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[4],'/regulon_info.rds'))
E18_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[5],'/regulon_info.rds'))
P0_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[6],'/regulon_info.rds'))
P2_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[7],'/regulon_info.rds'))
P5_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[8],'/regulon_info.rds'))
P8_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[9],'/regulon_info.rds'))
P14_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[10],'/regulon_info.rds'))

regulon_num_by_stage <- data.frame(stage = stage_names, 
								   regulon = c(nrow(E11_regulons_info$tf_info),
								   			   nrow(E12_regulons_info$tf_info),
								   			   nrow(E14_regulons_info$tf_info),
								   			   nrow(E16_regulons_info$tf_info),
								   			   nrow(E18_regulons_info$tf_info),
								   			   nrow(P0_regulons_info$tf_info),
								   			   nrow(P2_regulons_info$tf_info),
								   			   nrow(P5_regulons_info$tf_info),
								   			   nrow(P8_regulons_info$tf_info),
								   			   nrow(P14_regulons_info$tf_info)))

write.xlsx(regulon_num_by_stage,file = 'R_downstream/mouse/regulon_summary.xlsx',sheetName = 'summary',row.names = FALSE)
write.xlsx(regulon_num_by_stage,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_summary',row.names = FALSE)
E11_excel <- E11_regulons_info$tf_info[,-c(2,5)]
E11_excel <- E11_excel[,c(1,3,2)]
write.xlsx(E11_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_E11',append=TRUE, row.names = FALSE)

E12_excel <- E12_regulons_info$tf_info[,-c(2,5)]
E12_excel <- E12_excel[,c(1,3,2)]
write.xlsx(E12_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_E12',append=TRUE, row.names = FALSE)

E14_excel <- E14_regulons_info$tf_info[,-c(2,5)]
E14_excel <- E14_excel[,c(1,3,2)]
write.xlsx(E14_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_E14',append=TRUE, row.names = FALSE)

E16_excel <- E16_regulons_info$tf_info[,-c(2,5)]
E16_excel <- E16_excel[,c(1,3,2)]
write.xlsx(E16_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_E16',append=TRUE, row.names = FALSE)

E18_excel <- E18_regulons_info$tf_info[,-c(2,5)]
E18_excel <- E18_excel[,c(1,3,2)]
write.xlsx(E18_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_E18',append=TRUE, row.names = FALSE)

P0_excel <- P0_regulons_info$tf_info[,-c(2,5)]
P0_excel <- P0_excel[,c(1,3,2)]
write.xlsx(P0_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_P0',append=TRUE, row.names = FALSE)

P2_excel <- P2_regulons_info$tf_info[,-c(2,5)]
P2_excel <- P2_excel[,c(1,3,2)]
write.xlsx(P2_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_P2',append=TRUE, row.names = FALSE)

P5_excel <- P5_regulons_info$tf_info[,-c(2,5)]
P5_excel <- P5_excel[,c(1,3,2)]
write.xlsx(P5_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_P5',append=TRUE, row.names = FALSE)

P8_excel <- P8_regulons_info$tf_info[,-c(2,5)]
P8_excel <- P8_excel[,c(1,3,2)]
write.xlsx(P8_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_P8',append=TRUE, row.names = FALSE)

P14_excel <- P14_regulons_info$tf_info[,-c(2,5)]
P14_excel <- P14_excel[,c(1,3,2)]
write.xlsx(P14_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'mouse_P14',append=TRUE, row.names = FALSE)

#--------------------------------------------------E11-----------------------------------------------

E11_regulons_info$seurat <- subset(E11_regulons_info$seurat, SCT_snn_res.0.05 == '0')
DefaultAssay(E11_regulons_info$seurat) <- 'RNA'
E11_regulons_info$seurat <- E11_regulons_info$seurat %>%
		 SCTransform(variable.features.n = 1000, vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

p1 <- myPlotUMAP_clusters(E11_regulons_info$seura,groups = 'SCT_snn_res.0.2')
ggsave('R_downstream/mouse/E11/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)



E11_metadata <- E11_regulons_info$seurat@meta.data
E11_regulon_auc <- E11_regulons_info$filtered_regulons_auc
E11_regulon_auc <- E11_regulon_auc[,rownames(E11_metadata)]
rownames(E11_regulon_auc) <- E11_regulons_info$tf_info$new_names
top20_regulons <- names(sort(rowSums(E11_regulon_auc),decreasing = TRUE)[1:20])
top20_regulons_names <- str_split(top20_regulons, pattern=' ', simplify=TRUE)[,1]

for (x in 1:length(top20_regulons)){
	p1 <- plotBinarizedRegulonActivityAndTFExpression(E11_regulons_info$seurat,E11_regulon_auc,top20_regulons_names[x],top20_regulons[x])
	ggsave(paste0('R_downstream/mouse/E11/rss_info/',top20_regulons_names[x],'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
}

top20_regulons <- top20_regulons[-c(9,16)]
regulons_for_plot_df <- data.frame(Early.RPC = top20_regulons)
write.table(regulons_for_plot_df, 'R_downstream/mouse/E11/rss.txt', col.names = TRUE,sep = '\t', quote = FALSE, row.names = FALSE)


E11_filtered_auc <-  E11_regulon_auc[top20_regulons,]
E11_filtered_auc_scaled <- t(scale(t(E11_filtered_auc), center = T, scale=T))


pdf(file = paste0('R_downstream/mouse/E11/RegulonHeatmap.pdf'),width = 3,height = 5)
cys_heatmap <- 	Heatmap(E11_filtered_auc_scaled, 
                        name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        column_title = NULL,
                        heatmap_legend_param = list(title = 'Scaled Regulon Activity', direction = 'horizontal',title_position = 'lefttop'))

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()
#5X3
saveRDS(E11_regulons_info,'R_downstream/mouse/E11/regulon_info.rds')

E11_rpc_edge <- getTFNetworkDF(top20_regulons,E11_regulons_info, sample = 2000)
E11_rpc_node <- createNetworkNode(E11_rpc_edge)
pdf(file = paste0('R_downstream/mouse/E11/rpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E11_rpc_edge,E11_rpc_node)
dev.off()

#--------------------------------------------------E12-----------------------------------------------

E12_regulons_info$seurat <- subset(E12_regulons_info$seurat, SCT_snn_res.0.05 == '0')
DefaultAssay(E12_regulons_info$seurat) <- 'RNA'
E12_regulons_info$seurat <- E12_regulons_info$seurat %>%
		 SCTransform(variable.features.n = 1000, vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

E12_metadata <- E12_regulons_info$seurat@meta.data
E12_regulon_auc <- E12_regulons_info$filtered_regulons_auc
E12_regulon_auc <- E12_regulon_auc[,rownames(E12_metadata)]
rownames(E12_regulon_auc) <- E12_regulons_info$tf_info$new_names
E12_top20_regulons <- names(sort(rowSums(E12_regulon_auc),decreasing = TRUE)[1:20])
E12_top20_regulons_names <- str_split(E12_top20_regulons, pattern=' ', simplify=TRUE)[,1]

for (x in 1:length(E12_top20_regulons)){
	p1 <- plotBinarizedRegulonActivityAndTFExpression(E12_regulons_info$seurat,E12_regulon_auc,E12_top20_regulons_names[x],E12_top20_regulons[x])
	ggsave(paste0('R_downstream/mouse/E12/rss_info/',E12_top20_regulons_names[x],'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
}


E12_top20_regulons <- E12_top20_regulons[-9]
regulons_for_plot_df <- data.frame(Early.RPC = E12_top20_regulons)
write.table(regulons_for_plot_df, 'R_downstream/mouse/E12/rss.txt', col.names = TRUE,sep = '\t', quote = FALSE, row.names = FALSE)


E12_rpc_edge <- getTFNetworkDF(top20_regulons,E12_regulons_info, sample = 2000)
E12_rpc_node <- createNetworkNode(E12_rpc_edge)
pdf(file = paste0('R_downstream/mouse/E12/rpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E12_rpc_edge,E12_rpc_node)
dev.off()


E12_filtered_auc <-  E12_regulon_auc[E12_top20_regulons,]
E12_filtered_auc_scaled <- t(scale(t(E12_filtered_auc), center = T, scale=T))

pdf(file = paste0('R_downstream/mouse/E12/RegulonHeatmap.pdf'),width = 3,height = 5)
cys_heatmap <- 	Heatmap(E12_filtered_auc_scaled, 
                        name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        column_title = NULL,
                        heatmap_legend_param = list(title = 'Scaled Regulon Activity', direction = 'horizontal',title_position = 'lefttop'))

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()
#5X3
saveRDS(E12_regulons_info,'R_downstream/mouse/E12/regulon_info.rds')

#--------------------------------------------------E14-----------------------------------------------
E14_regulons_info$seurat <- subset(E14_regulons_info$seurat, cellType %in% c('Early RPC', 'Neurog.', 'RGC', 'HC','PR Precur.','Cone'))

DefaultAssay(E14_regulons_info$seurat) <- 'RNA'

E14_regulons_info$seurat <- E14_regulons_info$seurat %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

p1 <- myPlotUMAP_clusters(E14_regulons_info$seura,groups = 'cellType')
ggsave('R_downstream/mouse/E14/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)

E14_metadata <- E14_regulons_info$seurat@meta.data
E14_regulon_auc <- E14_regulons_info$filtered_regulons_auc
E14_regulon_auc <- E14_regulon_auc[,rownames(E14_metadata)]
rownames(E14_regulon_auc) <- E14_regulons_info$tf_info$new_names


E14_metadata$cross_species <- E14_metadata$cellType
E14_cellactivity_rpc <- calculateCellGroupsActivity(E14_regulon_auc, E14_metadata, 'cellType',scale=TRUE)
write.table(E14_cellactivity_rpc[,c('Early RPC', 'Neurog.')], 'R_downstream/mouse/E14/E14_rpc.txt', sep = '\t',quote = FALSE)

E14_metadata[E14_metadata$cross_species %in% c('Early RPC', 'Neurog.'),]$cross_species <- 'RPC'
E14_cellactivity <- calculateCellGroupsActivity(E14_regulon_auc, E14_metadata, 'cross_species',scale=TRUE)
write.table(E14_cellactivity, 'R_downstream/mouse/E14/E14_cellactivity.txt', sep = '\t',quote = FALSE)

E14_rss <- myCalRSS(E14_regulon_auc,E14_metadata, 'E14',groupby='cellType')

for (x in colnames(E14_rss)){
	if (!file.exists(paste0('R_downstream/mouse/E14/rss_info/',x))){
		dir.create((paste0('R_downstream/mouse/E14/rss_info/',x)))
	}
	for (regulons in E14_rss[,x]){

		regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
		p1 <- plotBinarizedRegulonActivityAndTFExpression(E14_regulons_info$seurat,
															 E14_regulon_auc,
															 regulon_name,
															 regulons,
															 'SCT')

		ggsave(paste0('R_downstream/mouse/E14/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
	} 

}

rss_for_plot <- c(E14_rss$Early.RPC[-c(10,16,17,19,20)],
				  E14_rss$Neurog.[c(1:4,6,7,9,18,19)],
				  E14_rss$RGC[-c(6,12,13,14,16,17,19,20)],
				  E14_rss$HC[c(1,2,5,6,12,15)],
				  E14_rss$PR.Precur.[c(3,4,7,8,17)],
				  E14_rss$Cone[c(3,4,7,19)])

rss_for_plot <- unique(rss_for_plot)

E14_metadata$barcode <- rownames(E14_metadata)
E14_metadata_sample <- E14_metadata %>% group_by(cellType) %>% slice_sample(n=100)


E14_filtered_auc <-  E14_regulon_auc[rss_for_plot,E14_metadata_sample$barcode]
E14_filtered_auc_scaled <- t(scale(t(E14_filtered_auc), center = T, scale=T))

cols <- c(plotColor[7], plotColor[1], plotColor[6], plotColor[2], plotColor[4], plotColor[3])
names(cols) <- names(table(E14_metadata$cellType))

cols <- c(cols['Early RPC'], cols['Neurog.'], cols['RGC'],cols['HC'],cols['PR Precur.'],cols['Cone'])

top_anno <- HeatmapAnnotation(
							  cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(E14_metadata_sample$cellType, levels = c('Early RPC', 'Neurog.', 'RGC', 'HC','PR Precur.','Cone'))

pdf(file = paste0('R_downstream/mouse/E14/RegulonHeatmap.pdf'),width = 4,height = 8)
cys_heatmap <- 	Heatmap(E14_filtered_auc_scaled, 
                        name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_title = NULL,
                        show_row_names = TRUE,
                        cluster_rows = FALSE,
                        top_annotation = top_anno,
                        column_split = column_splits,
                        use_raster = TRUE,
                        border_gp = gpar(col = "black", lwd = 0.3),
                        heatmap_legend_param = list(title = 'Scaled regulon activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm"),legend_gp = gpar(fontsize = 1)), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)
#direction = 'horizontal'
draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()
#figsize 6X4
saveRDS(E14_regulons_info,'R_downstream/mouse/E14/regulon_info.rds')

# only 4 networks are reconstructed
E14_erpc_edge <- getTFNetworkDF(E14_rss$Early.RPC[-c(10,16,17,19,20)],E14_regulons_info, sample = 1000)
E14_erpc_node <- createNetworkNode(E14_erpc_edge)
pdf(file = paste0('R_downstream/mouse/E14/erpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E14_erpc_edge,E14_erpc_node)
dev.off()

E14_neur_edge <- getTFNetworkDF(E14_rss$Neurog.[c(1:4,6,7,9,18,19)],E14_regulons_info, sample = 1000)
E14_neur_node <- createNetworkNode(E14_neur_edge)
pdf(file = paste0('R_downstream/mouse/E14/neur_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E14_neur_edge,E14_neur_node)
dev.off()

E14_rgc_edge <- getTFNetworkDF(E14_rss$RGC[-c(12,14,16,17,19,20)],E14_regulons_info, sample = 1000)
E14_rgc_node <- createNetworkNode(E14_rgc_edge)
pdf(file = paste0('R_downstream/mouse/E14/rgc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E14_rgc_edge,E14_rgc_node)
dev.off()

E14_hc_edge <- getTFNetworkDF(E14_rss$HC[c(1,2,5,6,12,15)],E14_regulons_info, sample = 100)
E14_hc_node <- createNetworkNode(E14_hc_edge)
pdf(file = paste0('R_downstream/mouse/E14/hc_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E14_hc_edge, E14_hc_node)
dev.off()

E14_prp_edge <- getTFNetworkDF(E14_rss$PR.Precur.[c(3,4,7,8,17)],E14_regulons_info, sample = 100)
E14_prp_node <- createNetworkNode(E14_prp_edge)
pdf(file = paste0('R_downstream/mouse/E14/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E14_prp_edge,E14_prp_node)
dev.off()


#--------------------------------------------------E16-----------------------------------------------
E16_regulons_info$seurat <- subset(E16_regulons_info$seurat, cellType != 'Late RPC')
DefaultAssay(E16_regulons_info$seurat) <- 'RNA'

E16_regulons_info$seurat <- E16_regulons_info$seurat %>%
		 SCTransform(vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

p1 <- myPlotUMAP_clusters(E16_regulons_info$seura,groups = 'cellType')
ggsave('R_downstream/mouse/E16/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)


E16_metadata <- E16_regulons_info$seurat@meta.data
E16_regulon_auc <- E16_regulons_info$filtered_regulons_auc
E16_regulon_auc <- E16_regulon_auc[,rownames(E16_metadata)]
rownames(E16_regulon_auc) <- E16_regulons_info$tf_info$new_names

E16_metadata$cross_species <- E16_metadata$cellType
E16_metadata[E16_metadata$cross_species %in% c('Early RPC', 'Neurog.'),]$cross_species <- 'RPC'
E16_cellactivity <- calculateCellGroupsActivity(E16_regulon_auc, E16_metadata, 'cross_species',scale=TRUE)
write.table(E16_cellactivity, 'R_downstream/mouse/E16/E16_cellactivity.txt', sep = '\t',quote = FALSE)

E16_rss <- myCalRSS(E16_regulon_auc,E16_metadata, 'E16',groupby='cellType')

for (x in colnames(E16_rss)){
    if (!file.exists(paste0('R_downstream/mouse/E16/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/E16/rss_info/',x)))
    }
    for (regulons in E16_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(E16_regulons_info$seurat,
                                                             E16_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/mouse/E16/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}


rss_for_plot <- c(E16_rss$Early.RPC[c(1:6,8,9)],
			      E16_rss$Neurog.[c(1,2,4,5,7,11,12)],
			      E16_rss$RGC[c(1:8,16)],
				  E16_rss$HC[c(5,10,20)],#HC
				  E16_rss$PR.Precur.[c(3,4,8,9,14)]) # PR

rss_for_plot <- unique(rss_for_plot)

E16_metadata$barcode <- rownames(E16_metadata)
E16_metadata_sample <- E16_metadata %>% group_by(cellType) %>% slice_sample(n=100)
E16_filtered_auc <-  E16_regulon_auc[rss_for_plot,E16_metadata_sample$barcode]
E16_filtered_auc_scaled <- t(scale(t(E16_filtered_auc), center = T, scale=T))

cols <- c(plotColor[1], plotColor[2], plotColor[3], plotColor[6], plotColor[4])
names(cols) <- c('Early RPC', 'Neurog.', 'RGC', 'HC', 'PR Precurs.')

top_anno <- HeatmapAnnotation(  cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(E16_metadata_sample$cellType, levels = c('Early RPC', 'Neurog.', 'RGC', 'HC','PR Precur.'))

pdf(file = paste0('R_downstream/mouse/E16/RegulonHeatmap.pdf'),width = 2,height = 5)
cys_heatmap <- 	Heatmap(E16_filtered_auc_scaled, 
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


# reconstruction of network
E16_erpc_edge <- getTFNetworkDF(E16_rss$Early.RPC[c(1:6,8,9)],E16_regulons_info, sample = 1000)
E16_erpc_node <- createNetworkNode(E16_erpc_edge)
pdf(file = paste0('R_downstream/mouse/E16/erpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E16_erpc_edge,E16_erpc_node)
dev.off()

E16_neur_edge <- getTFNetworkDF(E16_rss$Neurog.[c(1,2,4,5,7,11,12)],E16_regulons_info, sample = 1000)
E16_neur_node <- createNetworkNode(E16_neur_edge)
pdf(file = paste0('R_downstream/mouse/E16/neur_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E16_neur_edge,E16_neur_node)
dev.off()

E16_rgc_edge <- getTFNetworkDF(E16_rss$RGC[c(1:8,16)],E16_regulons_info, sample = 1000)
E16_rgc_node <- createNetworkNode(E16_rgc_edge)
pdf(file = paste0('R_downstream/mouse/E16/rgc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E16_rgc_edge,E16_rgc_node)
dev.off()

E16_hc_edge <- getTFNetworkDF(E16_rss$HC[c(5,10,20)],E16_regulons_info, sample = 100)
E16_hc_node <- createNetworkNode(E16_hc_edge)
pdf(file = paste0('R_downstream/mouse/E16/hc_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(E16_hc_edge,E16_hc_node,no_red = TRUE) 
dev.off()

E16_prp_edge <- getTFNetworkDF(E16_rss$PR.Precur.[c(3,4,8,9,14)],E16_regulons_info, sample = 1000)
E16_prp_node <- createNetworkNode(E16_prp_edge)
pdf(file = paste0('R_downstream/mouse/E16/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E16_prp_edge,E16_prp_node)
dev.off()

saveRDS(E16_regulons_info,'R_downstream/mouse/E16/regulon_info.rds')

#--------------------------------------------------E18-----------------------------------------------

DefaultAssay(E18_regulons_info$seurat) <- 'RNA'

E18_regulons_info$seurat <- E18_regulons_info$seurat %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

E18_regulons_info$seurat <- RunUMAP(E18_regulons_info$seurat, 
									n.neighbors = 80, 
									min.dist = 1.2, 
									dims = 1:20,
									reduction = 'canek', 
									verbose = FALSE)
# check the subtype of AC
#------
E18_AC <- subset(E18_regulons_info$seurat, cellType == 'AC')
E18_AC <- E18_AC %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(E18_AC) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    E18_AC <- FindClusters(E18_AC,resolution = i, verbose= FALSE)
}

# check the subtype of AC by the following markers
# GABA: Meis2, Gad1
# Glycinergic: Tcf4,GlyT1(Slc6a9)

GABA_AC <- colnames(subset(E18_AC, RNA_snn_res.0.05 %in% c('0', '2')))
Glycine_AC <- colnames(subset(E18_AC, RNA_snn_res.0.05 == '1'))


# reassign cell types
E18_regulons_info$seurat$newCellType <- E18_regulons_info$seurat$cellType
E18_regulons_info$seurat$newCellType[GABA_AC] <- 'GABA AC'
E18_regulons_info$seurat$newCellType[Glycine_AC] <- 'Glycine AC'

#------
p1 <- myPlotUMAP_clusters(E18_regulons_info$seura,groups = 'newCellType')
ggsave('R_downstream/mouse/E18/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)

E18_metadata <- E18_regulons_info$seurat@meta.data
E18_regulon_auc <- E18_regulons_info$filtered_regulons_auc
E18_regulon_auc <- E18_regulon_auc[,rownames(E18_metadata)]
rownames(E18_regulon_auc) <- E18_regulons_info$tf_info$new_names

E18_cellactivity_rpc <- calculateCellGroupsActivity(E18_regulon_auc, E18_metadata, 'cellType',scale=TRUE)
write.table(E18_cellactivity_rpc[,c('Late RPC', 'Neurog.')], 'R_downstream/mouse/E18/E18_rpc.txt', sep = '\t',quote = FALSE)

E18_metadata$cross_species <- E18_metadata$cellType
E18_metadata[E18_metadata$cross_species %in% c('Late RPC', 'Neurog.'),]$cross_species <- 'RPC'
E18_cellactivity <- calculateCellGroupsActivity(E18_regulon_auc, E18_metadata, 'cross_species', scale=TRUE)
write.table(E18_cellactivity, 'R_downstream/mouse/E18/E18_cellactivity.txt', sep = '\t',quote = FALSE)

E18_rss <- myCalRSS(E18_regulon_auc,E18_metadata, 'E18',groupby='newCellType')

for (x in colnames(E18_rss)){
    if (!file.exists(paste0('R_downstream/mouse/E18/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/E18/rss_info/',x)))
    }
    for (regulons in E18_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(E18_regulons_info$seurat,
                                                             E18_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/mouse/E18/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}

rss_for_plot <- c(E18_rss$Late.RPC[c(1:8,10:14,16,18,20)],
									E18_rss$Neurog.[c(1:3,5:9,12,16,17)],
									E18_rss$GABA.AC[c(1:5,9,11,14,20)],
				  				E18_rss$Glycine.AC[c(1:6,8:11,13:15,17)],#AC
				  				E18_rss$RGC[c(1,3:5,8,14,16,20)],
				  				E18_rss$PR.Precur.[c(1,2,4,6,7,13,14,16)],
				  				E18_rss$Cone[c(2,4:7,12,18)],#Cone
				  				E18_rss$Rod[c(1,2,11)]) # Rod

rss_for_plot <- unique(rss_for_plot)

E18_metadata$barcode <- rownames(E18_metadata)
E18_metadata_sample <- E18_metadata %>% group_by(newCellType) %>% slice_sample(n=100)
E18_filtered_auc <-  E18_regulon_auc[rss_for_plot,E18_metadata_sample$barcode]
E18_filtered_auc_scaled <- t(scale(t(E18_filtered_auc), center = T, scale=T))

cols <- c(plotColor[8], plotColor[2], plotColor[5],'#8a65d6', plotColor[3], plotColor[4],plotColor[7], plotColor[9])
names(cols) <- c('Late RPC', 'Neurog.', 'GABA AC','Glycine AC', 'RGC', 'PR Precur.', 'Cone', 'Rod')
top_anno <- HeatmapAnnotation(  cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

row_loc <- c(3,15:17,21,23,29,30, 37:40,42,45, 48,50,52:54)
row_label <- rownames(E18_filtered_auc_scaled)[row_loc]

row_anno <- rowAnnotation(
												foo = anno_mark(at = row_loc,
												 labels = row_label,
												 labels_gp = gpar(fontsize = 8)))



column_splits <- factor(E18_metadata_sample$newCellType, levels =c('Late RPC', 'Neurog.', 'GABA AC','Glycine AC', 'RGC', 'PR Precur.', 'Cone', 'Rod'))

pdf(file = paste0('R_downstream/mouse/E18/RegulonHeatmap.pdf'),,width = 2.6,height = 4)
cys_heatmap <- 	Heatmap(E18_filtered_auc_scaled, 
                        name="Regulon activity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_title = NULL,
                        show_row_names = FALSE,
                        top_annotation = top_anno,
                        right_annotation = row_anno,
                        column_split = column_splits,
                        use_raster = TRUE, cluster_rows = FALSE,
                        border_gp = gpar(col = "black", lwd = 0.3),
                        heatmap_legend_param = list(title = 'Scaled Regulon Activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()

# network
E18_gabaac_edge <- getTFNetworkDF(E18_rss$GABA.AC[c(1:5,9,11,14,20)],E18_regulons_info, sample = 1000)
E18_gabaac_node <- createNetworkNode(E18_gabaac_edge)
pdf(file = paste0('R_downstream/mouse/E18/gabaac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_gabaac_edge,E18_gabaac_node)
dev.off()

E18_glyac_edge <- getTFNetworkDF(E18_rss$Glycine.AC[c(1:6,8,10,11,13:15)],E18_regulons_info, sample = 500)
E18_glyac_node <- createNetworkNode(E18_glyac_edge)
pdf(file = paste0('R_downstream/mouse/E18/glyac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_glyac_edge,E18_glyac_node)
dev.off()

E18_cone_edge <- getTFNetworkDF(E18_rss$Cone[c(2,4:7,18)],E18_regulons_info, sample = 500)
E18_cone_node <- createNetworkNode(E18_cone_edge)
pdf(file = paste0('R_downstream/mouse/E18/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_cone_edge,E18_cone_node)
dev.off()

E18_neur_edge <- getTFNetworkDF(E18_rss$Neurog.[c(1:3,5:9,12,16,17)],E18_regulons_info, sample = 1000)
E18_neur_node <- createNetworkNode(E18_neur_edge)
pdf(file = paste0('R_downstream/mouse/E18/neur_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_neur_edge,E18_neur_node)
dev.off()

E18_prp_edge <- getTFNetworkDF(E18_rss$PR.Precur.[c(1,2,4,6,7,13,14,16)],E18_regulons_info, sample = 1000)
E18_prp_node <- createNetworkNode(E18_prp_edge)
pdf(file = paste0('R_downstream/mouse/E18/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_prp_edge,E18_prp_node)
dev.off()

E18_rgc_edge <- getTFNetworkDF(E18_rss$RGC[c(1,3:5,8,14,16,20)],E18_regulons_info, sample = 1000)
E18_rgc_node <- createNetworkNode(E18_rgc_edge)
pdf(file = paste0('R_downstream/mouse/E18/rgc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_rgc_edge,E18_rgc_node)
dev.off()

E18_rod_edge <- getTFNetworkDF(E18_rss$Rod[c(1,2,11)],E18_regulons_info, sample = 200)
E18_rod_node <- createNetworkNode(E18_rod_edge)
pdf(file = paste0('R_downstream/mouse/E18/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(E18_rod_edge,E18_rod_node)
dev.off()

saveRDS(E18_regulons_info,'R_downstream/mouse/E18/regulon_info.rds')


