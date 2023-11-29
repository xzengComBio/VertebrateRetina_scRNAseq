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


#
stage_names <- c('E11', 'E12', 'E14', 'E16', 'E18', 'P0', 'P2', 'P5', 'P8', 'P14')

#--------------------------------------------------P0-----------------------------------------------

P0_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[6],'/regulon_info.rds'))


DefaultAssay(P0_regulons_info$seurat) <- 'RNA'
P0_regulons_info$seurat <- P0_regulons_info$seurat %>%
		 SCTransform(variable.features.n = 1000, vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)
# check the subtype of AC
#------		 
P0_AC <- subset(P0_regulons_info$seurat, cellType == 'AC')
P0_AC <- P0_AC %>%
		 SCTransform(variable.features.n = 1000, vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    P0_AC <- FindClusters(P0_AC,resolution = i, verbose= FALSE)
}

# check the subtype of AC by the following markers
# GABA: Meis2, Gad1
# Glycinergic: Tcf4,GlyT1(Slc6a9)

GABA_AC <- colnames(subset(P0_AC, SCT_snn_res.0.05 == '0'))
Glycine_AC <- colnames(subset(P0_AC, SCT_snn_res.0.05 == '1'))

# reassign cell types
P0_regulons_info$seurat$newCellType <- P0_regulons_info$seurat$cellType
P0_regulons_info$seurat$newCellType[GABA_AC] <- 'GABA AC'
P0_regulons_info$seurat$newCellType[Glycine_AC] <- 'Glycine AC'

p1 <- myPlotUMAP_clusters(P0_regulons_info$seura,groups = 'newCellType')
ggsave('R_downstream/mouse/P0/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)
#------		

P0_metadata <- P0_regulons_info$seurat@meta.data
P0_regulon_auc <- P0_regulons_info$filtered_regulons_auc
P0_regulon_auc <- P0_regulon_auc[,rownames(P0_metadata)]
rownames(P0_regulon_auc) <- P0_regulons_info$tf_info$new_names

P0_metadata$cross_species <- P0_metadata$cellType
P0_metadata[P0_metadata$cross_species %in% c('Late RPC', 'Neurog.'),]$cross_species <- 'RPC'
P0_cellactivity <- calculateCellGroupsActivity(P0_regulon_auc, P0_metadata, 'cross_species',scale = TRUE)
write.table(P0_cellactivity, 'R_downstream/mouse/P0/P0_cellactivity.txt', sep = '\t',quote = FALSE)

P0_rss <- myCalRSS(P0_regulon_auc,P0_metadata, 'P0',groupby='newCellType')


for (x in colnames(P0_rss)){
    if (!file.exists(paste0('R_downstream/mouse/P0/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/P0/rss_info/',x)))
    }
    for (regulons in P0_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(P0_regulons_info$seurat,
                                                             P0_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/mouse/P0/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}

rss_for_plot <- c(P0_rss$Late.RPC[c(1:5,7,8,10,12,16,17,18)],
				  P0_rss$Neurog.[c(1:6,9)],
				  P0_rss$GABA.AC[c(2,3,4,11,14,16,17)], # AC
				  P0_rss$Glycine.AC[c(1:8,14)], 
				  P0_rss$PR.Precur[c(1,2,5,6,7,8,11,13,15,19)],
				  P0_rss$Cone[c(1,2,6,8,9,13,14,15,19)],
				  P0_rss$Rod[c(2,4,5,6,7,17,18)])  #RPC X2

rss_for_plot <- unique(rss_for_plot)
P0_metadata$barcode <- rownames(P0_metadata)
P0_metadata_sample <- P0_metadata %>% group_by(newCellType) %>% slice_sample(n=100)
P0_filtered_auc <-  P0_regulon_auc[rss_for_plot,P0_metadata_sample$barcode]
P0_filtered_auc_scaled <- t(scale(t(P0_filtered_auc), center = T, scale=T))

cols <- c(plotColor[8],plotColor[2],plotColor[5],'#8a65d6', plotColor[4], plotColor[7], plotColor[9])
names(cols) <- c('Late RPC', 'Neurog.', 'GABA AC','Glycine AC', 'PR Precur.', 'Cone', 'Rod')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(P0_metadata_sample$newCellType, levels = c('Late RPC', 'Neurog.', 'GABA AC','Glycine AC', 'PR Precur.', 'Cone', 'Rod'))

pdf(file = paste0('R_downstream/mouse/P0/RegulonHeatmap.pdf'),width = 2.7,height = 6)
cys_heatmap <- 	Heatmap(P0_filtered_auc_scaled, 
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


P0_lrpc_edge <- getTFNetworkDF(P0_rss$Late.RPC[c(1:5,7,8,10,12,16,17,18)],P0_regulons_info, sample = 1000)
P0_lrpc_node <- createNetworkNode(P0_lrpc_edge)
pdf(file = paste0('R_downstream/mouse/P0/lrpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P0_lrpc_edge,P0_lrpc_node)
dev.off()

P0_gabaac_edge <- getTFNetworkDF(P0_rss$GABA.AC[c(2,3,4,11,14,16,17)],P0_regulons_info, sample = 1000)
P0_gabaac_node <- createNetworkNode(P0_gabaac_edge)
pdf(file = paste0('R_downstream/mouse/P0/gabaac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P0_gabaac_edge,P0_gabaac_node)
dev.off()

P0_glyac_edge <- getTFNetworkDF(P0_rss$Late.RPC[c(1:5,7,8,10,12,16,17,18)],P0_regulons_info, sample = 1000)
P0_glyac_node <- createNetworkNode(P0_glyac_edge)
pdf(file = paste0('R_downstream/mouse/P0/glyac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P0_glyac_edge,P0_glyac_node)
dev.off()

P0_cone_edge <- getTFNetworkDF(P0_rss$Cone[c(1,2,6,8,9,13,14,15,19)],P0_regulons_info, sample = 1000)
P0_cone_node <- createNetworkNode(P0_cone_edge)
pdf(file = paste0('R_downstream/mouse/P0/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P0_cone_edge,P0_cone_node)
dev.off()

P0_neur_edge <- getTFNetworkDF(P0_rss$Neurog.[c(1:6,9)],P0_regulons_info, sample = 500)
P0_neur_node <- createNetworkNode(P0_neur_edge)
pdf(file = paste0('R_downstream/mouse/P0/neur_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P0_neur_edge,P0_neur_node)
dev.off()

P0_rod_edge <- getTFNetworkDF(P0_rss$Rod[c(2,4,5,6,7,17,18)],P0_regulons_info, sample = 500)
P0_rod_node <- createNetworkNode(P0_rod_edge)
pdf(file = paste0('R_downstream/mouse/P0/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P0_rod_edge,P0_rod_node)
dev.off()

saveRDS(P0_regulons_info,'R_downstream/mouse/P0/regulon_info.rds')

#--------------------------------------------------P2-----------------------------------------------
P2_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[7],'/regulon_info.rds'))

DefaultAssay(P2_regulons_info$seurat) <- 'RNA'
P2_regulons_info$seurat <- P2_regulons_info$seurat %>%
		 SCTransform(vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

# check the subtype of AC
#------
P2_AC <- subset(P2_regulons_info$seurat, cellType == 'AC')
P2_AC <- P2_AC %>%
		 SCTransform(vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    P2_AC <- FindClusters(P2_AC,resolution = i, verbose= FALSE)
}

# check the subtype of AC by the following markers
# GABA: Meis2, Gad1
# Glycinergic: Tcf4,GlyT1(Slc6a9)

GABA_AC <- colnames(subset(P2_AC, SCT_snn_res.0.1 == '0'))
Glycine_AC <- colnames(subset(P2_AC, SCT_snn_res.0.1 == '1'))

# reassign cell types
P2_regulons_info$seurat$newCellType <- P2_regulons_info$seurat$cellType
P2_regulons_info$seurat$newCellType[GABA_AC] <- 'GABA AC'
P2_regulons_info$seurat$newCellType[Glycine_AC] <- 'Glycine AC'

p1 <- myPlotUMAP_clusters(P2_regulons_info$seura,groups = 'newCellType')
ggsave('R_downstream/mouse/P2/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)
#------		

P2_metadata <- P2_regulons_info$seurat@meta.data
P2_regulon_auc <- P2_regulons_info$filtered_regulons_auc
P2_regulon_auc <- P2_regulon_auc[,rownames(P2_metadata)]
rownames(P2_regulon_auc) <- P2_regulons_info$tf_info$new_names

P2_rss <- myCalRSS(P2_regulon_auc,P2_metadata, 'P2',groupby='newCellType')

P2_metadata$cross_species <- P2_metadata$cellType
P2_metadata[P2_metadata$cross_species %in% c('Late RPC', 'Neurog.'),]$cross_species <- 'RPC'
P2_cellactivity <- calculateCellGroupsActivity(P2_regulon_auc, P2_metadata, 'cross_species',scale = TRUE)
write.table(P2_cellactivity, 'R_downstream/mouse/P2/P2_cellactivity.txt', sep = '\t',quote = FALSE)

for (x in colnames(P2_rss)){
    if (!file.exists(paste0('R_downstream/mouse/P2/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/P2/rss_info/',x)))
    }
    for (regulons in P2_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(P2_regulons_info$seurat,
                                                             P2_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/mouse/P2/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}


rss_for_plot <- c(P2_rss$Late.RPC[c(1:5,8,10,12:16,20)],
				  P2_rss$Neurog.[c(3:5,8,9,11,12,14,16,17)],
				  P2_rss$GABA.AC[c(2,3,4,6,14,20)],
				  P2_rss$Glycine.AC[c(1:10,13,14,19)], # AC#Cone
				  P2_rss$PR.Precur[c(1:3,8:13,16,17,20)],
				  P2_rss$Cone[c(5,7,10,16,18,20)], 
				  P2_rss$Rod[c(1,3:4,9,13,15,16,17)]) 


rss_for_plot <- unique(rss_for_plot)

P2_metadata$barcode <- rownames(P2_metadata)
P2_metadata_sample <- P2_metadata %>% group_by(newCellType) %>% slice_sample(n=100)
P2_filtered_auc <-  P2_regulon_auc[rss_for_plot,P2_metadata_sample$barcode]
P2_filtered_auc_scaled <- t(scale(t(P2_filtered_auc), center = T, scale=T))

cols <- c(plotColor[8],plotColor[2],plotColor[5],'#8a65d6', plotColor[4], plotColor[7], plotColor[9])
names(cols) <- c('Late RPC', 'Neurog.', 'GABA AC','Glycine AC', 'PR Precur.', 'Cone', 'Rod')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(P2_metadata_sample$newCellType, levels = c('Late RPC', 'Neurog.', 'GABA AC','Glycine AC', 'PR Precur.', 'Cone', 'Rod'))

pdf(file = paste0('R_downstream/mouse/P2/RegulonHeatmap.pdf'),width = 2.8,height = 6)
cys_heatmap <- 	Heatmap(P2_filtered_auc_scaled, 
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

P2_rod_edge <- getTFNetworkDF(P2_rss$Rod[c(1,3:4,9,13,15,16,17)],P2_regulons_info, sample = 500)
P2_rod_node <- createNetworkNode(P2_rod_edge)
pdf(file = paste0('R_downstream/mouse/P2/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_rod_edge,P2_rod_node)
dev.off()

P2_neur_edge <- getTFNetworkDF(P2_rss$Neurog.[c(3:5,8,9,11,12,14,16,17)],P2_regulons_info, sample = 500)
P2_neur_node <- createNetworkNode(P2_neur_edge)
pdf(file = paste0('R_downstream/mouse/P2/neur_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_neur_edge,P2_neur_node)
dev.off()

P2_prp_edge <- getTFNetworkDF(P2_rss$PR.Precur[c(1:3,8:13,16,17,20)],P2_regulons_info, sample = 500)
P2_prp_node <- createNetworkNode(P2_prp_edge)
pdf(file = paste0('R_downstream/mouse/P2/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_prp_edge,P2_prp_node)
dev.off()

P2_lrpc_edge <- getTFNetworkDF( P2_rss$Late.RPC[c(1:5,8,10,12:16,20)],P2_regulons_info, sample = 500)
P2_lrpc_node <- createNetworkNode(P2_lrpc_edge)
pdf(file = paste0('R_downstream/mouse/P2/lrpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_lrpc_edge,P2_lrpc_node)
dev.off()

P2_cone_edge <- getTFNetworkDF(P2_rss$Cone[c(5,7,10,16,18,20)],P2_regulons_info, sample = 500)
P2_cone_node <- createNetworkNode(P2_cone_edge)
pdf(file = paste0('R_downstream/mouse/P2/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_cone_edge,P2_cone_node)
dev.off()

P2_gly_edge <- getTFNetworkDF(P2_rss$Glycine.AC[c(1:10,13,14,19)],P2_regulons_info, sample = 500)
P2_gly_node <- createNetworkNode(P2_gly_edge)
pdf(file = paste0('R_downstream/mouse/P2/gly_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_gly_edge,P2_gly_node)
dev.off()

P2_gabaac_edge <- getTFNetworkDF(P2_rss$GABA.AC[c(2,3,4,6,14,20)],P2_regulons_info, sample = 500)
P2_gabaac_node <- createNetworkNode(P2_gabaac_edge)
pdf(file = paste0('R_downstream/mouse/P2/gabaac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P2_gabaac_edge,P2_gabaac_node)
dev.off()

saveRDS(P2_regulons_info,'R_downstream/mouse/P2/regulon_info.rds')

#--------------------------------------------------P5-----------------------------------------------
P5_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[8],'/regulon_info.rds'))

DefaultAssay(P5_regulons_info$seurat) <- 'RNA'
P5_regulons_info$seurat <- P5_regulons_info$seurat %>%
		 SCTransform(vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

p1 <- myPlotUMAP_clusters(P5_regulons_info$seura,groups = 'cellType')
ggsave('R_downstream/mouse/P5/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)

P5_metadata <- P5_regulons_info$seurat@meta.data
P5_regulon_auc <- P5_regulons_info$filtered_regulons_auc
P5_regulon_auc <- P5_regulon_auc[,rownames(P5_metadata)]
rownames(P5_regulon_auc) <- P5_regulons_info$tf_info$new_names

P5_metadata$cross_species <- P5_metadata$cellType
P5_metadata[P5_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
P5_cellactivity <- calculateCellGroupsActivity(P5_regulon_auc, P5_metadata, 'cross_species',scale = TRUE)
write.table(P5_cellactivity, 'R_downstream/mouse/P5/P5_cellactivity.txt', sep = '\t',quote = FALSE)

P5_rss <- myCalRSS(P5_regulon_auc,P5_metadata, 'P5',groupby='cellType')

for (x in colnames(P5_rss)){
    if (!file.exists(paste0('R_downstream/mouse/P5/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/P5/rss_info/',x)))
    }
    for (regulons in P5_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(P5_regulons_info$seurat,
                                                             P5_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/mouse/P5/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}

rss_for_plot <- c(P5_rss$Late.RPC[c(1,3,4,5,6,9,10,15,17,18)],
				  P5_rss$AC[c(1,4,5,9,11,12,13)], # AC
				  P5_rss$BC[c(7)], #BC
				  P5_rss$PR.Precur.[c(2:4)],
				  P5_rss$Rod[c(1:9,13,14,16)],
				  P5_rss$Cone[c(2:4,6,8:10,13,14,17)],
				  P5_rss$MG[c(1:8,12,19)]) 

rss_for_plot <- unique(rss_for_plot)
P5_metadata$barcode <- rownames(P5_metadata)
P5_metadata_sample <- P5_metadata %>% group_by(cellType) %>% slice_sample(n=100)
P5_filtered_auc <-  P5_regulon_auc[rss_for_plot,P5_metadata_sample$barcode]
P5_filtered_auc_scaled <- t(scale(t(P5_filtered_auc), center = T, scale=T))

cols <- c(plotColor[8],plotColor[5], plotColor[10],plotColor[4], plotColor[7], plotColor[9],plotColor[10])
names(cols) <- c('Late RPC', 'AC', 'BC','PR Precur.', 'Cone', 'Rod','MG')


top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(P5_metadata_sample$cellType, levels = c('Late RPC', 'AC', 'BC','PR Precur.', 'Cone', 'Rod','MG'))

pdf(file = paste0('R_downstream/mouse/P5/RegulonHeatmap.pdf'),width = 3,height = 6)
cys_heatmap <- 	Heatmap(P5_filtered_auc_scaled, 
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

P5_ac_edge <- getTFNetworkDF(P5_rss$AC[c(1,4,5,9,11,12,13)],P5_regulons_info, sample = 500)
P5_ac_node <- createNetworkNode(P5_ac_edge)
pdf(file = paste0('R_downstream/mouse/P5/ac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P5_ac_edge,P5_ac_node)
dev.off()


P5_lrpc_edge <- getTFNetworkDF(P5_rss$Late.RPC[c(1,3,4,5,6,9,10,15,17,18)],P5_regulons_info, sample = 500)
P5_lrpc_node <- createNetworkNode(P5_lrpc_edge)
pdf(file = paste0('R_downstream/mouse/P5/lrpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P5_lrpc_edge,P5_lrpc_node)
dev.off()

P5_mg_edge <- getTFNetworkDF(P5_rss$MG[c(1:8,10,12,19)],P5_regulons_info, sample = 500)
P5_mg_node <- createNetworkNode(P5_mg_edge)
pdf(file = paste0('R_downstream/mouse/P5/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P5_mg_edge,P5_mg_node)
dev.off()

P5_prp_edge <- getTFNetworkDF(P5_rss$PR.Precur.[c(2:4)],P5_regulons_info, sample = 500)
P5_prp_node <- createNetworkNode(P5_prp_edge)
pdf(file = paste0('R_downstream/mouse/P5/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P5_prp_edge,P5_prp_node)
dev.off()

P5_rod_edge <- getTFNetworkDF(P5_rss$Rod[c(1:9,13,14,16)],P5_regulons_info, sample = 500)
P5_rod_node <- createNetworkNode(P5_rod_edge)
pdf(file = paste0('R_downstream/mouse/P5/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P5_rod_edge,P5_rod_node)
dev.off()

P5_cone_edge <- getTFNetworkDF(P5_rss$Cone[c(2:4,6,8:10,13,14,17)],P5_regulons_info, sample = 500)
P5_cone_node <- createNetworkNode(P5_cone_edge)
pdf(file = paste0('R_downstream/mouse/P5/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P5_cone_edge,P5_cone_node)
dev.off()

saveRDS(P5_regulons_info,'R_downstream/mouse/P5/regulon_info.rds')
#--------------------------------------------------P8-----------------------------------------------

P8_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[9],'/regulon_info.rds'))

DefaultAssay(P8_regulons_info$seurat) <- 'RNA'
P8_regulons_info$seurat <- P8_regulons_info$seurat %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

# check the subtype of BC
#------		 
P8_bipolar <- subset(P8_regulons_info$seurat, cellType == 'BC')
DefaultAssay(P8_bipolar) <- 'RNA'
P8_bipolar <- P8_bipolar %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(P8_bipolar) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    P8_bipolar <- FindClusters(P8_bipolar,resolution = i, verbose= FALSE)
}

# BC marker
# Rod BC: Prkca,Car8, Sebox
# ON BC: Isl1, Grm6
# OFF BC: Scgn, Tacr3

Rod_BC <- colnames(subset(P8_bipolar, RNA_snn_res.0.1 %in% c('1','4')))
ON_BC <- colnames(subset(P8_bipolar, RNA_snn_res.0.1 %in% c('0','3')))
OFF_BC <- colnames(subset(P8_bipolar, RNA_snn_res.0.1 == '2'))

# reassign cell types
P8_regulons_info$seurat$newCellType <- P8_regulons_info$seurat$cellType
P8_regulons_info$seurat$newCellType[Rod_BC] <- 'Rod BC'
P8_regulons_info$seurat$newCellType[ON_BC] <- 'ON BC'
P8_regulons_info$seurat$newCellType[OFF_BC] <- 'OFF BC'

p1 <- myPlotUMAP_clusters(P8_regulons_info$seura,groups = 'newCellType')
ggsave('R_downstream/mouse/P8/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)


P8_metadata <- P8_regulons_info$seurat@meta.data
P8_regulon_auc <- P8_regulons_info$filtered_regulons_auc
P8_regulon_auc <- P8_regulon_auc[,rownames(P8_metadata)]
rownames(P8_regulon_auc) <- P8_regulons_info$tf_info$new_names

P8_metadata$cross_species <- P8_metadata$cellType
#P8_metadata[P8_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
P8_cellactivity <- calculateCellGroupsActivity(P8_regulon_auc, P8_metadata, 'cross_species',scale = TRUE)
write.table(P8_cellactivity, 'R_downstream/mouse/P8/P8_cellactivity.txt', sep = '\t',quote = FALSE)


P8_rss <- myCalRSS(P8_regulon_auc,P8_metadata, 'P8',groupby='newCellType')

for (x in colnames(P8_rss)){
    if (!file.exists(paste0('R_downstream/mouse/P8/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/P8/rss_info/',x)))
    }
    for (regulons in P8_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(P8_regulons_info$seurat,
                                                             P8_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/mouse/P8/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}

rss_for_plot <- c(P8_rss$ON.BC[c(2,3,12,13,17:19)],
				  P8_rss$OFF.BC[c(3)], 
				  P8_rss$Rod.BC[c(4,6,12:14)],  # BC
				  P8_rss$Cone[c(1,13,14,15)], #Cone
				  P8_rss$Rod[c(1,2,5,7,8,10)]) 


rss_for_plot <- unique(rss_for_plot)
P8_metadata$barcode <- rownames(P8_metadata)
P8_metadata_sample <- P8_metadata %>% group_by(cellType) %>% slice_sample(n=47)
P8_filtered_auc <-  P8_regulon_auc[rss_for_plot,P8_metadata_sample$barcode]
P8_filtered_auc_scaled <- t(scale(t(P8_filtered_auc), center = T, scale=T))

cols <- c(plotColor[11], plotColor[7], plotColor[9], plotColor[10])
names(cols) <- c('BC', 'Cone', 'Rod', 'MG')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))


column_splits <- factor(P8_metadata_sample$cellType, levels = c('BC', 'Cone', 'Rod', 'MG'))

pdf(file = paste0('R_downstream/mouse/P8/RegulonHeatmap.pdf'),width = 2,height = 4.8)
cys_heatmap <- 	Heatmap(P8_filtered_auc_scaled, 
                        name="Regulon rodtivity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_title = NULL,
                        top_annotation = top_anno,
                        column_split = column_splits,
                        use_raster = TRUE,  cluster_rows = FALSE,
                        border_gp = gpar(col = "black", lwd = 0.3),
                        heatmap_legend_param = list(title = 'Scaled Regulon activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()

P8_rod_edge <- getTFNetworkDF(P8_rss$Rod[c(1:9,13,14,16)],P8_regulons_info, sample = 500)
P8_rod_node <- createNetworkNode(P8_rod_edge)
pdf(file = paste0('R_downstream/mouse/P8/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P8_rod_edge,P8_rod_node)
dev.off()

bc_regulons <- unique(c(P8_rss$ON.BC[c(2,3,12,13,17:19)], P8_rss$OFF.BC[c(3)], P8_rss$Rod.BC[c(4,6,12:14)]))
P8_bc_edge <- getTFNetworkDF(bc_regulons,P8_regulons_info, sample = 500)
P8_bc_node <- createNetworkNode(P8_bc_edge)
pdf(file = paste0('R_downstream/mouse/P8/bc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P8_bc_edge,P8_bc_node)
dev.off()

P8_rodbc_edge <- getTFNetworkDF(P8_rss$Rod.BC[c(4,6,12:14)],P8_regulons_info, sample = 500)
P8_rodbc_node <- createNetworkNode(P8_rodbc_edge)
pdf(file = paste0('R_downstream/mouse/P8/rodbc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P8_rodbc_edge,P8_rodbc_node)
dev.off()

P8_cone_edge <- getTFNetworkDF(P8_rss$Cone[c(1,13,14,15)],P8_regulons_info, sample = 500)
P8_cone_node <- createNetworkNode(P8_cone_edge)
pdf(file = paste0('R_downstream/mouse/P8/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P8_cone_edge,P8_cone_node)
dev.off()

P8_mg_edge <- getTFNetworkDF(P8_rss$Rod[c(1,2,5,7,8,10)],P8_regulons_info, sample = 200)
P8_mg_node <- createNetworkNode(P8_mg_edge)
pdf(file = paste0('R_downstream/mouse/P8/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P8_mg_edge,P8_mg_node)
dev.off()

P8_rod_edge <- getTFNetworkDF(P8_rss$Rod[c(1,2,5,7,8,10)],P8_regulons_info, sample = 200)
P8_rod_node <- createNetworkNode(P8_rod_edge)
pdf(file = paste0('R_downstream/mouse/P8/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P8_rod_edge,P8_rod_node)
dev.off()


saveRDS(P8_regulons_info,'R_downstream/mouse/P8/regulon_info.rds')
#--------------------------------------------------P14-----------------------------------------------
P14_regulons_info <- readRDS(paste0('R_downstream/mouse/',stage_names[10],'/regulon_info.rds'))

DefaultAssay(P14_regulons_info$seurat) <- 'RNA'

P14_regulons_info$seurat <- P14_regulons_info$seurat %>%
		 SCTransform(vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

# check the subtype of BC
#------		 
P14_bipolar <- subset(P14_regulons_info$seurat, cellType == 'BC')
DefaultAssay(P14_bipolar) <- 'RNA'
P14_bipolar <- P14_bipolar %>%
		 SCTransform(vars.to.regress = c('percent.mt'),verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 FindNeighbors(reduction = 'pca', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20, verbose = FALSE)

resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    P14_bipolar <- FindClusters(P14_bipolar,resolution = i, verbose= FALSE)
}
# BC marker
# Rod BC: Prkca,Car8, Sebox
# ON BC: Isl1, Grm6
# OFF BC: Scgn, Tacr3

# reassign cell types
P14_regulons_info$seurat$newCellType <- P14_regulons_info$seurat$cellType
P14_regulons_info$seurat$newCellType[colnames(P14_bipolar)] <- 'Rod BC'

p1 <- myPlotUMAP_clusters(P14_regulons_info$seura,groups = 'newCellType')
ggsave('R_downstream/mouse/P14/umap.pdf', p1, dpi = 200, units = 'cm', width = 12, height = 10)

P14_metadata <- P14_regulons_info$seurat@meta.data
P14_regulon_auc <- P14_regulons_info$filtered_regulons_auc
P14_regulon_auc <- P14_regulon_auc[,rownames(P14_metadata)]
rownames(P14_regulon_auc) <- P14_regulons_info$tf_info$new_names

P14_metadata$cross_species <- P14_metadata$cellType
#P14_metadata[P14_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
P14_cellactivity <- calculateCellGroupsActivity(P14_regulon_auc, P14_metadata, 'cross_species',scale = TRUE)
write.table(P14_cellactivity, 'R_downstream/mouse/P14/P14_cellactivity.txt', sep = '\t',quote = FALSE)


P14_rss <- myCalRSS(P14_regulon_auc,P14_metadata, 'P14',groupby='newCellType')



#intersect(P14_markers[P14_markers$cluster=='MG',]$gene, str_split(P14_rss[,1], pattern = ' ', simplify = TRUE)[,1])
for (x in colnames(P14_rss)){
    if (!file.exists(paste0('R_downstream/mouse/P14/rss_info/',x))){
        dir.create((paste0('R_downstream/mouse/P14/rss_info/',x)))
    }
    for (regulons in P14_rss[,x]){

        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression(P14_regulons_info$seurat,
                                                             P14_regulon_auc,
                                                             regulon_name,
                                                             regulons,
                                                             'RNA')

        ggsave(paste0('R_downstream/mouse/P14/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 
}

rss_for_plot <- c(P14_rss$Rod.BC[c(2:4,7)], # BC
				  P14_rss$Cone[c(2,8,9,11,12,14,16,17,20)], #Cone # MG
				  P14_rss$Rod[c(1:3,6,7)],
				  P14_rss$MG[c(1:7,9:11,12,14,15,18)])  #Rod X2

rss_for_plot <- unique(rss_for_plot)
P14_metadata$barcode <- rownames(P14_metadata)
P14_metadata_sample <- P14_metadata %>% group_by(newCellType) %>% slice_sample(n=100)
P14_filtered_auc <-  P14_regulon_auc[rss_for_plot,P14_metadata_sample$barcode]
P14_filtered_auc_scaled <- t(scale(t(P14_filtered_auc), center = T, scale=T))

cols <- c(plotColor[11], plotColor[7], plotColor[9],plotColor[10])
names(cols) <-c('Rod BC', 'Cone', 'Rod', 'MG')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))


column_splits <- factor(P14_metadata_sample$newCellType, levels = c('Rod BC', 'Cone', 'Rod','MG'))

pdf(file = paste0('R_downstream/mouse/P14/RegulonHeatmap.pdf'),width = 2,height = 6)
cys_heatmap <- 	Heatmap(P14_filtered_auc_scaled, 
                        name="Regulon rodtivity",
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
                        heatmap_legend_param = list(title = 'Scaled Regulon rodtivity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()

P14_rod_edge <- getTFNetworkDF(P14_rss$Rod[c(1:3,6,7,11,15)],P14_regulons_info, sample = 200)
P14_rod_node <- createNetworkNode(P14_rod_edge)
pdf(file = paste0('R_downstream/mouse/P14/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P14_rod_edge,P14_rod_node)
dev.off()

P14_cone_edge <- getTFNetworkDF(P14_rss$Cone[c(2,8,9,11,12,14,16,17,20)],P14_regulons_info, sample = 200)
P14_cone_node <- createNetworkNode(P14_cone_edge)
pdf(file = paste0('R_downstream/mouse/P14/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P14_cone_edge,P14_cone_node)
dev.off()

P14_mg_edge <- getTFNetworkDF(P14_rss$MG[c(1:7,9:11,12,14,15,18)],P14_regulons_info, sample = 200)
P14_mg_node <- createNetworkNode(P14_mg_edge)
pdf(file = paste0('R_downstream/mouse/P14/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P14_mg_edge,P14_mg_node)
dev.off()

P14_rodbc_edge <- getTFNetworkDF(P14_rss$Rod.BC[c(2:4,7)],P14_regulons_info, sample = 200)
P14_rodbc_node <- createNetworkNode(P14_rodbc_edge)
pdf(file = paste0('R_downstream/mouse/P14/rodbc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(P14_rodbc_edge,P14_rodbc_node)
dev.off()

saveRDS(P14_regulons_info,'R_downstream/mouse/P14/regulon_info.rds')

