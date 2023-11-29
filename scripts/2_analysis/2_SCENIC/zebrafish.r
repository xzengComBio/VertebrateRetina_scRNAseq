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
library(xlsx)
library(circlize)
library(clusterProfiler)
library(Canek)
library(tidyverse)

library(igraph)
library(ggraph)
library(tidygraph)

setwd('~/Desktop/Research/ciona/SCENIC/results')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')

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
	#new_AUC_threshold <- AUCell_exploreThresholds(regulon_info$regulonAUC,plotHist = FALSE, assignCells = TRUE)
	#regulon_info$newRegulonAucThresholds <- getThresholdSelected(new_AUC_threshold)

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
	tf_info$gene_symbol <- zebrafish_gene_SymbolENSEMBL[tf_info$tf_name,]$Gene_symbol
	tf_info$new_names <- paste0(tf_info$gene_symbol, ' (', as.character(tf_info$target_gene_num), 'g)')
	rownames(tf_info) <- tf_info$tf_name

	# save the files
	regulon_info$seurat <- seurat_object
	regulon_info$filtered_regulons_auc <- regulon_auc_filtered
	regulon_info$tf_info <- tf_info

	print('saving regulon info....')
	saveRDS(regulon_info, paste0('R_downstream/zebrafish/',stage,'/regulon_info.rds'))

}

myCalRSS <- function(regulon_auc_filtered, metadata, stage,groupby='cellType'){

	#9. Calculate the RSS (regulon specific score) for each cell type.
	rss <- calcRSS(AUC=regulon_auc_filtered, cellAnnotation=metadata[, groupby])

	rss_list <- list()
	for (x in colnames(rss)){
		rss_list[[x]] <- names(sort(rss[,x],decreasing = TRUE)[1:30])
	}
	rss_df <- as.data.frame(rss_list)

	write.table(rss_df,paste0('R_downstream/zebrafish/',stage,'/rss.txt'), sep = '\t', quote = FALSE,row.names = FALSE)
	return(rss_df)

}

plotBinarizedRegulonActivityAndTFExpression_cz <- function(seuratObject,auc, tf_name,gene_symbol, new_names,assay_use='SCT'){

	
	DefaultAssay(seuratObject) <- assay_use
	seuratObject$TargetRegulon <- auc[new_names,]
	p1 <- FeaturePlot(seuratObject, features = 'TargetRegulon',reduction='umap') + 
			theme_bw() +
					theme(panel.grid.major=element_blank(),
						  panel.grid.minor=element_blank(),
						  axis.title = element_text(face = "bold",size = rel(1)),
						  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
					ggtitle(new_names) + NoLegend()

	if (gene_symbol %in% rownames(seuratObject)){

		p2 <- FeaturePlot(seuratObject, features = gene_symbol,reduction='umap') + 
				theme_bw() +
				theme(panel.grid.major=element_blank(),
							panel.grid.minor=element_blank(),
							axis.title = element_text(face = "bold",size = rel(1)),
							plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
				ggtitle(gene_symbol)
			}else if(tf_name %in% rownames(seuratObject)){

		p2 <- FeaturePlot(seuratObject, features = tf_name) + 
					theme_bw() +
							theme(panel.grid.major=element_blank(),
								  panel.grid.minor=element_blank(),
								  axis.title = element_text(face = "bold",size = rel(1)),
								  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
							ggtitle(tf_name)

			}else{
				x <- rnorm(100,mean = 2, sd = 3)
				y <- -1.5 + 2*x + rnorm(100)
				df <- data.frame(x = x, y = y)
				p2 <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point()
			}



	plot1 <- p1 | p2

	return(plot1)

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

zebrafish_gene_SymbolENSEMBL <- read.table('../../TF/zebrafishAndChicken/gtf/zebrafish/geneSymbolENSEMBL.txt', sep = '\t')
colnames(zebrafish_gene_SymbolENSEMBL) <- c('ENSEMBL', 'Gene_symbol')
zebrafish_gene_SymbolENSEMBL[zebrafish_gene_SymbolENSEMBL$Gene_symbol == '',]$Gene_symbol <- zebrafish_gene_SymbolENSEMBL[zebrafish_gene_SymbolENSEMBL$Gene_symbol == '',]$ENSEMBL
rownames(zebrafish_gene_SymbolENSEMBL) <- zebrafish_gene_SymbolENSEMBL$ENSEMBL

#--------------------------------------------------merge all-----------------------------------------------
merge_all <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/merge_all.rds')

loom_file = paste0('FromShirokane/Final/zebrafish_merge_all_auc_final.loom')
myPreprocessingRegulonActivity(loom_file,merge_all, 'merge_all')

#--------------------------------------------------14dpf-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
dpf14 <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/dpf14_filtered.rds')
#--------------------------------------------------72hpf-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
hpf72 <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/hpf72_filtered.rds')

#--------------------------------------------------48hpf-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
hpf48 <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/hpf48_filtered.rds')

#--------------------------------------------------36hpf-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
hpf36 <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/hpf36_filtered.rds')

#--------------------------------------------------24hpf-----------------------------------------------
# 1. Read the final loom file of regulon and the seurat files(for annotation)
hpf24 <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/hpf24_filtered.rds')

#-------------------------------------------------------------------------------------------------

stage_names <- c('24hpf', '36hpf', '48hpf', '72hpf', '14dpf')

seurat_list <- list(hpf24 = hpf24, hpf36 = hpf36, hpf48 = hpf48, hpf72 = hpf72, dpf14 = dpf14)

for (x in 1:5){

	print('-----------------------------------------------------')
	print(paste0('Processing ', stage_names[x],'...'))
	loom_file = paste0('FromShirokane/Final/', stage_names[x], '_auc_final.loom')
	myPreprocessingRegulonActivity(loom_file,seurat_list[[x]], stage_names[x])
}


hpf24_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[1],'/regulon_info.rds'))
hpf36_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[2],'/regulon_info.rds'))
hpf48_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[3],'/regulon_info.rds'))
hpf72_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[4],'/regulon_info.rds'))
dpf14_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[5],'/regulon_info.rds'))

zebrafish <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/merge_all.rds')

#------------ save the regulon informatrion to excel ----------
regulon_num_by_stage <- data.frame(stage = stage_names, 
								   regulon = c(nrow(hpf24_regulon_info$tf_info),
								   			   nrow(hpf36_regulon_info$tf_info),
								   			   nrow(hpf48_regulon_info$tf_info),
								   			   nrow(hpf72_regulon_info$tf_info),
								   			   nrow(dpf14_regulon_info$tf_info)))
write.xlsx(regulon_num_by_stage,file = 'R_downstream/regulon_summary.xlsx',sheetName = 'summary',,append=TRUE,row.names = FALSE)

hpf24_excel <- hpf24_regulon_info$tf_info[,-c(2,5)]
hpf24_excel <- hpf24_excel[,c(1,3,2)]
write.xlsx(hpf24_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = '24hpf',append=TRUE, row.names = FALSE)

hpf36_excel <- hpf36_regulon_info$tf_info[,-c(2,5)]
hpf36_excel <- hpf36_excel[,c(1,3,2)]
write.xlsx(hpf36_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = '36hpf',append=TRUE, row.names = FALSE)

hpf48_excel <- hpf48_regulon_info$tf_info[,-c(2,5)]
hpf48_excel <- hpf48_excel[,c(1,3,2)]
write.xlsx(hpf48_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = '48hpf',append=TRUE, row.names = FALSE)

hpf72_excel <- hpf72_regulon_info$tf_info[,-c(2,5)]
hpf72_excel <- hpf72_excel[,c(1,3,2)]
write.xlsx(hpf72_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = '72hpf',append=TRUE, row.names = FALSE)

dpf14_excel <- dpf14_regulon_info$tf_info[,-c(2,5)]
dpf14_excel <- dpf14_excel[,c(1,3,2)]
write.xlsx(dpf14_excel,file = 'R_downstream/regulon_summary.xlsx',sheetName = '14dpf',append=TRUE, row.names = FALSE)


#---------------------------------24hpf----------------------------------------
#myBinarizedandPlotRegulonActivity(hpf24_regulon_info,'24hpf')

hpf24_metadata <- hpf24_regulon_info$seurat@meta.data
hpf24_regulon_auc_filtered <- hpf24_regulon_info$filtered_regulons_auc
hpf24_regulon_auc_filtered <- hpf24_regulon_auc_filtered[,rownames(hpf24_metadata)]
rownames(hpf24_regulon_auc_filtered) <- hpf24_regulon_info$tf_info$new_names

p1 <- myPlotUMAP_clusters(hpf24_regulon_info$seurat,groups = 'module')
ggsave('R_downstream/zebrafish/24hpf/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)


# core_regulons_2 <- c('bhlhe40 (451g)', 'e2f3 (33g)', 'e2f4 (27g)', 'foxg1b (64g)', 'foxp4 (9g)', 'gsx2 (12g)', 'mitfa (1485g)',
# 					'msx1b (123g)', 'msx2b (421g)', 'msx3 (95g)', 'neurod1 (125g)', 'neurod4 (150g)', 'nfyal (124g)', 'nr2f1a (140g)',
# 					'nr2f1b (46g)', 'nr2f5 (206g)', 'otx1 (6g)', 'pax2a (160g)', 'rbpjb (122g)', 'sox19a (235g)', 'tfec (726g)', 'xbp1 (1653g)', 'zgc:101100 (102g)') 

hpf24_rss <- myCalRSS(hpf24_regulon_auc_filtered,hpf24_metadata, '24hpf',groupby='module')

for (x in colnames(hpf24_rss)){
	if (!file.exists(paste0('R_downstream/zebrafish/24hpf/rss_info/',x))){
		dir.create((paste0('R_downstream/zebrafish/24hpf/rss_info/',x)))
	}
	for (regulons in hpf24_rss[,x]){
		regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
		p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(hpf24_regulon_info$seurat,
															 hpf24_regulon_auc_filtered,
															 regulon_name,
															 regulon_name,
															 regulons,
															 'SCT')

		ggsave(paste0('R_downstream/zebrafish/24hpf/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
	} 

}

rss_for_plot <- unique(c(hpf24_rss$X1[c(1,2,4:9,10,12,13,14,25,30)], 
						  hpf24_rss$X2[c(1:8,10:13,19,20,30)], 
						  hpf24_rss$X3[c(2:7)]))


hpf24_filtered_auc <-  hpf24_regulon_auc_filtered[rss_for_plot,]
hpf24_filtered_auc_scaled <- t(scale(t(hpf24_filtered_auc), center = T, scale=T))

cols <- c(plotColor[1],'#d94865', '#d67e90')
names(cols) <- c('1', '2','3')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
			      labels = names(cols),
			      labels_gp = gpar(cex=0.5, col='black')))

column_splits <- factor(hpf24_metadata$module, levels =c('1', '2','3'))

pdf(file = paste0('R_downstream/zebrafish/24hpf/RegulonHeatmap.pdf'),width =5,height = 5)
cys_heatmap <- 	Heatmap(hpf24_filtered_auc_scaled, 
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
                        heatmap_legend_param = list(title = 'Scaled Regulon Activity', direction = 'horizontal',title_position = 'lefttop'), 
                        column_gap = unit('0.1', 'mm'))
draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()
saveRDS(hpf24_regulon_info,'R_downstream/zebrafish/24hpf/regulon_info.rds')

hpf24_module1_edge <- getTFNetworkDF(hpf24_rss$X1[c(1,2,4:9,10,12,13,14,25,30)],hpf24_regulon_info, sample = 500, is_mouse = FALSE)
hpf24_module1_node <- createNetworkNode(hpf24_module1_edge)
pdf(file = paste0('R_downstream/zebrafish/24hpf/module1_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf24_module1_edge,hpf24_module1_node)
dev.off()



#---------------------------------36hpf----------------------------------------

p1 <- myPlotUMAP_clusters(hpf36_regulon_info$seurat,groups = 'module')
ggsave('R_downstream/zebrafish/36hpf/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)

#myBinarizedandPlotRegulonActivity(hpf36_regulon_info,'36hpf')
hpf36_metadata <- hpf36_regulon_info$seurat@meta.data
hpf36_regulon_auc_filtered <- hpf36_regulon_info$filtered_regulons_auc
hpf36_regulon_auc_filtered <- hpf36_regulon_auc_filtered[,rownames(hpf36_metadata)]
rownames(hpf36_regulon_auc_filtered) <- hpf36_regulon_info$tf_info$new_names

hpf36_cellactivity_rpc <- calculateCellGroupsActivity(hpf36_regulon_auc_filtered, hpf36_metadata, 'module',scale = TRUE)
write.table(hpf36_cellactivity_rpc[,c(1,2,3)], 'R_downstream/zebrafish/36hpf/hpf36_cellactivity_rpc.txt', sep = '\t',quote = FALSE)


hpf36_metadata$cross_species <- hpf36_metadata$cellType

hpf36_cellactivity <- calculateCellGroupsActivity(hpf36_regulon_auc_filtered, hpf36_metadata, 'cross_species',scale = TRUE)
write.table(hpf36_cellactivity, 'R_downstream/zebrafish/36hpf/hpf36_cellactivity.txt', sep = '\t',quote = FALSE)


hpf36_rss <- myCalRSS(hpf36_regulon_auc_filtered,hpf36_metadata, '36hpf',groupby='module')

for (x in colnames(hpf36_rss)){
    if (!file.exists(paste0('R_downstream/zebrafish/36hpf/rss_info/',x))){
        dir.create((paste0('R_downstream/zebrafish/36hpf/rss_info/',x)))
    }
    for (regulons in hpf36_rss[,x]){
        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(hpf36_regulon_info$seurat,
                                                             hpf36_regulon_auc_filtered,
                                                             regulon_name,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/zebrafish/36hpf/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}

rss_for_plot <- c(hpf36_rss$X1[c(1:3,5,8,12,15,22)],
		  hpf36_rss$X2[-c(8,11,12,13,16,20:30)],#AC
		  hpf36_rss$X3[c(1:5,8:12,14,18)],#Cone
		  hpf36_rss$RGC[c(1:15,17,18,20,21,25,26)])
rss_for_plot <- unique(rss_for_plot)

hpf36_metadata$barcode <- rownames(hpf36_metadata)
hpf36_metadata_sample <- hpf36_metadata %>% group_by(module) %>% slice_sample(n=100)
hpf36_filtered_auc <-  hpf36_regulon_auc_filtered[rss_for_plot,hpf36_metadata_sample$barcode]
hpf36_filtered_auc_scaled <- t(scale(t(hpf36_filtered_auc), center = T, scale=T))

cols <- c(plotColor[1],'#d94865', '#d67e90', plotColor[3])
names(cols) <- c('1', '2','3', 'RGC')
top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(hpf36_metadata_sample$module, levels =c('1', '2','3', 'RGC'))

pdf(file = paste0('R_downstream/zebrafish/36hpf/RegulonHeatmap.pdf'),width =2,height = 6)
cys_heatmap <- 	Heatmap(hpf36_filtered_auc_scaled, 
                        name="Regulon rodbctivity",
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
                        heatmap_legend_param = list(title = 'Scaled Regulon rodbctivity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)


draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()

hpf36_module1_edge <- getTFNetworkDF(hpf36_rss$X1[c(1:3,5,8,12,15,22)],hpf36_regulon_info, sample = 500, is_mouse = FALSE)
hpf36_module1_node <- createNetworkNode(hpf36_module1_edge)
pdf(file = paste0('R_downstream/zebrafish/36hpf/module1_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf36_module1_edge,hpf36_module1_node)
dev.off()

hpf36_module2_edge <- getTFNetworkDF(hpf36_rss$X2[-c(8,11,12,13,16,20:30)],hpf36_regulon_info, sample = 1000, is_mouse = FALSE)
hpf36_module2_node <- createNetworkNode(hpf36_module2_edge)
pdf(file = paste0('R_downstream/zebrafish/36hpf/module2_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf36_module2_edge,hpf36_module2_node)
dev.off()

hpf36_module3_edge <- getTFNetworkDF(hpf36_rss$X3[c(1:5,8:12,14,18)],hpf36_regulon_info, sample = 1000, is_mouse = FALSE)
hpf36_module3_node <- createNetworkNode(hpf36_module3_edge)
pdf(file = paste0('R_downstream/zebrafish/36hpf/module3_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf36_module3_edge,hpf36_module3_node)
dev.off()

hpf36_rgc_edge <- getTFNetworkDF(hpf36_rss$RGC[c(1:15,17,18,20,21,25,26)],hpf36_regulon_info, sample = 1000, is_mouse = FALSE)
hpf36_rgc_node <- createNetworkNode(hpf36_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/36hpf/rgc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf36_rgc_edge,hpf36_rgc_node)
dev.off()

saveRDS(hpf36_regulon_info,'R_downstream/zebrafish/36hpf/regulon_info.rds')

for (x in colnames(rss)){

	hpf36_tf_info <- hpf36_regulon_info$tf_info[hpf36_regulon_info$tf_info$new_names %in% rss[,x],]
	hpf36_tf_target_list <- hpf36_regulon_info$regulons[paste0(hpf36_tf_info$tf_name, '_(+)')]

	names(hpf36_tf_target_list) <- hpf36_tf_info$gene_symbol

	hpf36_go <- compareCluster(hpf36_tf_target_list,
                             fun=enrichGO,
                             OrgDb = 'org.Dr.eg.db',
                             keyType = 'ENSEMBL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.005)

	saveRDS(hpf36_go,paste0('R_downstream/zebrafish/36hpf/',x,'_go.rds'))

	pdf(file = paste0('R_downstream/zebrafish/36hpf/',x,'_go.pdf'),
		width = 10,
		height = 12)

	dotplot(hpf36_go)
	dev.off()

}

#---------------------------------48hpf----------------------------------------
#myBinarizedandPlotRegulonActivity(hpf48_regulon_info,'48hpf')

hpf48_regulon_info$seurat$newCellType <- hpf48_regulon_info$seurat$cellType
hpf48_regulon_info$seurat$newCellType[hpf48_regulon_info$seurat$RNA_snn_res.0.6 == '9'] <- 'OFF BC'
hpf48_regulon_info$seurat$newCellType[hpf48_regulon_info$seurat$RNA_snn_res.0.6 == '10'] <- 'ON BC'
hpf48_regulon_info$seurat$newCellType[hpf48_regulon_info$seurat$newCellType == 'PR'] <- 'PR Precurs.'

p1 <- myPlotUMAP_clusters(hpf48_regulon_info$seurat,groups = 'newCellType')
ggsave('R_downstream/zebrafish/48hpf/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)


hpf48_metadata <- hpf48_regulon_info$seurat@meta.data
hpf48_regulon_auc_filtered <- hpf48_regulon_info$filtered_regulons_auc
hpf48_regulon_auc_filtered <- hpf48_regulon_auc_filtered[,rownames(hpf48_metadata)]
rownames(hpf48_regulon_auc_filtered) <- hpf48_regulon_info$tf_info$new_names
hpf48_cellactivity_rpc <- calculateCellGroupsActivity(hpf48_regulon_auc_filtered, hpf48_metadata, 'module',scale = TRUE)
write.table(hpf48_cellactivity_rpc[,c(1,2,3)], 'R_downstream/zebrafish/48hpf/hpf48_cellactivity_rpc.txt', sep = '\t',quote = FALSE)


hpf48_metadata$cross_species <- hpf48_metadata$cellType
hpf48_metadata[hpf48_metadata$cross_species %in% c('PR'),]$cross_species <- 'PR Precurs.'
hpf48_cellactivity <- calculateCellGroupsActivity(hpf48_regulon_auc_filtered, hpf48_metadata, 'cross_species',scale = TRUE)
write.table(hpf48_cellactivity, 'R_downstream/zebrafish/48hpf/hpf48_cellactivity.txt', sep = '\t',quote = FALSE)


hpf48_rss <- myCalRSS(hpf48_regulon_auc_filtered,hpf48_metadata, '48hpf',groupby='newCellType')

for (x in colnames(hpf48_rss)){
    if (!file.exists(paste0('R_downstream/zebrafish/48hpf/rss_info/',x))){
        dir.create((paste0('R_downstream/zebrafish/48hpf/rss_info/',x)))
    }
    for (regulons in hpf48_rss[,x]){
        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(hpf48_regulon_info$seurat,
                                                             hpf48_regulon_auc_filtered,
                                                             regulon_name,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/zebrafish/48hpf/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}


rss_for_plot <- c(hpf48_rss$RPC[c(1:13,15)],
			hpf48_rss$RGC[c(1:5,7,8,12,15:18,24)],
			hpf48_rss$AC[c(2:12)],
			hpf48_rss$HC[c(2:6,9,16)],
		  hpf48_rss$BC[c(1,4,6,7,11,14,22,27)],
		  hpf48_rss$OFF.BC[c(3,7,11)],
		  hpf48_rss$ON.BC[c(2)],
		  hpf48_rss$PR.Precurs.[c(1:3,6,7,9,12)],
		  hpf48_rss$MG[c(3,4,6,11,12,28)])
rss_for_plot <- unique(rss_for_plot)
hpf48_metadata$barcode <- rownames(hpf48_metadata)
hpf48_metadata_sample <- hpf48_metadata %>% group_by(cellType) %>% slice_sample(n=100)
hpf48_filtered_auc <-  hpf48_regulon_auc_filtered[rss_for_plot,hpf48_metadata_sample$barcode]
hpf48_filtered_auc_scaled <- t(scale(t(hpf48_filtered_auc), center = T, scale=T))

cols <- c(plotColor[1],plotColor[3],plotColor[5],plotColor[6],plotColor[11],plotColor[4],plotColor[10])
names(cols) <- c('RPC', 'RGC','AC', 'HC','BC', 'PR.Precurs.','MG')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

column_splits <- factor(hpf48_metadata_sample$cellType, levels =c('RPC','RGC','AC', 'HC','BC', 'PR Precurs.','MG'))
pdf(file = paste0('R_downstream/zebrafish/48hpf/RegulonHeatmap.pdf'),width =3,height = 6)
cys_heatmap <- 	Heatmap(hpf48_filtered_auc_scaled, 
                        name="Regulon rgctivity",
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
                        heatmap_legend_param = list(title = 'Scaled Regulon activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)


draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()

hpf48_rgc_edge <- getTFNetworkDF(hpf48_rss$AC[c(2:12)],hpf48_regulon_info, sample = 1000, is_mouse = FALSE)
hpf48_rgc_node <- createNetworkNode(hpf48_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/48hpf/ac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf48_rgc_edge,hpf48_rgc_node)
dev.off()

hpf48_rgc_edge <- getTFNetworkDF(hpf48_rss$BC[c(1,4,6,7,11,14,22,27)],hpf48_regulon_info, sample = 500, is_mouse = FALSE)
hpf48_rgc_node <- createNetworkNode(hpf48_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/48hpf/bc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf48_rgc_edge,hpf48_rgc_node)
dev.off()

hpf48_rgc_edge <- getTFNetworkDF(hpf48_rss$HC[c(2:6,9,16)],hpf48_regulon_info, sample = 500, is_mouse = FALSE)
hpf48_rgc_node <- createNetworkNode(hpf48_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/48hpf/hc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf48_rgc_edge,hpf48_rgc_node)
dev.off()

hpf48_rgc_edge <- getTFNetworkDF(hpf48_rss$MG[c(3,4,6,11,12,28)],hpf48_regulon_info, sample = 1000, is_mouse = FALSE)
hpf48_rgc_node <- createNetworkNode(hpf48_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/48hpf/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf48_rgc_edge,hpf48_rgc_node)
dev.off()

hpf48_rgc_edge <- getTFNetworkDF(hpf48_rss$OFF.BC[c(3,7,11)],hpf48_regulon_info, sample = 80, is_mouse = FALSE)
hpf48_rgc_node <- createNetworkNode(hpf48_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/48hpf/offbc_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(hpf48_rgc_edge,hpf48_rgc_node,no_red = TRUE)
dev.off()

saveRDS(hpf48_regulon_info,'R_downstream/zebrafish/48hpf/regulon_info.rds')


hpf48_regulonCellTypeActivity <- calculateCellGroupsActivity(hpf48_regulon_auc_filtered,
                                                             hpf48_metadata,
                                                             'cellType')

myPlotCellTypeHeatMap(hpf48_regulonCellTypeActivity, '48hpf', threshold = 0.03)
rss <- myCalRSS(hpf48_regulon_auc_filtered,hpf48_metadata, '48hpf',groupby='cellType')

for (x in colnames(rss)){

	hpf48_tf_info <- hpf48_regulon_info$tf_info[hpf48_regulon_info$tf_info$new_names %in% rss[,x],]
	hpf48_tf_target_list <- hpf48_regulon_info$regulons[paste0(hpf48_tf_info$tf_name, '_(+)')]

	names(hpf48_tf_target_list) <- hpf48_tf_info$gene_symbol

	hpf48_go <- compareCluster(hpf48_tf_target_list,
                             fun=enrichGO,
                             OrgDb = 'org.Dr.eg.db',
                             keyType = 'ENSEMBL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)

	saveRDS(hpf48_go,paste0('R_downstream/zebrafish/48hpf/',x,'_go.rds'))

	pdf(file = paste0('R_downstream/zebrafish/48hpf/',x,'_go.pdf'),
		width = 10,
		height = 12)

	dotplot(hpf48_go)
	dev.off()

}

#---------------------------------72hpf----------------------------------------
#myBinarizedandPlotRegulonActivity(hpf72_regulon_info,'72hpf')

hpf72_regulon_info$seurat$newCellType <- hpf72_regulon_info$seurat$cellType
hpf72_regulon_info$seurat$newCellType[hpf72_regulon_info$seurat$RNA_snn_res.0.2 %in% c('3','4')] <- 'ON BC'
hpf72_regulon_info$seurat$newCellType[hpf72_regulon_info$seurat$RNA_snn_res.0.2  == '0'] <- 'OFF BC'
hpf72_regulon_info$seurat$newCellType[hpf72_regulon_info$seurat$newCellType  == 'BC'] <- 'OFF BC'

#HC
hpf72_HC <- subset(hpf72_regulon_info$seurat, cellType == 'HC')
DefaultAssay(hpf72_HC) <- 'RNA'
hpf72_HC <-  hpf72_HC %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(hpf72_HC) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    hpf72_HC <- FindClusters(hpf72_HC,resolution = i, verbose= FALSE)
}

# T1 HC: lhx1a
# T2 HC: isl1
T1_HC <- colnames(subset(hpf72_HC, RNA_snn_res.0.05 == '1'))
T2_HC <- colnames(subset(hpf72_HC, RNA_snn_res.0.05 == '0'))

hpf72_regulon_info$seurat$newCellType[T1_HC] <- 'T1 HC'
hpf72_regulon_info$seurat$newCellType[T2_HC] <- 'T2 HC'

# Cone subtype
hpf72_Cone <- subset(hpf72_regulon_info$seurat, cellType == 'Cone')
DefaultAssay(hpf72_Cone) <- 'RNA'
hpf72_Cone <-  hpf72_Cone %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

DefaultAssay(hpf72_Cone) <- 'RNA'
resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
    hpf72_Cone <- FindClusters(hpf72_Cone,resolution = i, verbose= FALSE)
}
#

Red_Cone <- colnames(subset(hpf72_Cone, RNA_snn_res.0.4 == '0'))
Double_Cone <- colnames(subset(hpf72_Cone, RNA_snn_res.0.4 == '2'))
Green_Cone <- colnames(subset(hpf72_Cone, RNA_snn_res.0.4 == '1'))

hpf72_regulon_info$seurat$newCellType[Red_Cone] <- 'Red Cone'
hpf72_regulon_info$seurat$newCellType[Double_Cone] <- 'Double Cone'
hpf72_regulon_info$seurat$newCellType[Green_Cone] <- 'Green Cone'


hpf72_metadata <- hpf72_regulon_info$seurat@meta.data
hpf72_regulon_auc_filtered <- hpf72_regulon_info$filtered_regulons_auc
hpf72_regulon_auc_filtered <- hpf72_regulon_auc_filtered[,rownames(hpf72_metadata)]
rownames(hpf72_regulon_auc_filtered) <- hpf72_regulon_info$tf_info$new_names

hpf72_metadata$cross_species <- hpf72_metadata$cellType
hpf72_metadata[hpf72_metadata$cross_species %in% c('PR Precur.'),]$cross_species <- 'PR Precurs.'
hpf72_metadata[hpf72_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
hpf72_cellactivity <- calculateCellGroupsActivity(hpf72_regulon_auc_filtered, hpf72_metadata, 'cross_species',scale = TRUE)
write.table(hpf72_cellactivity, 'R_downstream/zebrafish/72hpf/hpf72_cellactivity.txt', sep = '\t',quote = FALSE)


p1 <- myPlotUMAP_clusters(hpf72_regulon_info$seurat,groups = 'newCellType')
ggsave('R_downstream/zebrafish/72hpf/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)
hpf72_rss <- myCalRSS(hpf72_regulon_auc_filtered,hpf72_metadata, '72hpf',groupby='newCellType')

for (x in colnames(hpf72_rss)){
    if (!file.exists(paste0('R_downstream/zebrafish/72hpf/rss_info/',x))){
        dir.create((paste0('R_downstream/zebrafish/72hpf/rss_info/',x)))
    }
    for (regulons in hpf72_rss[,x]){
        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(hpf72_regulon_info$seurat,
                                                             hpf72_regulon_auc_filtered,
                                                             regulon_name,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/zebrafish/72hpf/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}

rss_for_plot <- c(
			hpf72_rss$RPC[c(1:5,7:10,12)],
			hpf72_rss$Late.RPC[c(1:2,4:8)],
			hpf72_rss$MG[c(1:5,7,8,12)],
			hpf72_rss$RGC[c(1,4:8,10:16)],
	    hpf72_rss$AC[c(1:5,7,8,10,11,15:17,19)],
	    hpf72_rss$T1.HC[c(1:12)],
	    hpf72_rss$T2.HC[c(1:12)],
	    hpf72_rss$ON.BC[c(1:4,10,11,20)],
	    hpf72_rss$OFF.BC[c(1:3,7,14)],
		  hpf72_rss$PR.Precur.[c(1,2,5:9,12)],
		  hpf72_rss$Double.Cone[c(1,3:9)],#
		  hpf72_rss$Green.Cone[c(1,3,5:10)],#
		  hpf72_rss$Red.Cone[c(2,4,6,7,8,9,11:13)],
		  hpf72_rss$Rod[c(2:4,9,10)])

rss_for_plot <- unique(rss_for_plot)

hpf72_metadata$barcode <- rownames(hpf72_metadata)
#hpf72_metadata_sample <- hpf72_metadata %>% group_by(newCellType) %>% slice_sample(n=100)
hpf72_metadata$cellType[hpf72_metadata$cellType == 'Late RPC'] <- 'RPC'
hpf72_metadata_sample <- hpf72_metadata %>% group_by(cellType) %>% slice_sample(n=100)
hpf72_filtered_auc <-  hpf72_regulon_auc_filtered[rss_for_plot,hpf72_metadata_sample$barcode]
hpf72_filtered_auc_scaled <- t(scale(t(hpf72_filtered_auc), center = T, scale=T))

#cols <- c(plotColor[1], plotColor[8], plotColor[5], '#80800d','#999900',plotColor[6],'#87e087', plotColor[3],plotColor[4], 'Red',plotColor[7],'Green',plotColor[9],plotColor[10])
#names(cols)<- c('RPC', 'Late RPC', 'AC', 'OFF BC', 'ON BC', 'T1 HC', 'T2 HC', 'RGC', 'PR Precurs.', 'Red Cone', 'Double Cone', 'Green Cone', 'Rod', 'MG')

cols <- plotColor[c(1,10,3,5,6,11,4,7,9)]
names(cols)<- c('RPC', 'MG','RGC', 'AC', 'HC', 'BC', 'PR Precur.', 'Cone', 'Rod')
top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))

row_loc <- c(1,2,4,12,13,18:20,22,24:26,28,32,35,37,38,40,41,43,48,52,55,56,63,65,70,71)
row_label <- rownames(hpf72_filtered_auc_scaled)[row_loc]

row_anno <- rowAnnotation(
												foo = anno_mark(at = row_loc,
												 labels = row_label,
												 labels_gp = gpar(fontsize = 8)))


column_splits <- factor(hpf72_metadata_sample$cellType, levels =c('RPC','MG','RGC', 'AC', 'HC', 'BC', 'PR Precur.', 'Cone', 'Rod'))
pdf(file = paste0('R_downstream/zebrafish/72hpf/RegulonHeatmap.pdf'),width =3,height =6)
cys_heatmap <- 	Heatmap(hpf72_filtered_auc_scaled, 
                        name="Regulon rgctivity",
                        row_names_gp=grid::gpar(fontsize=8), 
                        show_column_names = FALSE,
                        show_row_names = FALSE, 
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_title = NULL,
                        right_annotation=row_anno,
                        top_annotation = top_anno,
                        column_split = column_splits,
                        use_raster = TRUE, cluster_rows = FALSE,
                        border_gp = gpar(col = "black", lwd = 0.3),
                        heatmap_legend_param = list(title = 'Scaled Regulon activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()



hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$Double.Cone[c(1,3:9)],hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/dcone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()


hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$Green.Cone[c(1,3,5:10)],hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/gcone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()

hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$Late.RPC[c(1:9)],hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/lrpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()

hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$MG[c(1:5,7,8,12,14,15)],hpf72_regulon_info, sample = 500, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()


hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$PR.Precur.[c(1,2,5:9,12)],hpf72_regulon_info, sample = 600, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()

hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$Red.Cone[c(2,4,6,7,8,9,11:13)],hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/rcone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()

hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$RGC[c(1,4:8,10:16)],hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/rgc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()

hpf72_rgc_edge <- getTFNetworkDF(hpf72_rss$Rod[c(2:4,9,10)],hpf72_regulon_info, sample = 200, is_mouse = FALSE)
hpf72_rgc_node <- createNetworkNode(hpf72_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_rgc_edge,hpf72_rgc_node)
dev.off()

hpf72_hc <- unique(c(hpf72_rss$T1.HC[c(1:5,7,9:12)],hpf72_rss$T2.HC[c(1:5,7,9,10,12)]))
hpf72_hc_edge <- getTFNetworkDF(hpf72_hc,hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_hc_node<- createNetworkNode(hpf72_hc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/hc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_hc_edge,hpf72_hc_node)
dev.off()


hpf72_bc <- unique(c(hpf72_rss$ON.BC[c(1:4,10,11,20)],hpf72_rss$OFF.BC[c(1:3,7,14)]))
hpf72_bc_edge <- getTFNetworkDF(hpf72_bc,hpf72_regulon_info, sample = 500, is_mouse = FALSE)
hpf72_bc_node<- createNetworkNode(hpf72_bc_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/bc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_bc_edge,hpf72_bc_node)
dev.off()

cone_regulons <-  c(hpf72_rss$Double.Cone[c(1,3:9)],#
		  hpf72_rss$Green.Cone[c(1,3,5:10)],#
		  hpf72_rss$Red.Cone[c(2,4,6,7,8,9,11:13)])
cone_regulons <- unique(cone_regulons)


hpf72_cone_edge <- getTFNetworkDF(cone_regulons,hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hpf72_cone_node <- createNetworkNode(hpf72_cone_edge)
pdf(file = paste0('R_downstream/zebrafish/72hpf/rcone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(hpf72_cone_edge,hpf72_cone_node)
dev.off()

saveRDS(hpf72_regulon_info,'R_downstream/zebrafish/72hpf/regulon_info.rds')

for (x in colnames(rss)){

    hpf72_tf_info <- hpf72_regulon_info$tf_info[hpf72_regulon_info$tf_info$new_names %in% rss[,x],]
    hpf72_tf_target_list <- hpf72_regulon_info$regulons[paste0(hpf72_tf_info$tf_name, '_(+)')]

    names(hpf72_tf_target_list) <- hpf72_tf_info$gene_symbol

    hpf72_go <- compareCluster(hpf72_tf_target_list,
                             fun=enrichGO,
                             OrgDb = 'org.Dr.eg.db',
                             keyType = 'ENSEMBL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)

    saveRDS(hpf72_go,paste0('R_downstream/zebrafish/72hpf/',x,'_go.rds'))

    pdf(file = paste0('R_downstream/zebrafish/72hpf/',x,'_go.pdf'),
        width = 10,
        height = 12)

    dotplot(hpf72_go)
    dev.off()

}

hc_edge <- getTFNetworkDF(hpf72_rss$T1.HC[c(1:12)], hpf72_regulon_info, sample = 1000, is_mouse = FALSE)
hc_node <- createNetworkNode(hc_edge)
p1 <- createAndPlotGraph(hc_edge,hc_node)


#---------------------------------14 dpf----------------------------------------


dpf14_metadata <- dpf14_regulon_info$seurat@meta.data
dpf14_regulon_auc_filtered <- dpf14_regulon_info$filtered_regulons_auc
dpf14_regulon_auc_filtered <- dpf14_regulon_auc_filtered[,intersect(colnames(dpf14_regulon_auc_filtered),rownames(dpf14_metadata))]
dpf14_metadata <- dpf14_metadata[intersect(colnames(dpf14_regulon_auc_filtered),rownames(dpf14_metadata)),]
rownames(dpf14_regulon_auc_filtered) <- dpf14_regulon_info$tf_info$new_names

dpf14_metadata$cross_species <- dpf14_metadata$cellType
dpf14_metadata[dpf14_metadata$cross_species %in% c('PR Precur.'),]$cross_species <- 'PR Precurs.'
#dpf14_metadata[dpf14_metadata$cross_species %in% c('Late RPC'),]$cross_species <- 'RPC'
dpf14_cellactivity <- calculateCellGroupsActivity(dpf14_regulon_auc_filtered, dpf14_metadata, 'cross_species',scale = TRUE)
write.table(dpf14_cellactivity, 'R_downstream/zebrafish/14dpf/dpf14_cellactivity.txt', sep = '\t',quote = FALSE)

dpf14_regulon_info$seurat <- subset(dpf14_regulon_info$seurat, cells=intersect(colnames(dpf14_regulon_auc_filtered),rownames(dpf14_metadata)))
p1 <- myPlotUMAP_clusters(dpf14_regulon_info$seurat,groups = 'cellType')
ggsave('R_downstream/zebrafish/14dpf/umap.pdf', p1, dpi = 200, units = 'cm', width = 14, height = 10)


#dpf14_regulonCellTypeActivity <- calculateCellGroupsActivity(dpf14_regulon_auc_filtered,
 #                                                            dpf14_metadata,
 #                                                            'cellType')

#myPlotCellTypeHeatMap(dpf14_regulonCellTypeActivity, '14dpf', threshold = 0.03)
dpf14_rss <- myCalRSS(dpf14_regulon_auc_filtered,dpf14_metadata, '14dpf',groupby='cellType')

for (x in colnames(dpf14_rss)){
    if (!file.exists(paste0('R_downstream/zebrafish/14dpf/rss_info/',x))){
        dir.create((paste0('R_downstream/zebrafish/14dpf/rss_info/',x)))
    }
    for (regulons in dpf14_rss[,x]){
        regulon_name <- str_split(regulons, pattern=' ', simplify = TRUE)[,1]
        p1 <- plotBinarizedRegulonActivityAndTFExpression_cz(dpf14_regulon_info$seurat,
                                                             dpf14_regulon_auc_filtered,
                                                             regulon_name,
                                                             regulon_name,
                                                             regulons,
                                                             'SCT')

        ggsave(paste0('R_downstream/zebrafish/14dpf/rss_info/',x,'/',regulon_name,'.pdf'), p1, dpi = 200, units = 'cm', width = 18, height = 8)
    } 

}


rss_for_plot <- c(
			dpf14_rss$RPC[c(1:9,11,14:16,19)],
			dpf14_rss$RGC[c(1,3:12,14:17,20)],
			dpf14_rss$AC[c(1:8,10,13:16,18,19)],
			dpf14_rss$HC[c(3,4,8,9,14,17)],
		  dpf14_rss$BC[c(1:7,10,13,14,16)],
		  dpf14_rss$PR.Precur.[c(1,2:9,11:14,16)],
		  dpf14_rss$Cone[c(1:3,5:9,11,14:16,18)],
		  dpf14_rss$Rod[c(1,5:8,13,14,17,19,20)],
		  dpf14_rss$MG[c(1:3,5,11)])

rss_for_plot <- unique(rss_for_plot)

dpf14_metadata$barcode <- rownames(dpf14_metadata)
dpf14_metadata_sample <- dpf14_metadata %>% group_by(cellType) %>% slice_sample(n=100)
dpf14_filtered_auc <-  dpf14_regulon_auc_filtered[rss_for_plot,dpf14_metadata_sample$barcode]
dpf14_filtered_auc_scaled <- t(scale(t(dpf14_filtered_auc), center = T, scale=T))

cols <- c(plotColor[1],plotColor[3], plotColor[5],plotColor[6],plotColor[11],plotColor[4],plotColor[7],plotColor[9],plotColor[10])
names(cols)<- c('RPC', 'RGC', 'AC','HC', 'BC', 'PR Precurs.', 'Cone', 'Rod', 'MG')

top_anno <- HeatmapAnnotation(cellType=anno_block(gp=gpar(fill=cols),
												  labels = names(cols),
												  labels_rot = 45,
												  labels_gp = gpar(cex=0.5, col='black'),
												  labels_just = 'center',
												   height = unit(2,'mm')))


column_splits <- factor(dpf14_metadata_sample$cellType, levels = c('RPC', 'RGC', 'AC','HC', 'BC', 'PR Precurs.', 'Cone', 'Rod', 'MG'))

pdf(file = paste0('R_downstream/zebrafish/14dpf/RegulonHeatmap.pdf'),width =3,height = 7.5)
cys_heatmap <- 	Heatmap(dpf14_filtered_auc_scaled, 
                        name="Regulon rgctivity",
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
                        heatmap_legend_param = list(title = 'Scaled Regulon activity', direction = 'horizontal',title_position = 'topcenter',
                        							at = c(-4, 0, 4),legend_width = unit(1, "cm")), 
                        column_gap = unit(0, 'mm'),
                        border = TRUE)

draw(cys_heatmap,  heatmap_legend_side = 'bottom')
dev.off()


dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$AC[c(1:8,10,13:16,18,19)],dpf14_regulon_info, sample = 1000, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/ac_network.pdf'),width = 5,height = 5)
createAndPlotGraph(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$BC[c(1:7,10,13,14,16)],dpf14_regulon_info, sample = 1000, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/bc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$Cone[c(1:3,5:9,11,14:16,18)],dpf14_regulon_info, sample = 1000, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/cone_network.pdf'),width = 5,height = 5)
createAndPlotGraph(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$HC[c(3,4,8,9,14,17)],dpf14_regulon_info, sample = 500, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/hc_network.pdf'),width = 5,height = 5)
createAndPlotGraph(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$MG[c(1:3,5,11,13,14,17,20)],dpf14_regulon_info, sample = 300, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/mg_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$PR.Precur.[c(1,2:9,11:14,16)],dpf14_regulon_info, sample = 300, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/prp_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$Rod[c(1,5:8,13,14,17,19,20)],dpf14_regulon_info, sample = 300, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/rod_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

dpf14_rgc_edge <- getTFNetworkDF(dpf14_rss$RPC[c(1:9,11,14:16,19)],dpf14_regulon_info, sample = 300, is_mouse = FALSE)
dpf14_rgc_node <- createNetworkNode(dpf14_rgc_edge)
pdf(file = paste0('R_downstream/zebrafish/14dpf/rpc_network.pdf'),width = 5,height = 5)
createAndPlotGraph_star(dpf14_rgc_edge,dpf14_rgc_node)
dev.off()

rss_for_plot <- c(dpf14_rss$AC[c(1:8,10,13:16,18,19)],
		  dpf14_rss$BC[c(1:7,10,13,14,16)],
		  dpf14_rss$Cone[c(1:3,5:9,11,14:16,18)],#
		  dpf14_rss$HC[c(3,4,8,9,14,17)],
		  dpf14_rss$MG[c(1:3,5,11,13,14,17,20)],
		  dpf14_rss$PR.Precurs.[c(1,2:9,11:14,16)],
		  dpf14_rss$RGC[c(1,3:12,14:17,20)],
		  dpf14_rss$Rod[c(1,5:8,13,14,17,19,20)],
		  dpf14_rss$RPC[c(1:9,11,14:16,19)])

saveRDS(dpf14_regulon_info,'R_downstream/zebrafish/14dpf/regulon_info.rds')


