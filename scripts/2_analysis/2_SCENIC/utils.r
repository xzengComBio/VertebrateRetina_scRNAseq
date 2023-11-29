library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(Seurat)

#  ----------------------------------------- basic regulon analysis --------------------------------------

# read the file loom file and save it to a list
# input: the file path of final loom
getRegulonInfoFromLoom <- function(scenicLoomFilePath){

	loom_object <- open_loom(scenicLoomFilePath)
	
	exprMat <- get_dgem(loom_object)
	regulons_incidMat <- get_regulons(loom_object, column.attr.name="Regulons")
	regulons <- regulonsToGeneLists(regulons_incidMat)
	regulonAUC <- get_regulons_AUC(loom_object,column.attr.name='RegulonsAUC')
	regulonAucThresholds <- get_regulon_thresholds(loom_object)
	# make sure to close the loom file
	close_loom(loom_object)

	regulon_info <- list(exprMat = exprMat,
						 regulons_incidMat = regulons_incidMat,
						 regulons = regulons,
						 regulonAUC = regulonAUC,
						 regulonAucThresholds = regulonAucThresholds)

	return(regulon_info)
}


# calculate the mean and sd of auc activity in all cells
# input: the auc matrix
calculateAUCInfo <- function(auc){

	# filter out regulons that are not active in all cells
	# otherwise the sd cannot be calculated.
	auc <- auc[rowSums(auc)> 0,]
	auc_t <- t(auc)

	auc_mean <- apply(auc_t,MARGIN = 2,mean)
	auc_sd <- apply(auc_t,MARGIN = 2,sd)
	auc_expressed <- colSums(auc_t > 0.01)
	auc_stat <- data.frame(Mean = auc_mean, SD = auc_sd, active_cell=auc_expressed)

	return(auc_stat)
}


# filter out regulons  with TF express in lower than 10(default) cells 
filterRegulonsByLowExpressedTF <- function(regulon_info,regulon_auc_stat, threshold=10){

	regulon_auc_stat$gene_name <- str_split(rownames(regulon_auc_stat), '_', simplify = TRUE)[,1]
	tf_count_matrix <- regulon_info$exprMat[regulon_auc_stat$gene_name,]
	regulon_auc_stat$expressed_cells <- rowSums(tf_count_matrix > 0)
	# the TF at least expressed in 10 cells
	saved_regulons <- regulon_auc_stat[regulon_auc_stat$expressed_cells > threshold, ]$gene_name
	regulon_auc_stat <- regulon_auc_stat[regulon_auc_stat$gene_name %in% saved_regulons, ]

	return(regulon_auc_stat)
}

# get the target genes by TF name
getTargetGenes <- function(tf_name, regulon_info){

	regulon_name <- paste0(tf_name, '_(+)')
	target_genes <- regulon_info$regulons[[regulon_name]]

	return(target_genes)

}

# get the target gene info by TF names
getTargetGeneInfo <- function(tf_names, regulon_info){

	target_genes_total <- c()
	target_gene_num_total <- c()

	for (x in tf_names){

		target_genes <- getTargetGenes(x, regulon_info)
		target_gene_num <- length(target_genes)

		target_genes_total <- c(target_genes_total,target_genes <- paste(target_genes, collapse=';'))
		target_gene_num_total <- c(target_gene_num_total, target_gene_num)

	}

	df <- data.frame(tf_name=tf_names, target_genes=target_genes_total, target_gene_num=target_gene_num_total)

	return(df)

}


# This function will be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds){

  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(auc[x,]>trh))
                                   }),names(thresholds))

  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}


plotBinarizedRegulonActivityAndTFExpression <- function(seuratObject,auc, tf_name, new_names,assay_use='SCT'){


	if (new_names %in% rownames(auc)){

		DefaultAssay(seuratObject) <- assay_use
		seuratObject$TargetRegulon <- auc[new_names,]
		p1 <- FeaturePlot(seuratObject, features = 'TargetRegulon',reduction='umap',raster = TRUE) + 
				theme_bw() +
						theme(panel.grid.major=element_blank(),
							  panel.grid.minor=element_blank(),
							  axis.title = element_text(face = "bold",size = rel(1)),
							  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
						ggtitle(new_names) + NoLegend()


		p2 <- FeaturePlot(seuratObject, features = tf_name,reduction='umap',raster = TRUE) + 
				theme_bw() +
						theme(panel.grid.major=element_blank(),
							  panel.grid.minor=element_blank(),
							  axis.title = element_text(face = "bold",size = rel(1)),
							  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
						ggtitle(tf_name)

		plot1 <- p1 | p2

	return(plot1)
}

}


plotBinarizedRegulonActivityAndTFExpression_cz <- function(seuratObject,auc, tf_name,gene_symbol, new_names,assay_use='RNA'){

	
	DefaultAssay(seuratObject) <- assay_use
	seuratObject$TargetRegulon <- auc[new_names,]
	p1 <- FeaturePlot(seuratObject, features = 'TargetRegulon',reduction='umap',raster = TRUE) + 
			theme_bw() +
					theme(panel.grid.major=element_blank(),
						  panel.grid.minor=element_blank(),
						  axis.title = element_text(face = "bold",size = rel(1)),
						  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
					ggtitle(new_names) + NoLegend()

	if (gene_symbol %in% rownames(seuratObject)){

		p2 <- FeaturePlot(seuratObject, features = gene_symbol,reduction='umap',raster = TRUE) + 
				theme_bw() +
				theme(panel.grid.major=element_blank(),
							panel.grid.minor=element_blank(),
							axis.title = element_text(face = "bold",size = rel(1)),
							plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
				ggtitle(gene_symbol)
			}else if(tf_name %in% rownames(seuratObject)){

		p2 <- FeaturePlot(seuratObject, features = tf_name,raster = TRUE) + 
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


#  -------------------------------------------------------------------------------------------------------

loadTFAnnotation <- function(path = '~/Desktop/Research/ciona/TF/GenomeInfo/ghost_tf.txt'){

	tf_info <- read.delim('~/Desktop/Research/ciona/TF/GenomeInfo/ghost_tf.txt')
	tf_info <- tf_info[!duplicated(tf_info$KH_name, fromLast = FALSE),]
	rownames(tf_info) <- tf_info$KH_name

	return(tf_info)
}


initializePathFromStage <- function(dev_stage){

	result_path =file.path(paste0("R_downstream/", dev_stage))

	motifEnrichmentFilePath <- file.path(paste0("Second_reg/", dev_stage, '_777_reg.csv'))
	scenicLoomFilePath <- file.path(paste0("Final/", dev_stage, '_777_auc_final.loom'))

	path_list <- list(result_path = result_path,
					  motifEnrichmentFilePath = motifEnrichmentFilePath,
					  scenicLoomFilePath = scenicLoomFilePath)

	return(path_list)
} 

replaceRegulonTFName <- function(regulonAUC, tf_info){

	kh_name <- str_extract(rownames(regulonAUC), '(?<=KH2013:).+?(?=\\_)')
	new_names <- tf_info[kh_name,]
	new_names$Notes <- str_split(new_names$Notes, pattern = ':', simplify = TRUE)[,2]
	new_names[new_names$Notes != "",]$GeneName <- new_names[new_names$Notes != "",]$Notes

	new_names$GeneName <- str_trim(new_names$GeneName, 'both')
	new_names$GeneName <- gsub('//','/', new_names$GeneName)

	rownames(regulonAUC) <- new_names$GeneName

	return(regulonAUC)
}



plotAUCStat <- function(auc_stat){
	p1 <- ggplot(auc_stat, aes(x = Mean, y = SD)) + 
  				geom_point() + 
  				geom_hline(aes(yintercept= 0.01)) +
  				geom_vline(aes(xintercept= 0.01)) + 
  				theme_bw() +
      			theme(panel.grid.major=element_blank(),
        				panel.grid.minor=element_blank(),
        				axis.title = element_text(face = "bold",size = rel(1)),
        				plot.title = element_blank(),
        				legend.position="none")
    return(p1)
}


filterRegulonsWithSD <- function(auc,auc_stat,sd_threshold){

  auc_filter_names <- rownames(subset(auc_stat,SD >= sd_threshold))
  auc_filter <- auc[auc_filter_names,]
  auc_scaled <- t(scale(t(auc_filter), center = T, scale=T))

  return(auc_scaled)
}


calculateCellGroupsActivity <- function(auc,cell_meta, group_by, scale=FALSE){

	regulonActivity_byTissue <- t(scale(t(auc), center = T, scale=T))
	cellsPerGroup <- split(rownames(cell_meta), cell_meta[,group_by])
	regulonActivity_byTissue <- sapply(cellsPerGroup,function(cells) rowMeans(regulonActivity_byTissue[,cells]))
	return(regulonActivity_byTissue)

}


plotRegulonActivityAndTFExpression <- function(seuratObject,auc, tf_name, assay_use='SCT'){

	regulon_name <- paste0(tf_name,'_(+)')
	DefaultAssay(seuratObject) <- assay_use
	seuratObject$TargetRegulon <- auc[regulon_name,]
	p1 <- FeaturePlot(seuratObject, features = 'TargetRegulon') + 
			theme_bw() +
					theme(panel.grid.major=element_blank(),
						  panel.grid.minor=element_blank(),
						  axis.title = element_text(face = "bold",size = rel(1)),
						  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
					ggtitle(regulon_name)


	p2 <- FeaturePlot(seuratObject, features = tf_name) + 
			theme_bw() +
					theme(panel.grid.major=element_blank(),
						  panel.grid.minor=element_blank(),
						  axis.title = element_text(face = "bold",size = rel(1)),
						  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
					ggtitle(tf_name)

	plot1 <- p1 | p2

	return(plot1)

}

newPlotRegulonActivityAndTFExpression <- function(seuratObject,auc, tf_name, assay_use='SCT'){

	DefaultAssay(seuratObject) <- assay_use
	seuratObject$TargetRegulon <- auc[tf_name,]
	p1 <- FeaturePlot(seuratObject, features = 'TargetRegulon') + 
			theme_bw() +
					theme(panel.grid.major=element_blank(),
						  panel.grid.minor=element_blank(),
						  axis.title = element_text(face = "bold",size = rel(1)),
						  plot.title = element_text(hjust = 0.5,colour = "black", face = "bold")) + 
					ggtitle(tf_name)

	return(p1)

}





newFilterOutRegulons <- function(file_path,gene_symbols_df,regulon_info, max_targets=1000, min_targets=10){

	regulon_names <- list.files(file_path)
	regulon_names <- str_split(regulon_names, pattern='\\.', simplify=TRUE)[,1]
	print(regulon_names)
	regulon_names_df <- gene_symbols_df[gene_symbols_df$V2 %in% regulon_names,]
	ensembl_names <- regulon_names[!regulon_names %in% regulon_names_df$V2]
	regulon_names_df <- rbind(regulon_names_df, data.frame(V1=ensembl_names, V2=rep("", length(ensembl_names))))
	regulon_names_df$V1 <- paste0(regulon_names_df$V1, '_(+)')

	num_target_genes <- c()
	for (x in regulon_names_df$V1){
		num_target_genes <- c(num_target_genes, length(regulon_info$regulons[[x]]))
	}

	saved_regulons <- regulon_names_df[(num_target_genes < max_targets) & (num_target_genes > min_targets),]
	num_target_genes <- num_target_genes[(num_target_genes < max_targets) & (num_target_genes > min_targets)]
	filter_out_regulons <- setdiff(regulon_names_df$V1, saved_regulons$V1)
	print("The following regluons are filtered out!!!!!")
	print("_____________________________________________")
	print(filter_out_regulons)
	print("_____________________________________________")

	#print(new_names)
	saved_regulons$new_names <- paste0(saved_regulons$V2, " (" , num_target_genes, 'g)')
	#print(length(new_names))

	return(saved_regulons)
}

FilterOutRegulons <- function(file_path,regulon_info, max_targets=1000, min_targets=10){

	regulon_names <- list.files(file_path)
	regulon_names <- str_split(regulon_names, pattern='\\.', simplify=TRUE)[,1]
	regulon_names <- paste0(regulon_names, '_(+)')
	print(length(regulon_names))
	num_target_genes <- c()
	for (x in regulon_names){
		num_target_genes <- c(num_target_genes, length(regulon_info$regulons[[x]]))
	}

	saved_regulons <- regulon_names[(num_target_genes < max_targets) & (num_target_genes > min_targets)]
	num_target_genes <- num_target_genes[(num_target_genes < max_targets) & (num_target_genes > min_targets)]
	print(saved_regulons)
	print(length(saved_regulons))
	filter_out_regulons <- setdiff(regulon_names, saved_regulons)
	print("The following regluons are filtered out!!!!!")
	print("_____________________________________________")
	print(filter_out_regulons)
	print("_____________________________________________")

	new_names <- str_split(saved_regulons, pattern='_', simplify=TRUE)[,1]
	#print(new_names)
	new_names <- paste0(new_names, " (" , num_target_genes, 'g)')
	#print(length(new_names))

	name_list <- list(regulon_name=saved_regulons, new_names=new_names,num_targets=num_target_genes)
	return(name_list)
}




getTargetGeneInfo <- function(tf_names, regulon_info){

	target_genes_total <- c()
	target_gene_num_total <- c()

	for (x in tf_names){

		target_genes <- getTargetGenes(x, regulon_info)
		target_gene_num <- length(target_genes)

		target_genes_total <- c(target_genes_total,target_genes <- paste(target_genes, collapse=';'))
		target_gene_num_total <- c(target_gene_num_total, target_gene_num)

	}

	df <- data.frame(tf_name=tf_names, target_genes=target_genes_total, target_gene_num=target_gene_num_total)

	return(df)

}



calOverlapOrthgroups <- function(orthogroups, species, geneSYMBOLENSEMBL){



	new_matrix <- matrix(nrow = dim(orthogroups[1]), ncol = dim(orthogroups)[2])
	colnames(new_matrix) <- colnames(orthogroups)
	rownames(new_matrix) <- rownames(orthogroups)

	# loop the orthogroups
	for (x in 1:dim(orthogroups)[1]){
		# loop the species
		for (y in 1:dim(orthogroups)[2]){
			gene_names <- str_split(orthogroups[x,y], pattern = ';')[[1]]
			genes_intersect <- intersect(gene_names,gene_list[[y]])

			if (length(genes_intersect) == 0){
				new_matrix[x,y] <- 0			
			}
			else if (length(genes_intersect) == 1){
				new_matrix[x,y] <- genes_intersect
			}
			else if (length(genes_intersect) > 1){
				new_matrix[x,y] <- paste(genes_intersect, collapse = ';')
			}
		}


	}

	new_matrix <- as.data.frame(new_matrix)

	return(new_matrix)

}

#human_gene_SymbolENSEMBL <- parseGeneSymbolEnsemblDF('../../TF/zebrafishAndChicken/gtf/human/geneSymbolENSEMBL.txt', make_unique = TRUE)
parseGeneSymbolEnsemblDF <- function(file_path, make_unique=FALSE){

	gene_symbol_ensembl_df <- read.table(file_path, sep = '\t')
	colnames(gene_symbol_ensembl_df) <- c('Ensembl', 'Gene_symbol')
	if (make_unique){
		gene_symbol_ensembl_df$Gene_symbol <- make.unique(gene_symbol_ensembl_df$Gene_symbol)
		rownames(gene_symbol_ensembl_df) <- gene_symbol_ensembl_df$Gene_symbol
	}

	return(gene_symbol_ensembl_df)

} 


# human_tf_info <- createGeneTFinfo(rownames(human_regulon_info$regulonAUC), human_gene_SymbolENSEMBL)
createGeneTFinfo <- function(regulon_names, SymbolENSEMBL_df){

	gene_symbols <- str_split(regulon_names, '_', simplify = TRUE)[,1]
	gene_names_df <- data.frame(gene_symbol=gene_symbols)
	gene_names_df$Ensembl <- SymbolENSEMBL_df[gene_symbols,]$Ensembl

	return(gene_names_df)
}

AddorthgroupInfo <- function(tf_info_df,orthogroup_df){

	rownames(orthogroup_df) <- 	orthogroup_df$Ensembl
	tf_info_df$Orthogroups <- orthogroup_df[tf_info_df$Ensembl,]$Orthogroups

	return(tf_info_df)
}

newAddorthgroupInfo <- function(genes,orthogroup_df){

	rownames(orthogroup_df) <- 	orthogroup_df$Ensembl

	orthogroups <- orthogroup_df[genes,]$Orthogroups
	return(orthogroups)
}


regulonAUCVlnplot <- function(seuratObject,regulonFeatureName){

	p1 <- VlnPlot(seuratObject, features = regulonFeatureName, pt.size = FALSE) + 
				geom_boxplot(width=0.1,fill="white",outlier.size = 0) + 
				theme(axis.title.x=element_blank(),
							axis.title.y = element_text(face="bold"),
							plot.title = element_blank(),
							axis.text.x=element_text(angle=0),
							axis.ticks.x=element_blank()) +
				 NoLegend()

	return(p1)

}


getOrthologous <- function(genes,out_species, orthogroups){

	index_list <- list(human = 'ENSG00',
										 mouse = 'ENSMUSG',
										 chicken = 'ENSGALG',
										 zebrafish = 'ENSDARG',
										 ciona = 'KH')
	genes_orth <- orthogroups[orthogroups$Ensembl %in% genes,]
	orth_genes <- orthogroups[(startsWith(orthogroups$Ensembl, index_list[[out_species]])) & (orthogroups$Orthogroups %in% genes_orth$Orthogroups),]

	output_df <- merge(genes_orth,orth_genes,by ='Orthogroups', all=FALSE)

	return(output_df)

}

changeGeneNames <- function(genes, geneSymbol_df,in_type='gene_symbol'){

	if (in_type == 'gene_symbol'){
		out_genes <- geneSymbol_df[geneSymbol_df$Gene_symbol %in% genes, ]$Ensembl
	} else if (in_type == 'ensembl'){

		out_genes <- geneSymbol_df[geneSymbol_df$Ensembl %in% genes, ]$Gene_symbol
	}
	
	return(out_genes)

}




