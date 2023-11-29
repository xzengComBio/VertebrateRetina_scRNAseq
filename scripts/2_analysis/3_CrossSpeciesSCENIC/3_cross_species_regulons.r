library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(homologene)
library(ComplexHeatmap)
library(tidyverse)
library(biomaRt)
library(Canek)
library(xlsx)
library(networkD3)

setwd('~/Desktop/Research/ciona/SCENIC/results')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')

#------------------------------------------Function--------------------------------------
myCalRSS <- function(regulon_auc_filtered, metadata, groupby){

	#9. Calculate the RSS (regulon specific score) for each cell type.
	rss <- calcRSS(AUC=regulon_auc_filtered, cellAnnotation=metadata[, groupby])

	rss_list <- list()
	for (x in colnames(rss)){
		rss_list[[x]] <- names(sort(rss[,x],decreasing = TRUE)[1:20])
	}
	rss_df <- as.data.frame(rss_list)

	return(rss_df)

}

myCrossSPeciesPreprocessing <- function(regulon_info,select_regulons){

	regulon_genes <- regulon_info$regulons
	regulon_genes_filtered <- regulon_genes[rownames(regulon_info$filtered_regulons_auc)]
	names(regulon_genes_filtered) <- regulon_info$tf_info$new_names

	regulon_info$filter_genes <- regulon_genes_filtered[select_regulons]


	regulon_info$cell_rankings <- AUCell_buildRankings(regulon_info$exprMat[,colnames(regulon_info$seurat)])

	return(regulon_info)
}

Mouse2Others <- function(temp_regulon,out_species){

	txid <- c(zebrafish=7955,chicken=9031)
	homo_genes <- homologene(temp_regulon, inTax = 10090, outTax = txid[out_species])
	homo_genes <- changeGeneNames(homo_genes[,2],SymbolENSEMBL[[out_species]],in_type = 'gene_symbol')

	orth_genes <- getLDS(attributes = 'mgi_symbol', 
				 filters = 'mgi_symbol', 
				 mart = ensembl_db$mouse,
				 values = temp_regulon, 
				 attributesL ='ensembl_gene_id', 
				 martL = ensembl_db[[out_species]], 
				 uniqueRows = TRUE)
	out_regulon <- union(orth_genes$Gene.stable.ID, homo_genes)

	return(out_regulon)
} 

Others2Mouse <- function(gene_names,in_species){

    txid <- c(zebrafish=7955,chicken=9031)
    attr_list <- c(mouse = 'mgi_symbol', zebrafish = 'zfin_id_symbol', chicken = 'external_gene_name')
    homo_genes <- homologene(gene_names, inTax = txid[in_species], outTax = 10090)
    homo_genes <- homo_genes[,c(1,2)]
    colnames(homo_genes) <- c(in_species, 'mouse')
    orth_genes <- getLDS(attributes =attr_list[[in_species]], 
                 filters = 'zfin_id_symbol', 
                 mart = ensembl_db[[in_species]],
                 values = gene_names, 
                 attributesL ='mgi_symbol', 
                 martL = ensembl_db$mouse, 
                 uniqueRows = TRUE)
    colnames(orth_genes) <- c(in_species, 'mouse')
    out_regulon <- rbind(homo_genes,orth_genes)
    out_regulon <- out_regulon[!duplicated(out_regulon$mouse),]
    return(out_regulon)
} 



Others2Others <- function(temp_regulon, in_species, out_species){

	txid <- c(zebrafish=7955,chicken=9031)
	homo_genes <- changeGeneNames(temp_regulon,SymbolENSEMBL[[in_species]],in_type = 'ensembl')
	homo_genes <- homologene(homo_genes, inTax = txid[in_species], outTax = txid[out_species])

	homo_genes_ensembl <- changeGeneNames(homo_genes[,2],SymbolENSEMBL[[out_species]],in_type = 'gene_symbol')
	orth_genes <- getLDS(attributes = 'ensembl_gene_id', 
				 filters = 'ensembl_gene_id', 
				 mart = ensembl_db[[in_species]],
				 values = temp_regulon, 
				 attributesL ='ensembl_gene_id', 
				 martL = ensembl_db[[out_species]], 
				 uniqueRows = TRUE)
	out_regulon <- union(orth_genes[,2], homo_genes_ensembl)

}

getGeneListOrthologousGenes <- function(gene_list,in_species,out_species){

	out_gene_list <- list()
	for (x in 1:length(gene_list)){
		temp_regulon <- gene_list[[x]]
		if (in_species == 'mouse'){

			temp_regulon_orth <- Mouse2Others(temp_regulon, out_species)
		}else{

			if(out_species == 'mouse'){
				temp_regulon_orth <- Others2Mouse(temp_regulon, in_species)

			}else{

				temp_regulon_orth <- Others2Others(temp_regulon, in_species, out_species)
			}
		}
		out_gene_list[[x]] <- temp_regulon_orth
	}
	names(out_gene_list) <- names(gene_list)

	return(out_gene_list) 
}

calGeneSetActivity <- function(cells_rankings, gene_sets){

	regulon_info <- list()
	geneset_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank = nrow(cells_rankings)*0.05, nCores = 4)
	geneset_AUC_thresholds <- AUCell_exploreThresholds(geneset_AUC,plotHist = FALSE, assignCells = TRUE, nCores = 4)

	regulon_info$AUC <- geneset_AUC
	regulon_info$regulonAucThresholds <- getThresholdSelected(geneset_AUC_thresholds)

	return(regulon_info)
}

myCalRSS <- function(regulon_auc_filtered, metadata, stage,groupby='cellType', top=20){

	#9. Calculate the RSS (regulon specific score) for each cell type.
	rss <- calcRSS(AUC=regulon_auc_filtered, cellAnnotation=metadata[, groupby])

	rss_list <- list()
	for (x in colnames(rss)){
		rss_list[[x]] <- names(sort(rss[,x],decreasing = TRUE)[1:top])
	}
	rss_df <- as.data.frame(rss_list)

	write.table(rss_df,paste0('R_downstream/cross_species/',stage,'/rss.txt'), sep = '\t', quote = FALSE,row.names = FALSE)
	return(rss_df)

}

symbolChange <- function(gene_list,in_species, out_species){

	txid <- c(mouse=10090,zebrafish=7955,chicken=9031)
	homo_genes <- homologene(gene_list, inTax = txid[in_species], outTax = txid[out_species])[,c(1,2)]
	colnames(homo_genes) <- c(in_species,out_species)
	attr_list <- c(mouse = 'mgi_symbol', zebrafish = 'zfin_id_symbol', chicken = 'external_gene_name')
	orth_genes <- getLDS(attributes = attr_list[in_species], 
				 filters = attr_list[in_species], 
				 mart = ensembl_db[[in_species]],
				 values = gene_list, 
				 attributesL =attr_list[out_species], 
				 martL = ensembl_db[[out_species]], 
				 uniqueRows = TRUE)
	colnames(orth_genes) <- c(in_species,out_species)

	out_gene_all <- rbind(orth_genes, homo_genes)
	out_gene_all <- out_gene_all[!duplicated(out_gene_all),]

	return(out_gene_all)
} 

crossSpeciesRegulonIntersect <- function(tf_info,rss, species){

	intersect_regulons <- as.data.frame(matrix(ncol=2))
	colnames(intersect_regulons) <- c('cellType', 'tf_name')
	for (x in colnames(rss)){
		if ( x %in% tf_info$newCellType){

			tf_name <- paste(tf_info[tf_info$newCellType == x,]$stage,tf_info[tf_info$newCellType == x,]$tf_name, sep = '_') 
			temp_intersect <- intersect(rss[,x], tf_name)

			temp_df <- data.frame(cellType = rep(x, length(temp_intersect)),
							      tf_name = temp_intersect)

			intersect_regulons <- rbind(intersect_regulons, temp_df)

		}else{
			print(paste0(x, ' is not in the tf_info file, please check!'))
		}
	}

	intersect_regulons <- intersect_regulons[-1,]
	intersect_regulons$species <- species
	intersect_regulons <- intersect_regulons[,c(3,1,2)]
	return(intersect_regulons)

}

plotRegulonActivity  <- function(intersec_names, regulon_info, save_name){

	regulon_auc <- getAUC(regulon_info$AUC)
	for (regulon_name in intersec_names){

		regulon_info$seurat@meta.data['regulon_plot'] <- regulon_auc[regulon_name,]
		p1 <- FeaturePlot(regulon_info$seurat, features = 'regulon_plot', raster = TRUE, reduction='umap') + 
				theme_bw() +
				theme(panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				axis.title = element_text(face = "bold",size = rel(1)),
				plot.title = element_text(hjust = 0.5),
				legend.margin=margin(t = 0, unit='cm'),
				legend.key.size = unit(0.2, 'cm')) +
				ggtitle(regulon_name)

		ggsave(paste0('R_downstream/cross_species/',save_name, '/rss_activity/',regulon_name ,'.pdf'), 
			   p1,
			   dpi = 200, 
			   units = 'cm', 
			   width = 12,
			   height = 10)

	}	
}


#--------------------------------------------Preprocessing-----------------------------------------
# read symbol ENSEMBL info
SymbolENSEMBL <- list()
SymbolENSEMBL$mouse <- parseGeneSymbolEnsemblDF('../../TF/zebrafishAndChicken/gtf/mouse/geneSymbolENSEMBL.txt', make_unique = TRUE)
SymbolENSEMBL$chicken <- read.table('../../TF/zebrafishAndChicken/gtf/chicken/geneSymbolENSEMBL.txt', sep = '\t')
SymbolENSEMBL$zebrafish <- read.table('../../TF/zebrafishAndChicken/gtf/zebrafish/geneSymbolENSEMBL.txt', sep = '\t')

colnames(SymbolENSEMBL$chicken) <- colnames(SymbolENSEMBL$mouse)
colnames(SymbolENSEMBL$zebrafish) <- colnames(SymbolENSEMBL$mouse)

SymbolENSEMBL$chicken$Gene_symbol[SymbolENSEMBL$chicken$Gene_symbol == ''] <- SymbolENSEMBL$chicken[SymbolENSEMBL$chicken$Gene_symbol == '',]$Ensembl
SymbolENSEMBL$chicken <- SymbolENSEMBL$chicken[!duplicated(SymbolENSEMBL$chicke$Gene_symbol),]
rownames(SymbolENSEMBL$chicken) <- SymbolENSEMBL$chicken$Ensembl

SymbolENSEMBL$zebrafish$Gene_symbol[SymbolENSEMBL$zebrafish$Gene_symbol == ''] <- SymbolENSEMBL$zebrafish[SymbolENSEMBL$zebrafish$Gene_symbol == '',]$Ensembl
rownames(SymbolENSEMBL$zebrafish) <- SymbolENSEMBL$zebrafish$Ensembl

# 
db_host <- 'https://dec2021.archive.ensembl.org'
ensembl_db <- list()
ensembl_db$mouse <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host = db_host, mirror = 'asia')
ensembl_db$chicken <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'ggallus_gene_ensembl', host = db_host, mirror = 'asia')
ensembl_db$zebrafish <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl', host = db_host, mirror = 'asia')
# host = 'https://dec2021.archive.ensembl.org/'

tf_info <- readRDS('R_downstream/cross_species/tf_info.rds')
newCellType <- str_split(tf_info$mouse$cellTyp, pattern = ' ', simplify = TRUE)
newCellType[newCellType[,2] == '',][,2] <- newCellType[newCellType[,2] == '',][,1]
tf_info$mouse$newCellType <- newCellType[,2]
tf_info$mouse[tf_info$mouse$newCellType == 'Neurog.',]$newCellType <- 'RPC' 
tf_info$mouse[tf_info$mouse$newCellType == 'Precur.',]$newCellType <- 'PR Precur.' 


newCellType <- str_split(tf_info$chicken$cellTyp, pattern = ' ', simplify = TRUE)
newCellType[newCellType[,3] == '',][,3] <- newCellType[newCellType[,3] == '',][,2]
newCellType[newCellType[,3] == '',][,3] <- newCellType[newCellType[,3] == '',][,1]
tf_info$chicken$newCellType <- newCellType[,3]
tf_info$chicken[tf_info$chicken$newCellType == 'Precurs.',]$newCellType <- 'PR Precur.' 

newCellType <- str_split(tf_info$zebrafish$cellTyp, pattern = ' ', simplify = TRUE)
newCellType[newCellType[,2] == '',][,2] <- newCellType[newCellType[,2] == '',][,1]
tf_info$zebrafish$newCellType <- newCellType[,2]
tf_info$zebrafish[tf_info$zebrafish$newCellType == 'Precurs.',]$newCellType <- 'PR Precur.' 

tf_info$mouse$stage_regulon <- paste(tf_info$mouse$stage, tf_info$mouse$tf_name, sep = '_')
tf_info$chicken$stage_regulon <- paste(tf_info$chicken$stage, tf_info$chicken$tf_name, sep = '_')
tf_info$zebrafish$stage_regulon <- paste(tf_info$zebrafish$stage, tf_info$zebrafish$tf_name, sep = '_')



geneSymbolList <- list()
geneSymbolList$c2z <- symbolChange(unique(tf_info$chicken$tf_name),'chicken','zebrafish')
geneSymbolList$z2c <- symbolChange(unique(tf_info$zebrafish$tf_name),'zebrafish','chicken')

geneSymbolList$m2c <- symbolChange(unique(tf_info$mouse$tf_name),'mouse','chicken')
geneSymbolList$m2z <- symbolChange(unique(tf_info$mouse$tf_name),'mouse','zebrafish')



# get the regulon_gene_list
regulon_list <- list()
regulon_list$mouse <- readRDS('R_downstream/mouse/all_tf_info.rds')
regulon_list$chicken <- readRDS('R_downstream/chicken/all_tf_info.rds')
regulon_list$zebrafish <- readRDS('R_downstream/zebrafish/all_tf_info.rds')

species <- c('mouse', 'chicken', 'zebrafish')
regulon_select_list <- list()
for (x in 1:3){
	regulon_list[[x]]$stage_regulon <- str_replace(regulon_list[[x]]$stage_regulon, 'ENSGALG00000027907', 'NR2F1')
	rownames(regulon_list[[x]]) <- regulon_list[[x]]$stage_regulon
	
	temp_regulon <- regulon_list[[x]][unique(tf_info[[x]]$stage_regulon), ]$target_genes
	names(temp_regulon) <- unique(tf_info[[x]]$stage_regulon)

	regulon_select_list[[species[x]]] <- list()
	for (t in 1: length(temp_regulon)){
		regulon_select_list[[species[x]]][[names(temp_regulon[t])]] <- str_split(temp_regulon[t], pattern = ';')[[1]]

	}

}

# convert the gene names in the regulon to other two species
m2z_regulons <- getGeneListOrthologousGenes(regulon_select_list$mouse, 'mouse', 'zebrafish')
m2c_regulons <- getGeneListOrthologousGenes(regulon_select_list$mouse, 'mouse', 'chicken')

c2m_regulons <- getGeneListOrthologousGenes(regulon_select_list$chicken, 'chicken', 'mouse')
c2z_regulons <- getGeneListOrthologousGenes(regulon_select_list$chicken, 'chicken', 'zebrafish')

z2m_regulons <- getGeneListOrthologousGenes(regulon_select_list$zebrafish, 'zebrafish', 'mouse')
z2c_regulons <- getGeneListOrthologousGenes(regulon_select_list$zebrafish, 'zebrafish', 'chicken')

saveRDS(m2z_regulons,'R_downstream/cross_species/homologous/m2z_regulons.rds')
saveRDS(m2c_regulons,'R_downstream/cross_species/homologous/m2c_regulons.rds')
saveRDS(c2m_regulons,'R_downstream/cross_species/homologous/c2m_regulons.rds')
saveRDS(c2z_regulons,'R_downstream/cross_species/homologous/c2z_regulons.rds')
saveRDS(z2m_regulons,'R_downstream/cross_species/homologous/z2m_regulons.rds')
saveRDS(z2c_regulons,'R_downstream/cross_species/homologous/z2c_regulons.rds')

#----evaluate the activity of regulons  generated from mouse and chicken in zebrafish ---------------------------------------

# Load Zebrafish count matrix
sample_paths <- paste0(list.dirs('../../singleCell/other_species/Xu2020Development/count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
names(sample_paths) <- c('14dpf1', '14dpf2', '24hpf', '36hpf', '48hpf1',
						 '48hpf2', '72hpf1', '72hpf2')
zebrafish_retina_counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE, gene.column = 1)
sce <- CreateSeuratObject(zebrafish_retina_counts,names.field = 1, min.cells = 3, min.features = 200, project = 'zebrafish retina')

zebrafish_seurat <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/merge_all.rds')

sce <- subset(sce, cells = colnames(zebrafish_seurat))

zebrafish_count <- sce@assays$RNA@counts
zebrafish_count <- as.matrix(zebrafish_count)
zebrafish_count_ranking <- AUCell_buildRankings(zebrafish_count)

rm(zebrafish_retina_counts)
rm(sce)
gc()


c2z_info <- calGeneSetActivity(zebrafish_count_ranking, c2z_regulons)
m2z_info <- calGeneSetActivity(zebrafish_count_ranking, m2z_regulons)


saveRDS(c2z_info,'R_downstream/cross_species/homologous/c2z_info.rds')
saveRDS(m2z_info,'R_downstream/cross_species/homologous/m2z_info.rds')

#----evaluate the activity of regulons generated from mouse and zebrafish in chicken ---------------------------------------

# load chicken count matrix
sample_paths <- paste0(list.dirs('../../singleCell/other_species/Yamagata2021eLife/count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
names(sample_paths) <- c('E18A', 'E18B', 'E18C', 'E18D', 
						'dRGC1','dRGC2','E12A', 'E12B',
						'E12C', 'E12D', 'vRGC1', 'vRGC2')

chicken_retina_seurat <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_filtered_processed.rds')
chicken_counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE, gene.column = 1)
sce <- CreateSeuratObject(chicken_counts,names.field = 1, min.cells = 3, min.features = 200, project = 'chick retina')
sce <- subset(sce, cells = colnames(chicken_retina_seurat))

E18_names <- names(chicken_retina_seurat$stage[chicken_retina_seurat$stage == 'E18'])
E12_names <- names(chicken_retina_seurat$stage[chicken_retina_seurat$stage == 'E12'])
E16_names <- names(chicken_retina_seurat$stage[chicken_retina_seurat$stage == 'E16'])
E12_name_10000 <- sample(E12_names, size = 10000, replace = FALSE)

chichen_all_cell <- c(E12_name_10000,E16_names,E18_names)
chicken_all <- subset(sce, cells = chichen_all_cell)

chicken_count <- as.matrix(chicken_all@assays$RNA@counts)
chicken_count_ranking <- AUCell_buildRankings(chicken_count)

rm(chicken_counts)
rm(sce)
gc()

z2c_info <- calGeneSetActivity(chicken_count_ranking, z2c_regulons)
m2c_info <- calGeneSetActivity(chicken_count_ranking, m2c_regulons)


saveRDS(z2c_info,'R_downstream/cross_species/homologous/z2c_info.rds')
saveRDS(m2c_info,'R_downstream/cross_species/homologous/m2c_info.rds')

#----evaluate the activity of regulons generated from chicken and zebrafish in mouse ---------------------------------------

# load mouse matrix
sce_finalKeep <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')
Idents(sce_finalKeep) <- sce_finalKeep$cellType

new_cellType_ID <- c('Early RPC', 'Neurog.', 'RGC', 'Photo. Precurs.', 'AC', 'HC', 'Cone', 'Late RPC', 'Rod', 'MG', 'BC')
names(new_cellType_ID) <- levels(sce_finalKeep)
sce_finalKeep <- RenameIdents(sce_finalKeep, new_cellType_ID)
sce_finalKeep <- subset(sce_finalKeep, stage %in% c('E14', 'E18', 'P8', 'P14'))

mouse_count <- as.matrix(sce_finalKeep@assays$RNA@counts)
mouse_count_ranking <- AUCell_buildRankings(mouse_count)

rm(sce_finalKeep)
rm(mouse_count)
gc()

z2m_info <- calGeneSetActivity(mouse_count_ranking, z2m_regulons)
c2m_info <- calGeneSetActivity(mouse_count_ranking, c2m_regulons)

saveRDS(z2m_info,'R_downstream/cross_species/homologous/z2m_info.rds')
saveRDS(c2m_info,'R_downstream/cross_species/homologous/c2m_info.rds')


#---------------------------add seurat object to the auc file of mouse-------------------------------------------------
## mouse
sce_finalKeep <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')
Idents(sce_finalKeep) <- sce_finalKeep$cellType

new_cellType_ID <- c('Early RPC', 'Neurog.', 'RGC', 'Photo. Precurs.', 'AC', 'HC', 'Cone', 'Late RPC', 'Rod', 'MG', 'BC')
names(new_cellType_ID) <- levels(sce_finalKeep)
sce_finalKeep <- RenameIdents(sce_finalKeep, new_cellType_ID)
sce_finalKeep <- subset(sce_finalKeep, stage %in% c('E14', 'E18', 'P8', 'P14'))
sce_finalKeep$newCellType <- Idents(sce_finalKeep)

mouse_metadata <- sce_finalKeep@meta.data

#---

c2m_info <- readRDS('R_downstream/cross_species/homologous/c2m_info.rds')
c2m_info$seurat <- sce_finalKeep
saveRDS(c2m_info,'R_downstream/cross_species/c2m/regulon_info.rds')

#------
z2m_info <- readRDS('R_downstream/cross_species/homologous/z2m_info.rds')
z2m_info$seurat <- sce_finalKeep
saveRDS(z2m_info,'R_downstream/cross_species/z2m/regulon_info.rds')


#-------------------------------------add seurat object to the auc file of chicken------------------------------------------
chicken <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_filtered_processed_final.rds')

m2c_info <- readRDS('R_downstream/cross_species/homologous/m2c_info.rds')
chicken <- subset(chicken, cells = intersect(colnames(chicken), colnames(m2c_info$AUC)))

DefaultAssay(chicken) <- 'RNA'
chicken <- chicken %>% 
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("orig.ident",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

chicken_metadata <- chicken@meta.data
m2c_info$filtered_auc <- getAUC(m2c_info$AUC)[,rownames(chicken_metadata)]
m2c_info$seurat <- chicken
saveRDS(m2c_info,'R_downstream/cross_species/m2c/regulon_info.rds')

#---

z2c_info <- readRDS('R_downstream/cross_species/homologous/z2c_info.rds')
z2c_info$filtered_auc <- getAUC(z2c_info$AUC)[,rownames(chicken_metadata)]
z2c_info$seurat <- chicken
saveRDS(z2c_info,'R_downstream/cross_species/z2c/regulon_info.rds')


#----------------------------------------------------------------------------------------------------
#zebrafish
zebrafish_seurat <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/merge_all.rds')
zebrafish_metadata <- zebrafish_seurat@meta.data

m2z_info <- readRDS('R_downstream/cross_species/homologous/m2z_info.rds')
m2z_info$seurat <- zebrafish_seurat

saveRDS(m2z_info,'R_downstream/cross_species/m2z/regulon_info.rds')

c2z_info <- readRDS('R_downstream/cross_species/homologous/c2z_info.rds')
c2z_info$seurat <- zebrafish_seurat
saveRDS(c2z_info,'R_downstream/cross_species/c2z/regulon_info.rds')


#----------------------------------------------------------------------------------------------------




##### ------------------analysis of regulon from mouse
mouse_scaled_regulons <- list()
stage_names <- c('E14', 'E16', 'E18', 'P0', 'P2', 'P5', 'P8', 'P14')
for (stage in stage_names){

	file_name <- paste0('R_downstream/mouse/',stage, '/', stage,'_cellactivity.txt')
	temp_file <- read.table(file_name, sep = '\t')

	regulon_names <- str_split(rownames(temp_file), pattern = ' ',simplify = TRUE)[,1]
	rownames(temp_file) <- paste0(stage,'_',regulon_names)
	colnames(temp_file)[colnames(temp_file) == 'PR.Precur.'] <- 'PR Precur.'

	mouse_scaled_regulons[[stage]] <- temp_file
}

regulon_from_mouse <- tf_info$mouse[tf_info$mouse$stage %in%stage_names, ]
regulon_from_mouse$act_mouse <- ''

for (line in 1:dim(regulon_from_mouse)[1]){

	regulon_from_mouse[line,'act_mouse'] <- mouse_scaled_regulons[[regulon_from_mouse[line,]$stage]][regulon_from_mouse[line,]$stage_regulon,regulon_from_mouse[line,]$newCellType]
}

regulon_from_mouse <- regulon_from_mouse[regulon_from_mouse$act_mouse > 0,]

m2z_info <- readRDS('R_downstream/cross_species/m2z/regulon_info.rds')
m2z_auc <- getAUC(m2z_info$AUC)
m2z_metadata <- m2z_info$seurat@meta.data

m2z_metadata[m2z_metadata$cellType == 'PR',]$cellType <- 'PR Precur.'
m2z_scaled_activity <- calculateCellGroupsActivity(m2z_auc,m2z_metadata, 'cellType',scale = TRUE)

regulon_from_mouse <- regulon_from_mouse[regulon_from_mouse$stage_regulon %in% intersect(rownames(m2z_scaled_activity),regulon_from_mouse$stage_regulon),]
regulon_from_mouse$act_zebrafish <- ''

for (line in 1:dim(regulon_from_mouse)[1]){
	regulon_from_mouse[line,'act_zebrafish'] <- m2z_scaled_activity[regulon_from_mouse[line,]$stage_regulon,regulon_from_mouse[line,]$newCellType]
}


m2c_info <- readRDS('R_downstream/cross_species/m2c/regulon_info.rds')
m2c_auc <- getAUC(m2c_info$AUC)
m2c_metadata <- m2c_info$seurat@meta.data
m2c_metadata[m2c_metadata$cellType == 'PR Precurs.',]$cellType <- 'PR Precur.'
m2c_scaled_activity <- calculateCellGroupsActivity(m2c_auc,m2c_metadata, 'cellType',scale = TRUE)

regulon_from_mouse <- regulon_from_mouse[regulon_from_mouse$stage_regulon %in% intersect(rownames(m2c_scaled_activity),regulon_from_mouse$stage_regulon),]
regulon_from_mouse$act_chicken <- ''
regulon_from_mouse <- regulon_from_mouse[regulon_from_mouse$newCellType != 'RPC',]

for (line in 1:dim(regulon_from_mouse)[1]){
	regulon_from_mouse[line,'act_chicken'] <- m2c_scaled_activity[regulon_from_mouse[line,]$stage_regulon,regulon_from_mouse[line,]$newCellType]
}

write.table(regulon_from_mouse, 'R_downstream/cross_species/regulons_from_mouse.txt', sep = '\t', quote = FALSE)


##### ------------------analysis of regulon from chicken
chicken_scaled_regulons <- list()
stage_names <- c('E12', 'E18')
for (stage in stage_names){

	file_name <- paste0('R_downstream/chicken/',stage, '/', stage,'_cellactivity.txt')
	temp_file <- read.table(file_name, sep = '\t')

	regulon_names <- str_split(rownames(temp_file), pattern = ' ',simplify = TRUE)[,1]
	rownames(temp_file) <- paste0(stage,'_',regulon_names)
	colnames(temp_file)[colnames(temp_file) == 'PR.Precurs.'] <- 'PR Precur.'

	chicken_scaled_regulons[[stage]] <- temp_file
}


regulon_from_chicken <- tf_info$chicken[tf_info$chicken$stage %in%stage_names, ]
regulon_from_chicken$act_chicken <- ''

for (line in 1:dim(regulon_from_chicken)[1]){

	regulon_from_chicken[line,'act_chicken'] <- chicken_scaled_regulons[[regulon_from_chicken[line,]$stage]][regulon_from_chicken[line,]$stage_regulon,regulon_from_chicken[line,]$newCellType]
}

regulon_from_chicken <- regulon_from_chicken[regulon_from_chicken$act_chicken > 0,]


c2z_info <- readRDS('R_downstream/cross_species/c2z/regulon_info.rds')
c2z_auc <- getAUC(c2z_info$AUC)
c2z_metadata <- c2z_info$seurat@meta.data

c2z_metadata[c2z_metadata$cellType == 'PR',]$cellType <- 'PR Precur.'
c2z_scaled_activity <- calculateCellGroupsActivity(c2z_auc,c2z_metadata, 'cellType',scale = TRUE)

regulon_from_chicken <- regulon_from_chicken[regulon_from_chicken$stage_regulon %in% intersect(rownames(c2z_scaled_activity),regulon_from_chicken$stage_regulon),]
regulon_from_chicken$act_zebrafish <- ''

for (line in 1:dim(regulon_from_chicken)[1]){
	regulon_from_chicken[line,'act_zebrafish'] <- c2z_scaled_activity[regulon_from_chicken[line,]$stage_regulon,regulon_from_chicken[line,]$newCellType]
}

c2m_info <- readRDS('R_downstream/cross_species/c2m/regulon_info.rds')
c2m_auc <- getAUC(c2m_info$AUC)
c2m_metadata <- c2m_info$seurat@meta.data
c2m_metadata$newCellType <-  as.character(c2m_metadata$newCellType)
c2m_metadata[c2m_metadata$newCellType == 'Photo. Precurs.',]$newCellType <- 'PR Precur.'
c2m_scaled_activity <- calculateCellGroupsActivity(c2m_auc,c2m_metadata, 'newCellType',scale = TRUE)

regulon_from_chicken <- regulon_from_chicken[regulon_from_chicken$stage_regulon %in% intersect(rownames(c2m_scaled_activity),regulon_from_chicken$stage_regulon),]
regulon_from_chicken$act_mouse <- ''
for (line in 1:dim(regulon_from_chicken)[1]){
	regulon_from_chicken[line,'act_mouse'] <- c2m_scaled_activity[regulon_from_chicken[line,]$stage_regulon,regulon_from_chicken[line,]$newCellType]
}

write.table(regulon_from_chicken, 'R_downstream/cross_species/regulon_from_chicken.txt', sep = '\t', quote = FALSE)


##### ------------------analysis of regulon from zebrafish
zebrafish_scaled_regulons <- list()
stage_names <- c('36hpf', '48hpf', '72hpf', '14dpf')
stage_files <- c('hpf36', 'hpf48', 'hpf72', 'dpf14')
names(stage_files) <- stage_names
for (x in 1:4){

	file_name <- paste0('R_downstream/zebrafish/',stage_names[x], '/', stage_files[x],'_cellactivity.txt')
	temp_file <- read.table(file_name, sep = '\t')

	regulon_names <- str_split(rownames(temp_file), pattern = ' ',simplify = TRUE)[,1]
	rownames(temp_file) <- paste0(stage_names[x],'_',regulon_names)
	colnames(temp_file)[colnames(temp_file) == 'PR.Precurs.'] <- 'PR Precur.'

	zebrafish_scaled_regulons[[stage_files[x]]] <- temp_file
}

regulon_from_zebrafish <- tf_info$zebrafish[tf_info$zebrafish$stage %in%stage_names, ]
regulon_from_zebrafish$act_zebrafish <- ''

for (line in 1:dim(regulon_from_zebrafish)[1]){

	regulon_from_zebrafish[line,'act_zebrafish'] <- zebrafish_scaled_regulons[[stage_files[regulon_from_zebrafish[line,]$stage]]][regulon_from_zebrafish[line,]$stage_regulon,regulon_from_zebrafish[line,]$newCellType]
}

regulon_from_zebrafish <- regulon_from_zebrafish[regulon_from_zebrafish$act_zebrafish > 0,]

z2c_info <- readRDS('R_downstream/cross_species/z2c/regulon_info.rds')
z2c_auc <- getAUC(z2c_info$AUC)
z2c_metadata <- z2c_info$seurat@meta.data

z2c_metadata[z2c_metadata$cellType == 'PR Precurs.',]$cellType <- 'PR Precur.'
z2c_scaled_activity <- calculateCellGroupsActivity(z2c_auc,z2c_metadata, 'cellType',scale = TRUE)

regulon_from_zebrafish <- regulon_from_zebrafish[regulon_from_zebrafish$stage_regulon %in% intersect(rownames(z2c_scaled_activity),regulon_from_zebrafish$stage_regulon),]
regulon_from_zebrafish$act_chicken <- ''

regulon_from_zebrafish <- regulon_from_zebrafish[regulon_from_zebrafish$newCellType != 'RPC',]

for (line in 1:dim(regulon_from_zebrafish)[1]){
	regulon_from_zebrafish[line,'act_chicken'] <- z2c_scaled_activity[regulon_from_zebrafish[line,]$stage_regulon,regulon_from_zebrafish[line,]$newCellType]
}

z2m_info <- readRDS('R_downstream/cross_species/z2m/regulon_info.rds')
z2m_auc <- getAUC(z2m_info$AUC)
z2m_metadata <- z2m_info$seurat@meta.data
z2m_metadata$newCellType <-  as.character(z2m_metadata$newCellType)
z2m_metadata[z2m_metadata$newCellType == 'Photo. Precurs.',]$newCellType <- 'PR Precur.'
z2m_scaled_activity <- calculateCellGroupsActivity(z2m_auc,z2m_metadata, 'newCellType',scale = TRUE)

regulon_from_zebrafish <- regulon_from_zebrafish[regulon_from_zebrafish$stage_regulon %in% intersect(rownames(z2m_scaled_activity),regulon_from_zebrafish$stage_regulon),]
regulon_from_zebrafish$act_mouse <- ''

for (line in 1:dim(regulon_from_zebrafish)[1]){
	regulon_from_zebrafish[line,'act_mouse'] <- z2m_scaled_activity[regulon_from_zebrafish[line,]$stage_regulon,regulon_from_zebrafish[line,]$newCellType]
}

write.table(regulon_from_zebrafish, 'R_downstream/cross_species/regulon_from_zebrafish.txt', sep = '\t', quote = FALSE)

#--------------------------------------------------scatter plot -------------------------------------
regulon_from_mouse_high <- regulon_from_mouse[(regulon_from_mouse$act_zebrafish > 1.6) & (regulon_from_mouse$act_chicken > 1.6),]
p1 <-	ggplot(regulon_from_mouse,aes(x = act_chicken, y = act_zebrafish, color = newCellType)) + 
		geom_point(size = 1) + 
#		scale_x_continuous(limits = c(-1, 3)) + 
#		scale_y_continuous(limits = c(-1, 3)) + 
		theme_bw() +
		ggtitle("Regulons from Mouse") +
		labs(x = 'Scaled regulon activity in Chicken',
			   y = 'Scaled Regulon activity in Zebrafish',
			   title = 'Regulons from Mouse') +
		theme(panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(),
					axis.title = element_text(size = rel(1)),
					legend.margin=margin(t = 0, unit='cm'),
					legend.title=element_blank(),
					plot.title = element_text(face = "bold",hjust = 0.5)) +
		scale_color_manual(values = plotColor) + 
		geom_hline(aes(yintercept=1.6), linetype="dashed") + 
		geom_vline(aes(xintercept=1.6), linetype="dashed")

p2 <- p1 + geom_text_repel(aes(label = stage_regulon),regulon_from_mouse_high ,size = 4, min.segment.length = unit(0.1, "lines"))

ggsave('R_downstream/cross_species/regulons_from_mouse.pdf',p2,dpi = 200, units = 'cm', width = 15, height = 12)



regulon_from_chicken_high <- regulon_from_chicken[(regulon_from_chicken$act_mouse > 1.6) & (regulon_from_chicken$act_zebrafish > 1.6),]
regulon_from_chicken_high <- regulon_from_chicken_high[!duplicated(regulon_from_chicken_high$stage_regulon),]

p1 <-	ggplot(regulon_from_chicken,aes(x = act_mouse, y = act_zebrafish, color = newCellType)) + 
		geom_point(size = 1) + 
#		scale_x_continuous(limits = c(-1.5, 6)) + 
#		scale_y_continuous(limits = c(-1.5, 6)) + 
		theme_bw() +
		ggtitle("Regulons from Chicken") +
		labs(x = 'Scaled Regulon activity in Mouse',
			   y = 'Scaled Regulon activity in Zebrafish',
			   title = 'Regulons from Chicken') +
		theme(panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(),
					axis.title = element_text(size = rel(1)),
					legend.margin=margin(t = 0, unit='cm'),
					legend.title=element_blank(),
					plot.title = element_text(face = "bold",hjust = 0.5)) +
		scale_color_manual(values = plotColor) + 
		geom_hline(aes(yintercept=1.6), linetype="dashed") + 
		geom_vline(aes(xintercept=1.6), linetype="dashed")

p2 <- p1 + geom_text_repel(aes(label = stage_regulon),regulon_from_chicken_high ,size = 4, min.segment.length = unit(0.1, "lines"))

ggsave('R_downstream/cross_species/regulons_from_chicken.pdf',p2,dpi = 200, units = 'cm', width = 15, height = 12)



regulon_from_zebrafish_high <- regulon_from_zebrafish[(regulon_from_zebrafish$act_mouse > 1.6) & (regulon_from_zebrafish$act_chicken > 1.6),]

p1 <-	ggplot(regulon_from_zebrafish,aes(x = act_mouse, y = act_chicken, color = newCellType)) + 
    geom_point(size = 1) + 
    # scale_x_continuous(limits = c(-1.5, 5)) + 
    # scale_y_continuous(limits = c(-1.5, 5)) + 
    theme_bw() +
    labs(x = 'Scaled Regulon activity in Mouse',
         y = 'Scaled Regulon activity in Chicken',
         title = 'Regulons from Zebrafish') +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title = element_text(size = rel(1)),
          legend.margin=margin(t = 0, unit='cm'),
          legend.title=element_blank(),
          plot.title = element_text(face = "bold",hjust = 0.5)) +
    scale_color_manual(values = plotColor) + 
    geom_hline(aes(yintercept=1.6), linetype="dashed") + 
    geom_vline(aes(xintercept=1.6), linetype="dashed")



p2 <- p1 + geom_text_repel(aes(label = stage_regulon),regulon_from_zebrafish_high ,box.padding = 0.5,point.padding = 0.1, min.segment.length = 0.1, segment.color='black', max.overlaps = 50)

ggsave('R_downstream/cross_species/regulons_from_zebrafish.pdf',p2,dpi = 200, units = 'cm', width = 15, height = 12)

#--------------------------------------------------plot conserved regulons ------------------------------------- 

rownames(SymbolENSEMBL$chicken) <- SymbolENSEMBL$chicken$Ensembl
sce_finalKeep$newCellType <- Idents(sce_finalKeep)

#-----------E14_Lhx4
E14_lhx4_target_genes <- str_split(mouse_all_tf_df[mouse_all_tf_df$stage_regulon == 'E14_Lhx4',]$target_genes, pattern = ';')[[1]]
E14_lhx4_target_counts <- sce_finalKeep@assays$SCT@data[E14_lhx4_target_genes,names(sce_finalKeep$newCellType[sce_finalKeep$newCellType == 'Cone'])]
E14_lhx4_target_genes_top50 <- names(sort(rowSums(lhx4_target_counts), decreasing = TRUE)[1:50])

E14_lhx4_target_genes_c <- unique(na.omit(SymbolENSEMBL$chicken[m2c_regulons$E14_Lhx4,]$Gene_symbol))
E14_lhx4_target_count_c <- chicken@assays$SCT@data[intersect(E14_lhx4_target_genes_c, rownames(chicken@assays$SCT@data)),names(chicken$cellType[chicken$cellType == 'Cone'])]
E14_lhx4_target_genes_top50_c <- names(sort(rowSums(E14_lhx4_target_count_c), decreasing = TRUE)[1:50])
E14_lhx4_target_genes_top50_c <- Others2Mouse(E14_lhx4_target_genes_top50_c, 'chicken')

E14_lhx4_target_genes_z <- unique(na.omit(SymbolENSEMBL$zebrafish[m2z_regulons$E14_Lhx4,]$Gene_symbol))
E14_lhx4_target_count_z <- zebrafish@assays$SCT@data[intersect(E14_lhx4_target_genes_z, rownames(zebrafish@assays$SCT@data)),names(zebrafish$newCellType[zebrafish$newCellType == 'Cone'])]
E14_lhx4_target_genes_top50_z <- names(sort(rowSums(E14_lhx4_target_count_z), decreasing = TRUE)[1:50])
E14_lhx4_target_genes_top50_z <- Others2Mouse(E14_lhx4_target_genes_top50_z, 'zebrafish')

E14_lhx4_target_genes_intersect <- intersect(intersect(E14_lhx4_target_genes_top50_c$mouse, E14_lhx4_target_genes_top50_z$mouse), E14_lhx4_target_genes_top50)

p1 <- DotPlot(sce_finalKeep, features = E14_lhx4_target_genes_intersect) +  
	theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E14_Lhx4_m.pdf', p1, dpi = 200, units = 'cm', width = 21, height = 8)

p2 <- DotPlot(chicken, features = E14_lhx4_target_genes_top50_c[E14_lhx4_target_genes_top50_c$mouse %in% E14_lhx4_target_genes_intersect,]$chicken) +
	theme(axis.title.x= element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E14_Lhx4_c.pdf', p2, dpi = 200, units = 'cm', width = 21, height = 8)

p3 <- DotPlot(zebrafish, features = E14_lhx4_target_genes_top50_z[E14_lhx4_target_genes_top50_z$mouse %in% E14_lhx4_target_genes_intersect,]$zebrafish) + 
		theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E14_Lhx4_z.pdf', p3, dpi = 200, units = 'cm', width = 21, height = 8)


#-----------Ehpf72_mafba
hpf72_mafba_target_genes <- changeGeneNames(str_split(zebrafish_all_tf_df[zebrafish_all_tf_df$stage_regulon == '72hpf_mafba',]$target_genes, pattern = ';')[[1]],SymbolENSEMBL$zebrafish, in_type = 'ensembl')
p3 <- DotPlot(zebrafish, features = hpf72_mafba_target_genes) + 
		theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))
ggsave('R_downstream/cross_species/regulons/hpf72_mafba_z.pdf', p3, dpi = 200, units = 'cm', width = 21, height = 8)

hpf72_mafba_target_genes_c <- changeGeneNames(z2c_regulons$`72hpf_mafba`, SymbolENSEMBL$chicken,in_type = 'ensembl' )
p2 <- DotPlot(chicken, features = hpf72_mafba_target_genes_c) +
	theme(axis.title.x= element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))
ggsave('R_downstream/cross_species/regulons/hpf72_mafba_c.pdf', p2, dpi = 200, units = 'cm', width = 21, height = 8)

p1 <- DotPlot(sce_finalKeep, features = z2m_regulons$`72hpf_mafba`) +  
	theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/72hpf_mafba_m.pdf', p1, dpi = 200, units = 'cm', width = 23, height = 8)


#-----------E18_VAX1
E18_VAX1_target_genes_c <- str_split(chicken_all_tf_df[chicken_all_tf_df$stage_regulon == 'E18_VAX1',]$target_genes, pattern = ';')[[1]]
E18_VAX1_target_genes_c <- changeGeneNames(E18_VAX1_target_genes_c, SymbolENSEMBL$chicken, in_type = 'ensembl')
E18_VAX1_c_target_counts <- chicken@assays$RNA@counts[E18_VAX1_target_genes_c,names(chicken$cellType[chicken$cellType == 'Rod'])]
E18_VAX1_c_target_genes_top12 <- names(sort(rowSums(E18_VAX1_c_target_counts), decreasing = TRUE)[1:12])

p2 <- DotPlot(chicken, features = E18_VAX1_c_target_genes_top12) +
	theme(axis.title.x= element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E18_VAX1_c.pdf', p2, dpi = 200, units = 'cm', width = 21, height = 8)

E18_VAX1_target_genes_z <- changeGeneNames(c2z_regulons$E18_VAX1, SymbolENSEMBL$zebrafish, in_type = 'ensembl')
E18_VAX1_z_target_counts <- zebrafish@assays$SCT@data[intersect(E18_VAX1_target_genes_z, rownames(zebrafish@assays$SCT@data)),names(zebrafish$newCellType[zebrafish$newCellType == 'Rod'])]
E18_VAX1_z_target_genes_top <- names(sort(rowSums(E18_VAX1_z_target_counts), decreasing = TRUE)[1:10])

p3 <- DotPlot(zebrafish, features = E18_VAX1_z_target_genes_top) + 
		theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))
ggsave('R_downstream/cross_species/regulons/E18_VAX1_z.pdf', p3, dpi = 200, units = 'cm', width = 21, height = 8)


E18_VAX1_m_target_counts <- sce_finalKeep@assays$SCT@data[intersect(c2m_regulons$E18_VAX1,rownames(sce_finalKeep)),names(sce_finalKeep$newCellType[sce_finalKeep$newCellType == 'Rod'])]
E18_VAX1_m_target_genes_top <- names(sort(rowSums(E18_VAX1_m_target_counts), decreasing = TRUE)[1:12])

p1 <- DotPlot(sce_finalKeep, features = E18_VAX1_m_target_genes_top) +  
	theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E18_VAX1_m.pdf', p1, dpi = 200, units = 'cm', width = 23, height = 8)

#-----------E12_PBX1
E12_PBX1_target_genes <- changeGeneNames(str_split(chicken_all_tf_df[chicken_all_tf_df$stage_regulon == 'E12_PBX1',]$target_genes, pattern = ';')[[1]],SymbolENSEMBL$chicken, in_type = 'ensembl')

p2 <- DotPlot(chicken, features = E12_PBX1_target_genes) +
	theme(axis.title.x= element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E12_PBX1_c.pdf', p2, dpi = 200, units = 'cm', width = 21, height = 8)

p1 <- DotPlot(sce_finalKeep, features = c2m_regulons$E12_PBX1) +  
	theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/E12_PBX1_m.pdf', p1, dpi = 200, units = 'cm', width = 23, height = 8)

E12_PBX1_target_gene_z <- changeGeneNames(c2z_regulons$E12_PBX1,SymbolENSEMBL$zebrafish, in_type = 'ensembl')
p3 <- DotPlot(zebrafish, features = unique(E12_PBX1_target_gene_z)) + 
		theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))
ggsave('R_downstream/cross_species/regulons/E12_PBX1_z.pdf', p3, dpi = 200, units = 'cm', width = 21, height = 8)

## P14_Sox9
P14_Sox9_target_genes <- str_split(mouse_all_tf_df[mouse_all_tf_df$stage_regulon == 'P14_Sox9',]$target_genes, pattern = ';')[[1]]
P14_Sox9_target_counts <- sce_finalKeep@assays$SCT@data[P14_Sox9_target_genes,names(sce_finalKeep$newCellType[sce_finalKeep$newCellType == 'MG'])]
P14_Sox9_target_genes_top30 <- names(sort(rowSums(P14_Sox9_target_counts), decreasing = TRUE)[1:30])


P14_Sox9_target_genes_c <- unique(na.omit(SymbolENSEMBL$chicken[m2c_regulons$P14_Sox9,]$Gene_symbol))
P14_Sox9_target_count_c <- chicken@assays$SCT@data[intersect(P14_Sox9_target_genes_c, rownames(chicken@assays$SCT@data)),names(chicken$cellType[chicken$cellType == 'MG'])]
P14_Sox9_target_genes_top30_c <- names(sort(rowSums(P14_Sox9_target_count_c), decreasing = TRUE)[1:30])
P14_Sox9_target_genes_top30_c <- Others2Mouse(P14_Sox9_target_genes_top30_c, 'chicken')

P14_Sox9_target_genes_z <- unique(na.omit(SymbolENSEMBL$zebrafish[m2z_regulons$P14_Sox9,]$Gene_symbol))
P14_Sox9_target_count_z <- zebrafish@assays$SCT@data[intersect(P14_Sox9_target_genes_z, rownames(zebrafish@assays$SCT@data)),names(zebrafish$newCellType[zebrafish$newCellType == 'MG'])]
P14_Sox9_target_genes_top30_z <- names(sort(rowSums(P14_Sox9_target_count_z), decreasing = TRUE)[1:30])
P14_Sox9_target_genes_top30_z <- Others2Mouse(P14_Sox9_target_genes_top30_z, 'zebrafish')

P14_Sox9_target_genes_intersect <- intersect(intersect(P14_Sox9_target_genes_top30_c$mouse, P14_Sox9_target_genes_top30_z$mouse), P14_Sox9_target_genes_top30)

p1 <- DotPlot(sce_finalKeep, features = P14_Sox9_target_genes_intersect) +  
	theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/P14_Sox9_m.pdf', p1, dpi = 200, units = 'cm', width = 21, height = 8)


p2 <- DotPlot(chicken, features = P14_Sox9_target_genes_top30_c[P14_Sox9_target_genes_top30_c$mouse %in% P14_Sox9_target_genes_intersect,]$chicken) +
	theme(axis.title.x= element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/P14_Sox9_c.pdf', p2, dpi = 200, units = 'cm', width = 21, height = 8)

p3 <- DotPlot(zebrafish, features = P14_Sox9_target_genes_top30_z[P14_Sox9_target_genes_top30_z$mouse %in% P14_Sox9_target_genes_intersect,]$zebrafish) + 
		theme(axis.title.x= element_blank(), axis.title.y= element_blank(),axis.text.x = element_text(angle = 20))

ggsave('R_downstream/cross_species/regulons/P14_Sox9_z.pdf', p3, dpi = 200, units = 'cm', width = 21, height = 8)
