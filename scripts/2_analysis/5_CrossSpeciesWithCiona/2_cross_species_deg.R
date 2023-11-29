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


#function
calGeneSetActivity <- function(cells_rankings, gene_sets){

	regulon_info <- list()
	geneset_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank = nrow(cells_rankings)*0.05, nCores = 4)

	regulon_info$AUC <- geneset_AUC

	return(regulon_info)
}



retina_celltype <- c( 'RGC', 'AC', 'HC', 'Rod', 'MG', 'BC', 'Cone')
#------------------ Preprocessing of orth----------------------- 
orthogroups <- readRDS('~/Desktop/Research/ciona/orthofinder/Orthogroups/fourSpecies/ciona_orth_one.rds')
orthogenes <- read.table('~/Desktop/Research/ciona/orthofinder/Orthogroups/fourSpecies/gene_orthogroup_info.txt', sep = '\t', header = TRUE)


# read symbol ENSEMBL info
SymbolENSEMBL <- list()
SymbolENSEMBL$mouse <- parseGeneSymbolEnsemblDF('../../TF/zebrafishAndChicken/gtf/mouse/geneSymbolENSEMBL.txt', make_unique = TRUE)
SymbolENSEMBL$chicken <- read.table('../../TF/zebrafishAndChicken/gtf/chicken/geneSymbolENSEMBL.txt', sep = '\t')
SymbolENSEMBL$zebrafish <- read.table('../../TF/zebrafishAndChicken/gtf/zebrafish/geneSymbolENSEMBL.txt', sep = '\t')

colnames(SymbolENSEMBL$chicken) <- colnames(SymbolENSEMBL$mouse)
colnames(SymbolENSEMBL$zebrafish) <- colnames(SymbolENSEMBL$mouse)
#----------------------------------------- 


#------------------------------ Identify the marker genes for Ciona photosensing related cells-------------------------- 
larva_ner <- readRDS('~/Desktop/Research/ciona/singleCell/result/larve/larva_ner/larva_ner_processed.rds')
#ciona_count_ranking <- AUCell_buildRankings(ciona_count)

cellType_selected <- c('Rx+ aSV', 'Opsin1+ PTPRB+ aSV', 'Opsin1+ STUM+ aSV', 'Ci-VP+ pSV', 'GLGB+ pSV', 'coronet cells', 'pigment cells', 'eminence')
photo_related_cells <-  subset(larva_ner, cellType %in% cellType_selected)

photo_related_cells  <- photo_related_cells %>%
		 SCTransform(verbose = FALSE) %>% 
		 RunPCA(verbose = FALSE) %>%
		 RunCanek("sample",correctEmbeddings = TRUE, assay = 'SCT') %>%
		 FindNeighbors(reduction = 'canek', dims = 1:20, verbose = FALSE) %>%
		 RunUMAP(dims = 1:20,reduction = 'canek', verbose = FALSE)

ciona_metadata <- photo_related_cells@meta.data
ciona_cellGroup <- split(rownames(ciona_metadata), ciona_metadata[,'cellType'])

Idents(photo_related_cells) <- photo_related_cells$cellType
photo_related_markers <- FindAllMarkers(photo_related_cells, only.pos = TRUE)

photo_related_markers <- photo_related_markers[photo_related_markers$p_val_adj < 0.01,]
photo_related_markers <- photo_related_markers %>% group_by(cluster) %>% top_n(n=100, wt = avg_log2FC)

p1 <- myPlotUMAP_clusters(photo_related_cells, 'cellType') + 
			theme(axis.ticks.x = element_blank(),  
						axis.ticks.y = element_blank(), 
						axis.text.x  = element_blank(),  
						axis.text.y  = element_blank(),
						panel.border = element_blank(),
						axis.title= element_blank())

saveRDS(photo_related_cells, 'R_downstream/ciona/larva_ner/photo_cell.rds')
ggsave('R_downstream/ciona/larva_ner/photo_cellTypes.pdf',p1,width = 12.5, height = 8, units = 'cm', dpi = 200)

#------------------------------ Identify the marker genes for mouse -------------------------- 

sce_finalKeep <- readRDS('../../singleCell/other_species/Clark2019neuron/results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')

Idents(sce_finalKeep) <- sce_finalKeep$cellType

new_cellType_ID <- c('Early RPC', 'Neurog.', 'RGC', 'PR Precurs.', 'AC', 'HC', 'Cone', 'Late RPC', 'Rod', 'MG', 'BC')
names(new_cellType_ID) <- levels(sce_finalKeep)
sce_finalKeep <- RenameIdents(sce_finalKeep, new_cellType_ID)
sce_finalKeep <- PrepSCTFindMarkers(sce_finalKeep, assay = 'SCT')

mouse_markers <- FindAllMarkers(sce_finalKeep, only.pos = TRUE)
write.table(mouse_markers, 'R_downstream/cross_species/deg/mouse_all_deg.txt',sep = '\t', row.names = FALSE, quote = FALSE)

#------------------------------ Identify the marker genes for chicken -------------------------- 

chicken <- readRDS('../../singleCell/other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_filtered_processed_final.rds')

Idents(chicken) <- chicken$cellType
chicken_markers <-  FindAllMarkers(chicken, only.pos = TRUE)

write.table(chicken_markers, 'R_downstream/cross_species/deg/chicken_all_deg.txt',sep = '\t', row.names = FALSE, quote = FALSE)

#------------------------------ Identify the marker genes for zebrafish -------------------------- 

zebrafish <- readRDS('../../singleCell/other_species/Xu2020Development/results/rds/merge_all.rds')
zebrafish_markers <-  FindAllMarkers(zebrafish, only.pos = TRUE)

write.table(zebrafish_markers, 'R_downstream/cross_species/deg/zebrafish_all_deg.txt',sep = '\t', row.names = FALSE, quote = FALSE)

#-------------------------------- Identify the marker genes for pseudo-cell -------------------------------- 
merge_sce <- readRDS('../../singleCell/Integration/result/merge_sce.rds')
Idents(merge_sce) <- merge_sce$RNA_snn_res.0.001

merge_all_markers <- FindAllMarkers(merge_sce, only.pos = TRUE)
write.table(merge_all_markers,'R_downstream/cross_species/deg/all_deg.txt',sep = '\t', row.names = FALSE, quote = FALSE)


#-------------read the orth_df generated by ENSEMBL and NCBI------------------ 
merge_all_markers <- merge_all_markers[merge_all_markers$cluster != 0,]

orth_df <- read.table('../../singleCell/Integration/result/ensembl_orth.txt', sep = '\t',header = TRUE)
rownames(orth_df) <- orth_df$mouse_ensembl
rownames(SymbolENSEMBL$mouse) <- SymbolENSEMBL$mouse$Gene_symbol
mouse_markers$mouse_ensembl <- SymbolENSEMBL$mouse[mouse_markers$gene, ]$Ensembl
mouse_markers$orth <- orth_df[mouse_markers$mouse_ensembl,]$orth
mouse_markers <- na.omit(mouse_markers)

rownames(orth_df) <- orth_df$chicken_ensembl
SymbolENSEMBL$chicken <- SymbolENSEMBL$chicken[!duplicated(SymbolENSEMBL$chicken$Gene_symbol),]
rownames(SymbolENSEMBL$chicken) <- SymbolENSEMBL$chicken$Gene_symbol
chicken_markers$chicken_ensembl <- SymbolENSEMBL$chicken[chicken_markers$gene, ]$Ensembl
chicken_markers$orth <- orth_df[chicken_markers$chicken_ensembl,]$orth
chicken_markers <- na.omit(chicken_markers)

rownames(orth_df) <- orth_df$zebrafish_ensembl
SymbolENSEMBL$zebrafish <- SymbolENSEMBL$zebrafish[!duplicated(SymbolENSEMBL$zebrafish$Gene_symbol),]
rownames(SymbolENSEMBL$zebrafish) <- SymbolENSEMBL$zebrafish$Gene_symbol
zebrafish_markers$zebrafish_ensembl <- SymbolENSEMBL$zebrafish[zebrafish_markers$gene, ]$Ensembl
zebrafish_markers$orth <- orth_df[zebrafish_markers$zebrafish_ensembl,]$orth
zebrafish_markers <- na.omit(zebrafish_markers)


#----------------------take the intersection by homologous genes ---------------
cell_type <- c('AC', 'BC', 'HC', 'RGC', 'Cone', 'Rod')

for (x in cell_type){
	mouse_marker_temp <- mouse_markers[mouse_markers$cluster == x,c(7,8,9)]
	colnames(mouse_marker_temp) <- c('mouse_symbol', 'mouse_ensembl', 'orth')
	chicken_marker_temp <- chicken_markers[chicken_markers$cluster == x,c(7,8,9)]
	colnames(chicken_marker_temp) <- c('chicken_symbol', 'chicken_ensembl', 'orth')
	zebrafish_marker_temp <- zebrafish_markers[zebrafish_markers$cluster == x,c(7,8,9)]
	colnames(zebrafish_marker_temp) <- c('zebrafish_symbol', 'zebrafish_ensembl', 'orth')

	mouse_chicken_temp <- merge(mouse_marker_temp, chicken_marker_temp, by='orth', all = FALSE)
	all_temp <- merge(mouse_chicken_temp, zebrafish_marker_temp, by= 'orth', all = FALSE)
	all_temp$cellType <- x
	if(!(exists("all_markers"))){
		all_markers <- all_temp
	}else{
		all_markers <- rbind(all_markers, all_temp)
	}

}

write.table(all_markers, 'R_downstream/cross_species/deg/all_deg_intersect.txt',sep = '\t', row.names = FALSE, quote = FALSE)


#----------------------keep the degs ---------------
merge_all_markers <- merge_all_markers[merge_all_markers$gene %in% union(merge_var.genes$hvg,all_markers$orth),]

rownames(orth_df) <- orth_df$orth
merge_all_markers$mouse_ensembl <- orth_df[merge_all_markers$gene, ]$mouse_ensembl
merge_all_markers$chicken_ensembl <- orth_df[merge_all_markers$gene, ]$chicken_ensembl
merge_all_markers$zebrafish_ensembl <- orth_df[merge_all_markers$gene, ]$zebrafish_ensembl

rownames(SymbolENSEMBL$mouse) <- SymbolENSEMBL$mouse$Ensembl
rownames(SymbolENSEMBL$chicken) <- SymbolENSEMBL$chicken$Ensembl
rownames(SymbolENSEMBL$zebrafish) <- SymbolENSEMBL$zebrafish$Ensembl

merge_all_markers$mouse_symbol <- SymbolENSEMBL$mouse[merge_all_markers$mouse_ensembl,]$Gene_symbol
merge_all_markers$chicken_symbol <- SymbolENSEMBL$chicken[merge_all_markers$chicken_ensembl,]$Gene_symbol
merge_all_markers$zebrafish_symbol <- SymbolENSEMBL$zebrafish[merge_all_markers$zebrafish_ensembl,]$Gene_symbol
write.table(merge_all_markers, 'R_downstream/cross_species/deg/all_deg_filtered.txt',sep = '\t', row.names = FALSE, quote = FALSE)

#----------------------read the ortholog information generated by orthoFinder ---------------
orth_genes <- read.table('../../orthofinder/Orthogroups/fourSpecies/gene_orthogroup_info.txt', sep = '\t', header = TRUE)
OrthogroupsGeneCount <- read.table('../../orthofinder/Orthogroups/fourSpecies/ciona_orth_less4.txt', sep = '\t', header = TRUE)

orth_genes <- orth_genes[orth_genes$Orthogroups %in% rownames(OrthogroupsGeneCount),]

orth_genes$species <- 'Ciona'
orth_genes[startsWith(orth_genes$Ensembl, 'ENSMUSG'),]$species <- 'Mouse'
orth_genes[startsWith(orth_genes$Ensembl, 'ENSGALG'),]$species <- 'Chicken'
orth_genes[startsWith(orth_genes$Ensembl, 'ENSDARG'),]$species <- 'Zebrafish'


merge_all_markers_same_orth <- merge_all_markers %>%
  left_join(orth_genes %>% filter(species == "Mouse"), by = c("mouse_ensembl" = "Ensembl")) %>%
  left_join(orth_genes %>% filter(species == "Chicken"), by = c("chicken_ensembl" = "Ensembl"), suffix = c("_mouse", "_chicken")) %>%
  left_join(orth_genes %>% filter(species == "Zebrafish"), by = c("zebrafish_ensembl" = "Ensembl"), suffix = c("", "_zebrafish"))

merge_all_markers_same_orth <- merge_all_markers_same_orth %>%
    filter(Orthogroups == Orthogroups_mouse & Orthogroups == Orthogroups_chicken)


orth_genes_ciona <- orth_genes[(orth_genes$species == 'Ciona') & (orth_genes$Orthogroups %in% (unique(merge_all_markers_same_orth$Orthogroups))),]
orth_genes_ciona$Ensembl <- paste0('KH2013:',orth_genes_ciona$Ensembl)
orth_genes_ciona_keep <- orth_genes_ciona[orth_genes_ciona$Ensembl %in% photo_related_markers$gene, ]

rownames(orth_genes_ciona_keep) <- orth_genes_ciona_keep$Orthogroups

merge_all_markers_same_orth <- merge_all_markers_same_orth[merge_all_markers_same_orth$Orthogroups %in% orth_genes_ciona_keep$Orthogroups, ]

merge_all_markers_same_orth$ciona <- orth_genes_ciona_keep[merge_all_markers_same_orth$Orthogroups,]$Ensembl


merge_all_markers_same_orth$ciona_cluster <- ''
for (x in 1: nrow(merge_all_markers_same_orth)){
  
  ciona_gene <- merge_all_markers_same_orth[x,'ciona']
  ciona_clusters <- photo_related_markers[photo_related_markers$gene %in% ciona_gene, ]$cluster
  ciona_clusters <- paste0(ciona_clusters, collapse = ';')
  merge_all_markers_same_orth[x,'ciona_cluster'] <- ciona_clusters

}

merge_all_markers_same_orth <- merge_all_markers_same_orth[,c(6,7,8,9,10,11,12,13,17,20,21)]


# for plotting
merge_sce <- readRDS('../../singleCell/Integration/result/merge_sce.rds')
Idents(merge_sce) <- merge_sce$RNA_snn_res.0.001
new_cellType <- c('RPC', 'PR', 'AC', 'RGC', 'BC', 'HC')
names(new_cellType) <- levels(merge_sce)
merge_sce <- RenameIdents(merge_sce, new_cellType)

merge_sce$all_cellType <- Idents(merge_sce)
merge_all_markers <- merge_all_markers[merge_all_markers$pct.1 > 0.8,]
merge_all_markers <- merge_all_markers[merge_all_markers$cluster != 'RPC',]
merge_all_markers_top10 <- merge_all_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)

merge_sce_onrpc <- subset(merge_sce, all_cellType != 'RPC')
merge_sce_matrix <- merge_sce_onrpc@assays$RNA@data[merge_all_markers_top10$gene,]
merge_sce_matrix <- as.matrix(merge_sce_matrix)

rownames(merge_sce_matrix) <- paste0(merge_all_markers_top10$mouse_symbol, '-homolog')

merge_sce_onrpc_metadata <- merge_sce_onrpc@meta.data
merge_sce_onrpc_metadata$all_cellType <- droplevels(merge_sce_onrpc_metadata$all_cellType)
merge_sce_onrpc_metadata$species <- merge_sce_onrpc_metadata$orig.ident

cellType_cols <- plotColor[c(4,5,3,11,6)]
names(cellType_cols) <-new_cellType[-1]
species_col <- plotColor[c(38,39,40)]
names(species_col) <- c('mouse', 'chicken', 'zebrafish')


top_anno <- HeatmapAnnotation(df = merge_sce_onrpc_metadata[,c( 'all_cellType','species')],
																	col = list(species = species_col,
																							all_cellType = cellType_cols))



column_splits <- factor(merge_sce_onrpc_metadata$all_cellType, levels = new_cellType[-1])

pdf('R_downstream/cross_species/deg/top10_deg.pdf', width = 14,height = 8)
cys_heatmap <- 	Heatmap(merge_sce_matrix,
												name="Expression Level",
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
                        row_names_side='left',
                        border = TRUE)

draw(cys_heatmap)
dev.off()

vertebrate_cellType <- c('PR', 'AC', 'RGC', 'BC', 'HC')

sankey_df <- data.frame(vertebrate = character(0),
												ciona = character(0))

for (x in 1: nrow(merge_all_markers_same_orth)){
  
  ciona_celltype <- merge_all_markers_same_orth[x,'ciona_cluster']
  ciona_celltype <- str_split(ciona_celltype, pattern = ';')[[1]]

  vertebrate <- vertebrate_cellType[merge_all_markers_same_orth[x,'cluster']]
  vertebrate <- rep(vertebrate, length(ciona_celltype))

  temp_df <- data.frame(vertebrate = vertebrate,
  											ciona = ciona_celltype)

  sankey_df <- rbind(sankey_df, temp_df)

}

sankey_df <- sankey_df %>% make_long(ciona, vertebrate)

plot_sankey_color <- plotColor[c(4,5,3,11,6,18:26)]
names(plot_sankey_color) <- c(new_cellType[-1],cellType_selected) 

p1 <- ggplot(sankey_df, aes(x = x, next_x = next_x,
                      node = node, next_node = next_node,
                      fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = .6) +
    geom_sankey_text(size = 3, color = "black", hjust = 0) +
    scale_fill_manual(values = plot_sankey_color)+
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5))

ggsave('R_downstream/cross_species//photo_cellTypes.pdf',p1,width = 12, height = 8, units = 'cm', dpi = 200)


my_gene_plot <- function(sce, feature, raster = TRUE, order = TRUE, pt.size = 1){
    p1 <- FeaturePlot(sce, features = feature, reduction = 'umap', raster = raster, order = order, pt.size = pt.size) + 
        theme_bw() +
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x  = element_blank(),
              axis.text.y  = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              legend.margin=margin(t = 0, unit='cm')) + 
        NoLegend()
    
    return(p1)
    
}

for (x in 1:nrow(merge_all_markers_same_orth)){
	p1 <- my_gene_plot(sce_finalKeep, merge_all_markers_same_orth[x,'mouse_symbol'])
	p2 <- my_gene_plot(chicken, merge_all_markers_same_orth[x,'chicken_symbol'])
	p3 <- my_gene_plot(zebrafish, merge_all_markers_same_orth[x,'zebrafish_symbol'])
	p4 <- my_gene_plot(photo_related_cells, merge_all_markers_same_orth[x,'ciona'], raster = FALSE, order = FALSE, pt.size = 0.3)

	plt_figure <- p1 | p2 | p3 | p4
	existing_files <- list.files('R_downstream/cross_species/comparedCIona/')
	file_name = paste0(merge_all_markers_same_orth[x,'mouse_symbol'], '.pdf')
	if (file_name %in% existing_files) {
		timestamp <- Sys.time() %>% as.character() %>% gsub("[ :]", "_", .)
		file_name <- paste0("unique_", timestamp, "_", file_name)
		print(paste("File already exists. Unique file name will be:", file_name))
	}else{
		print("File does not exist. You can proceed.")
	}
	ggsave(paste0('R_downstream/cross_species/comparedCIona/', file_name),plt_figure,width = 22, height = 6, units = 'cm', dpi = 200)
}
