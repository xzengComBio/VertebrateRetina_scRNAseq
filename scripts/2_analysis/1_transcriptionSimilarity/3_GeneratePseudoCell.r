
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))

setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/')
set.seed(123)
#--------------------------------------------Function------------------------------------------------------------------


# Define the myPseudocell function with several arguments:
# count_matrix: Original single-cell expression matrix.
# metadata: Metadata that includes cell type information.
# species: A string, probably to annotate the pseudocells with the source species.
# pseudocell.size: The number of cells you want in each pseudocell.
# rep_times: The number of times to replicate each cell for sampling.


myPseudocell <- function(count_matrix, metadata, species, pseudocell.size = 10, rep_times=5){

	metadata$CellType <- droplevels(metadata$CellType)
	metadata$CellType <- as.factor(metadata$CellType)

	new_ids_list = list()

	# Loop through each unique cell type in the metadata
	for (i in 1:length(levels(metadata$CellType))) {

        # Identify the cells that belong to the current cluster/cell type
		cluster_id = levels(metadata$CellType)[i]
		cluster_cells <- rownames(metadata[metadata$CellType == cluster_id,])
		cluster_size <- length(cluster_cells)

		# Generate pseudocell IDs. Replicate the cells, then define pseudocell groups.		
		pseudo_ids <- floor(seq_along(rep(cluster_cells,rep_times))/pseudocell.size) 
		pseudo_ids <- paste0(cluster_id, "_Cell_", pseudo_ids)
		names(pseudo_ids) <- sample(cluster_cells,size=length(pseudo_ids), replace = TRUE)	
		new_ids_list[[i]] <- pseudo_ids		
	}

	# Combine all the new IDs across cell types
	new_ids <- unlist(new_ids_list)

    # Create a data frame mapping the new pseudocell IDs to the original cell IDs
	new_ids_df <- data.frame(new_id = new_ids, old_id = names(new_ids))
	new_ids_length <- table(new_ids_df$new_id)

	new_colnames <- new_ids_df$old_id
	new_count_matrix<-count_matrix[,as.character(new_colnames)]

	new_count_matrix <- t(new_count_matrix)
	new.data<-aggregate(list(new_count_matrix[,1:length(new_count_matrix[1,])]),
		list(name=new_ids_df[,1]),FUN=mean)
	rownames(new.data) <- new.data$name

	new.data<-new.data[,-1]
	new_ids_length <- as.matrix(new_ids_length)

    # Check if any pseudocells represent fewer cells than expected and handle those
	if (sum(new_ids_length<5) != 0){

		short <- which(new_ids_length<5)
		new_good_ids <- as.matrix(new_ids_length[-short,])
		result <- t(new.data)[,rownames(new_good_ids)]
	}else{
		result <- t(new.data)[,rownames(new_ids_length)]
	}

	# Add the species prefix to the pseudocell names
	colnames(result)<-paste(species,colnames(result),sep="_")
	rownames(result)<- rownames(count_matrix)

	return_result <- list(psuedocell_matrix=result, ids=new_ids_df)

    # Return a list containing the pseudocell expression matrix and the ID mapping
	return(return_result)
	
}


#--------------------------------------------Mouse------------------------------------------------------------------

mouse_matrix <- readRDS('other_species/Clark2019neuron/results/rds/Clark2019neuron_forIntegration.rds')

mouse_seurat <- readRDS('other_species/Clark2019neuron/results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')
mouse_seurat <- RenameCells(mouse_seurat, new.names = paste0("mouse-", colnames(mouse_seurat)))

Idents(mouse_seurat) <- mouse_seurat$cellType
new_cellType_ID <- c('ERPC', 'Neurog', 'RGC', 'PRPrecurs', 'AC', 'HC', 'Cone', 'LRPC', 'Rod', 'MG', 'BC')
names(new_cellType_ID) <- levels(mouse_seurat)
mouse_seurat <- RenameIdents(mouse_seurat, new_cellType_ID)

mouse_seurat$cellType <- Idents(mouse_seurat)
mouse_metadata <- mouse_seurat@meta.data

# Create a new data frame that contains the cell IDs and their associated cell types
mouse_cell_type_info <- data.frame(CellID = rownames(mouse_metadata),
								  CellType = mouse_metadata$cellType) 
rownames(mouse_cell_type_info) <- mouse_cell_type_info$CellID

keep_cells <- c()
# Loop through each unique cell type
for (x in 1:length(table(mouse_cell_type_info$CellType))){
	if (table(mouse_cell_type_info$CellType)[x] > 3000){
		temp_cells <- rownames(mouse_cell_type_info[mouse_cell_type_info$CellType == names(table(mouse_cell_type_info$CellType)[x]),])
		temp_select_cells <- sample(temp_cells, 3000)
	}
	else {

		# If the number of cells of a specific type is 3000 or less, keep all of them
		temp_select_cells <- rownames(mouse_cell_type_info[mouse_cell_type_info$CellType == names(table(mouse_cell_type_info$CellType)[x]),])
	}

	keep_cells <- c(keep_cells, temp_select_cells)
}

mouse_cell_type_info_selected <- mouse_cell_type_info[keep_cells,]

mouse_matrix_filtered <- mouse_matrix[,mouse_cell_type_info_selected$CellID]

mouse_matrix_filtered <- as.matrix(mouse_matrix_filtered)
mouse_cell_type_info_selected$CellType <- as.factor(mouse_cell_type_info_selected$CellType)

mouse_all_pseudoCell <- myPseudocell(mouse_matrix_filtered, mouse_cell_type_info_selected, species = 'mouse',pseudocell.size = 10, rep_times=2)

saveRDS(mouse_all_pseudoCell$psuedocell_matrix, 'Integration/result/pseudoCell/mouse_all_pseudoCelll10_2.rds')
write.table(mouse_all_pseudoCell$ids,'Integration/result/pseudoCell/mouse_all_pseudoCell10_2_cellid.txt',row.names = FALSE,quote = FALSE, col.names = FALSE, sep='\t')

#--------------------------------------------Chicken------------------------------------------------------------------
chicken <- readRDS('other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_forIntegration.rds')

chicken_all_seurat <- readRDS('other_species/Yamagata2021eLife/results/rds/Yamagata2021eLife_filtered_processed_final.rds')
chicken_all_seurat <- RenameCells(chicken_all_seurat, new.names = paste0("chicken-", colnames(chicken_all_seurat)))
chicken <- chicken[,colnames(chicken_all_seurat)]

chicken_all_metadata <- chicken_all_seurat@meta.data

# Create a new data frame that contains the cell IDs and their associated cell types
chicken_all_cell_type_info <- data.frame(CellID = rownames(chicken_all_metadata),
								  CellType = chicken_all_metadata$cellType) 
chicken_all_cell_type_info$CellType[chicken_all_cell_type_info$CellType == 'PR Precurs.'] <- 'PRPrecurs'
rownames(chicken_all_cell_type_info) <- chicken_all_cell_type_info$CellID

keep_cells <- c()
for (x in 1:length(table(chicken_all_cell_type_info$CellType))){
	if (table(chicken_all_cell_type_info$CellType)[x] > 3000){
		temp_cells <- rownames(chicken_all_cell_type_info[chicken_all_cell_type_info$CellType == names(table(chicken_all_cell_type_info$CellType)[x]),])
		temp_select_cells <- sample(temp_cells, 3000)
	}
	else {
		temp_select_cells <- rownames(chicken_all_cell_type_info[chicken_all_cell_type_info$CellType == names(table(chicken_all_cell_type_info$CellType)[x]),])
	}

	keep_cells <- c(keep_cells, temp_select_cells)
}

chicken_cell_type_info_selected <- chicken_all_cell_type_info[keep_cells,]

chicken_cell_type_info_selected$CellType <- as.factor(chicken_cell_type_info_selected$CellType)

chicken <- chicken[,chicken_cell_type_info_selected$CellID]
chicken <- as.matrix(chicken)

# Run myPseudo Cell
chicken_all_pseudoCell <- myPseudocell(chicken, chicken_cell_type_info_selected, species = 'chicken', pseudocell.size = 10, rep_times=2)

saveRDS(chicken_all_pseudoCell$psuedocell_matrix, 'Integration/result/pseudoCell/chicken_all_pseudoCell10_2.rds')
write.table(chicken_all_pseudoCell$ids,'Integration/result/pseudoCell/chicken_all_pseudoCell10_2_cellid.txt',row.names = FALSE,quote = FALSE, col.names = FALSE, sep='\t')


#--------------------------------------------Zebrafish------------------------------------------------------------------
zebrafish <- readRDS('other_species/Xu2020Development/results/rds/Xu2020Development_forIntegration.rds')

zebrafish_all_seurat <- readRDS('other_species/Xu2020Development/results/rds/merge_all.rds')
zebrafish_all_seurat <- RenameCells(zebrafish_all_seurat, new.names = paste0("zebrafish-", colnames(zebrafish_all_seurat)))


zebrafish_all_metadata <- zebrafish_all_seurat@meta.data

# Create a new data frame that contains the cell IDs and their associated cell types
zebrafish_all_cell_type_info <- data.frame(CellID = rownames(zebrafish_all_metadata),
								  CellType = zebrafish_all_metadata$newCellType)
zebrafish_all_cell_type_info$CellType[zebrafish_all_cell_type_info$CellType == 'PR'] <- 'PRPrecurs'

keep_cells <- c()
for (x in 1:length(table(zebrafish_all_cell_type_info$CellType))){
	if (table(zebrafish_all_cell_type_info$CellType)[x] > 3000){
		temp_cells <- rownames(zebrafish_all_cell_type_info[zebrafish_all_cell_type_info$CellType == names(table(zebrafish_all_cell_type_info$CellType)[x]),])
		temp_select_cells <- sample(temp_cells, 3000)
	}
	else {
		temp_select_cells <- rownames(zebrafish_all_cell_type_info[zebrafish_all_cell_type_info$CellType == names(table(zebrafish_all_cell_type_info$CellType)[x]),])
	}

	keep_cells <- c(keep_cells, temp_select_cells)
}

zebrafish_cell_type_info_selected <- zebrafish_all_cell_type_info[keep_cells,]

zebrafish_cell_type_info_selected$CellType <- as.factor(zebrafish_cell_type_info_selected$CellType)
rownames(zebrafish_cell_type_info_selected) <- zebrafish_cell_type_info_selected$CellID

zebrafish <- zebrafish[,intersect(zebrafish_cell_type_info_selected$CellID, colnames(zebrafish))]
zebrafish <- as.matrix(zebrafish)

zebrafish_cell_type_info_selected <- zebrafish_cell_type_info_selected[intersect(zebrafish_cell_type_info_selected$CellID, colnames(zebrafish)),]
zebrafish_all_pseudoCell <- myPseudocell(zebrafish, zebrafish_cell_type_info_selected, species = 'zebrafish',pseudocell.size = 10, rep_times=5)

saveRDS(zebrafish_all_pseudoCell$psuedocell_matrix, 'Integration/result/pseudoCell/zebrafish_all_pseudoCell10_5.rds')
write.table(zebrafish_all_pseudoCell$ids,'Integration/result/pseudoCell/zebrafish_all_pseudoCell10_5_cellid.txt',row.names = FALSE,quote = FALSE, col.names = FALSE, sep='\t')
