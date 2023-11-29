library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(ggvenn)
library(homologene)
library(biomaRt)
library(xlsx)

setwd('~/Desktop/Research/ciona/SCENIC/results')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')

#------------------------------------------function------------------------------------------------

plot_three_features <- function(mouse_feature, chicken_feature, zebrafish_feature, folder, file_name){

	p1 <- FeaturePlot(sce_finalKeep, features = mouse_feature,reduction='umap',raster = TRUE, order = TRUE) + 
    	  theme_bw() +
    	  theme(panel.grid.major=element_blank(),
          		panel.grid.minor=element_blank(), 
          		axis.text.x=element_blank(),
          		axis.ticks.x=element_blank(),
          		axis.title.x=element_blank(),
          		axis.ticks.y=element_blank(),
          		axis.text.y=element_blank(),
          		axis.title.y=element_blank(),  
          		plot.title = element_text(hjust = 0.5,colour = "black")) +
          		NoLegend() +
    ggtitle(mouse_feature)

    p2 <- FeaturePlot(chicken, features = chicken_feature,reduction='umap',raster = TRUE, order = TRUE) + 
    	  theme_bw() +
    	  theme(panel.grid.major=element_blank(),
          		panel.grid.minor=element_blank(), 
          		axis.text.x=element_blank(),
          		axis.title.x=element_blank(),
          		axis.ticks.x=element_blank(),
          		axis.ticks.y=element_blank(),
          		axis.text.y=element_blank(),
          		axis.title.y=element_blank(), 
          		plot.title = element_text(hjust = 0.5,colour = "black")) + 
    	  NoLegend()  +
    ggtitle(chicken_feature)

    p3 <- FeaturePlot(zebrafish, features = zebrafish_feature,reduction='umap',raster = TRUE, order = TRUE) + 
    	  theme_bw() +
    	  theme(panel.grid.major=element_blank(),
          		panel.grid.minor=element_blank(), 
          		axis.text.x=element_blank(),
          		axis.ticks.x=element_blank(),
          		axis.title.x=element_blank(), 
          		axis.ticks.y=element_blank(),
          		axis.text.y=element_blank(), 
          		axis.title.y=element_blank(), 
          		plot.title = element_text(hjust = 0.5,colour = "black")) + 
    	  NoLegend() + 
    ggtitle(zebrafish_feature)

    all_plot <- p1 | p2 | p3
    ggsave(paste0('R_downstream/cross_species/three_features/',folder, '/',file_name,'.pdf'), all_plot, dpi = 200, units = 'cm', width = 16, height = 6)
    return(all_plot)
}


#------------------------------------------get the regulon tf names------------------------------------------------

# Load the regulon info for mouse
mouse_rss <- list()
mouse_stages <- c('E11', 'E12', 'E14', 'E16', 'E18', 'P0', 'P2', 'P5', 'P8', 'P14')

for (x in mouse_stages){
	file_names <- 'rss.txt'
	rss <-  read.table(paste0('R_downstream/mouse/',x,'/rss.txt'), header = TRUE,sep = '\t')
	mouse_rss[[x]] <- rss
} 


mouse_E11_tf <- mouse_rss$E11$Early.RPC
mouse_E11_cellType <- rep('RPC', length(mouse_E11_tf))


mouse_E12_tf <- mouse_rss$E12$Early.RPC
mouse_E12_cellType <- rep('RPC', length(mouse_E12_tf))


mouse_E14_tf  <- c(mouse_rss$E14$Early.RPC[-c(10,16,17,19,20)], # RPC
				  mouse_rss$E14$RGC[-c(6,12,13,14,16,17,19,20)],
				  mouse_rss$E14$Neurog.[c(1:4,6,7,9,18,19)],
				  mouse_rss$E14$PR.Precur.[c(3,4,7,8,17)],
				  mouse_rss$E14$HC[c(1,2,5,12,15)],
				  mouse_rss$E14$Cone[c(2,3,4,7,19)])

mouse_E14_cellType <- c(rep('Early RPC', 15), rep('RGC', 14),rep('Neurog.',9),rep('PR Precur.', 5), rep('HC',5), rep('Cone', 5))


mouse_E16_tf <- c(mouse_rss$E16$HC[c(5,10,20)],#HC
				  mouse_rss$E16$Early.RPC[c(1:6,8,9)],
				  mouse_rss$E16$RGC[c(1:8,16)],
				  mouse_rss$E16$Neurog.[c(1,2,4,5,7,11,12)],
				  mouse_rss$E16$PR.Precur.[c(3,4,8,9,14)]) # PR

mouse_E16_cellType <- c(rep('HC', 3),rep('Early RPC', 8),rep('RGC', 9),rep('Neurog.', 7), rep('PR Precur.', 5))

mouse_E18_tf <- c(mouse_rss$E18$GABA.AC[c(1:5,9,11,14,20)],
				  mouse_rss$E18$Glycine.AC[c(1:6,8:11,13:15,17)],#AC
				  mouse_rss$E18$Cone[c(2,4:7,12,18)],#Cone
				  mouse_rss$E18$Late.RPC[c(1:8,10:14,16,18,20)],
				  mouse_rss$E18$Neurog.[c(1:3,5:9,12,16,17)],
				  mouse_rss$E18$PR.Precur.[c(1,2,4,6,7,13,14,16)],
				  mouse_rss$E18$RGC[c(1,3:5,8,14,16,20)],
				  mouse_rss$E18$Rod[c(1,2,11)]) # PR

mouse_E18_cellType <- c(rep('GABA AC', 9),rep('Glycine AC', 14),rep('Cone', 7),rep('Late RPC', 16), rep('Neurog.', 11),
						rep('PR Precur.', 8),rep('RGC', 8),rep('Rod', 3))


mouse_P0_tf <-  c(mouse_rss$P0$GABA.AC[c(2,3,4,11,14,16,17)], # AC
				  mouse_rss$P0$Glycine.AC[c(1:8,14)], #Cone
				  mouse_rss$P0$Late.RPC[c(1:5,7,8,10,12,16,17,18)], # R
				  mouse_rss$P0$PR.Precur[c(1,2,5,6,7,8,11,13,15,19)],
				  mouse_rss$P0$Cone[c(1,2,6,8,9,13,14,15,19)],
				  mouse_rss$P0$Neurog.[c(1:6,9)],
				  mouse_rss$P0$Rod[c(2,4,5,6,7,17,18)])

mouse_P0_cellType <- c(rep('GABA AC', 7),rep('Glycine AC', 9),rep('Late RPC', 12),
						rep('PR Precur.', 10), rep('Cone', 9), rep('Neurog.', 7),rep('Rod', 7))


mouse_P2_tf <-  c(mouse_rss$P2$GABA.AC[c(2,3,4,6,14,20)],
				  mouse_rss$P2$Glycine.AC[c(1:10,13,14,19)], # AC
				  mouse_rss$P2$Cone[c(5,7,10,16,18,20)], #Cone
				  mouse_rss$P2$Late.RPC[c(1:5,8,10,12:16,20)], # R
				  mouse_rss$P2$PR.Precur[c(1:3,8:13,16,17,20)],
				  mouse_rss$P2$Neurog.[c(3:5,8,9,11,12,14,16,17)],
				  mouse_rss$P2$Rod[c(1,3:4,9,13,15,16,17)]) 

mouse_P2_cellType <- c(rep('GABA AC', 6),rep('Glycine AC', 13),rep('Cone', 6),rep('Late RPC', 13),
						rep('PR Precur.', 12), rep('Neurog.', 10),rep('Rod', 8))

mouse_P5_tf <-  c(mouse_rss$P5$AC[c(1,4,5,9,11,12,13,19)], # AC
				  mouse_rss$P5$BC[c(7)], #BC
				  mouse_rss$P5$Cone[c(2:4,6,8:10,13,14,17)], # R
				  mouse_rss$P5$Late.RPC[c(1,3,4,5,6,9,10,15,17,18)],
				  mouse_rss$P5$MG[c(1:8,12,19)],
				  mouse_rss$P5$PR.Precur.[c(2:4)],
				  mouse_rss$P5$Rod[c(1:9,13,14,16)]) 

mouse_P5_cellType <- c(rep('AC', 8),rep('BC', 1),rep('Cone', 10),rep('Late RPC', 10),rep('MG', 11),
						rep('PR Precur.', 3),rep('Rod', 12))


mouse_P8_tf <-  c(mouse_rss$P8$ON.BC[c(2,3,12,13,17:19)],
				  mouse_rss$P8$OFF.BC[c(3)], 
				  mouse_rss$P8$Rod.BC[c(4,6,12:14)],  # BC
				  mouse_rss$P8$Cone[c(1,13,14,15)], #Cone
				  mouse_rss$P8$Rod[c(1,2,5,7,8,10)]) 

mouse_P8_cellType <- c(rep('ON BC', 7),rep('OFF BC', 1),rep('Rod BC', 5),rep('Cone', 4),rep('Rod', 6))

mouse_P14_tf <-  c(mouse_rss$P14$Rod.BC[c(2:4,7)], # BC
				  mouse_rss$P14$Cone[c(2,8,9,11,12,14,16,17,20)], #Cone
				  mouse_rss$P14$MG[c(1:7,9:11,12,14,15,18)], # MG
				  mouse_rss$P14$Rod[c(1:3,6,7)])

mouse_P14_cellType <- c(rep('BC', 4),rep('Cone', 9),
					 rep('MG', 14),rep('Rod', 5))

mouse_stage_rep <- c(rep('E11',length(mouse_E11_tf)),
				 rep('E12',length(mouse_E12_tf)),
				 rep('E14',length(mouse_E14_tf)),
				 rep('E16',length(mouse_E16_tf)),
				 rep('E18',length(mouse_E18_tf)),
				 rep('P0',length(mouse_P0_tf)),
				 rep('P2',length(mouse_P2_tf)),
				 rep('P5',length(mouse_P5_tf)),
				 rep('P8',length(mouse_P8_tf)),
				 rep('P14',length(mouse_P14_tf)))

mouse_all_tf <- c(mouse_E11_tf, mouse_E12_tf, mouse_E14_tf, mouse_E16_tf,mouse_E18_tf,
				  mouse_P0_tf, mouse_P2_tf, mouse_P5_tf, mouse_P8_tf, mouse_P14_tf)

mouse_all_cellType <- c(mouse_E11_cellType, mouse_E12_cellType,mouse_E14_cellType,mouse_E16_cellType,mouse_E18_cellType,
						mouse_P0_cellType, mouse_P2_cellType, mouse_P5_cellType, mouse_P8_cellType, mouse_P14_cellType)
mouse_tf_df <- data.frame(stage = mouse_stage_rep,
						  regulon = mouse_all_tf,
						  tf_name = str_split(mouse_all_tf, pattern = ' ', simplify = TRUE)[,1],
						  cellType = mouse_all_cellType)


# Load the regulon info for chicken
chicken_rss <- list()
chicken_stages <- c('E12', 'E16', 'E18')

for (x in chicken_stages){
	file_names <- 'rss.txt'
	rss <-  read.table(paste0('R_downstream/chicken/',x,'/rss.txt'), header = TRUE,sep = '\t')
	chicken_rss[[x]] <- rss
} 

chicken_E12_tf <-  c(chicken_rss$E12$AC[c(1:5,8)], # RPC X 0
				  chicken_rss$E12$RGC[c(1,4,5,7)], #RPC X1
				  chicken_rss$E12$BC[c(6,8)],
				  chicken_rss$E12$HC[c(1,2)],
				  chicken_rss$E12$MG[c(1,2,4,6,7,10)],
				  chicken_rss$E12$PR.Precurs.[c(3)])

chicken_E12_cellType <- c(rep('AC', 6),rep('RGC', 4),rep('BC', 2),rep('HC', 2),
					 rep('MG', 6),rep('PR Precurs.', 1))

chicken_E16_tf <- chicken_rss$E16$RGC
chicken_E16_cellType <- rep('RGC', length(chicken_E16_tf))

chicken_E18_tf <- c(chicken_rss$E18$Blue.Cone[c(1,2,5)],
				  chicken_rss$E18$Dev..Cone[c(2,3)],
				  chicken_rss$E18$Double.Cone[c(1,3,4)],
				  chicken_rss$E18$GABA.AC[c(1:4,7)],
				  chicken_rss$E18$Glycine.AC[c(1:4,9)],
				  chicken_rss$E18$Green.Cone[c(1,3,4,5)],
				  chicken_rss$E18$MG[c(3,4,5,6,9)],
				  chicken_rss$E18$OFF.BC[c(1,5)],
				  chicken_rss$E18$ON.BC[c(5,7)],
				  chicken_rss$E18$Red.Cone[c(1,2,3)],
				  chicken_rss$E18$Rod[c(1,4,6)],
				  chicken_rss$E18$Rod.BC[c(1)],
				  chicken_rss$E18$T1.HC[c(1,2)],
				  chicken_rss$E18$T2.HC[c(1,2)],
				  chicken_rss$E18$Violet.Cone[c(1,3,4)])

chicken_E18_cellType <- rep(str_replace_all(colnames(chicken_rss$E18)[c(15,16,1,8,5,12,4,2,3,6,11,10,7,13,9,14)], '\\.', ' '), c(3,2,3,5,5,4,5,2,2,3,3,1,2,2,0,3))
#[1] "IKZF3 (7g)"   "IRF8 (52g)"   "E2F4 (5g)"    "PRRX2 (11g)"  "NR1H4 (118g)"

chicken_all_cellType <- c(chicken_E12_cellType, chicken_E16_cellType, chicken_E18_cellType)

chicken_stage_rep <- c(rep('E12',length(chicken_E12_tf)),
				 rep('E16',length(chicken_E16_tf)),
				 rep('E18',length(chicken_E18_tf)))

chicken_all_tf <- c(chicken_E12_tf, chicken_E16_tf,chicken_E18_tf)
chicken_tf_df <- data.frame(stage = chicken_stage_rep,
						  regulon = chicken_all_tf,
						  tf_name = str_split(chicken_all_tf, pattern = ' ', simplify = TRUE)[,1],
						  cellType = chicken_all_cellType)


# Load the regulon info for zebrafish

zebrafish_rss <- list()
zebrafish_stages <- c('24hpf', '36hpf', '48hpf', '72hpf', '14dpf')

for (x in zebrafish_stages){
	file_names <- 'rss.txt'
	rss <-  read.table(paste0('R_downstream/zebrafish/',x,'/rss.txt'), header = TRUE,sep = '\t')
	zebrafish_rss[[x]] <- rss
} 

zebrafish_24hpf_tf  <- c(zebrafish_rss$`24hpf`$X1[c(1,2,4:9,10,12,13,14,25,30)], 
						  zebrafish_rss$`24hpf`$X2[c(1:8,10:13,19,20,30)], 
						  zebrafish_rss$`24hpf`$X3[c(2:7)])

zebrafish_24hpf_cellType <- rep('RPC', length(zebrafish_24hpf_tf))


zebrafish_36hpf_tf <- c(zebrafish_rss$`36hpf`$X1[c(1:3,5,8,12,15,22)],
		  zebrafish_rss$`36hpf`$X2[-c(8,11,12,13,16,20:30)],#AC
		  zebrafish_rss$`36hpf`$X3[c(1:5,8:12,14,18)],#Cone
		  zebrafish_rss$`36hpf`$RGC[c(1:15,17,18,20,21,25,26)])
zebrafish_36hpf_cellType <- c(rep('RPC', 34), rep('RGC', 21))

zebrafish_48_hpf_tf <- c(zebrafish_rss$`48hpf`$AC[c(2:12)],
		  zebrafish_rss$`48hpf`$BC[c(1,4,6,7,11,14,22,27)],#AC
		  zebrafish_rss$`48hpf`$HC[c(2:6,9,16)],#Cone
		  zebrafish_rss$`48hpf`$MG[c(3,4,6,11,12,28)],
		  zebrafish_rss$`48hpf`$OFF.BC[c(3,7,11)],
		  zebrafish_rss$`48hpf`$ON.BC[c(2)],
		  zebrafish_rss$`48hpf`$PR.Precurs.[c(1:3,6,7,9,12)],
		  zebrafish_rss$`48hpf`$RGC[c(1:5,7,8,12,15:18,24)],
		  zebrafish_rss$`48hpf`$RPC[c(1:13,15)])

zebrafish_48hpf_cellType <- c(rep('AC', 11), rep('BC', 8),rep('HC', 7), rep('MG', 6),rep('OFF BC', 3), rep('ON BC', 1),
							  rep('PR Precurs.', 7), rep('RGC', 13),rep('RPC', 14))


zebrafish_72hpf_tf <- c(zebrafish_rss$`72hpf`$AC[c(1:5,7,8,10,11,15:19,22)],
		  zebrafish_rss$`72hpf`$Double.Cone[c(1,3:9)],#
		  zebrafish_rss$`72hpf`$Green.Cone[c(1,3,5:10)],#
		  zebrafish_rss$`72hpf`$Late.RPC[c(1,2,4:8)],
		  zebrafish_rss$`72hpf`$MG[c(1:5,7,8,12,14,15)],
		  zebrafish_rss$`72hpf`$OFF.BC[c(1:3,7,14)],
		  zebrafish_rss$`72hpf`$ON.BC[c(1:4,10,11,20)],
		  zebrafish_rss$`72hpf`$PR.Precur.[c(1,2,5:9,12)],
		  zebrafish_rss$`72hpf`$Red.Cone[c(2,4,6,7,8,9,11:13)],
		  zebrafish_rss$`72hpf`$RGC[c(1,4:8,10:16)],
		  zebrafish_rss$`72hpf`$Rod[c(2:4,9,10)],
		  zebrafish_rss$`72hpf`$RPC[c(1:5,7:10,12)],
		  zebrafish_rss$`72hpf`$T1.HC[c(1:12)],
		  zebrafish_rss$`72hpf`$T2.HC[c(1:12)])

zebrafish_72hpf_cellType <- c(rep('AC', 15), rep('Double Cone', 8),rep('Green Cone', 8), rep('Late RPC', 7),rep('MG', 10), rep('OFF BC', 5),
							  rep('ON BC', 7), rep('PR Precurs.', 8),rep('Red Cone', 9),rep('RGC', 13), rep('Rod', 5),rep('RPC', 10),
							  rep('T1 HC', 12), rep('T2 HC', 12))

zebrafish_14dpf_tf <- c(zebrafish_rss$`14dpf`$AC[c(1:8,10,13:16,18,19)],
		  zebrafish_rss$`14dpf`$BC[c(1:7,10,13,14,16)],
		  zebrafish_rss$`14dpf`$Cone[c(1:3,5:9,11,14:16,18)],#
		  zebrafish_rss$`14dpf`$HC[c(3,4,8,9,14,17)],
		  zebrafish_rss$`14dpf`$MG[c(1:3,5,11)],
		  zebrafish_rss$`14dpf`$PR.Precur.[c(1,2:9,11:14,16)],
		  zebrafish_rss$`14dpf`$RGC[c(1,3:12,14:17,20)],
		  zebrafish_rss$`14dpf`$Rod[c(1,5:8,13,14,17,19,20)],
		  zebrafish_rss$`14dpf`$RPC[c(1:9,11,14:16,19)])

zebrafish_14dpf_cellType <- c(rep('AC', 15), rep('BC', 11),rep('Cone', 13),rep('HC', 6), rep('MG', 5),
							  rep('PR Precurs.', 14), rep('RGC', 16),rep('Rod', 10),rep('RPC', 14))


zebrafish_stage_rep <- c(rep('24hpf',length(zebrafish_24hpf_tf)),
				 rep('36hpf',length(zebrafish_36hpf_tf)),
				 rep('48hpf',length(zebrafish_48_hpf_tf)),
				 rep('72hpf',length(zebrafish_72hpf_tf)),
				 rep('14dpf',length(zebrafish_14dpf_tf)))

zebrafish_all_cellType <- c(zebrafish_24hpf_cellType, zebrafish_36hpf_cellType, zebrafish_48hpf_cellType,
							zebrafish_72hpf_cellType, zebrafish_14dpf_cellType)
zebrafish_all_tf <- c(zebrafish_24hpf_tf, zebrafish_36hpf_tf,zebrafish_48_hpf_tf,zebrafish_72hpf_tf,zebrafish_14dpf_tf)

zebrafish_tf_df <- data.frame(stage = zebrafish_stage_rep,
						  regulon = zebrafish_all_tf,
						  tf_name = str_split(zebrafish_all_tf, pattern = ' ', simplify = TRUE)[,1],
						  cellType = zebrafish_all_cellType)




write.table(mouse_tf_df,'R_downstream/mouse/regulon_for_analysis.txt', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(chicken_tf_df,'R_downstream/chicken/regulon_for_analysis.txt', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(zebrafish_tf_df,'R_downstream/zebrafish/regulon_for_analysis.txt', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

mouse_tf_info <- table(str_split(mouse_all_tf, pattern = ' ', simplify = TRUE)[,1])
mouse_tf_df_unique <- data.frame(mouse_tf_info)
colnames(mouse_tf_df_unique) <- c('TF', 'Freq')
write.table(mouse_tf_df_unique,'R_downstream/mouse/tf_df.txt', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

chicken_tf_info <- table(str_split(chicken_all_tf, pattern = ' ', simplify = TRUE)[,1])
chicken_tf_df_unique <- data.frame(mouse_tf_info)
colnames(chicken_tf_df_unique) <- c('TF', 'Freq')
write.xlsx(chicken_tf_df_unique,file = 'R_downstream/chicken/tf_df.xlsx',sheetName = 'summary',row.names = FALSE)

zebrafish_tf_info <- table(str_split(zebrafish_all_tf, pattern = ' ', simplify = TRUE)[,1])
zebrafish_tf_df_unique <- data.frame(zebrafish_tf_info)
colnames(zebrafish_tf_df_unique) <- c('TF', 'Freq')
write.xlsx(zebrafish_tf_df_unique,file = 'R_downstream/zebrafish/tf_df.xlsx',sheetName = 'summary',row.names = FALSE)

#---------------------------------------regulon tf annotation------------------------------------------------

# read symbol ENSEMBL info
SymbolENSEMBL <- list()
SymbolENSEMBL$mouse <- parseGeneSymbolEnsemblDF('../../TF/zebrafishAndChicken/gtf/mouse/geneSymbolENSEMBL.txt', make_unique = TRUE)
SymbolENSEMBL$chicken <- read.table('../../TF/zebrafishAndChicken/gtf/chicken/geneSymbolENSEMBL.txt', sep = '\t')
SymbolENSEMBL$zebrafish <- read.table('../../TF/zebrafishAndChicken/gtf/zebrafish/geneSymbolENSEMBL.txt', sep = '\t')

# parse mouse regulon tf gene annotation
mouse_tf_df$ensembl <- SymbolENSEMBL$mouse[mouse_tf_df$tf_name,]$Ensembl

# parse chicken regulon tf gene annotation
colnames(SymbolENSEMBL$chicken) <- c('Ensembl', 'Gene_symbol')
SymbolENSEMBL$chicken[SymbolENSEMBL$chicken$Gene_symbol == '',]$Gene_symbol <- SymbolENSEMBL$chicken[SymbolENSEMBL$chicken$Gene_symbol == '',]$Ensembl
SymbolENSEMBL$chicken <- SymbolENSEMBL$chicken[!duplicated(SymbolENSEMBL$chicken[,2]),]
rownames(SymbolENSEMBL$chicken) <- SymbolENSEMBL$chicken$Gene_symbol
chicken_tf_df$ensembl <- SymbolENSEMBL$chicken[chicken_tf_df$tf_name,]$Ensembl

# parse zebrafish regulon tf gene annotation
colnames(SymbolENSEMBL$zebrafish) <- c('Ensembl', 'Gene_symbol')
SymbolENSEMBL$zebrafish[SymbolENSEMBL$zebrafish$Gene_symbol == '',]$Gene_symbol <- SymbolENSEMBL$zebrafish[SymbolENSEMBL$zebrafish$Gene_symbol == '',]$Ensembl
SymbolENSEMBL$zebrafish <- SymbolENSEMBL$zebrafish[!duplicated(SymbolENSEMBL$zebrafish[,2]),]
rownames(SymbolENSEMBL$zebrafish) <- SymbolENSEMBL$zebrafish$Gene_symbol

zebrafish_tf_df$ensembl <- SymbolENSEMBL$zebrafish[zebrafish_tf_df$tf_name,]$Ensembl

# update the gene symbol
chicken_tf_df[chicken_tf_df$tf_name == 'ENSGALG00000027907',]$tf_name <- 'NR2F1'

#---------------------------------------regulon tf union------------------------------------------------
mouse_ensembl <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', mirror = 'asia', host = 'https://dec2021.archive.ensembl.org')
zebrafish_ensembl <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl', mirror = 'asia',host = 'https://dec2021.archive.ensembl.org')
chicken_ensembl <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'ggallus_gene_ensembl', mirror = 'asia',host = 'https://dec2021.archive.ensembl.org')

#convert zebrafish and chicken tf to mouse symbol for comparison.
z2m_tf <- getLDS(attributes = 'ensembl_gene_id', 
				 filters = 'ensembl_gene_id', 
				 mart = zebrafish_ensembl,
				 values = unique(zebrafish_tf_df$ensembl), 
				 attributesL ='mgi_symbol', 
				 martL = mouse_ensembl, 
				 uniqueRows = TRUE)
z2m_tf_homo <- homologene(unique(zebrafish_tf_df$tf_name), inTax = 7955, outTax = 10090)

#add to tf_df
rownames(z2m_tf_homo) <- z2m_tf_homo$`7955`

z2m_tf_new <- as.data.frame(matrix(ncol = 2))
colnames(z2m_tf_new) <- c('ensembl', 'mouse_homo')
for (x in unique(z2m_tf$Gene.stable.ID)){

	temp <- paste(z2m_tf[z2m_tf$Gene.stable.ID == x,]$MGI.symbol, collapse =  ';')
    z2m_tf_new[nrow(z2m_tf_new)+1, ] <- c(x,temp)   
}
z2m_tf_new <- z2m_tf_new[-1,]
rownames(z2m_tf_new) <- z2m_tf_new$ensembl
zebrafish_tf_df$mouse_symbol <- z2m_tf_new[zebrafish_tf_df$ensembl,]$mouse_homo

diff_homo <-z2m_tf_homo[z2m_tf_homo$`10090` %in% setdiff(z2m_tf_homo$`10090`,z2m_tf$MGI.symbol),]
zebrafish_tf_df[!is.na(diff_homo[zebrafish_tf_df$tf_name,]$`10090`),]$mouse_symbol <- na.omit(diff_homo[zebrafish_tf_df$tf_name,]$`10090`)

z2m_tf_total <- unique(union(z2m_tf$MGI.symbol,z2m_tf_homo$`10090`))

#chicken
c2m_tf <- getLDS(attributes = 'ensembl_gene_id', 
				 filters = 'ensembl_gene_id', 
				 mart = chicken_ensembl,
				 values = unique(chicken_tf_df$ensembl), 
				 attributesL ='mgi_symbol', 
				 martL = mouse_ensembl, 
				 uniqueRows = TRUE)
c2m_tf_homo <- homologene(chicken_tf_df$tf_name, inTax = 9031, outTax = 10090)

rownames(c2m_tf_homo) <- c2m_tf_homo$`9031`
c2m_tf_new <- as.data.frame(matrix(ncol = 2))
colnames(c2m_tf_new) <- c('ensembl', 'mouse_homo')
for (x in unique(c2m_tf$Gene.stable.ID)){

	temp <- paste(c2m_tf[c2m_tf$Gene.stable.ID == x,]$MGI.symbol, collapse =  ';')
    c2m_tf_new[nrow(c2m_tf_new)+1, ] <- c(x,temp)   
}
c2m_tf_new <- c2m_tf_new[-1,]
rownames(c2m_tf_new) <- c2m_tf_new$ensembl
chicken_tf_df$mouse_symbol <- c2m_tf_new[chicken_tf_df$ensembl,]$mouse_homo

diff_homo <-c2m_tf_homo[c2m_tf_homo$`10090` %in% setdiff(c2m_tf_homo$`10090`,c2m_tf$MGI.symbol),]
chicken_tf_df[!is.na(diff_homo[chicken_tf_df$tf_name,]$`10090`),]$mouse_symbol <- na.omit(diff_homo[chicken_tf_df$tf_name,]$`10090`)

c2m_tf_total <- unique(union(c2m_tf$MGI.symbol,c2m_tf_homo$`10090`))


tf_info <- list(mouse=mouse_tf_df,chicken=chicken_tf_df, zebrafish=zebrafish_tf_df)
saveRDS(tf_info,'R_downstream/cross_species/tf_info.rds')


# Venn plot
orth_list <- list(Mouse = mouse_tf_df$tf_name,
                  Chicken = c2m_tf_total,
                  Zebrafish = z2m_tf_total)

orthogroup_overlap_plot <- ggvenn(orth_list, 
						fill_color =  c("#0073C2FF", "#EFC000FF", "#868686FF"), 
						fill_alpha = 0.5, 
						show_percentage = FALSE)




ggsave('R_downstream/cross_species/intersect_regulon_venn.pdf',orthogroup_overlap_plot,width = 15, height = 8, units = 'cm', dpi = 200)




#------------------------------------------Analysis of regulators-----------------------------------------------

# Analysis of regulators

# Analysis of the intersect part
three_intersect  <- Reduce(intersect,orth_list) 
# [1] "Nr2f1" "Jun"   "Maf"   "Pax6"  "Sox9"  "Fos"   "Tfdp2" "Thrb"  "Lhx1"  "Klf7"  "Gbx2"  "Etv5"  "Pbx1"  "Nr2f2" "Lhx9"  "Sox8" 

mouse_three <- mouse_tf_df[mouse_tf_df$tf_name %in% three_intersect,]
chicken_three <- chicken_tf_df[chicken_tf_df$mouse_symbol %in%three_intersect,]
zebrafish_three <- zebrafish_tf_df[zebrafish_tf_df$mouse_symbol %in%three_intersect,]

# mouse_chicken
mouse_zebrafish_intersect  <- setdiff(intersect(orth_list$Mouse,orth_list$Chicken),orth_list$Zebrafish)


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

tf_info$mouse$celltype_regulon <- paste0(tf_info$mouse$newCellType, '_', tf_info$mouse$tf_name)
tf_info$chicken$celltype_regulon <- paste0(tf_info$chicken$newCellType, '_', tf_info$chicken$mouse_symbol)
tf_info$zebrafish$celltype_regulon <- paste0(tf_info$zebrafish$newCellType, '_', tf_info$zebrafish$mouse_symbol)

#--- regulators that are active in all the three species

three_intersect  <- Reduce(intersect,orth_list) 


for (x in three_intersect){
	
	mouse_name <- tf_info$mouse[tf_info$mouse$celltype_regulon == x,]$tf_name
	mouse_name <- unique(mouse_name)
	chicken_name <- unique(chicken_name)
	chicken_name <- tf_info$chicken[tf_info$chicken$celltype_regulon == x,]$tf_name
	chicken_name <- unique(chicken_name)
	zebrafish_name <- tf_info$zebrafish[tf_info$zebrafish$celltype_regulon == x,]$tf_name
	zebrafish_name <- unique(zebrafish_name)
	plot_three_features(mouse_name,chicken_name,zebrafish_name, 'three_intersect', x)

}
plot_three_features(Klf7,ENSGALG00000008501,klf7b, 'three_intersect', ’RGC_Klf7‘)


#--- regulators that are only active in mouse and zebrafish

mouse_zebrafish <- sort(setdiff(intersect(tf_info$mouse$celltype_regulon,tf_info$zebrafish$celltype_regulon), three_intersect))

for (x in mouse_zebrafish[23:length(mouse_zebrafish)]){

	mouse_feature <- str_split(x, pattern = '_', simplify = TRUE)[,2]
	p1 <- FeaturePlot(sce_finalKeep, features = mouse_feature,reduction='umap',raster = TRUE, order = TRUE) + 
    	  theme_bw() +
    	  theme(panel.grid.major=element_blank(),
          		panel.grid.minor=element_blank(), 
          		axis.text.x=element_blank(),
          		axis.ticks.x=element_blank(),
          		axis.title.x=element_blank(),
          		axis.ticks.y=element_blank(),
          		axis.text.y=element_blank(),
          		axis.title.y=element_blank(),  
          		plot.title = element_text(hjust = 0.5,colour = "black")) +
          		NoLegend() + 
    ggtitle(mouse_feature)

    zebrafish_feature <- tf_info$zebrafish[tf_info$zebrafish$celltype_regulon == x,]$tf_name
    zebrafish_feature <- unique(zebrafish_feature)
    p2 <- FeaturePlot(zebrafish, features = zebrafish_feature,reduction='umap',raster = TRUE, order = TRUE) + 
    	  theme_bw() +
    	  theme(panel.grid.major=element_blank(),
          		panel.grid.minor=element_blank(), 
          		axis.text.x=element_blank(),
          		axis.ticks.x=element_blank(),
          		axis.title.x=element_blank(), 
          		axis.ticks.y=element_blank(),
          		axis.text.y=element_blank(), 
          		axis.title.y=element_blank(), 
          		plot.title = element_text(hjust = 0.5,colour = "black")) + 
    	  NoLegend() + 
    ggtitle(zebrafish_feature)
 	all_plot <- p1 | p2
    ggsave(paste0('R_downstream/cross_species/three_features/mouse_zebrafish/',x,'.pdf'), all_plot, dpi = 200, units = 'cm', width = 12, height = 6)


}


plot_three_features('Barhl2', 'ENSGALG00000040069', 'barhl2',folder = 'mouse_zebrafish', 'AC_Barhl2')
plot_three_features('Nr2f2', 'NR2F2', 'nr2f2',folder = 'mouse_zebrafish', 'AC_Nr2f2')
plot_three_features('Sox4', 'ENSGALG00000036497', 'sox4a',folder = 'mouse_zebrafish', 'AC_Sox4')
plot_three_features('Crx', 'OTX5', 'otx5',folder = 'mouse_zebrafish', 'Cone_Crx')
plot_three_features('Lhx4', 'ENSGALG00000036587', 'lhx4',folder = 'mouse_zebrafish', 'Cone_Lhx4')
plot_three_features('Nr3c1', 'NR3C1', 'nr3c1',folder = 'mouse_zebrafish', 'Cone_Nr3c1')
plot_three_features('Rxrg', 'RXRG', 'rxrgb',folder = 'mouse_zebrafish', 'Cone_Rxrg')
plot_three_features('Tcf7l1', 'ENSGALG00000032656', 'tcf7l1a',folder = 'mouse_zebrafish', 'MG_Tcf7l1')
plot_three_features('Ebf1', 'EBF1', 'ebf3a',folder = 'mouse_zebrafish', 'RGC_Ebf1')
plot_three_features('Isl2', 'ISL2', 'isl2b',folder = 'mouse_zebrafish', 'RGC_Isl2')


#--- regulators that are only active in mouse and chicken

mouse_chicken <- sort(setdiff(intersect(tf_info$mouse$celltype_regulon,tf_info$chicken$celltype_regulon), three_intersect))
plot_three_features('Tfap2c', 'ENSGALG00000007690', 'tfap2c',folder = 'mouse_chicken', 'AC_Barhl2')
plot_three_features('E2f1', 'E2F1', 'e2f1',folder = 'mouse_chicken', 'MG_E2f1')
plot_three_features('Nr2e1', 'NR2E1', 'nr2e1',folder = 'mouse_chicken', 'MG_Nr2e1')
plot_three_features('Eomes', 'EOMES', 'eomesa',folder = 'mouse_chicken', 'RGC_Eomes')

#--- regulators that are only active in zebrafish and chicken
chicken_zebrafish <- sort(setdiff(intersect(tf_info$zebrafish$celltype_regulon,tf_info$chicken$celltype_regulon), three_intersect))
plot_three_features('Pax6', 'PAX6', 'pax6a',folder = 'chicken_zebrafish', 'AC_Pax6')
plot_three_features('Dmbx1', 'DMBX1', 'dmbx1a',folder = 'chicken_zebrafish', 'BC_Dmbx1')
plot_three_features('Hmx3', 'HMX3', 'hmx3a',folder = 'chicken_zebrafish', 'HC_Hmx3')
plot_three_features('Onecut3', 'ENSGALG00000028081', 'onecut3b',folder = 'chicken_zebrafish', 'HC_Onecut3')
plot_three_features('Sox6', 'SOX6', 'sox6',folder = 'chicken_zebrafish', 'RGC_Sox6')




