library(SCopeLoomR)
library(SCENIC)
library(Seurat)

setwd('~/Desktop/Research/ciona/SCENIC/results')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')
set.seed(123)


# mouse regulons
stage_names <- c('E11', 'E12', 'E14', 'E16', 'E18', 'P0', 'P2', 'P5', 'P8', 'P14')

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

mouse_all_tf_info <- list(E11 = E11_regulons_info$tf_info,
					E12 = E12_regulons_info$tf_info,
					E14 = E14_regulons_info$tf_info,
					E16 = E16_regulons_info$tf_info,
					E18 = E18_regulons_info$tf_info,
					P0 = P0_regulons_info$tf_info,
					P2 = P2_regulons_info$tf_info,
					P5 = P5_regulons_info$tf_info,
					P8 = P8_regulons_info$tf_info,
					P14 = P14_regulons_info$tf_info)


stages <- c('E11','E12','E14','E16','E18','P0','P2','P5','P8','P14')

for (x in 1:10){
	mouse_all_tf_info[[x]]$stage_regulon <- paste(stages[x],mouse_all_tf_info[[x]]$tf_name, sep = '_')
}


mouse_all_tf_df <- Reduce(rbind,mouse_all_tf_info)

saveRDS(mouse_all_tf_df, 'R_downstream/mouse/all_tf_info.rds')

# chicken regulons
stage_names <- c('E12', 'E16', 'E18')

E12_regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[1],'/regulon_info.rds'))
E16_regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[2],'/regulon_info.rds'))
E18_regulon_info <- readRDS(paste0('R_downstream/chicken/',stage_names[3],'/regulon_info.rds'))



chicken_all_tf_info <- list(E12 = E12_regulon_info$tf_info,
							E16 = E16_regulon_info$tf_info,
							E18 = E18_regulon_info$tf_info)

for (x in 1:3){
	chicken_all_tf_info[[x]]$stage_regulon <- paste(stage_names[x],chicken_all_tf_info[[x]]$gene_symbol, sep = '_')
}

chicken_all_tf_df <- Reduce(rbind,chicken_all_tf_info)

saveRDS(chicken_all_tf_df, 'R_downstream/chicken/all_tf_info.rds')


# zebrafish regulons

stage_names <- c('24hpf', '36hpf', '48hpf', '72hpf', '14dpf')
hpf24_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[1],'/regulon_info.rds'))
hpf36_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[2],'/regulon_info.rds'))
hpf48_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[3],'/regulon_info.rds'))
hpf72_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[4],'/regulon_info.rds'))
dpf14_regulon_info <- readRDS(paste0('R_downstream/zebrafish/',stage_names[5],'/regulon_info.rds'))

zebrafish_all_tf_info <- list(hpf24 = hpf24_regulon_info$tf_info,
							  hpf36 = hpf36_regulon_info$tf_info,
							  hpf48 = hpf48_regulon_info$tf_info,
							  hpf72 = hpf72_regulon_info$tf_info,
							  dpf14 = dpf14_regulon_info$tf_info)

for (x in 1:5){
	zebrafish_all_tf_info[[x]]$stage_regulon <- paste(stage_names[x],zebrafish_all_tf_info[[x]]$gene_symbol, sep = '_')
}

zebrafish_all_tf_df <- Reduce(rbind,zebrafish_all_tf_info)


saveRDS(zebrafish_all_tf_df, 'R_downstream/zebrafish/all_tf_info.rds')




