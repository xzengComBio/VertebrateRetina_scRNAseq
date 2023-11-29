library(stringr)
library(ggplot2)


library(homologene)
library(biomaRt)

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Dr.eg.db)

setwd('~/Desktop/Research/ciona/singleCell/velocyto/')
source('~/Desktop/Research/ciona/singleCell/codes/util.r')
source('/Users/xzeng/Desktop/Research/ciona/SCENIC/code/downstream/utils.r')

ensembl_db <- list()
db_host <- 'https://dec2021.archive.ensembl.org'
ensembl_db$mouse <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host = db_host, mirror = 'asia')
ensembl_db$zebrafish <- useEnsembl('ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl', host = db_host, mirror = 'asia')

SymbolENSEMBL <- list()
SymbolENSEMBL$mouse <- parseGeneSymbolEnsemblDF('../../TF/zebrafishAndChicken/gtf/mouse/geneSymbolENSEMBL.txt', make_unique = TRUE)
SymbolENSEMBL$zebrafish <- read.table('../../TF/zebrafishAndChicken/gtf/zebrafish/geneSymbolENSEMBL.txt', sep = '\t')

colnames(SymbolENSEMBL$zebrafish) <- colnames(SymbolENSEMBL$mouse)


### Function
Others2Mouse <- function(gene_names,in_species){

    txid <- c(zebrafish=7955,chicken=9031)
    attr_list <- c(mouse = 'mgi_symbol', zebrafish = 'zfin_id_symbol', chicken = 'external_gene_name')
    homo_genes <- homologene(gene_names, inTax = txid[in_species], outTax = 10090)
    homo_genes <- homo_genes[,c(1,2)]
    colnames(homo_genes) <- c(in_species, 'mouse')
    orth_genes <- getLDS(attributes =attr_list[[in_species]], 
                 filters = attr_list[[in_species]], 
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



###
## Mouse
mouse_driver <- read.table('results/mouse/merge/lineage_drivers.txt', header = TRUE, sep = '\t')

mouse_top100_list <- list()
mouse_cellTypes <- unique(str_split(colnames(mouse_driver)[-1], pattern = '_', simplify = TRUE)[,1])

mouse_driver_celltype_genes <- c()
for (cellType in mouse_cellTypes){
    cellType_genes <- paste0(cellType, '_', mouse_driver[order(mouse_driver[,paste0(cellType,'_corr')], decreasing = TRUE),][1:100,]$Gene)
    mouse_driver_celltype_genes <- c(mouse_driver_celltype_genes,cellType_genes)
	mouse_top100_list[[cellType]] <- mouse_driver[order(mouse_driver[,paste0(cellType,'_corr')], decreasing = TRUE),][1:100,]$Gene
}


tf_info <- readRDS('../../SCENIC/results/R_downstream/cross_species/tf_info.rds')

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



mouse_driver_go_BP <- compareCluster(mouse_top100_list,
                             fun=enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             keyType = 'SYMBOL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)
write.table(as.data.frame(mouse_driver_go_BP), 'results/mouse/merge/GO/mouse_driver_go_BP.csv',sep = '\t', row.names = FALSE, quote = F)

p1 <- dotplot(mouse_driver_go_BP, showCategory=5) + 
		scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + 
		theme(axis.title.x = element_blank())

ggsave('results/mouse/merge/GO/mouse_driver_go_BP.pdf', p1, dpi = 200, units = 'cm', width = 25, height = 18)

# Zebrafish
zebrafish_driver <- read.table('results/zebrafish/merge/merge_lineage_drivers.txt', header = TRUE, sep = '\t')

zebrafish_top100_list <- list()
zebrafish_cellTypes <- unique(str_split(colnames(zebrafish_driver)[-1], pattern = '_', simplify = TRUE)[,1])

zebrafish_driver_celltype_genes <- c()
for (cellType in zebrafish_cellTypes){
    cellType_genes <- paste0(cellType, '_', zebrafish_driver[order(zebrafish_driver[,paste0(cellType,'_corr')], decreasing = TRUE),][1:100,]$Gene)
    zebrafish_driver_celltype_genes <- c(zebrafish_driver_celltype_genes,cellType_genes)
	zebrafish_top100_list[[cellType]] <- zebrafish_driver[order(zebrafish_driver[,paste0(cellType,'_corr')], decreasing = TRUE),][1:100,]$Gene

}

tf_info$zebrafish$celltype_regulon <- paste0(tf_info$zebrafish$newCellType, '_', tf_info$zebrafish$tf_name)



zebrafish_top100_list2m <- list()
for (x in names(zebrafish_top100_list)){
    print(x)
     zebrafish_top100_list2m[[x]] <- Others2Mouse(zebrafish_top100_list[[x]], 'zebrafish')
}

zebrafish_driver_go_BP <- compareCluster(zebrafish_top100_list,
                             fun=enrichGO,
                             OrgDb = 'org.Dr.eg.db',
                             keyType = 'SYMBOL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)
write.table(as.data.frame(zebrafish_driver_go_BP), 'results/zebrafish/merge/GO/zebrafish_driver_go_BP.csv',sep = '\t', row.names = FALSE, quote = F)

p1 <- dotplot(zebrafish_driver_go_BP, showCategory=5) + 
		scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + 
		theme(axis.title.x = element_blank())

ggsave('results/zebrafish/merge/GO/zebrafish_driver_go_BP.pdf', p1, dpi = 200, units = 'cm', width = 23, height = 18)


cellType <- intersect(names(zebrafish_top100_list2m), names(mouse_top100_list))

intersect_num <- c()
intersect_genes <- list()
for (x in cellType){
    print(x)
    print(length(intersect(zebrafish_top100_list2m[[x]][,2],mouse_top100_list[[x]])))
    intersect_genes[[x]] <- intersect(zebrafish_top100_list2m[[x]][,2],mouse_top100_list[[x]])
    intersect_num <- c(intersect_num, length(intersect(zebrafish_top100_list2m[[x]][,2],mouse_top100_list[[x]])))
}

intersect_num_df <- data.frame(cellType = cellType, intersect_num = intersect_num)

ggplot(intersect_num_df, aes(x=reorder(cellType, -intersect_num), y=intersect_num)) + 
    geom_bar(stat="identity") +
    xlab("Cell Type") +
    ylab("Shared Driver Genes in Mouse and Zebrafish") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())

intersect_genes_z <- list()
for (x in cellType){
    print(x)
    intersect_genes_z[[x]] <- zebrafish_top100_list2m[[x]][zebrafish_top100_list2m[[x]][,2] %in% intersect_genes[[x]],][,1]
}







