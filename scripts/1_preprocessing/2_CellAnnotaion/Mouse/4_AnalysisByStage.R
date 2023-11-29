suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(clustree))
suppressMessages(library(readxl))
suppressMessages(library(do))
suppressMessages(library(igraph))
suppressMessages(library(SeuratDisk))
suppressMessages(library(Canek))
suppressMessages(library(homologene))

source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Clark2019neuron/')
set.seed(123)

options(future.globals.maxSize = 15000 * 1024^2)

sce_finalKeep <- readRDS('results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')

# We have data from 10 different stages,
#
#

#----------------------------------------P14------------------------------
mouse_P14 <- subset(sce_finalKeep, stage == 'P14')
DefaultAssay(mouse_P14) <- 'RNA'

mouse_P14 <- SCTransform(mouse_P14, method = "glmGamPoi")
mouse_P14 <- RunPCA(mouse_P14, verbose = FALSE)
mouse_P14 <- RunUMAP(mouse_P14, dims = 1:20)

mouse_P14 <- FindNeighbors(mouse_P14, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.2,0.4, 0.6, 0.8,1)
for (i in resolutions){
     mouse_P14 <- FindClusters(mouse_P14,resolution = i, verbose= FALSE)
}

# filter out rare cell types
mouse_P14$newCellType <- 'Others'
mouse_P14$newCellType[mouse_P14$SCT_snn_res.0.2 == "0"] <- "Rod_PR_1"
mouse_P14$newCellType[mouse_P14$SCT_snn_res.0.2 == "1"] <- "Muller_glia_1"
mouse_P14$newCellType[mouse_P14$SCT_snn_res.0.2 == "2"] <- "Biopolar_cell_1"
mouse_P14$newCellType[(mouse_P14$SCT_snn_res.0.2 == "3") & (mouse_P14$cellType == "Rods")] <- "Rod_PR_2"
mouse_P14$newCellType[(mouse_P14$SCT_snn_res.0.2 == "3") & (mouse_P14$cellType == "Cones")] <- "Cone_PR_1"
mouse_P14$newCellType[(mouse_P14$SCT_snn_res.0.2 == "4") & (mouse_P14$cellType == "Rods")] <- "Rod_PR_3"
mouse_P14$newCellType[(mouse_P14$SCT_snn_res.0.2 == "4") & (mouse_P14$cellType == "Cones")] <- "Cone_PR_2"
mouse_P14$newCellType[(mouse_P14$SCT_snn_res.0.2 == "4") & (mouse_P14$cellType == "Bipolar Cells")] <- "Biopolar_cell_2"
mouse_P14$newCellType[(mouse_P14$SCT_snn_res.0.2 == "5") | (mouse_P14$SCT_snn_res.0.2 == "6")] <- "Biopolar_cell_3"
mouse_P14$newCellType[mouse_P14$SCT_snn_res.0.2 == "8"] <- "Muller_glia_2"

mouse_P14 <- subset(mouse_P14, newCellType != 'Others')

saveRDS(mouse_P14,'results/rds/P14_filtered_lowcells.rds')

#----------------------------------------P8------------------------------
#sce_finalKeep <- readRDS('results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')
mouse_P8 <- subset(sce_finalKeep, stage == 'P8')
DefaultAssay(mouse_P8) <- 'RNA'

mouse_P8 <- SCTransform(mouse_P8, method = "glmGamPoi")
mouse_P8 <- RunPCA(mouse_P8, verbose = FALSE)

# run integration
mouse_P8 <- RunCanek(mouse_P8, "orig.ident")
mouse_P8<- ScaleData(mouse_P8)
mouse_P8 <- RunPCA(mouse_P8)

mouse_P8 <- RunUMAP(mouse_P8, dims = 1:20)

mouse_P8 <- FindNeighbors(mouse_P8, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.2,0.4, 0.6, 0.8,1)
for (i in resolutions){
     mouse_P8 <- FindClusters(mouse_P8,resolution = i, verbose= FALSE)
}

# filter out rare cell types
mouse_P8$newCellType <- 'Others'
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "0"] <- "Rod_PR_1"
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "1"] <- "Rod_PR_2"
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "2"] <- "Biopolar_cell_1"
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "3"] <- "Rod_PR_3"
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "4"] <- "Cone_PR"
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "5"] <- "Biopolar_cell_2"
mouse_P8$newCellType[mouse_P8$Canek_snn_res.0.2 == "6"] <- "Biopolar_cell_3"
mouse_P8$newCellType[(mouse_P8$Canek_snn_res.0.2 == "7") & (mouse_P8$cellType == "Muller Glia")] <- "Muller_glia"

mouse_P8 <- subset(mouse_P8, newCellType != 'Others')
saveRDS(mouse_P8,'results/rds/P8_filtered_lowcells.rds')

#----------------------------------------P5------------------------------
#sce_finalKeep <- readRDS('results/rds/Clark2019neuron_raw_integrated_finalKeep.rds')

mouse_P5 <- subset(sce_finalKeep, stage == 'P5')
DefaultAssay(mouse_P5) <- 'RNA'

mouse_P5 <- SCTransform(mouse_P5, method = "glmGamPoi")
mouse_P5 <- RunPCA(mouse_P5, verbose = FALSE)
mouse_P5 <- RunUMAP(mouse_P5, dims = 1:20)

mouse_P5 <- FindNeighbors(mouse_P5, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.2,0.4, 0.6, 0.8,1)
for (i in resolutions){
     mouse_P5 <- FindClusters(mouse_P5,resolution = i, verbose= FALSE)
}

mouse_P5$newCellType <- 'Others'
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "0") & (mouse_P5$cellType == "Muller Glia")] <- "Muller_glia"
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "1") & (mouse_P5$cellType == "Rods")] <- "Rod_PR_1"
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "2") & (mouse_P5$cellType == "Rods")] <- "Rod_PR_2"
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "3") & (mouse_P5$cellType == "Late RPCs")] <- "Late_RPCs"
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "4") & (mouse_P5$cellType == "Photoreceptor Precursors")] <- "PR_Precursors"
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "7") & (mouse_P5$cellType == "Amacrine Cells")] <- "Amacrine_cell"
mouse_P5$newCellType[(mouse_P5$SCT_snn_res.0.4 == "9") & (mouse_P5$cellType == "Cones")] <- "Cone_PR"
mouse_P5$newCellType[mouse_P5$cellType == "Bipolar Cells"] <- "Biopolar_cell"

mouse_P5 <- subset(mouse_P5, newCellType != 'Others')
saveRDS(mouse_P5,'results/rds/P5_filtered_lowcells.rds')

#----------------------------------------P2------------------------------
mouse_P2 <- subset(sce_finalKeep, stage == 'P2')
DefaultAssay(mouse_P2) <- 'RNA'

mouse_P2 <- SCTransform(mouse_P2, method = "glmGamPoi")
mouse_P2 <- RunPCA(mouse_P2, verbose = FALSE)
mouse_P2 <- RunUMAP(mouse_P2, dims = 1:20)

mouse_P2 <- FindNeighbors(mouse_P2, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.2,0.4, 0.6, 0.8,1)
for (i in resolutions){
     mouse_P2 <- FindClusters(mouse_P2,resolution = i, verbose= FALSE)
}

mouse_P2$newCellType <- 'Others'
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "0") & (mouse_P2$cellType == "Late RPCs")] <- "Late_RPCs_1"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "1") & (mouse_P2$cellType == "Late RPCs")] <- "Late_RPCs_2"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "2") & (mouse_P2$cellType == "Rods")] <- "Rod_PR"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "2") & (mouse_P2$cellType == "Photoreceptor Precursors")] <- "PR_Precursors_1"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "3") & (mouse_P2$cellType == "Neurogenic Cells")] <- "Neurogenic_cell_1"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "3") & (mouse_P2$cellType == "Photoreceptor Precursors")] <- "PR_Precursors_2"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "4") & (mouse_P2$cellType == "Amacrine Cells")] <- "Amacrine_cell"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "5") & (mouse_P2$cellType == "Cones")] <- "Cone_PR"
mouse_P2$newCellType[(mouse_P2$SCT_snn_res.0.2 == "7") & (mouse_P2$cellType == "Neurogenic Cells")] <- "Neurogenic_cell_2"

mouse_P2 <- subset(mouse_P2, newCellType != 'Others')

saveRDS(mouse_P2,'results/rds/P2_filtered_lowcells.rds')

#----------------------------------------P0------------------------------
mouse_P0 <- subset(sce_finalKeep, stage == 'P0')
DefaultAssay(mouse_P0) <- 'RNA'

mouse_P0 <- SCTransform(mouse_P0, method = "glmGamPoi")
mouse_P0 <- RunPCA(mouse_P0, verbose = FALSE)
mouse_P0 <- RunUMAP(mouse_P0, dims = 1:20)

mouse_P0 <- FindNeighbors(mouse_P0, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.1,0.2,0.4)
for (i in resolutions){
     mouse_P0 <- FindClusters(mouse_P0,resolution = i, verbose= FALSE)
}

mouse_P0$newCellType <- 'Others'
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "0") & (mouse_P0$cellType == "Late RPCs")] <- "Late_RPCs_1"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "1") & (mouse_P0$cellType == "Late RPCs")] <- "Late_RPCs_2"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "2") & (mouse_P0$cellType == "Neurogenic Cells")] <- "Neurogenic_cell_1"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "2") & (mouse_P0$cellType == "Photoreceptor Precursors")] <- "PR_Precursors_1"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "3") & (mouse_P0$cellType == "Amacrine Cells")] <- "Amacrine_cell"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "4") & (mouse_P0$cellType == "Cones")] <- "Cone_PR"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "4") & (mouse_P0$cellType == "Rods")] <- "Rod_PR"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "4") & (mouse_P0$cellType == "Photoreceptor Precursors")] <- "PR_Precursors_2"

mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "5") & (mouse_P0$cellType == "Cones")] <- "Cone_PR"
mouse_P0$newCellType[(mouse_P0$SCT_snn_res.0.1 == "7") & (mouse_P0$cellType == "Muller Glia")] <- "Muller_glia"

mouse_P0 <- subset(mouse_P0, newCellType != 'Others')

saveRDS(mouse_P0,'results/rds/P0_filtered_lowcells.rds')

#----------------------------------------E18------------------------------
mouse_E18 <- subset(sce_finalKeep, stage == 'E18')
DefaultAssay(mouse_E18) <- 'RNA'

mouse_E18 <- SCTransform(mouse_E18, method = "glmGamPoi")
mouse_E18 <- RunPCA(mouse_E18, verbose = FALSE)

mouse_E18 <- RunCanek(mouse_E18, "orig.ident")
mouse_E18<- ScaleData(mouse_E18)
mouse_E18 <- RunPCA(mouse_E18)

mouse_E18 <- RunUMAP(mouse_E18, dims = 1:20)

mouse_E18 <- FindNeighbors(mouse_E18, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     mouse_E18 <- FindClusters(mouse_E18,resolution = i, verbose= FALSE)
}

mouse_E18$newCellType <- 'Others'
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "0") & (mouse_E18$cellType == "Late RPCs")] <- "Late_RPCs_1"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "1") & (mouse_E18$cellType == "Late RPCs")] <- "Late_RPCs_2"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "2") & (mouse_E18$cellType == "Amacrine Cells")] <- "Amacrine_cell_1"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "3") & (mouse_E18$cellType == "Amacrine Cells")] <- "Amacrine_cell_2"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "4") & (mouse_E18$cellType == "Photoreceptor Precursors")] <- "PR_Precursors"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "4") & (mouse_E18$cellType == "Rods")] <- "Rod_PR"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "5") & (mouse_E18$cellType == "Cones")] <- "Cone_PR"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "6") & (mouse_E18$cellType == "Neurogenic Cells")] <- "Neurogenic_cell_1"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "7") & (mouse_E18$cellType == "Neurogenic Cells")] <- "Neurogenic_cell_2"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "8") & (mouse_E18$cellType == "Retinal Ganglion Cells")] <- "Retinal_ganglion"
mouse_E18$newCellType[(mouse_E18$Canek_snn_res.0.2 == "9") & (mouse_E18$cellType == "Amacrine Cells")] <- "Amacrine_cell_3"

mouse_E18 <- subset(mouse_E18, newCellType != 'Others')

saveRDS(mouse_E18,'results/rds/E18_filtered_lowcells.rds')

#----------------------------------------E16------------------------------
mouse_E16 <- subset(sce_finalKeep, stage == 'E16')
DefaultAssay(mouse_E16) <- 'RNA'

mouse_E16 <- SCTransform(mouse_E16, method = "glmGamPoi")
mouse_E16 <- RunPCA(mouse_E16, verbose = FALSE)
mouse_E16 <- RunUMAP(mouse_E16, dims = 1:20)
mouse_E16 <- FindNeighbors(mouse_E16, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     mouse_E16 <- FindClusters(mouse_E16,resolution = i, verbose= FALSE)
}

mouse_E16$newCellType <- 'Others'
mouse_E16$newCellType[(mouse_E16$cellType == "Amacrine Cells")] <- "Amacrine_cell"
mouse_E16$newCellType[(mouse_E16$cellType == "Late RPCs")] <- "Late_RPCs"
mouse_E16$newCellType[(mouse_E16$cellType == "Photoreceptor Precursors")] <- "PR_Precursors"
mouse_E16$newCellType[(mouse_E16$cellType == "Rods")] <- "Rod_PR"
mouse_E16$newCellType[(mouse_E16$cellType == "Cones")] <- "Cone_PR"
mouse_E16$newCellType[(mouse_E16$cellType == "Neurogenic Cells")] <- "Neurogenic_cell"
mouse_E16$newCellType[(mouse_E16$cellType == "Retinal Ganglion Cells")] <- "Retinal_ganglion"
mouse_E16$newCellType[(mouse_E16$cellType == "Early RPCs") & (mouse_E16$SCT_snn_res.0.2 == "0")] <- "Early_RPC_1"
mouse_E16$newCellType[(mouse_E16$cellType == "Early RPCs") & (mouse_E16$SCT_snn_res.0.2 == "3")] <- "Early_RPC_2"
mouse_E16$newCellType[(mouse_E16$cellType == "Early RPCs") & (mouse_E16$SCT_snn_res.0.2 == "5")] <- "Early_RPC_3"
mouse_E16$newCellType[(mouse_E16$cellType == "Horizontal Cells")] <- "Horizontal_cell"

mouse_E16 <- subset(mouse_E16, newCellType != 'Others')

saveRDS(mouse_E16,'results/rds/E16_filtered_lowcells.rds')

#----------------------------------------E14------------------------------
mouse_E14 <- subset(sce_finalKeep, stage == 'E14')
DefaultAssay(mouse_E14) <- 'RNA'

mouse_E14 <- SCTransform(mouse_E14, method = "glmGamPoi")
mouse_E14 <- RunPCA(mouse_E14, verbose = FALSE)
mouse_E14 <- RunCanek(mouse_E14, "orig.ident")
mouse_E14<- ScaleData(mouse_E14)
mouse_E14 <- RunPCA(mouse_E14)

mouse_E14 <- RunUMAP(mouse_E14, dims = 1:20)

mouse_E14 <- FindNeighbors(mouse_E14, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     mouse_E14 <- FindClusters(mouse_E14,resolution = i, verbose= FALSE)
}

mouse_E14$newCellType <- 'Others'
mouse_E14$newCellType[(mouse_E14$cellType == "Amacrine Cells")] <- "Amacrine_cell"
mouse_E14$newCellType[(mouse_E14$cellType == "Late RPCs")] <- "Late_RPCs"
mouse_E14$newCellType[(mouse_E14$cellType == "Photoreceptor Precursors")] <- "PR_Precursors"
mouse_E14$newCellType[(mouse_E14$cellType == "Rods")] <- "Rod_PR"
mouse_E14$newCellType[(mouse_E14$cellType == "Cones")] <- "Cone_PR"
mouse_E14$newCellType[(mouse_E14$cellType == "Neurogenic Cells")] <- "Neurogenic_cell"
mouse_E14$newCellType[(mouse_E14$cellType == "Retinal Ganglion Cells")] <- "Retinal_ganglion"
mouse_E14$newCellType[(mouse_E14$cellType == "Early RPCs") & (mouse_E14$Canek_snn_res.0.2 == "0")] <- "Early_RPC_1"
mouse_E14$newCellType[(mouse_E14$cellType == "Early RPCs") & (mouse_E14$Canek_snn_res.0.2 == "1")] <- "Early_RPC_2"
mouse_E14$newCellType[(mouse_E14$cellType == "Early RPCs") & (mouse_E14$Canek_snn_res.0.2 == "5")] <- "Early_RPC_3"
mouse_E14$newCellType[(mouse_E14$cellType == "Early RPCs") & (mouse_E14$Canek_snn_res.0.2 == "8")] <- "Early_RPC_4"
mouse_E14$newCellType[(mouse_E14$cellType == "Horizontal Cells")] <- "Horizontal_cell"

mouse_E14 <- subset(mouse_E14, newCellType != 'Others')

saveRDS(mouse_E14,'results/rds/E14_filtered_lowcells.rds')

#----------------------------------------E12------------------------------
mouse_E12 <- subset(sce_finalKeep, stage == 'E12')
DefaultAssay(mouse_E12) <- 'RNA'

mouse_E12 <- SCTransform(mouse_E12, method = "glmGamPoi")
mouse_E12 <- RunPCA(mouse_E12, verbose = FALSE)
mouse_E12 <- RunUMAP(mouse_E12, dims = 1:20)
mouse_E12 <- FindNeighbors(mouse_E12, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     mouse_E12 <- FindClusters(mouse_E12,resolution = i, verbose= FALSE)
}

mouse_E12$newCellType <- 'Others'
mouse_E12$newCellType[(mouse_E12$cellType == "Neurogenic Cells")] <- "Neurogenic_cell"
mouse_E12$newCellType[(mouse_E12$cellType == "Retinal Ganglion Cells")] <- "Retinal_ganglion"
mouse_E12$newCellType[(mouse_E12$cellType == "Early RPCs") & (mouse_E12$SCT_snn_res.0.05 == "0")] <- "Early_RPC_1"
mouse_E12$newCellType[(mouse_E12$cellType == "Early RPCs") & (mouse_E12$SCT_snn_res.0.05 == "1")] <- "Early_RPC_2"
mouse_E12$newCellType[(mouse_E12$cellType == "Early RPCs") & (mouse_E12$SCT_snn_res.0.05 == "2")] <- "Early_RPC_3"
mouse_E12$newCellType[(mouse_E12$cellType == "Early RPCs") & (mouse_E12$SCT_snn_res.0.05 == "3")] <- "Early_RPC_4"

mouse_E12 <- subset(mouse_E12, newCellType != 'Others')

saveRDS(mouse_E12,'results/rds/E12_filtered_lowcells.rds')

#----------------------------------------E11------------------------------
mouse_E11 <- subset(sce_finalKeep, stage == 'E11')
DefaultAssay(mouse_E11) <- 'RNA'

mouse_E11 <- SCTransform(mouse_E11, method = "glmGamPoi")
mouse_E11 <- RunPCA(mouse_E11, verbose = FALSE)
mouse_E11 <- RunUMAP(mouse_E11, dims = 1:20)

mouse_E11 <- FindNeighbors(mouse_E11, reduction = 'pca', dims = 1:20, verbose = FALSE)
resolutions <- c(0.05,0.1,0.2,0.4, 0.6)
for (i in resolutions){
     mouse_E11 <- FindClusters(mouse_E11,resolution = i, verbose= FALSE)
}

mouse_E11$newCellType <- 'Others'
mouse_E11$newCellType[(mouse_E11$SCT_snn_res.0.05 == "0")] <- "Early_RPC_1"
mouse_E11$newCellType[(mouse_E11$SCT_snn_res.0.05 == "1")] <- "Early_RPC_2"
mouse_E11$newCellType[(mouse_E11$SCT_snn_res.0.05 == "2")] <- "Early_RPC_3"
mouse_E11$newCellType[(mouse_E11$SCT_snn_res.0.05 == "3")] <- "Early_RPC_4"
mouse_E11$newCellType[(mouse_E11$SCT_snn_res.0.05 == "4")] <- "Early_RPC_5"
mouse_E11$newCellType[(mouse_E11$SCT_snn_res.0.05 == "5")] <- "Early_RPC_6"

mouse_E11 <- subset(mouse_E11, cellType=='Early RPCs')

saveRDS(mouse_E11,'results/rds/E11_filtered_lowcells.rds')