suppressMessages(library(Seurat))
suppressMessages(library(homologene))

set.seed(123)
options(future.globals.maxSize = 15000 * 1024^2)
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Clark2019neuron/')

#-----------------------------------------sample info-------------------------------------------#
sample_paths <- paste0(list.dirs('count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
sample_paths <- sample_paths[-11]
names(sample_paths) <- c('E11', 'E12', 'E14R1', 'E14R2', 'E16', 'E18R2', 'E18R3',
						 'P0', 'P14', 'P2R2', 'P5', 'P8R1', 'P8R2')

mcc.s <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090)$`10090`
mcc.g2m <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090)$`10090`


counts <- Read10X(data.dir = sample_paths, strip.suffix = FALSE)
sce <- CreateSeuratObject(counts,names.field = 1, min.cells = 3, min.features = 200, project = 'mouse retina')

sce[['percent.mt']] <- PercentageFeatureSet(sce, pattern = '^mt-')

sce <- subset(sce, subset = nFeature_RNA > 200 & percent.mt < 10)

sce$stage <- 'P8'
sce$stage[sce$orig.ident == "E11"] <- "E11"
sce$stage[sce$orig.ident == "E12"] <- "E12"
sce$stage[sce$orig.ident == "E14R1"] <- "E14"
sce$stage[sce$orig.ident == "E14R2"] <- "E14"
sce$stage[sce$orig.ident == "E16"] <- "E16"
sce$stage[sce$orig.ident == "E18R2"] <- "E18"
sce$stage[sce$orig.ident == "E18R3"] <- "E18"

sce$stage[sce$orig.ident == "P0"] <- "P0"
sce$stage[sce$orig.ident == "P14"] <- "P14"
sce$stage[sce$orig.ident == "P2R2"] <- "P2"
sce$stage[sce$orig.ident == "P5"] <- "P5"


sce.list <- SplitObject(sce, split.by = 'orig.ident')
sce.list$E11 <- subset(sce.list$E11, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
#sce.list$E12 <- subset(sce.list$E12, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
sce.list$E14R1 <- subset(sce.list$E14R1, subset = nCount_RNA < 10000 & nFeature_RNA < 4000)
sce.list$E14R2 <- subset(sce.list$E14R2, subset = nCount_RNA < 13000 & nFeature_RNA < 5000)
sce.list$E16 <- subset(sce.list$E16, subset = nCount_RNA < 10000 & nFeature_RNA < 4000)
sce.list$E18R2 <- subset(sce.list$E18R2, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
sce.list$E18R3 <- subset(sce.list$E18R3, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
sce.list$P0 <- subset(sce.list$P0, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
#sce.list$P14 <- subset(sce.list$P14, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
sce.list$P2R2 <- subset(sce.list$P2R2, subset = nCount_RNA < 15000 & nFeature_RNA < 5000)
sce.list$P5 <- subset(sce.list$P5, subset = nCount_RNA < 10000 & nFeature_RNA < 4000)
sce.list$P8R1 <- subset(sce.list$P8R1, subset = nCount_RNA < 10000 & nFeature_RNA < 4000)
sce.list$P8R2 <- subset(sce.list$P8R2, subset = nCount_RNA < 10000 & nFeature_RNA < 4000)

saveRDS('results/rds/Clark2019neuron_raw.rds')
