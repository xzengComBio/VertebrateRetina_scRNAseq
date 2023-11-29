suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(homologene))

# set the working directory
source('/Users/xzeng/Desktop/Research/ciona/singleCell/codes/util.r')
setwd('/Users/xzeng/Desktop/Research/ciona/singleCell/other_species/Xu2020Development/')
set.seed(123)

# read the samples
sample_paths <- paste0(list.dirs('count_matrix', recursive = FALSE), '/filtered_feature_bc_matrix/')
names(sample_paths) <- c('14dpf1', '14dpf2', '24hpf', '36hpf', '48hpf1',
						 '48hpf2', '72hpf1', '72hpf2')
counts <- Read10X(data.dir = sample_paths, strip.suffix = TRUE)

# create seurat object and add stage data
sce <- CreateSeuratObject(counts,names.field = 1, min.cells = 3, min.features = 200, project = 'zebrafish retina')

# add metadata stage
sce$stage <- '14dpf'
sce$stage[sce$orig.ident == "24hpf"] <- "24hpf"
sce$stage[sce$orig.ident == "36hpf"] <- "36hpf"
sce$stage[sce$orig.ident == "48hpf1"] <- "48hpf"
sce$stage[sce$orig.ident == "48hpf2"] <- "48hpf"
sce$stage[sce$orig.ident == "72hpf1"] <- "72hpf"
sce$stage[sce$orig.ident == "72hpf2"] <- "72hpf"

#calcuate the mitochondria rate as well as cell cycle scores
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")

zfcc.g2m <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 7955)
zfcc.g2m  <- zfcc.g2m$`7955`
zfcc.s <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 7955)
zfcc.s  <- zfcc.s$`7955`

sce <- CellCycleScoring(sce, s.features = zfcc.s, g2m.features = zfcc.g2m)


#---------------------------------------------------------------


# filter out cells
sce <- subset(sce, subset = percent.mt < 10)

sce.list <- SplitObject(sce, split.by = 'orig.ident')
sce.list$`14dpf1` <- subset(sce.list$`14dpf1`, subset = nCount_RNA < 19713 & nFeature_RNA < 4000)
sce.list$`14dpf2` <- subset(sce.list$`14dpf2`, subset = nCount_RNA < 27021 & nFeature_RNA < 4000)
sce.list$`24hpf` <- subset(sce.list$`24hpf`, subset = nCount_RNA < 44089 & nFeature_RNA < 4500)
sce.list$`36hpf` <- subset(sce.list$`36hpf`, subset = nCount_RNA < 23070 & nFeature_RNA < 3900)
sce.list$`48hpf1` <- subset(sce.list$`48hpf1`, subset = nCount_RNA < 27296 & nFeature_RNA < 4200)
sce.list$`48hpf2` <- subset(sce.list$`48hpf2`, subset = nCount_RNA < 20114 & nFeature_RNA < 4200)
sce.list$`72hpf1` <- subset(sce.list$`72hpf1`, subset = nCount_RNA < 18039 & nFeature_RNA < 4000)
sce.list$`72hpf2` <- subset(sce.list$`72hpf2`, subset = nCount_RNA < 21424 & nFeature_RNA < 3500)

saveRDS(sce.list,'results/rds/Xu2020Development_raw.rds')


