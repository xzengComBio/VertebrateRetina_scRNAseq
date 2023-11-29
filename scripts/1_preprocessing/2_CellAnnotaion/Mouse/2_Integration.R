suppressMessages(library(Seurat))
suppressMessages(library(Canek))
suppressMessages(library(ggplot2))
set.seed(123)

options(future.globals.maxSize = 1000000 * 1024^2)

args <- commandArgs()
sce_merge <- readRDS(args[6])

plotColor <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
		"#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
		"#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
		"#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


print("Runing SCTransform for samples....")
sce_merge <- SCTransform(sce_merge, method = "glmGamPoi")
sce_merge <- RunPCA(sce_merge)

p1 <- DimPlot(sce_merge, reduction = 'pca', group.by = 'orig.ident', cols = plotColor) +
			theme_bw() +
			theme(panel.grid.major=element_blank(),
				  panel.grid.minor=element_blank(),
				  axis.title = element_text(face = "bold",size = rel(1))) + 
			ggtitle('Before integration')

print("Runing Canek for samples....")
sce_merge <- RunCanek(sce_merge, "orig.ident")
sce_merge<- ScaleData(sce_merge)
sce_merge <- RunPCA(sce_merge)

p2<-  DimPlot(sce_merge, reduction = 'pca', group.by = 'orig.ident', cols = plotColor) +
			theme_bw() +
			theme(panel.grid.major=element_blank(),
				  panel.grid.minor=element_blank(),
				  axis.title = element_text(face = "bold",size = rel(1))) + 
			ggtitle('After integration')


Integration_plot <- p1 | p2

ggsave('/home/xzeng/project/Ciona/output/integration/integ_compare_pca.png',
		Integration_plot,
		dpi = 200,
		units = 'cm',
		width = 45,
		height = 15)

print("Saveing integrated rds")
saveRDS(sce_merge,args[7])
