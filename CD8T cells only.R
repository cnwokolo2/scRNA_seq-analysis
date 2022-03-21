DefaultAssay(data.combined) <- "integrated"
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster2.rds")

levels(data.combined)
cd8t.cells <- subset(data.combined, idents = "CD8_T")
cd8t.cells$celltype <- "CD8_T"
plot1 <- FeatureScatter(cd8t.cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cd8t.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
head(table(Idents(cd8t.cells))) #604

cd8t.cells <- NormalizeData(cd8t.cells)
cd8t.cells <- FindVariableFeatures(cd8t.cells, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top30_cd8t.cells <- head(VariableFeatures(cd8t.cells), 30)
##write.csv(top30_cd8t.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/top30variantgenes.csv")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cd8t.cells)
plot2 <- LabelPoints(plot = plot1, points = top30_cd8t.cells, repel = TRUE, max.overlaps = 40)
plot1 + plot2

all.genes <- rownames(cd8t.cells)
cd8t.cells <- ScaleData(cd8t.cells, features = all.genes)
cd8t.cells <- RunPCA(cd8t.cells, features = VariableFeatures(cd8t.cells))
ElbowPlot(cd8t.cells, ndims = 30)
cd8t.cells <- RunUMAP(cd8t.cells, dims = 1:10)
cd8t.cells <- FindNeighbors(cd8t.cells, dims = 1:10)
cd8t.cells <- FindClusters(cd8t.cells, resolution = 0.7)

DimPlot(cd8t.cells, reduction = "umap", label = TRUE)
DimPlot(cd8t.cells, reduction = "umap", split.by = "KO", label = TRUE)
p1 <- DimPlot(cd8t.cells, reduction = "umap", group.by = "KO")
p2 <- DimPlot(cd8t.cells, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DefaultAssay(cd8t.cells) <- "RNA"
cd8t.markers <- FindAllMarkers(cd8t.cells, only.pos = TRUE, logfc.threshold = 0.8)
##write.csv(meta.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/Markers.csv")
saveRDS(cd8t.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/AllClusters1.rds")
readRDS(cd8t.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/AllClusters1.rds")

table(Idents(cd8t.cells))
#0   1   2   3   4   5 
#188 113 112  88  58  45 
prop.table(table(Idents(cd8t.cells)))
#  0            1          2          3          4          5 
#0.31125828 0.18708609 0.18543046 0.14569536 0.09602649 0.07450331

VlnPlot(cd8t.cells, features = c("SLC11A1", "PSAP", "NFKBID", "SOCS6"), ncol = 1, pt.size = 0.1)
library(metap)
plot_paper<-VlnPlot(cd8t.cells, features=c("SLC11A1", "PSAP", "NFKBID", "SOCS6", "GSN", "GRN", "CD3"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue", "red"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=3, legend = 'right')
plot_paper<-VlnPlot(cd8t.cells, features=c("CD8A", "CD8B","NKG7", "GZMB", "CD52"), pt.size = 0.1, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue", "red"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=2, legend = 'right')
DotPlot(cd8t.cells, features=c("CD8A", "CD8B", "LYZ","NKG7", "GZMB"), split.by="KO", cols = c("yellow","blue", "red"))
DotPlot(cd8t.cells, features=c("CD8A", "CD8B", "LYZ","NKG7", "GZMB"),  cols = c("yellow","blue", "red"))

library(dplyr)
DefaultAssay(cd8t.cells) <- "integrated"
cd8t.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top20
DoHeatmap(subset(cd8t.cells, downsample = 100), features = top20$gene)

cluster0.markers <- FindMarkers(cd8t.cells, ident.1 = 0, logfc.threshold = 0.8)
cluster1.markers <- FindMarkers(cd8t.cells, ident.1 = 1, logfc.threshold = 0.8)
cluster2.markers <- FindMarkers(cd8t.cells, ident.1 = 2, logfc.threshold = 0.8)
cluster3.markers <- FindMarkers(cd8t.cells, ident.1 = 3, logfc.threshold = 0.8)
cluster4.markers <- FindMarkers(cd8t.cells, ident.1 = 4, logfc.threshold = 0.8)
cluster5.markers <- FindMarkers(cd8t.cells, ident.1 = 5, logfc.threshold = 0.8)
#write.csv(cluster0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/cluster0.markers.csv")
#write.csv(cluster1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/cluster1.markers.csv")
#write.csv(cluster2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/cluster2.markers.csv")
#write.csv(cluster3.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/cluster3.markers.csv")
#write.csv(cluster4.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/cluster4.markers.csv")
#write.csv(cluster5.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/cluster5.markers.csv")

Zero.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 0, grouping.var = "KO", verbose = FALSE)
One.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 1, grouping.var = "KO", verbose = FALSE)
Two.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 2, grouping.var = "KO", verbose = FALSE)
Three.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 3, grouping.var = "KO", verbose = FALSE)
Four.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 4, grouping.var = "KO", verbose = FALSE)
Five.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 5, grouping.var = "KO", verbose = FALSE)
Six.markers <- FindConservedMarkers(cd8t.cells, ident.1 = 6, grouping.var = "KO", verbose = FALSE)
#write.csv(Zero.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C0Markers.csv")
#write.csv(One.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C1Markers.csv")
#write.csv(Two.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C2Markers.csv")
#write.csv(Three.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C3Markers.csv")
#write.csv(Four.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C4Markers.csv")
#write.csv(Five.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C5Markers.csv")
#write.csv(Six.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/csv/C6Markers.csv")

saveRDS(cd8t.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/AllClusters2.rds")
readRDS(cd8t.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/AllClusters2.rds")


#CLUSTER 3&5
cd8t_c3_c5 <- subset(cd8t.cells, idents = c(3,5))
DefaultAssay(cd8t_c3_c5) <- "RNA"

#DE genes expression
plot_paper<-VlnPlot(cd8t_c3_c5, features=c("CCR7", "LEF1", "TCF7", "SELL", "CD52", "CCL4", "CCL3", "CCL2", "CD14", "IL1B"),
                    pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("blue","yellow"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=5, legend = 'right')
plot_paper<-VlnPlot(cd8t_c3_c5, features=c("IL4R", "IL1B", "TNF", "CCR7", "AIF1", "CCL2", "SERPINF1","CCL4", "HLA-DRB1","JUN"),
                    pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("blue","yellow"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=5, legend = 'right')

library(dplyr)
Idents(cd8t_c3_c5) <- "KO"
cd8t_c3_c5.markers <- FindAllMarkers(cd8t_c3_c5, only.pos = TRUE, logfc.threshold = 0.5)
DefaultAssay(cd8t_c3_c5) <- "integrated"
cd8t_c3_c5.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top20
DoHeatmap(subset(cd8t_c3_c5, downsample = 50), features = top20$gene)

cd8t_c3.cells <- subset(cd8t.cells, idents = 3)
cd8t_c5.cells <- subset(cd8t.cells, idents = 5)
DefaultAssay(cd8t_c3.cells) <- "RNA"
DefaultAssay(cd8t_c5.cells) <- "RNA"
cd8t_c3.cells$celltype.KO <- paste(Idents(cd8t_c3.cells), cd8t_c3.cells$KO, sep = "_")
cd8t_c5.cells$celltype.KO <- paste(Idents(cd8t_c5.cells), cd8t_c5.cells$KO, sep = "_")
cd8t_c3.cells$celltype <- Idents(cd8t_c3.cells)
cd8t_c5.cells$celltype <- Idents(cd8t_c5.cells)
Idents(cd8t_c3.cells) <- "celltype.KO"
Idents(cd8t_c5.cells) <- "celltype.KO"

CD8T.c3.cd.markers <- FindMarkers(cd8t_c3.cells, ident.1 = "3_cd", verbose = FALSE)
#write.csv(CD8T.c3.cd.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/CSV/CD8T.c3.cd.markers.csv")
CD8T.c3.mts.markers <- FindMarkers(cd8t_c3.cells, ident.1 = "3_mts", verbose = FALSE)
#write.csv(CD8T.c3.mts.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/CSV/CD8T.c3.mts.markers.csv")
CD8T.c5.cd.markers <- FindMarkers(cd8t_c5.cells, ident.1 = "5_cd", verbose = FALSE)
#write.csv(CD8T.c5.cd.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/CSV/CD8T.c5.cd.markers.csv")
CD8T.c5.mts.markers <- FindMarkers(cd8t_c5.cells, ident.1 = "5_mts", verbose = FALSE)
#write.csv(CD8T.c5.mts.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/CSV/CD8T.c5.mts.markers.csv")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
Idents(cd8t.cells) <- "KO"
cd8t.cd.markers <- FindMarkers(cd8t.cells, ident.1 = "cd", verbose = FALSE)
#write.csv(cd8t.cd.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/CSV/cd8t.cd.markers.csv")
cd8t.mts.markers <- FindMarkers(cd8t.cells, ident.1 = "mts", verbose = FALSE)
#write.csv(cd8t.mts.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8T_only/CSV/cd8t.mts.markers.csv")


