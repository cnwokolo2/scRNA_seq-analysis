DefaultAssay(data.combined) <- "integrated"
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster2.rds")

levels(data.combined)
nk.cells <- subset(data.combined, idents = "NK cells")
nk.cells$celltype <- "NK cells"
plot1 <- FeatureScatter(nk.cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(nk.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
head(table(Idents(nk.cells))) #159

nk.cells <- NormalizeData(nk.cells)
nk.cells <- FindVariableFeatures(nk.cells, selection.method = "vst", nfeatures = 4000)

# Identify the 10 most highly variable genes
top30_nk.cells <- head(VariableFeatures(nk.cells), 30)
#write.csv(top30_nk.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/top30variantgenes.csv")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nk.cells)
plot2 <- LabelPoints(plot = plot1, points = top30_nk.cells, repel = TRUE, max.overlaps = 40)
plot1 + plot2

all.genes <- rownames(nk.cells)
nk.cells <- ScaleData(nk.cells, features = all.genes)
nk.cells <- RunPCA(nk.cells, features = VariableFeatures(nk.cells))
ElbowPlot(nk.cells, ndims = 50)
nk.cells <- RunUMAP(nk.cells, dims = 1:20)
nk.cells <- FindNeighbors(nk.cells, dims = 1:20)
nk.cells <- FindClusters(nk.cells, resolution = 0.7)

DimPlot(nk.cells, reduction = "umap", label = TRUE)
DimPlot(nk.cells, reduction = "umap", split.by = "KO", label = TRUE)
p1 <- DimPlot(nk.cells, reduction = "umap", group.by = "KO")
p2 <- DimPlot(nk.cells, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DefaultAssay(nk.cells) <- "RNA"
nk.markers <- FindAllMarkers(nk.cells, only.pos = TRUE, logfc.threshold = 0.8)
write.csv(meta.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/Markers.csv")
saveRDS(nk.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/AllClusters1.rds")
readRDS(nk.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/AllClusters1.rds")

table(Idents(nk.cells))
#0   1   2 
#76 66 17 
#mts   cd pbmc 
#17    1  141 
prop.table(table(Idents(nk.cells)))

VlnPlot(nk.cells, features = c("SLC11A1", "PSAP", "NFKBID", "SOCS6"), ncol = 1, pt.size = 0.1)
library(metap)
plot_paper<-VlnPlot(nk.cells, features=c("SLC11A1", "PSAP", "NFKBID", "SOCS6", "GSN", "GRN", "CD3"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("blue","yellow", "red"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')
plot_paper<-VlnPlot(nk.cells, features=c("CD8A", "CD8B","NKG7", "GZMB"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue", "red"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=2, legend = 'right')
DotPlot(nk.cells, features=c("CD8A", "CD8B", "LYZ","NKG7", "GZMB"), split.by="KO", cols = c("yellow","blue", "red"))
DotPlot(nk.cells, features=c("CD8A", "CD8B", "LYZ","NKG7", "GZMB"),  cols = c("yellow","blue", "red"))

library(dplyr)
DefaultAssay(nk.cells) <- "integrated"
nk.markers %>%
  group_by(cluster) %>%
  top_n(n = 40, wt = avg_log2FC) -> top20
DoHeatmap(subset(nk.cells, downsample = 50), features = top20$gene)

cluster0.markers <- FindMarkers(nk.cells, ident.1 = 0, logfc.threshold = 0.8)
cluster1.markers <- FindMarkers(nk.cells, ident.1 = 1, logfc.threshold = 0.8)
cluster2.markers <- FindMarkers(nk.cells, ident.1 = 2, logfc.threshold = 0.8)
write.csv(cluster0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/cluster0.markers.csv")
write.csv(cluster1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/cluster1.markers.csv")
write.csv(cluster2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/cluster2.markers.csv")

Zero.markers <- FindConservedMarkers(nk.cells, ident.1 = 0, grouping.var = "KO", verbose = FALSE)
One.markers <- FindConservedMarkers(nk.cells, ident.1 = 1, grouping.var = "KO", verbose = FALSE)
Two.markers <- FindConservedMarkers(nk.cells, ident.1 = 2, grouping.var = "KO", verbose = FALSE)
write.csv(Zero.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/C0Markers.csv")
write.csv(One.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/C1Markers.csv")
write.csv(Two.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/csv/C2Markers.csv")

saveRDS(nk.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/AllClusters2.rds")
readRDS(nk.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/AllClusters2.rds")


#CLUSTER 4 & 6
nk_c2.cells <- subset(nk.cells, idents = 2)
DefaultAssay(nk_c2.cells) <- "RNA"
VlnPlot(nk_c2.cells, features = c("GZMB","GNLY", "GZMM", "PRF1", "NKG7","CCL4",  "IFITM2", "JAK1"), ncol = 4, pt.size = 0.1, cols = "sky blue")
VlnPlot(nk.cells, features = c("GZMB","GNLY", "GZMM", "PRF1", "NKG7","CCL4",  "IFITM2", "JAK1"), ncol = 4, pt.size = 0.1, cols = c("yellow", "blue", "magenta"))
plot_paper<-VlnPlot(nk.cells, features= c("GZMB","GNLY", "GZMM", "PRF1", "NKG7","CCL4",  "IFITM2", "JAK1"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue", "red"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')

nk.c2.mts.markers <- FindMarkers(nk_c2.cells, ident.1 = "2_mts")
write.csv(nk.c2.mts.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/CSV/nk.c2.mts.markers.csv")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
Idents(nk.cells) <- "KO"
nk.cd.markers <- FindMarkers(nk.cells, ident.1 = "cd", verbose = FALSE)
write.csv(nk.cd.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/CSV/nk.cd.markers.csv")
nk.mts.markers <- FindMarkers(nk.cells, ident.1 = "mts", verbose = FALSE)
write.csv(nk.mts.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/nk_only/CSV/nk.mts.markers.csv")


