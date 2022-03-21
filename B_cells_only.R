DefaultAssay(data.combined) <- "integrated"
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster6.rds")

levels(data.combined)
b.cells <- subset(data.combined, idents = "B cells")
b.cells$celltype <- "b cells"
plot1 <- FeatureScatter(b.cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(b.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
head(table(Idents(b.cells))) #1128

b.cells <- NormalizeData(b.cells)
b.cells <- FindVariableFeatures(b.cells, selection.method = "vst", nfeatures = 4000)

# Identify the 10 most highly variable genes
top30_b.cells <- head(VariableFeatures(b.cells), 30)
#write.csv(top30_b.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/top30variantgenes.csv")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(b.cells)
plot2 <- LabelPoints(plot = plot1, points = top30_b.cells, repel = TRUE, max.overlaps = 40)
plot1 + plot2

all.genes <- rownames(b.cells)
b.cells <- ScaleData(b.cells, features = all.genes)
b.cells <- RunPCA(b.cells, features = VariableFeatures(b.cells))
ElbowPlot(b.cells, ndims = 50)
b.cells <- RunUMAP(b.cells, dims = 1:20)
b.cells <- FindNeighbors(b.cells, dims = 1:20)
b.cells <- FindClusters(b.cells, resolution = 0.7)

DimPlot(b.cells, reduction = "umap", label = TRUE)
DimPlot(b.cells, reduction = "umap", split.by = "KO", label = TRUE)

#Idents(b.cells) <- 'KO'
#head(table(Idents(b.cells)))
#mts   cd pbmc 
#692  111  325 
DefaultAssay(b.cells) <- "RNA"
b.markers <- FindAllMarkers(b.cells, logfc.threshold = 0.8, only.pos = TRUE )
#write.csv(b.markers, file = '/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/all_b_cell_markers.csv')
unique_c3.marker <- FindMarkers(b.cells, ident.1 = 3)
#write.csv(unique_c3.marker, file = '/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/unique_c3.marker.csv')

#visualization
library(dplyr)
DefaultAssay(b.cells) <- "integrated"
b.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(subset(b.cells, downsample = 50), features = top20$gene)

VlnPlot(b.cells, features = c("TNF", "CXCL8", "PER3", "CCL2"), ncol = 2, pt.size = 0.1)
cluster3 <- subset(b.cells, idents = 3)
plot_paper<-VlnPlot(cluster3, features=c("TNF", "NFKBIA", "FOS", "CCL2", "MT-CO1", "MT-CO2"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("blue","yellow"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=3, legend = 'right')

saveRDS(b.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/bcells_only.rds")
readRDS(b.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/bcells_only.rds")


