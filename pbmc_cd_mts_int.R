#Plan: To integrated PBMC, CD, MTS dataset and look for transcriptomic differences.

library(Seurat)
library(patchwork)
library(ggplot2)
library(metap)
library(dplyr)
library(Matrix)

mts <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/Temporal_Sclerosis_Biopsy/")
cd <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/Cortical_Dysplasia_Biopsy/")
pbmc <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21_22/PBMCs/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc, project = "pbmc", min.cells = 4, min.features = 300)
pbmc$KO <- "pbmc"
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(table(Idents(pbmc))) #2685

mts <- CreateSeuratObject(counts = mts, project = "mts", min.cells = 4, min.features = 300)
mts$KO <- "mts"
mts[["percent.mt"]] <- PercentageFeatureSet(mts, pattern = "^MT-")
head(table(Idents(mts))) #8412

cd <- CreateSeuratObject(counts = cd, project = "cd", min.cells = 4, min.features = 300)
cd$KO <- "cd"
cd[["percent.mt"]] <- PercentageFeatureSet(cd, pattern = "^MT-")
head(table(Idents(cd))) #1612

VlnPlot(mts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(cd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(mts, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot4 <- FeatureScatter(cd, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(cd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 + plot5

plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 + plot5

mts <- subset(mts, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 15000 & percent.mt <=15)
plot7 <- FeatureScatter(mts, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot8 <- FeatureScatter(mts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot7 + plot8
head(table(Idents(mts))) #8049

cd <- subset(cd, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 50000 & percent.mt <=10)
plot11 <- FeatureScatter(cd, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot12 <- FeatureScatter(cd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot11 + plot12
head(table(Idents(cd))) #1265

pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 10000 & percent.mt <=5)
plot11 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot12 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot11 + plot12
head(table(Idents(pbmc))) #2450

VlnPlot(mts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(cd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data.list <- list(mts, cd, pbmc)
data.list <- lapply(X= data.list, FUN = function(x) { 
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)

diseased.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = diseased.anchors)
DefaultAssay(data.combined) <- "integrated"

all.genes <- rownames(data.combined)
# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE, features = all.genes)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:10)
ElbowPlot(data.combined, ndims = 30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:10)
data.combined <- FindClusters(data.combined, resolution = 0.5)
DimPlot(data.combined, reduction = "umap", label = TRUE)

table(Idents(data.combined))
DefaultAssay(data.combined) <- "RNA"
data.markers <- FindAllMarkers(data.combined, logfc.threshold = 0.8, only.pos = TRUE )
#write.csv(data.markers, file = '/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/data_markers.csv')

p1 <- DimPlot(data.combined, reduction = "umap", group.by = "KO")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(data.combined, reduction = "umap", split.by = "KO", label = TRUE)
DimPlot(data.combined, reduction = "umap", label = T)

VlnPlot(data.combined, features = c("CCR7", "IL7R", "S100A4", "MS4A7", "FCGR3A"), ncol = 1, pt.size = 0.1)
VlnPlot(data.combined, features = c("CD14", "LYZ", "CD68", "FCER1G"), ncol = 1, pt.size = 0.1)    
VlnPlot(data.combined, features = c("CD79B", "CD3G", "CD3E", "CD8A", "CD8B", "NKG7"), ncol = 3, pt.size = 0.1)
VlnPlot(data.combined, features = c("CD4", "NCR1", "TBX21", "NKG7", "GZMB"), ncol = 1, pt.size = 0.1)
VlnPlot(data.combined, features = c("CCR7","IL7R", "CST3","FCER1G"), ncol = 2, pt.size = 0.1)
VlnPlot(data.combined, features = c("PTPRC","TREM2", "TMEM119","P2RY12", "GPR34", "CX3CR1", "SPI1", "CEBPA", "IRF9", "MEF2A", "IRF8", "EGR2", "CD3E"), ncol = 6, pt.size = 0) #Microglia genes
VlnPlot(data.combined, features = c("CD3E", "CD3G", "CD8A", "CD8B", "NKG7", "CD4", "IL7R", "SIGLEC1",  "CST3", "CD79B"), ncol = 5, pt.size = 0)
#no hdc, ppbp, SIGLECH, CD141

saveRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/Trial1Cluster.rds")
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/Trial1Cluster.rds")

#Visualization
library(metap)
plot_paper<-VlnPlot(data.combined, features=c("CD3G", "CD3E", "CCR7", "CD8A", "CD8B", "NKG7", "GZMB", "CD79B", "FCER1G","CST3"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("blue","yellow", "red"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=5, legend = 'right')
DotPlot(data.combined, features=c("CCR7", "CD3G", "CD8A", "CD8B", "NKG7", "GZMB", "CD79B", "FCER1G","CST3"), split.by="KO", cols = c("yellow","blue", "red"))
DotPlot(data.combined, features=c("CCR7", "CD3G", "CD8A", "CD8B", "NKG7", "GZMB", "CD79B", "FCER1G","CST3"))


cluster0.markers <- FindMarkers(data.combined, ident.1 = 0, logfc.threshold = 0.8)
cluster1.markers <- FindMarkers(data.combined, ident.1 = 1, logfc.threshold = 0.8)
cluster2.markers <- FindMarkers(data.combined, ident.1 = 2, logfc.threshold = 0.8)
cluster3.markers <- FindMarkers(data.combined, ident.1 = 3, logfc.threshold = 0.8)
cluster4.markers <- FindMarkers(data.combined, ident.1 = 4, logfc.threshold = 0.8)
cluster5.markers <- FindMarkers(data.combined, ident.1 = "UNK", logfc.threshold = 0.5)
cluster6.markers <- FindMarkers(data.combined, ident.1 = 6, logfc.threshold = 0.8)
cluster7.markers <- FindMarkers(data.combined, ident.1 = 7, logfc.threshold = 0.8)
cluster8.markers <- FindMarkers(data.combined, ident.1 = 8, logfc.threshold = 0.8)
cluster9.markers <- FindMarkers(data.combined, ident.1 = 9, logfc.threshold = 0.8)
cluster10.markers <- FindMarkers(data.combined, ident.1 = 10, logfc.threshold = 0.8)
cluster11.markers <- FindMarkers(data.combined, ident.1 = 11, logfc.threshold = 0.8)
cluster12.markers <- FindMarkers(data.combined, ident.1 = 12, logfc.threshold = 0.8)
cluster13.markers <- FindMarkers(data.combined, ident.1 = 13, logfc.threshold = 0.8)

write.csv(cluster0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster0.markers.csv")
write.csv(cluster1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster1.markers.csv")
write.csv(cluster2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster2.markers.csv")
write.csv(cluster3.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster3.markers.csv")
write.csv(cluster4.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/PBMC_CD_MTS_Int/cluster4.markers.csv")
write.csv(cluster5.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/PBMC_CD_MTS_Int/cluster5.markers.csv")
write.csv(cluster6.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster6.markers.csv")
write.csv(cluster7.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster7.markers.csv")
write.csv(cluster8.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster8.markers.csv")
write.csv(cluster9.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster9.markers.csv")
write.csv(cluster10.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster10.markers.csv")
write.csv(cluster11.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster11.markers.csv")
write.csv(cluster12.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/cluster12.markers.csv")
write.csv(cluster13.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/PBMC_CD_MTS_Int/CSV/cluster13.markers.csv")

Zero.markers <- FindConservedMarkers(data.combined, ident.1 = 0, grouping.var = "KO", verbose = FALSE)
One.markers <- FindConservedMarkers(data.combined, ident.1 = 1, grouping.var = "KO", verbose = FALSE)
Two.markers <- FindConservedMarkers(data.combined, ident.1 = 2, grouping.var = "KO", verbose = FALSE)
Three.markers <- FindConservedMarkers(data.combined, ident.1 = 3, grouping.var = "KO", verbose = FALSE)
Four.markers <- FindConservedMarkers(data.combined, ident.1 = 4, grouping.var = "KO", verbose = FALSE)
Five.markers <- FindConservedMarkers(data.combined, ident.1 = 5, grouping.var = "KO", verbose = FALSE)
Six.markers <- FindConservedMarkers(data.combined, ident.1 = 6, grouping.var = "KO", verbose = FALSE)
Seven.markers <- FindConservedMarkers(data.combined, ident.1 = 7, grouping.var = "KO", verbose = FALSE)
Eight.markers <- FindConservedMarkers(data.combined, ident.1 = 8, grouping.var = "KO", verbose = FALSE)
Nine.markers <- FindConservedMarkers(data.combined, ident.1 = 9, grouping.var = "KO", verbose = FALSE)
Ten.markers <- FindConservedMarkers(data.combined, ident.1 = 10, grouping.var = "KO", verbose = FALSE)
Eleven.markers <- FindConservedMarkers(data.combined, ident.1 = 11, grouping.var = "KO", verbose = FALSE)
Twelve.markers <- FindConservedMarkers(data.combined, ident.1 = 12, grouping.var = "KO", verbose = FALSE)
Thirteen.markers <- FindConservedMarkers(data.combined, ident.1 = 13, grouping.var = "KO", verbose = FALSE)
CD8T.markers <- FindMarkers(data.combined, ident.1 = "CD8+T", grouping.var = "KO", verbose = FALSE)

#no pbmc in 5, 11-13

write.csv(Zero.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C0Markers.csv")
write.csv(One.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C1Markers.csv")
write.csv(Two.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C2Markers.csv")
write.csv(Three.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C3Markers.csv")
write.csv(Four.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C4Markers.csv")
write.csv(Five.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C5Markers.csv")
write.csv(Six.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C6Markers.csv")
write.csv(Seven.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C7Markers.csv")
write.csv(Eight.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C8Markers.csv")
write.csv(Nine.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C9Markers.csv")
write.csv(Ten.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C10Markers.csv")
write.csv(Eleven.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C11Markers.csv")
write.csv(Twelve.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C12Markers.csv")
write.csv(Thirteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/C13Markers.csv")
write.csv(CD8T.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CD8TMarkers.csv")

saveRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster.rds")
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster.rds")

data.combined <- RenameIdents(data.combined, `0` = "Microglia", `1` = "Microglia", `2` = "Microglia", `3` = "B cells", 
                                  `4` = "Microglia", `5` = "UNK", `6` = "CD8_T", `7` = "CD4_T", `8` = "CD4_T", `9` = "NK cells", 
                                  `10` = "DC", `11` = "Brain cells1", `12` = "Brain cells2", `13` = "Microglia")

DimPlot(data.combined, label = TRUE)
DefaultAssay(data.combined) <- "integrated"
saveRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster6.rds")
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster6.rds")

p1 <- DimPlot(data.combined, reduction = "umap", group.by = "KO")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(data.combined, reduction = "umap", split.by = "KO", label = F)
DotPlot(data.combined, features=c("S100A4", "IL7R", "CD8A", "CD8B", "LYZ", "CD14", "FCGR3A", "FCER1G","CST3", "NKG7", "GZMB"), split.by="KO", cols = c("yellow","blue", "red"))
DotPlot(data.combined, features=c("S100A4", "IL7R", "CD8A", "CD8B", "LYZ", "CD14", "FCGR3A", "FCER1G","CST3", "NKG7", "GZMB"))


table(Idents(data.combined))
#Microglia  UNK     CD8+T   CD4+T       DC 
#8909       649     763      1226       217 
prop.table(table(Idents(data.combined)))
#Microglia        UNK      CD8+T       CD4+T         DC 
#0.75731044   0.05516831  0.06485889 0.10421625  0.01844611 

#Split by KO (CD8+T, CD4+T Genes)
plots <- VlnPlot(data.combined, features = c("CD8A", "CD8B", "NKG7", "IL7R", "S100A4"), split.by = "KO",  
                 pt.size = 0, combine = FALSE, ncol = 1, cols = c("yellow","blue", "red"))
wrap_plots(plots = plots, ncol = 1)


#DEG with Heatmap
library(dplyr)
DefaultAssay(data.combined) <- "integrated"
data.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(subset(data.combined, downsample = 50), features = top20$gene)

saveRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster2.rds")
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster2.rds")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
cd8t <- subset(data.combined, idents = "CD8+T")
Idents(cd8t) <- "KO"
avg.cd8t <- as.data.frame(log1p(AverageExpression(cd8t, verbose = FALSE)$RNA))
avg.cd8t$gene <- rownames(avg.cd8t)
#write.csv(avg.cd8t, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/avgCD8TExp.csv")

cd4t <- subset(data.combined, idents = "CD4+T")
Idents(cd4t) <- "KO"
avg.cd4t <- as.data.frame(log1p(AverageExpression(cd4t, verbose = FALSE)$RNA))
avg.cd4t$gene <- rownames(avg.cd4t)
#write.csv(avg.cd4t, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/avgCD4TExp.csv")

#Genes with similar expression levels reflect conserved pathway in both cell types
#DEG with scatterplots

#CD4T: CD vs MTS
genes.to.label = c("CD8A", "CD8B", "NKG7", "IL7R", "S100A4")
p1 <- ggplot(avg.cd8t, aes(cd, mts)) + geom_point() + ggtitle("CD8T") + geom_text(aes(label=ifelse(cd>3,as.character(gene),'')),hjust=0,vjust=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p1
#CD vs PBMC
genes.to.label = c("CD8A", "CD8B", "NKG7", "IL7R", "S100A4")
p1 <- ggplot(avg.cd8t, aes(cd, pbmc)) + geom_point() + ggtitle("CD8T") + geom_text(aes(label=ifelse(cd>3,as.character(gene),'')),hjust=0,vjust=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p1
#mst vs pbmc
genes.to.label = c("CD8A", "CD8B", "NKG7", "IL7R", "S100A4")
p1 <- ggplot(avg.cd8t, aes(mts, pbmc)) + geom_point() + ggtitle("CD8T") + geom_text(aes(label=ifelse(mts>3,as.character(gene),'')),hjust=0,vjust=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p1

saveRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster3.rds")
readRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/AfterConservedCluster3.rds")

#DEGs of CD8T in each condition
DefaultAssay(data.combined) <- "RNA"
data.combined$celltype.KO <- paste(Idents(data.combined), data.combined$KO, sep = "_")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.KO"

CD8T.cd.mts.response <- FindMarkers(data.combined, ident.1 = "CD8+T_cd", ident.2 = "CD8+T_mts", verbose = FALSE)
#write.csv(CD8T.cd.mts.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/CD8T.cd.mts.response.csv")
CD8T.cd.pbmc.response <- FindMarkers(data.combined, ident.1 = "CD8+T_cd", ident.2 = "CD8+T_pbmc", verbose = FALSE)
#write.csv(CD8T.cd.pbmc.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/CD8T.cd.pbmc.response.csv")
CD8T.mts.pbmc.response <- FindMarkers(data.combined, ident.1 = "CD8+T_mts", ident.2 = "CD8+T_pbmc", verbose = FALSE)
#write.csv(CD8T.mts.pbmc.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_Int/CSV/CD8T.mts.pbmc.response.csv")

saveRDS(data.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/microglia_data/PBMC_CD_MTS_IntAfterConservedCluster4.rds")

