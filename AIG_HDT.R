#Plan: To integrated AIG and Infection dataset and analyze the SPEM only clusters to discover unique DEGs, understand its pseudotime
#and understand the differences between the spem clusters

library(Seurat)
library(patchwork)
library(ggplot2)

TxA23 <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21/TestData/TxA23")
HDT <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21/TestData/HDT")

TxA23 <- CreateSeuratObject(counts = TxA23, project = "TxA23", min.cells = 4, min.features = 300)
TxA23$KO <- "TxA23"
TxA23[["percent.mt"]] <- PercentageFeatureSet(TxA23, pattern = "^mt-")
head(table(Idents(TxA23))) #12721

HDT <- CreateSeuratObject(counts = HDT, project = "HDT", min.cells = 4, min.features = 300)
HDT$KO <- "HDT"
HDT[["percent.mt"]] <- PercentageFeatureSet(HDT, pattern = "^mt-")
head(table(Idents(HDT))) #10262

VlnPlot(TxA23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(HDT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot4 <- FeatureScatter(HDT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(HDT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 + plot5

TxA23 <- subset(TxA23, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 40000 & percent.mt <=15 & Epcam > 0 & Pecam1 <= 0 & Ptprc <=0)
plot7 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot8 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot7 + plot8
head(table(Idents(TxA23))) #5066

HDT <- subset(HDT, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 50000 & percent.mt <=15 & Epcam > 0 & Pecam1 <= 0 & Ptprc <=0)
plot11 <- FeatureScatter(HDT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot12 <- FeatureScatter(HDT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot11 + plot12
head(table(Idents(HDT))) #7185

VlnPlot(TxA23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(HDT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data.list <- list(TxA23, HDT)
data.list <- lapply(X= data.list, FUN = function(x) { 
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)

diseased.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
diseased.combined <- IntegrateData(anchorset = diseased.anchors)
DefaultAssay(diseased.combined) <- "integrated"

all.genes <- rownames(diseased.combined)
# Run the standard workflow for visualization and clustering
diseased.combined <- ScaleData(diseased.combined, verbose = FALSE, features = all.genes)
diseased.combined <- RunPCA(diseased.combined, npcs = 30, verbose = FALSE)
diseased.combined <- RunUMAP(diseased.combined, reduction = "pca", dims = 1:20)
ElbowPlot(diseased.combined, ndims = 30)
#write.csv(all.genes, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/Allgenes.csv")
diseased.combined <- FindNeighbors(diseased.combined, reduction = "pca", dims = 1:20)
diseased.combined <- FindClusters(diseased.combined, resolution = 0.6)
DimPlot(diseased.combined, reduction = "umap", label = TRUE)

table(Idents(diseased.combined))
p1 <- DimPlot(diseased.combined, reduction = "umap", group.by = "KO")
p2 <- DimPlot(diseased.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1
p1 + p2
DimPlot(diseased.combined, reduction = "umap", split.by = "KO", label = TRUE)

saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/Trial1Cluster.rds")
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/Trial1Cluster.rds")

#Visualization
library(metap)
plot_paper<-VlnPlot(diseased.combined, features=c("Gif", "Muc6", "Muc5ac", "Tff2", "Gkn3", "Stmn1", "Birc5", "Mki67"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')

FeaturePlot(diseased.combined, features = c("Tff2", "Stmn1", "Gkn3", "Birc5"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey","red"), label = TRUE)
FeaturePlot(diseased.combined, features = c("Gif","Tff2","Muc5ac","Stmn1","Muc6","Gkn3","Mki67","Birc5"), max.cutoff = 3, 
            cols = c("grey","red"), label = TRUE, ncol = 4)
DotPlot(diseased.combined, features=c("Gif", "Muc6", "Muc5ac", "Tff2", "Gkn3", "Stmn1", "Birc5", "Mki67"), split.by="KO")
DotPlot(diseased.combined, features=c("Gif", "Muc6", "Muc5ac", "Tff2", "Gkn3", "Stmn1", "Birc5", "Mki67"))

#The following requested variables were not found:  Wfdc4, Popdc1 
#red represents high expression level

DefaultAssay(diseased.combined) <- "RNA"
diseased.markers <- FindAllMarkers(diseased.combined, logfc.threshold = 1.1, only.pos = TRUE )
#write.csv(diseased.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/AllMarkers.csv")

cluster0.markers <- FindMarkers(diseased.combined, ident.1 = 0, logfc.threshold = 1.1, test.use = "roc")
cluster1.markers <- FindMarkers(diseased.combined, ident.1 = 1, logfc.threshold = 1.1, test.use = "roc")
cluster2.markers <- FindMarkers(diseased.combined, ident.1 = 2, logfc.threshold = 1.1, test.use = "roc")
cluster3.markers <- FindMarkers(diseased.combined, ident.1 = 3, logfc.threshold = 1.1, test.use = "roc")
cluster4.markers <- FindMarkers(diseased.combined, ident.1 = "Neck", logfc.threshold = 1.1)
cluster5.markers <- FindMarkers(diseased.combined, ident.1 = 5, logfc.threshold = 1.1, test.use = "roc")
cluster6.markers <- FindMarkers(diseased.combined, ident.1 = 6, logfc.threshold = 1.1, test.use = "roc")
cluster7.markers <- FindMarkers(diseased.combined, ident.1 = 7, logfc.threshold = 1.1, test.use = "roc")
cluster8.markers <- FindMarkers(diseased.combined, ident.1 = 8, logfc.threshold = 1.1, test.use = "roc")
cluster9.markers <- FindMarkers(diseased.combined, ident.1 = 9, logfc.threshold = 1.1, test.use = "roc")
cluster10.markers <- FindMarkers(diseased.combined, ident.1 = 10, logfc.threshold = 1.1, test.use = "roc")
cluster11.markers <- FindMarkers(diseased.combined, ident.1 = 11, logfc.threshold = 1.1, test.use = "roc")
cluster12.markers <- FindMarkers(diseased.combined, ident.1 = 12, logfc.threshold = 1.1, test.use = "roc")
cluster13.markers <- FindMarkers(diseased.combined, ident.1 = 13, logfc.threshold = 1.1, test.use = "roc")
cluster14.markers <- FindMarkers(diseased.combined, ident.1 = 14, logfc.threshold = 1.1, test.use = "roc")
cluster15.markers <- FindMarkers(diseased.combined, ident.1 = 15, logfc.threshold = 1.1, test.use = "roc")
cluster16.markers <- FindMarkers(diseased.combined, ident.1 = 16, logfc.threshold = 1.1, test.use = "roc")
cluster17.markers <- FindMarkers(diseased.combined, ident.1 = 17, logfc.threshold = 1.1, test.use = "roc")
cluster18.markers <- FindMarkers(diseased.combined, ident.1 = 18, logfc.threshold = 1.1, test.use = "roc")

write.csv(cluster0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster0.markers.csv")
write.csv(cluster1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster1.markers.csv")
write.csv(cluster2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster2.markers.csv")
write.csv(cluster3.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster3.markers.csv")
write.csv(cluster4.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster4.markers.csv")
write.csv(cluster5.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster5.markers.csv")
write.csv(cluster6.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster6.markers.csv")
write.csv(cluster7.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster7.markers.csv")
write.csv(cluster8.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster8.markers.csv")
write.csv(cluster9.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster9.markers.csv")
write.csv(cluster10.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster10.markers.csv")
write.csv(cluster11.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster11.markers.csv")
write.csv(cluster12.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster12.markers.csv")
write.csv(cluster13.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster13.markers.csv")
write.csv(cluster14.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster14.markers.csv")
write.csv(cluster15.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster15.markers.csv")
write.csv(cluster16.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster16.markers.csv")
write.csv(cluster17.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster17.markers.csv")
write.csv(cluster18.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/cluster18.markers.csv")

Zero.markers <- FindConservedMarkers(diseased.combined, ident.1 = 0, grouping.var = "KO", verbose = FALSE)
One.markers <- FindConservedMarkers(diseased.combined, ident.1 = 1, grouping.var = "KO", verbose = FALSE)
Two.markers <- FindConservedMarkers(diseased.combined, ident.1 = 2, grouping.var = "KO", verbose = FALSE)
Three.markers <- FindConservedMarkers(diseased.combined, ident.1 = 3, grouping.var = "KO", verbose = FALSE)
Four.markers <- FindConservedMarkers(diseased.combined, ident.1 = 4, grouping.var = "KO", verbose = FALSE)
Five.markers <- FindConservedMarkers(diseased.combined, ident.1 = 5, grouping.var = "KO", verbose = FALSE)
Six.markers <- FindConservedMarkers(diseased.combined, ident.1 = 6, grouping.var = "KO", verbose = FALSE)
Seven.markers <- FindConservedMarkers(diseased.combined, ident.1 = 7, grouping.var = "KO", verbose = FALSE)
Eight.markers <- FindConservedMarkers(diseased.combined, ident.1 = 8, grouping.var = "KO", verbose = FALSE)
Nine.markers <- FindConservedMarkers(diseased.combined, ident.1 = 9, grouping.var = "KO", verbose = FALSE)
Ten.markers <- FindConservedMarkers(diseased.combined, ident.1 = 10, grouping.var = "KO", verbose = FALSE)
Eleven.markers <- FindConservedMarkers(diseased.combined, ident.1 = 11, grouping.var = "KO", verbose = FALSE)
Twelve.markers <- FindConservedMarkers(diseased.combined, ident.1 = 12, grouping.var = "KO", verbose = FALSE)
Thirteen.markers <- FindConservedMarkers(diseased.combined, ident.1 = 13, grouping.var = "KO", verbose = FALSE)
Fourteen.markers <- FindConservedMarkers(diseased.combined, ident.1 = 14, grouping.var = "KO", verbose = FALSE)
Fifteen.markers <- FindConservedMarkers(diseased.combined, ident.1 = 15, grouping.var = "KO", verbose = FALSE)
Sixteen.markers <- FindConservedMarkers(diseased.combined, ident.1 = 16, grouping.var = "KO", verbose = FALSE)
Seventeen.markers <- FindConservedMarkers(diseased.combined, ident.1 = 17, grouping.var = "KO", verbose = FALSE)
Eighteen.markers <- FindConservedMarkers(diseased.combined, ident.1 = 18, grouping.var = "KO", verbose = FALSE)

write.csv(Zero.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C0Markers.csv")
write.csv(One.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C1Markers.csv")
write.csv(Two.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C2Markers.csv")
write.csv(Three.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C3Markers.csv")
write.csv(Four.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C4Markers.csv")
write.csv(Five.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C5Markers.csv")
write.csv(Six.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C6Markers.csv")
write.csv(Seven.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C7Markers.csv")
write.csv(Eight.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C8Markers.csv")
write.csv(Nine.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C9Markers.csv")
write.csv(Ten.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C10Markers.csv")
write.csv(Eleven.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C11Markers.csv")
write.csv(Twelve.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C12Markers.csv")
write.csv(Thirteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C13Markers.csv")
write.csv(Fourteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C14Markers.csv")
write.csv(Fifteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C15Markers.csv")
write.csv(Sixteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C16Markers.csv")
write.csv(Seventeen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C17Markers.csv")
write.csv(Eighteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/C18Markers.csv")

saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster.rds")

readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster.rds")
diseased.combined <- RenameIdents(diseased.combined, `0` = "Pit", `1` = "Chief", `2` = "Stem", `3` = "Parietal", 
                                  `4` = "Neck", `5` = "Pit", `6` = "Stem", `7` = "Pit", `8` = "SPEM", `9` = "NE", 
                                  `10` = "Parietal", `11` = "SPEM", `12` = "Parietal", `13` = "NE", `14` = "FS", 
                                  `15` = "Fibroblast", `16` ="Tuft", `17` = "UNK", `18` = "HPC")
#NE = Neuroendocrine, HPC = Hematopoietic cells, SM = Smooth Muscle, FS = Forestomach
DimPlot(diseased.combined, label = TRUE, label.size = )

DefaultAssay(diseased.combined) <- "integrated"
saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster6.rds")
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster6.rds")

p1 <- DimPlot(diseased.combined, reduction = "umap", group.by = "KO", cols = c("blue", "red"))
p2 <- DimPlot(diseased.combined, reduction = "umap", label = TRUE, repel = FALSE)
p1 + p2
DimPlot(diseased.combined, reduction = "umap", split.by = "KO")

#visualization
DefaultAssay(diseased.combined) <- "RNA"
plot_paper<-VlnPlot(diseased.combined, features=c("Gif", "Muc6", "Muc5ac", "Tff2", "Gkn3", "Stmn1", "Birc5", "Mki67"), pt.size = 0, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')

FeaturePlot(diseased.combined, features = c("Tff2”, “Stmn1”, “Gkn3”, “Birc5"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey","red"), label = TRUE)
FeaturePlot(diseased.combined, features = c("Muc6","Gkn3","Mki67","Birc5"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey","red"), label = TRUE)
FeaturePlot(diseased.combined, features = c("Gif","Tff2","Muc5ac","Stmn1","Muc6","Gkn3","Mki67","Birc5"), max.cutoff = 3, 
            cols = c("grey","red"), label = TRUE, ncol = 4)
DotPlot(diseased.combined, features=c("Gif", "Muc6", "Muc5ac", "Tff2", "Gkn3", "Stmn1", "Birc5", "Mki67"), split.by="KO")
DotPlot(diseased.combined, features=c("Gif", "Muc6", "Muc5ac", "Tff2", "Gkn3", "Stmn1", "Birc5", "Mki67"))


table(Idents(diseased.combined))
prop.table(table(Idents(diseased.combined)))
SPEM <-  subset(diseased.combined, idents = "SPEM")
RidgePlot(diseased.combined, features = c("Tff2", "Gkn3"))
RidgePlot(SPEM, features = c("Tff2", "Gkn3"), group.by = "KO", cols = c("yellow", "blue"))
VlnPlot(diseased.combined, features = c("Gif","Tff2", "Gkn3", "Muc5ac", "Mki67", "Stmn1"), ncol = 3, pt.size = 0, group.by = "KO")
VlnPlot(diseased.combined, features = c("Iqgap3", "Mki67", "Stmn1", "Birc5"), ncol = 1, pt.size = 0.1)
VlnPlot(diseased.combined, features = c("Hdc", "Dclk1", "Dcn", "Vim", "Myl9", "Ceacam10"), ncol = 1, pt.size = 0)
VlnPlot(diseased.combined, features = c("Cd44", "Wfdc2", "Wfdc18", "Wfdc4", "Klk1", "Aqp5"), ncol = 1, pt.size = 0)
VlnPlot(diseased.combined, features = c("Bmi1", "Tnfrsf19", "Lrig1", "Bhlha15"), ncol = 1, pt.size = 0)
VlnPlot(diseased.combined, features = c("Tff3", "Tacstd2", "Muc2"), ncol = 1, pt.size = 0)
VlnPlot(diseased.combined, features = c("Barx1", "Sox2", "Popdc1", "Dmbt1", "Trp63"), ncol = 1, pt.size = 0)
#The following requested variables were not found:  Wfdc4, Popdc1 

#Split by KO
plots <- VlnPlot(diseased.combined, features = c("Gif", "Tff2", "Gkn3", "Muc5ac", "Stmn1", "Mki67"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 3)
wrap_plots(plots = plots, ncol = 3)
plots <- VlnPlot(diseased.combined, features = c("Atp4b", "Atp4a", "Muc5ac", "Mki67", "Stmn1", "Birc5"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(diseased.combined, features = c("Hdc", "Dclk1", "Dcn", "Vim", "Myl9", "Ceacam10"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(diseased.combined, features = c("Cd44", "Wfdc2", "Wfdc18", "Wfdc4", "Klk1", "Aqp5"), split.by = "KO",  
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(diseased.combined, features = c("Bmi1", "Tnfrsf19", "Lrig1", "Bhlha15"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(diseased.combined, features = c("Tff3", "Cd44", "Tacstd2", "Muc2"), split.by = "KO",  
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(diseased.combined, features = c("Barx1", "Sox2", "Popdc1", "Popdc3", "Trp63"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(diseased.combined, features = c("Lgr5", "Bmi1", "Dmbt1", "Popdc3", "Tff1"), split.by = "KO",
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
#The following requested variables were not found:  Wfdc4, Popdc1 


FeaturePlot(diseased.combined, features = c("Gkn3", "Tff2", "Stmn1"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(diseased.combined, features = c("Mki67", "Sox4", "Lgr5", "Agr2"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(diseased.combined, features = c("Birc5", "Dmbt1", "Bhlha15"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)

#DEG with Heatmap

library(dplyr)
DefaultAssay(diseased.combined) <- "integrated"
diseased.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10
DoHeatmap(subset(diseased.combined, downsample = 100), features = top10$gene)

#SUBSETFORPAPER
DefaultAssay(diseased.combined) <- "RNA"
dotPlot.cells <- subset(diseased.combined, idents = c("Stem","SPEM"))
levels(dotPlot.cells) <- c("SPEM", "Stem")
DotPlot(dotPlot.cells, features=c("Stmn1", "Birc5", "Mki67", "Gif", "Muc6", "Tff2", "Gkn3"))
levels(dotPlot.cells)
Idents(dotPlot.cells)

saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/AIG_HDT_Int/AfterConservedCluster2.rds")
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/AIG_HDT_Int/AfterConservedCluster2.rds")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
SPEM.cells <- subset(diseased.combined, idents = "SPEM")
Idents(SPEM.cells) <- "KO"
avg.SPEM.cells <- as.data.frame(log1p(AverageExpression(SPEM.cells, verbose = FALSE)$RNA))
avg.SPEM.cells$gene <- rownames(avg.SPEM.cells)
#write.csv(avg.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/avgSpemExp.csv")

stem.cells <- subset(diseased.combined, idents = "Stem")
Idents(stem.cells) <- "KO"
avg.stem.cells <- as.data.frame(log1p(AverageExpression(stem.cells, verbose = FALSE)$RNA))
avg.stem.cells$gene <- rownames(avg.stem.cells)
#write.csv(avg.stem.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/avgStemExp.csv")

#Genes with similar expression levels reflect conserved pathway in both cell types
#DEG with scatterplots
#AIGvsHDT
genes.to.label = c("Gkn3", "Tff2", "Stmn1", "Mki67", "Birc5", "Spink4", "Sox4", "Sox2", "Lgr5", "Chil4", "Dmbt1", "Bhlha15")
p1 <- ggplot(avg.SPEM.cells, aes(HDT, TxA23)) + geom_point() + ggtitle("SPEM") + geom_text(aes(label=ifelse(TxA23>3,as.character(gene),'')),hjust=0,vjust=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p2 <- ggplot(avg.stem.cells, aes(HDT, AIG)) + geom_point() + ggtitle("Stem") + geom_text(aes(label=ifelse(AIG>3,as.character(gene),'')),hjust=0,vjust=0)
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p1 + p2


#Stem vs spem comparison
stem.spem.combined <- subset(diseased.combined, idents = c("SPEM", "Stem"))
VlnPlot(stem.spem.combined, features = c("Gkn3", "Tff2", "Stmn1"), pt.size = 1, ncol = 1)
VlnPlot(stem.spem.combined, features = c("Aqp5", "Klk1", "Birc5", "Mki67"), pt.size = 1, ncol = 1)

stem.vs.spem.markers <- FindMarkers(diseased.combined, ident.1 = "Stem", ident.2 = "SPEM", verbose = FALSE)
#write.csv(stem.vs.spem.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/Stem.vs.Spem.markers.csv")

FindMarkers(diseased.combined, ident.1 = "Stem", ident.2 = "SPEM", verbose = FALSE, 
            features = c("Stmn1", "Gkn3", "Tff2", "Aqp5", "Klk1", "Birc5", "Mki67"), logfc.threshold = 0, min.pct = 0)

diseased.combined$celltype <- Idents(diseased.combined)
Idents(diseased.combined) <-"celltype"
avg.meta <- log1p(AverageExpression(diseased.combined, verbose = FALSE)$RNA)
avg.meta<-as.data.frame(avg.meta)
avg.meta$gene <- rownames(avg.meta)
meta_comp<-ggplot(avg.meta, aes(Stem, SPEM)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Stem vs SPEM Average Gene Expression") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5))+ xlab("Stem") + ylab("SPEM")
meta_comp<-LabelPoints(plot = meta_comp, points = c("Stmn1", "Gkn3", "Tff2", "Aqp5", "Klk1", "Birc5", "Mki67"), color="red", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(meta_comp)

saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster3.rds")
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster3.rds")

DefaultAssay(diseased.combined) <- "RNA"
diseased.combined$celltype.KO <- paste(Idents(diseased.combined), diseased.combined$KO, sep = "_")
diseased.combined$celltype <- Idents(diseased.combined)
Idents(diseased.combined) <- "celltype.KO"
SPEM.response <- FindMarkers(diseased.combined, ident.1 = "SPEM_HDT", ident.2 = "SPEM_AIG", verbose = FALSE)
write.csv(SPEM.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/SPEMResponse.csv")

stem.response <- FindMarkers(diseased.combined, ident.1 = "Stem_HDT", ident.2 = "Stem_AIG", verbose = FALSE)
write.csv(stem.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/CSV/Stemresponse.csv")

saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster4.rds")

#TrajectoryAnalysis
library(patchwork)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
DefaultAssay(diseased.combined) <- "RNA"
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterConservedCluster4.rds")
DimPlot(diseased.combined, reduction = "umap", split.by = "KO")

diseased.combined$celltype <- Idents(diseased.combined)
AIG <- diseased.combined[, diseased.combined$celltype.KO %in% c("Neck_AIG", "Chief_AIG", "Stem_AIG", "SPEM_AIG")]
HDT <- diseased.combined[,  diseased.combined$celltype.KO %in% c("Neck_HDT", "Chief_HDT", "Stem_HDT", "SPEM_HDT")]
AIG.HDT <- diseased.combined[, diseased.combined$celltype.KO %in% c("Neck_AIG", "Chief_AIG", "Stem_AIG", "SPEM_AIG", "Neck_HDT", "Chief_HDT", "Stem_HDT", "SPEM_HDT")]

AIG.cds <- as.cell_data_set(AIG)
## Calculate size factors using built-in function in monocle3
AIG.cds <- estimate_size_factors(AIG.cds)
## Add gene names into CDS
AIG.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(AIG.cells[["RNA"]])

AIG.cds <- preprocess_cds(AIG.cds, reduction_method = "PCA")
AIG.cds <- reduce_dimension(AIG.cds, reduction_method = "UMAP") #also preprocess_cds with PCA by default
AIG.cds <- cluster_cells(cd = AIG.cds, reduction_method = "UMAP")
AIG.cds <- learn_graph(AIG.cds, use_partition = TRUE)
colData(AIG.cds)


HDT.cds <- as.cell_data_set(HDT)
## Calculate size factors using built-in function in monocle3
HDT.cds <- estimate_size_factors(HDT.cds)
## Add gene names into CDS
HDT.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(HDT.cells[["RNA"]])

HDT.cds <- preprocess_cds(HDT.cds, reduction_method = "PCA")
HDT.cds <- reduce_dimension(HDT.cds, reduction_method = "UMAP")
HDT.cds <- cluster_cells(cds = HDT.cds, reduction_method = "UMAP")
HDT.cds <- learn_graph(HDT.cds, use_partition = TRUE)
colData(HDT.cds)

AIG.HDT.cds <- as.cell_data_set(AIG.HDT)
## Calculate size factors using built-in function in monocle3
AIG.HDT.cds <- estimate_size_factors(AIG.HDT.cds)
## Add gene names into CDS
AIG.HDT.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(AIG.HDT.cells[["RNA"]])

AIG.HDT.cds <- preprocess_cds(AIG.HDT.cds, reduction_method = "PCA")
AIG.HDT.cds <- reduce_dimension(AIG.HDT.cds, reduction_method = "UMAP")
AIG.HDT.cds <- cluster_cells(cds = AIG.HDT.cds, reduction_method = "UMAP")
AIG.HDT.cds <- learn_graph(AIG.HDT.cds, use_partition = TRUE)
colData(AIG.HDT.cds)

plot_cells(
  cds = AIG.cds,
  color_cells_by = "celltype.KO", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE)


plot_cells(
  cds = HDT.cds,
  color_cells_by = "celltype.KO",label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE)

plot_cells(
  cds = AIG.HDT.cds,
  color_cells_by = "celltype.KO", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE)

AIG.cds <- order_cells(cds = AIG.cds, reduction = "UMAP")
HDT.cds <- order_cells(cds = HDT.cds, reduction = "UMAP")
AIG.HDT.cds <- order_cells(cds = AIG.HDT.cds, reduction = "UMAP")

plot_cells(
  cds = AIG.cds,
  color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = TRUE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = HDT.cds,
  color_cells_by = "pseudotime",label_cell_groups = FALSE,
  label_leaves = TRUE, label_branch_points = TRUE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = AIG.HDT.cds,
  color_cells_by = "pseudotime",label_cell_groups = FALSE,
  label_leaves = TRUE, label_branch_points = TRUE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = HDT.cds,
  genes = c("Gkn3", "Tff2", "Chil4", "Spp1", "Stmn1", "Birc5", "Mki67"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

plot_cells(
  cds = AIG.cds,
  genes = c("Gkn3", "Tff2", "Chil4", "Spp1", "Stmn1", "Birc5", "Mki67"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

plot_cells(
  cds = AIG.HDT.cds,
  genes = c("Gkn3", "Tff2", "Chil4", "Spp1", "Stmn1", "Birc5", "Mki67"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

diseased.combined <- AddMetaData(
  object = diseased.combined,
  metadata = AIG.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "AIG"
)

diseased.combined <- AddMetaData(
  object = diseased.combined,
  metadata = HDT.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "HDT"
)
diseased.combined <- AddMetaData(
  object = diseased.combined,
  metadata = HDT.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "AIG.HDT"
)
FeaturePlot(diseased.combined, c("AIG", "HDT", "AIG.HDT"), pt.size = 0.1) & scale_color_viridis_c()
saveRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterPseudotime1.rds")
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/AfterPseudotime1.rds")
FeaturePlot(diseased.combined, c("AIG", "HDT", "AIG.HDT"), label = TRUE, pt.size = 0.1)




#SPEMSTEMCLUSTER
DefaultAssay(int.SPEM.cells) <- "RNA"
SPEM.Stem.cells <- subset(int.SPEM.cells, idents = "Stem-SPEM")
SPEM.Stem.cells$celltype <- "SPEM-Stem"
plot1 <- FeatureScatter(SPEM.Stem.cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SPEM.Stem.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
head(table(Idents(SPEM.Stem.cells)))

VlnPlot(int.Stem.cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(SPEM.Stem.cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data.list <- list(int.Stem.cells, SPEM.Stem.cells)
data.list <- lapply(X= data.list, FUN = function(x) { 
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)

STEM.SPEM.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
STEM.SPEM.combined <- IntegrateData(anchorset = diseased.anchors)
DefaultAssay(STEM.SPEM.combined) <- "integrated"

all.genes <- rownames(STEM.SPEM.combined)
# Run the standard workflow for visualization and clustering
diseased.combined <- ScaleData(STEM.SPEM.combined, verbose = FALSE, features = all.genes)
diseased.combined <- RunPCA(STEM.SPEM.combined, npcs = 30, verbose = FALSE)
diseased.combined <- RunUMAP(STEM.SPEM.combined, reduction = "pca", dims = 1:20)
ElbowPlot(STEM.SPEM.combined, ndims = 30)
write.csv(all.genes, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_HDT_Int/STEMSPEM/CSV/Allgenes.csv")
diseased.combined <- FindNeighbors(STEM.SPEM.combined, reduction = "pca", dims = 1:20)
diseased.combined <- FindClusters(STEM.SPEM.combined, resolution = 0.6)

p1 <- DimPlot(STEM.SPEM.combined, reduction = "umap", group.by = "KO")
p2 <- DimPlot(STEM.SPEM.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(STEM.SPEM.combined, reduction = "umap", split.by = "KO", label = TRUE)
DimPlot(STEM.SPEM.combined, reduction = "umap", label = TRUE)



