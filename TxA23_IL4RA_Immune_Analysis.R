#Plan: To integrated TxA23 and TxA23 IL4rA knockout dataset, and anlyze the immune cells for DEGs in celltypes

library(Seurat)
library(patchwork)
library(ggplot2)
library(metap)

TxA23 <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21/TestData/TxA23")
IL4RA <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21/TestData/IL4RA")
Tx23_Immune <- Read10X("/Users/dipaololabmacbook/Desktop/scRNA21/TestData/TxA23_Immune")

TxA23 <- CreateSeuratObject(counts = TxA23, project = "TxA23", min.cells = 4, min.features = 300)
TxA23$KO <- "TxA23"
TxA23[["percent.mt"]] <- PercentageFeatureSet(TxA23, pattern = "^mt-")
head(table(Idents(TxA23))) #12721

IL4RA <- CreateSeuratObject(counts = IL4RA, project = "IL4RA", min.cells = 4, min.features = 300)
IL4RA$KO <- "KO"
IL4RA[["percent.mt"]] <- PercentageFeatureSet(IL4RA, pattern = "^mt-")
head(table(Idents(IL4RA)))#9042

Tx23_Immune <- CreateSeuratObject(counts = Tx23_Immune, project = "Tx23_Immune", min.cells = 4, min.features = 300)
Tx23_Immune$KO <- "TxA23"
Tx23_Immune[["percent.mt"]] <- PercentageFeatureSet(Tx23_Immune, pattern = "^mt-")
head(table(Idents(Tx23_Immune))) #6502

VlnPlot(TxA23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(IL4RA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Tx23_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot4 <- FeatureScatter(IL4RA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(IL4RA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 + plot5

plot4 <- FeatureScatter(Tx23_Immune, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(Tx23_Immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 + plot5

TxA23 <- subset(TxA23, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 40000 & percent.mt <=15 & Epcam <= 0 & Pecam1 <= 0 & Ptprc >0 & Dcn <= 0)
plot7 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot8 <- FeatureScatter(TxA23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot7 + plot8
head(table(Idents(TxA23))) #689

IL4RA <- subset(IL4RA, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 50000 & percent.mt <=15 & Epcam <= 0 & Pecam1 <= 0 & Ptprc >0 & Dcn <= 0)
plot11 <- FeatureScatter(IL4RA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot12 <- FeatureScatter(IL4RA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot11 + plot12
head(table(Idents(IL4RA))) #1028

Tx23_Immune <- subset(Tx23_Immune, subset = nFeature_RNA > 500 & nCount_RNA >= 100 & nCount_RNA < 20000 & percent.mt <=15 & Epcam <= 0 & Pecam1 <= 0 & Ptprc >0 & Dcn <= 0)
plot11 <- FeatureScatter(Tx23_Immune, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot12 <- FeatureScatter(Tx23_Immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot11 + plot12
head(table(Idents(Tx23_Immune))) #2859

VlnPlot(TxA23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(IL4RA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Tx23_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data.list <- list(TxA23, IL4RA, Tx23_Immune)
data.list <- lapply(X= data.list, FUN = function(x) { 
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)

diseased.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = diseased.anchors)
DefaultAssay(immune.combined) <- "integrated"

all.genes <- rownames(immune.combined)
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE, features = all.genes)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
ElbowPlot(immune.combined, ndims = 30)
#write.csv(all.genes, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Allgenes.csv")
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.6)
DimPlot(immune.combined, reduction = "umap", label = TRUE)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "KO")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(immune.combined, reduction = "umap", split.by = "KO", label = TRUE)

saveRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/Trial1Cluster.rds")
readRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/Trial1Cluster.rds")

VlnPlot(immune.combined, features = c("Cd3e", "Cd3g", "Cd79a", "Cd79b", "Cd19"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Tbx21", "Ifng", "Gata3"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Foxp3", "Cd4", "Rorc", "Il17a", "Il2rg", "Aighe"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Cd4", "Cd8a", "Cd8b1", "Nkg7", "Trdc", "Tcrg-C1"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Ncr1", "Nkg7", "Gzmb", "Siglech"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Cd68", "Lyz2", "C1qb", "Cd14"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Mcpt1", "Mcpt2", "Fcer1a", "Fcer1g", "Kit", "Hdc"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Il13", "Il13ra1", "Il4", "Il4ra" ,"Il5", "Il1rl1"), ncol = 1, pt.size = 1)

DefaultAssay(immune.combined) <- "RNA"
diseased.markers <- FindAllMarkers(immune.combined, logfc.threshold = 1.1, only.pos = TRUE )
#write.csv(diseased.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/AllMarkers.csv")

cluster0.markers <- FindMarkers(immune.combined, ident.1 = 0, logfc.threshold = 1.1)
cluster1.markers <- FindMarkers(immune.combined, ident.1 = 1, logfc.threshold = 1.1)
cluster2.markers <- FindMarkers(immune.combined, ident.1 = 2, logfc.threshold = 1.1)
cluster3.markers <- FindMarkers(immune.combined, ident.1 = 3, logfc.threshold = 1.1)
cluster4.markers <- FindMarkers(immune.combined, ident.1 = 4, logfc.threshold = 1.1)
cluster5.markers <- FindMarkers(immune.combined, ident.1 = 5, logfc.threshold = 1.1)
cluster6.markers <- FindMarkers(immune.combined, ident.1 = 6, logfc.threshold = 1.1)
cluster7.markers <- FindMarkers(immune.combined, ident.1 = 7, logfc.threshold = 1.1)
cluster8.markers <- FindMarkers(immune.combined, ident.1 = 8, logfc.threshold = 1.1 )
cluster9.markers <- FindMarkers(immune.combined, ident.1 = 9, logfc.threshold = 1.1)
cluster10.markers <- FindMarkers(immune.combined, ident.1 = 10, logfc.threshold = 1.1 )
cluster11.markers <- FindMarkers(immune.combined, ident.1 = 11, logfc.threshold = 1.1)
cluster12.markers <- FindMarkers(immune.combined, ident.1 = 12, logfc.threshold = 1.1)
cluster13.markers <- FindMarkers(immune.combined, ident.1 = 13, logfc.threshold = 1.1)
cluster14.markers <- FindMarkers(immune.combined, ident.1 = 14, logfc.threshold = 1.1)
cluster15.markers <- FindMarkers(immune.combined, ident.1 = 15, logfc.threshold = 1.1)
cluster16.markers <- FindMarkers(immune.combined, ident.1 = 16, logfc.threshold = 1.1)
cluster17.markers <- FindMarkers(immune.combined, ident.1 = 17, logfc.threshold = 1.1)
cluster18.markers <- FindMarkers(immune.combined, ident.1 = 18, logfc.threshold = 1.1)

write.csv(cluster0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster0.markers.csv")
write.csv(cluster1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster1.markers.csv")
write.csv(cluster2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster2.markers.csv")
write.csv(cluster3.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster3.markers.csv")
write.csv(cluster4.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster4.markers.csv")
write.csv(cluster5.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster5.markers.csv")
write.csv(cluster6.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster6.markers.csv")
write.csv(cluster7.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster7.markers.csv")
write.csv(cluster8.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster8.markers.csv")
write.csv(cluster9.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster9.markers.csv")
write.csv(cluster10.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster10.markers.csv")
write.csv(cluster11.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster11.markers.csv")
write.csv(cluster12.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster12.markers.csv")
write.csv(cluster13.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster13.markers.csv")
write.csv(cluster14.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster14.markers.csv")
write.csv(cluster15.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster15.markers.csv")
write.csv(cluster16.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster16.markers.csv")
write.csv(cluster17.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster17.markers.csv")
write.csv(cluster18.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cluster18.markers.csv")


Zero.markers <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "KO", verbose = FALSE)
One.markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "KO", verbose = FALSE)
Two.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "KO", verbose = FALSE)
Three.markers <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "KO", verbose = FALSE)
Four.markers <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "KO", verbose = FALSE)
Five.markers <- FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = "KO", verbose = FALSE)
Six.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "KO", verbose = FALSE)
Seven.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "KO", verbose = FALSE)
Eight.markers <- FindConservedMarkers(immune.combined, ident.1 = 8, grouping.var = "KO", verbose = FALSE)
Nine.markers <- FindConservedMarkers(immune.combined, ident.1 = 9, grouping.var = "KO", verbose = FALSE)
Ten.markers <- FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = "KO", verbose = FALSE)
Eleven.markers <- FindConservedMarkers(immune.combined, ident.1 = 11, grouping.var = "KO", verbose = FALSE)
Twelve.markers <- FindConservedMarkers(immune.combined, ident.1 = 12, grouping.var = "KO", verbose = FALSE)
Thirteen.markers <- FindConservedMarkers(immune.combined, ident.1 = 13, grouping.var = "KO", verbose = FALSE)
Fourteen.markers <- FindConservedMarkers(immune.combined, ident.1 = 14, grouping.var = "KO", verbose = FALSE)
Fifteen.markers <- FindConservedMarkers(immune.combined, ident.1 = 15, grouping.var = "KO", verbose = FALSE)
Sixteen.markers <- FindConservedMarkers(immune.combined, ident.1 = 16, grouping.var = "KO", verbose = FALSE)
Seventeen.markers <- FindConservedMarkers(immune.combined, ident.1 = 17, grouping.var = "KO", verbose = FALSE)
Eighteen.markers <- FindConservedMarkers(immune.combined, ident.1 = 18, grouping.var = "KO", verbose = FALSE)

write.csv(Zero.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C0Markers.csv")
write.csv(One.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C1Markers.csv")
write.csv(Two.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C2Markers.csv")
write.csv(Three.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C3Markers.csv")
write.csv(Four.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C4Markers.csv")
write.csv(Five.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C5Markers.csv")
write.csv(Six.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C6Markers.csv")
write.csv(Seven.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C7Markers.csv")
write.csv(Eight.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C8Markers.csv")
write.csv(Nine.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C9Markers.csv")
write.csv(Ten.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C10Markers.csv")
write.csv(Eleven.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C11Markers.csv")
write.csv(Twelve.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C12Markers.csv")
write.csv(Thirteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C13Markers.csv")
write.csv(Fourteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C14Markers.csv")
write.csv(Fifteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C15Markers.csv")
write.csv(Sixteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C16Markers.csv")
write.csv(Seventeen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C17Markers.csv")
write.csv(Eighteen.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/C18Markers.csv")

saveRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster.rds")

readRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster.rds")

immune.combined <- RenameIdents(immune.combined, `0` = "B_cells", `1` = "B_cells", `2` = "Treg", `3` = "CD4_T", 
                                `4` = "CD8_T", `5` = "Stem", `6` = "Macrophages", `7` = "B_cells", `8` = "CD4_T", `9` = "CD8_T", 
                                `10` = "B_cells", `11` = "CD4_T", `12` = "NK_cells", `13` = "Mast", `14` = "Epithelial_junk", 
                                `15` = "Macrophages", `16` = "ILC2/GD_T", `17` = "DC", `18` = "B_cells")
#immune.combined <- RenameIdents(immune.combined, `0` = "B_cells1", `1` = "B_cells2", `2` = "Treg", `3` = "CD4_T1", 
#                                  `4` = "CD8_T1", `5` = "Stem", `6` = "Macrophages", `7` = "B_cells3", `8` = "CD4_T*", `9` = "CD8_T2", 
#                                  `10` = "B_cells4", `11` = "CD4_T2", `12` = "NK_cells", `13` = "Mast", `14` = "Epithelial_junk", 
#                                `15` = "Macrophages", `16` = "ILC2/GD_T", `17` = "DC", `18` = "B_cells5")
#CD4_T*  has unconventional DEGs in csv file
DimPlot(immune.combined, label = FALSE)
DefaultAssay(immune.combined) <- "integrated"
saveRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster6.rds")
readRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster6.rds")

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "KO")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(immune.combined, reduction = "umap", split.by = "KO", label = FALSE)

table(Idents(immune.combined))
#B_cells          Treg           CD4_T           CD8_T            Stem     Macrophages          CD4_T*        NK_cells 
#2351             456             497             404             246             237             115              94 
#Mast      Epithelial_junk       ILC2/GD_T       DC 
#66              53              32              25 
prop.table(table(Idents(immune.combined)))
#B_cells1        B_cells2            Treg          CD4_T1          CD8_T1            Stem     Macrophages        B_cells3 
#0.305288462     0.140734266     0.099650350     0.088068182     0.064029720     0.053758741     0.051791958     0.039991259 
#CD4_T*          CD8_T2        B_cells4          CD4_T2        NK_cells            Mast Epithelial_junk       ILC2/GD_T 
#0.025131119     0.024256993     0.023819930     0.020541958     0.020541958     0.014423077     0.011582168     0.006993007 
#DC        B_cells5 
#0.005463287     0.003933566 
#DEGAnalysis
Bcells2vsBcells4.response <- FindMarkers(immune.combined, ident.1 = "B_cells2", ident.2 = "B_cells4", verbose = FALSE, logfc.threshold = 1.1)
#write.csv(Bcells2vsBcells4.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Bcells2vsBcells4.csv")

#DefaultAssay(immune.combined) <- "RNA"
#Idents(immune.combined) <- "KO"
#IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "IL4RA", ident.2 = "TxA23", verbose = FALSE, logfc.threshold = 1.1)
#write.csv(IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/IL4RAvsAIG.csv")

VlnPlot(immune.combined, features = c("Il4ra", "Il4"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Cd3e", "Cd3g", "Cd79a", "Cd79b", "Cd19"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Tbx21", "Gata3", "Foxp3", "Cd4", "Rorc"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Cd4", "Cd8a", "Cd8b1", "Nkg7", "Trdc", "Tcrg-C1"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Ncr1", "Nkg7", "Gzmb", "Cd4", "Cd8a", "Siglech"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Cd68", "Lyz2", "C1qb", "Cd14"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Mcpt1", "Mcpt2", "Fcer1a", "Fcer1g", "Hdc"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5", "Il18"), ncol = 1, pt.size = 1)
VlnPlot(immune.combined, features = c("Il4","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il13ra2","Il13","Il13ra1"), ncol = 1, pt.size = 1)

#Cytokines: "Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5"
#Receptors: "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1"

#Split by KO
plots <- VlnPlot(immune.combined, features = c("Il4ra", "Il4"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(immune.combined, features = c("Ncr1", "Nkg7", "Gzmb", "Cd4", "Cd8a", "Siglech"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)

#ApoptoticExpression
plots <- VlnPlot(immune.combined, features = c("Chia1", "Jun", "Nr4a1", "Dusp1"), group.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(immune.combined, features = c("Cd3g", "Ppp1r15a", "Fos", "Lgals1"), group.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
#FeaturePlot(immune.combined, features = c("Chia1", "Jun", "Nr4a1", "Dusp1"), split.by = "KO", max.cutoff = 3, 
#            cols = c("grey", "red"), label = FALSE)
#Antibody production
B.cells <- subset(immune.combined, idents = "B_cells")
plots <- VlnPlot(B.cells, features = c("Ighd", "Ighm", "Igkc", "Iglc3"), group.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 2)
wrap_plots(plots = plots, ncol = 2, label_value(labels = label, multi_line = FALSE))
plots <- VlnPlot(Bcells, features = c("Tipin", "Bnip3l", "Klrc3", "Ncoa3"), group.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)

#bcellsgenes
plots <- VlnPlot(immune.combined, features = c("H2-DMb2", "Cd74", "Cd79a", "Cd79b"), SPLIT.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(immune.combined, features = c("Ighm", "Ighd", "Igkc"), split.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)



#variantgenes
stressed.cells <- subset(immune.combined, idents = c("B_cells1","B_cells3", "B_cells5", "B_cells2", "B_cells4", "CD4_T*", "ILC2/GD_T", "DC", "Mast"))

plots <- VlnPlot(stressed.cells, features = c("Nr4a1", "Chia1", "Ccl5", "Jun", "Fos", "Lgals1", "Dusp1", "Tsc22d3"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 4)
wrap_plots(plots = plots, ncol = 4)
plots <- VlnPlot(stressed.cells, features = c( "Jun", "Fos", "Chia1", "Dusp1", "Ighm", "Ighd", "Hbb-bs", "H2-DMb2"), split.by = "KO", 
                 pt.size = 0.1, combine = FALSE, ncol = 4)
wrap_plots(plots = plots, ncol = 4)
plots <- VlnPlot(stressed.cells, features = c("Tnf", "Top2a", "Timp1", "Tap1"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(stressed.cells, features = c("Nr4a1", "Ighm", "Cd3g", "Cd3e"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(stressed.cells, features = c("Lars2", "Rps21", "Rps15"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(stressed.cells, features = c("Junb", "Pnrc1", "Ier2"), split.by = "KO", 
                 pt.size = 1, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(stressed.cells, features = c("Nr4a1", "Chia1", "Ccl5", "Jun", "Fos"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(immune.combined, features = c("Il4ra", "Il4", "Il13"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(immune.combined, features = c("Ncr1", "Nkg7", "Cd4", "Cd8a", "Siglech"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(immune.combined, features = c("Tbx21", "Gata3", "Foxp3", "Rorc", "Tcrg-C1"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(immune.combined, features = c("Cd3e", "Cd3g", "Cd79a", "Cd79b", "Cd19"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)

#DEG with Heatmap
library(dplyr)
DefaultAssay(immune.combined) <- "integrated"
diseased.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(immune.combined, downsample = 100), features = top10$gene)

Bcells <- subset(stressed.cells, idents = c("B_cells1", "B_cells2", "B_cells3", "B_cells4", "B_cells5"))
saveRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster2.rds")
readRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster2.rds")



library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
macro.cells <- subset(immune.combined, idents = "Macrophages")
Idents(macro.cells) <- "KO"
avg.macro.cells <- as.data.frame(log1p(AverageExpression(macro.cells, verbose = FALSE)$RNA))
avg.macro.cells$gene <- rownames(avg.macro.cells)
#write.csv(avg.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/avgMacroExp.csv")

cd8t1.cells <- subset(immune.combined, idents = "CD8_T1")
Idents(cd8t1.cells) <- "KO"
avg.cd8t1.cells <- as.data.frame(log1p(AverageExpression(cd8t1.cells, verbose = FALSE)$RNA))
avg.cd8t1.cells$gene <- rownames(avg.cd8t1.cells)
#write.csv(avg.stem.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/avgCD8tExp.csv")

cd8t2.cells <- subset(immune.combined, idents = "CD8_T2")
Idents(cd8t2.cells) <- "KO"
avg.cd8t2.cells <- as.data.frame(log1p(AverageExpression(cd8t2.cells, verbose = FALSE)$RNA))
avg.cd8t2.cells$gene <- rownames(avg.cd8t2.cells)

cd4t1.cells <- subset(immune.combined, idents = "CD4_T1")
Idents(cd4t1.cells) <- "KO"
avg.cd4t1.cells <- as.data.frame(log1p(AverageExpression(cd4t1.cells, verbose = FALSE)$RNA))
avg.cd4t1.cells$gene <- rownames(avg.cd4t1.cells)

cd4t2.cells <- subset(immune.combined, idents = "CD4_T2")
Idents(cd4t2.cells) <- "KO"
avg.cd4t2.cells <- as.data.frame(log1p(AverageExpression(cd4t2.cells, verbose = FALSE)$RNA))
avg.cd4t2.cells$gene <- rownames(avg.cd4t2.cells)

b.cells1 <- subset(immune.combined, idents = "B_cells1")
Idents(b.cells1) <- "KO"
avg.b.cells1 <- as.data.frame(log1p(AverageExpression(b.cells1, verbose = FALSE)$RNA))
avg.b.cells1$gene <- rownames(avg.b.cells1)

b.cells2 <- subset(immune.combined, idents = "B_cells2")
Idents(b.cells2) <- "KO"
avg.b.cells2 <- as.data.frame(log1p(AverageExpression(b.cells2, verbose = FALSE)$RNA))
avg.b.cells2$gene <- rownames(avg.b.cells2)

b.cells3 <- subset(immune.combined, idents = "B_cells3")
Idents(b.cells3) <- "KO"
avg.b.cells3 <- as.data.frame(log1p(AverageExpression(b.cells3, verbose = FALSE)$RNA))
avg.b.cells$gene <- rownames(avg.b.cells3)

b.cells4 <- subset(immune.combined, idents = "B_cells4")
Idents(b.cells4) <- "KO"
avg.b.cells4 <- as.data.frame(log1p(AverageExpression(b.cells4, verbose = FALSE)$RNA))
avg.b.cells4$gene <- rownames(avg.b.cells4)

b.cells5 <- subset(immune.combined, idents = "B_cells5")
Idents(b.cells5) <- "KO"
avg.b.cells5 <- as.data.frame(log1p(AverageExpression(b.cells5, verbose = FALSE)$RNA))
avg.b.cells5$gene <- rownames(avg.b.cells5)
#too few cells to be analyzed
cd4tspecial.cells <- subset(immune.combined, idents = "CD4_T*")
Idents(cd4tspecial.cells) <- "KO"
avg.cd4tspecial.cells <- as.data.frame(log1p(AverageExpression(cd4tspecial.cells, verbose = FALSE)$RNA))
avg.cd4tspecial.cells$gene <- rownames(avg.cd4tspecial.cells)

treg.cells <- subset(immune.combined, idents = "Treg")
Idents(treg.cells) <- "KO"
avg.treg.cells <- as.data.frame(log1p(AverageExpression(treg.cells, verbose = FALSE)$RNA))
avg.treg.cells$gene <- rownames(avg.treg.cells)

stem.cells <- subset(immune.combined, idents = "Stem")
Idents(stem.cells) <- "KO"
avg.stem.cells <- as.data.frame(log1p(AverageExpression(stem.cells, verbose = FALSE)$RNA))
avg.stem.cells$gene <- rownames(avg.stem.cells)

nk.cells <- subset(immune.combined, idents = "NK_cells")
Idents(nk.cells) <- "KO"
avg.nk.cells <- as.data.frame(log1p(AverageExpression(nk.cells, verbose = FALSE)$RNA))
avg.nk.cells$gene <- rownames(avg.nk.cells)

mast.cells <- subset(immune.combined, idents = "Mast")
Idents(mast.cells) <- "KO"
avg.mast.cells <- as.data.frame(log1p(AverageExpression(mast.cells, verbose = FALSE)$RNA))
avg.mast.cells$gene <- rownames(avg.mast.cells)

epithelial.cells <- subset(immune.combined, idents = "Epithelial_junk")
Idents(epithelial.cells) <- "KO"
avg.epithelial.cells <- as.data.frame(log1p(AverageExpression(epithelial.cells, verbose = FALSE)$RNA))
avg.epithelial.cells$gene <- rownames(avg.epithelial.cells)

ilc2.gdt.cells <- subset(immune.combined, idents = "ILC2/GD_T")
Idents(ilc2.gdt.cells) <- "KO"
avg.ilc2.gdt.cells <- as.data.frame(log1p(AverageExpression(ilc2.gdt.cells, verbose = FALSE)$RNA))
avg.ilc2.gdt.cells$gene <- rownames(avg.ilc2.gdt.cells)

dc.cells <- subset(immune.combined, idents = "DC")
Idents(dc.cells) <- "KO"
avg.dc.cells <- as.data.frame(log1p(AverageExpression(dc.cells, verbose = FALSE)$RNA))
avg.dc.cells$gene <- rownames(avg.dc.cells)

plots <- VlnPlot(immune.combined, features = c("Ccl5", "Nkg7", "Il2rb", "Ccl3", "Ccl4", "Cxcr6"), split.by = "KO", split.plot = TRUE,
                 pt.size = 1, combine = FALSE, ncol = 3)
CombinePlots(plots = plots, ncol = 1)

plot_paper<-VlnPlot(immune.combined, features=c("Ccl5", "Nkg7", "Il2rb", "Ccl3", "Ccl4", "Cxcr6"), pt.size = 0.1, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=3, legend = 'right')

#Genes with similar expression levels reflect conserved pathway in both cell types
#DEG with scatterplots
#AIGvsIL4RA: Macrophages and CD8_T cells
genes.to.label = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5","Il18", "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1", "Il13ra2")
p1 <- ggplot(avg.macro.cells, aes(IL4RA, AIG)) + geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Macrophages") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p1
p2 <- ggplot(avg.cd8t1.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD8_T1") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p2
p2 <- ggplot(avg.cd8t2.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD8_T2") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p2
p3 <- ggplot(avg.b.cells1, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("B_cells1") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p3
p3 <- ggplot(avg.b.cells2, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("B_cells2") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p3
#p3 <- ggplot(avg.b.cells3, aes(IL4RA, AIG)) + geom_point() + ggtitle("B_cells3") + geom_text(aes(label=ifelse(IL4RA>1,as.character(gene),'')),hjust=0,vjust=0)
#p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
#p3 genes not found
p3 <- ggplot(avg.b.cells4, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("B_cells4") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p3
#p3 <- ggplot(avg.b.cells5, aes(IL4RA, AIG)) + geom_point() + ggtitle("B_cells5") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
#p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
#p3 has too few cells
p4 <- ggplot(avg.cd4tspecial.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD4_T*") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p4 <- LabelPoints(plot = p4, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p4
p5 <- ggplot(avg.treg.cells, aes(IL4RA, AIG)) + geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Treg") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p5
p6 <- ggplot(avg.stem.cells, aes(IL4RA, AIG)) + geom_smooth(method = "lm", se = FALSE, color = "gray") + geom_point() + ggtitle("Stem") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p6
p7 <- ggplot(avg.nk.cells, aes(IL4RA, AIG)) + geom_smooth(method = "lm", se = FALSE, color = "gray") + geom_point() + ggtitle("NK_cells") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p7
p8 <- ggplot(avg.mast.cells, aes(IL4RA, AIG)) + geom_smooth(method = "lm", se = FALSE, color = "gray") + geom_point() + ggtitle("Mast") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p8 <- LabelPoints(plot = p8, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p8
p9 <- ggplot(avg.epithelial.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Epithelial_junk") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p9 <- LabelPoints(plot = p9, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p9
p10 <- ggplot(avg.ilc2.gdt.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("ILC2/GD_T") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p10 <- LabelPoints(plot = p10, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p10
p11 <- ggplot(avg.dc.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("DC") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p11 <- LabelPoints(plot = p11, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p11
p12 <- ggplot(avg.cd4t1.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD4_T1") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p12 <- LabelPoints(plot = p12, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p12
p12 <- ggplot(avg.cd4t2.cells, aes(IL4RA, AIG)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD4_T2") + geom_text(aes(label=ifelse(IL4RA>3,as.character(gene),'')),hjust=0,vjust=0)
p12 <- LabelPoints(plot = p12, points = genes.to.label, repel = TRUE, color = "red", xnudge = 0, ynudge = 0)
p12
#MacrophagesvsCD8Tcells
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <-"celltype"
avg.meta <- log1p(AverageExpression(immune.combined, verbose = FALSE)$RNA)
avg.meta<-as.data.frame(avg.meta)
avg.meta$gene <- rownames(avg.meta)
meta_comp<-ggplot(avg.meta, aes(Macrophages, CD8_T)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Macrophages vs CD8_T Average Gene Expression") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5))+ xlab("Macrophage") + ylab("CD8_T")
meta_comp<-LabelPoints(plot = meta_comp, points = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5", "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1"), color="red", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(meta_comp)

meta_comp<-ggplot(avg.meta, aes(CD4_T1, CD4_T2)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD4_T1 vs CD4_T2 Average Gene Expression") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5))+ xlab("CD4_T1") + ylab("CD4_T2")
meta_comp<-LabelPoints(plot = meta_comp, points = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5", "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1", "Il13ra2"), color="red", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(meta_comp)

meta_comp<-ggplot(avg.meta, aes(CD8_T1, CD8_T2)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("CD8_T1 vs CD8_T2 Average Gene Expression") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5))+ xlab("CD8_T1") + ylab("CD8_T2")
meta_comp<-LabelPoints(plot = meta_comp, points = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5", "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1", "Il13ra2"), color="red", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(meta_comp)

meta_comp<-ggplot(avg.meta, aes(B_cells1, B_cells2)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("B_cells1 vs B_cells2 Average Gene Expression") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5))+ xlab("B_cells1") + ylab("B_cells2")
meta_comp<-LabelPoints(plot = meta_comp, points = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5", "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1", "Il13ra2"), color="red", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(meta_comp)

meta_comp<-ggplot(avg.meta, aes(B_cells3, B_cells4)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("B_cells3 vs B_cells4 Average Gene Expression") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5))+ xlab("B_cells3") + ylab("B_cells4")
meta_comp<-LabelPoints(plot = meta_comp, points = c("Il1a","Il6","Il7","Il10","Ifng","Il17a","Il13","Il4","Il5", "Kit","Il3ra","Il1rl1","Il2ra","Il2rg","Il4ra","Il6ra","Il6st","Il13ra1", "Il13ra2"), color="red", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(meta_comp)

saveRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster3.rds")
immune.combined <- readRDS(file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster3.rds")



DefaultAssay(immune.combined) <- "RNA"
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

macrophagesIL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "Macrophages_IL4RA", ident.2 = "Macrophages_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(macrophagesIL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/MacrophagesResponseIL4RAvsAIG.csv")

Bcells1IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "B_cells1_IL4RA", ident.2 = "B_cells1_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(Bcells1IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Bcells1ResponseIL4RAvsAIG.csv")

Bcells2IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "B_cells2_IL4RA", ident.2 = "B_cells2_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(Bcells2IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Bcells2ResponseIL4RAvsAIG.csv")

Bcells3IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "B_cells3_IL4RA", ident.2 = "B_cells3_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(Bcells3IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Bcells3ResponseIL4RAvsAIG.csv")

Bcells4IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "B_cells4_IL4RA", ident.2 = "B_cells4_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(Bcells4IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Bcells4ResponseIL4RAvsAIG.csv")

Bcells5IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "B_cells5_IL4RA", ident.2 = "B_cells5_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(Bcells1IL5RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/Bcells5ResponseIL4RAvsAIG.csv")

cd4t1IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "CD4_T1_IL4RA", ident.2 = "CD4_T1_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(cd4t1IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cd4t1ResponseIL4RAvsAIG.csv")

cd4t2IL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "CD4_T2_IL4RA", ident.2 = "CD4_T2_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(cd4t2IL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cd4t2ResponseIL4RAvsAIG.csv")

tregIL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "Treg_IL4RA", ident.2 = "Treg_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(tregIL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/tregResponseIL4RAvsAIG.csv")

nkIL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "NK_cells_IL4RA", ident.2 = "NK_cells_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(nkIL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/nkResponseIL4RAvsAIG.csv")

dcIL4RAvsAIG.response <- FindMarkers(immune.combined, ident.1 = "DC_IL4RA", ident.2 = "DC_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(dcIL4RAvsAIG.response, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/dcResponseIL4RAvsAIG.csv")

cd8t1.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "CD8_T1_IL4RA", ident.2 = "CD8_T1_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(cd8t1.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cd8t1responseIL4RAvsAIG.csv")

cd8t2.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "CD8_T2_IL4RA", ident.2 = "CD8_T2_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(cd8t2.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cd8t2responseIL4RAvsAIG.csv")

stem.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "Stem_IL4RA", ident.2 = "Stem_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(stem.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/stemresponseIL4RAvsAIG.csv")

mast.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "Mast_IL4RA", ident.2 = "Mast_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(mast.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/mastresponseIL4RAvsAIG.csv")

cd4tspecial.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "CD4_T*_IL4RA", ident.2 = "CD4_T*_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(cd4tspecial.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/cd4tspecialresponseIL4RAvsAIG.csv")

epithelial.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "Epithelial_junk_IL4RA", ident.2 = "Epithelial_junk_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(epithelial.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/epithelialresponseIL4RAvsAIG.csv")

ilc2.gd.t.responseIL4RAvsAIG <- FindMarkers(immune.combined, ident.1 = "ILC2/GD_T_IL4RA", ident.2 = "ILC2/GD_T_AIG", verbose = FALSE, logfc.threshold = 1.1)
write.csv(ilc2.gd.t.responseIL4RAvsAIG, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/CSV/ilc2.gd.tresponseIL4RAvsAIG.csv")

saveRDS(immune.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/AIG_Immune_Int/AfterConservedCluster4.rds")

plots <- VlnPlot(immune.combined, features = c("Il4ra"), split.by = "celltype", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(immune.combined, features = c("Nkg7", "Ccl12", "Cxcl10", "Ccl4"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(immune.combined, features = c("Cxcl2", "Cxcl13", "Ccl2"), split.by = "KO", 
                 pt.size = 0, combine = FALSE, ncol = 1)
wrap_plots(plots = plots, ncol = 1)


