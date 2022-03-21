#SPEM ONLY CLUSTER
#Analysis on SPEM only clusters identifying ckuster 8, 11 as SPEM clusters from int data
DefaultAssay(diseased.combined) <- "integrated"
readRDS(diseased.combined, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/AfterConservedCluster6.rds")

levels(diseased.combined)
int.SPEM.cells <- subset(diseased.combined, idents = "SPEM")
int.SPEM.cells$celltype <- "SPEM"
plot1 <- FeatureScatter(int.SPEM.cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(int.SPEM.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
head(table(Idents(int.SPEM.cells))) #659 

int.STEM.cells <- subset(diseased.combined, idents = "Stem")
int.STEM.cells$celltype <- "STEM"
plot1 <- FeatureScatter(int.STEM.cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(int.STEM.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
head(table(Idents(int.STEM.cells))) #1073

int.SPEM.cells <- NormalizeData(int.SPEM.cells)
int.SPEM.cells <- FindVariableFeatures(int.SPEM.cells, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top30_int.SPEM.cells <- head(VariableFeatures(int.SPEM.cells), 30)
#write.csv(top30_int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/top30variantgenes.csv")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(int.SPEM.cells)
plot2 <- LabelPoints(plot = plot1, points = top30_int.SPEM.cells, repel = TRUE, max.overlaps = 40)
plot1 + plot2

all.genes <- rownames(int.SPEM.cells)
int.SPEM.cells <- ScaleData(int.SPEM.cells, features = all.genes)
int.SPEM.cells <- RunPCA(int.SPEM.cells, features = VariableFeatures(int.SPEM.cells))
int.SPEM.cells <- RunUMAP(int.SPEM.cells, dims = 1:25)
DimPlot(int.SPEM.cells, reduction = "pca")
ElbowPlot(int.SPEM.cells, ndims = 50)

int.SPEM.cells <- FindNeighbors(int.SPEM.cells, dims = 1:25)
int.SPEM.cells <- FindClusters(int.SPEM.cells, resolution = 0.6)
DimPlot(int.SPEM.cells, reduction = "umap", label = TRUE)
DimPlot(int.SPEM.cells, reduction = "umap", split.by = "KO")
DimPlot(int.SPEM.cells, reduction = "umap", group.by = "KO", label = TRUE)

saveRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/AllClusters1.rds")
readRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/AllClusters1.rds")

DefaultAssay(int.SPEM.cells) <- "RNA"
meta.markers <- FindAllMarkers(int.SPEM.cells, only.pos = TRUE, logfc.threshold = 0.8)
#write.csv(meta.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/Markers.csv")

#Visualization
VlnPlot(int.SPEM.cells, features = c("Gif","Tff2","Muc5ac","Stmn1","Muc6","Gkn3","Mki67","Birc5"), ncol = 4, pt.size = 0.1)
VlnPlot(int.SPEM.cells, features = c("Wfdc2","Gkn3", "Lgr5", "Klk1", "Wfdc18"), ncol = 1, pt.size = 0.1)
VlnPlot(int.SPEM.cells, features = c("Aqp5", "Agr2", "Birc5", "Mki67", "Bhlha15", "Gif"), ncol = 1, pt.size = 0.1)

plot_paper<-VlnPlot(int.SPEM.cells, features=c("Tff2", "Gkn3", "Gif", "Muc6", "Rps2", "Stmn1", "Birc5", "Mki67"), pt.size = 0.1, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')
plot_paper<-VlnPlot(int.SPEM.cells, features=c("Gkn3", "Tff2", "Stmn1", "Birc5"), pt.size = 0.1, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')

#VlnPlot(int.SPEM.cells, features = c("Tff3", "Cd44", "Tacstd2", "Muc2", "Reg3g"), ncol = 1, pt.size = 0)
#VlnPlot(int.SPEM.cells, features = c("Barx1", "Sox2", "Sox4", "Bmi1", "Trp63"), ncol = 1, pt.size = 0)
#The following requested variables were not found: Defa5, Popdc1. All cells have the same value of rna_Popdc3.

Scluster0.markers <- FindMarkers(int.SPEM.cells, ident.1 = 0, logfc.threshold = 1.1, test.use = "roc")
Scluster1.markers <- FindMarkers(int.SPEM.cells, ident.1 = 1, logfc.threshold = 1.1, test.use = "roc")
Scluster2.markers <- FindMarkers(int.SPEM.cells, ident.1 = 2, logfc.threshold = 1.1, test.use = "roc")
Scluster3.markers <- FindMarkers(int.SPEM.cells, ident.1 = "Proliferative-SPEM", logfc.threshold = 1.1)

write.csv(Scluster0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/Cluster0.csv")
write.csv(Scluster1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/Cluster1.csv")
write.csv(Scluster2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/Cluster2.csv")
write.csv(Scluster3.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21_22/AIG_HDT_Int/SPEM/CSV/Proliferative.csv")
#groupedbyKO
SConservedC0.markers <- FindConservedMarkers(int.SPEM.cells, ident.1 = 0, grouping.var = "KO", verbose = FALSE)
SConservedC1.markers <- FindConservedMarkers(int.SPEM.cells, ident.1 = 1, grouping.var = "KO", verbose = FALSE)
SConservedC2.markers <- FindConservedMarkers(int.SPEM.cells, ident.1 = 2, grouping.var = "KO", verbose = FALSE)
SConservedC3.markers <- FindConservedMarkers(int.SPEM.cells, ident.1 = 3, grouping.var = "KO", verbose = FALSE)

write.csv(SConservedC0.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/ConserverdC0.csv")
write.csv(SConservedC1.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/ConserverdC1.csv")
write.csv(SConservedC2.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/ConserverdC2.csv")
write.csv(SConservedC3.markers, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/ConserverdC3.csv")

saveRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/Final.rds")
readRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/Final.rds")

p1 <- DimPlot(int.SPEM.cells, reduction = "umap", group.by = "KO", label = TRUE)
p2 <- DimPlot(int.SPEM.cells, reduction = "umap", label = TRUE, repel = TRUE)
p1
p1 + p2
DimPlot(int.SPEM.cells, label = TRUE)
DimPlot(int.SPEM.cells, reduction = "umap", split.by = "KO")

int.SPEM.cells <- RenameIdents(int.SPEM.cells, `0` = "Translationally activeSPEM", `1` = "Mucous SPEM",  `2` = "Zymogenic SPEM", `3` = "Proliferative SPEM")

DimPlot(int.SPEM.cells, reduction = "umap", label = TRUE)
DimPlot(int.SPEM.cells, reduction = "umap", split.by = "KO")
DimPlot(int.SPEM.cells, reduction = "umap", group.by = "KO", label = TRUE)

VlnPlot(int.SPEM.cells, features = c("Gif","Tff2","Muc5ac","Stmn1","Muc6","Gkn3","Mki67","Birc5"), ncol = 4, pt.size = 0.1)
plot_paper<-VlnPlot(int.SPEM.cells, features=c("Tff2", "Gkn3", "Muc6", "Gif", "Stmn1", "Birc5", "Mki67", "Ube2c"), pt.size = 0.1, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')
plot_paper<-VlnPlot(int.SPEM.cells, features=c("Gkn3", "Tff2", "Stmn1", "Birc5"), pt.size = 0.1, combine=FALSE, split.plot = TRUE, cols = c("yellow","blue"), split.by="KO")
CombinePlots(plots=plot_paper, ncol=4, legend = 'right')

DefaultAssay(int.SPEM.cells) <- "integrated"

table(Idents(int.SPEM.cells))
#Neck-SPEM Super-SPEM  Stem-SPEM   Pit-SPEM   UNK-SPEM 
#177        139        125        122         96 
prop.table(table(Idents(int.SPEM.cells)))
#Neck-SPEM Super-SPEM  Stem-SPEM   Pit-SPEM   UNK-SPEM 
#0.2685888  0.2109256  0.1896813  0.1851290  0.1456753 

FeaturePlot(int.SPEM.cells, features = c("Gkn3", "Tff2", "Stmn1"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(int.SPEM.cells, features = c("Mki67", "Sox4", "Lgr5", "Agr2"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)
FeaturePlot(int.SPEM.cells, features = c("Birc5", "Dmbt1", "Bhlha15"), split.by = "KO", max.cutoff = 3, 
            cols = c("grey", "red"), label = TRUE)

VlnPlot(int.SPEM.cells, features = c("Gif", "Muc6", "Gkn3", "Tff2", "Stmn1", "Muc5ac"), ncol = 1, pt.size = 0)
VlnPlot(int.SPEM.cells, features = c("Wfdc2","Gkn3", "Lgr5", "Klk1", "Wfdc18"), ncol = 1, pt.size = 0)
VlnPlot(int.SPEM.cells, features = c("Aqp5", "Agr2", "Birc5", "Mki67", "Bhlha15", "Gif"), ncol = 1, pt.size = 0)
VlnPlot(int.SPEM.cells, features = c("Tff3", "Cd44", "Tacstd2", "Muc2"), ncol = 1, pt.size = 0)
VlnPlot(int.SPEM.cells, features = c("Barx1", "Sox2", "Sox4", "Bmi1", "Trp63"), ncol = 1, pt.size = 0)

plots <- VlnPlot(int.SPEM.cells, features = c("Gif", "Muc6", "Gkn3", "Tff2", "Stmn1", "Muc5ac"), ncol = 1, pt.size = 0.1, split.by = "KO", combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(int.SPEM.cells, features = c("Iqgap3","Stmn1", "Birc5"), ncol = 1, pt.size = 0.1, split.by = "KO", combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(int.SPEM.cells, features = c("Aqp5", "Agr2", "Birc5", "Mki67", "Bhlha15", "Gif"), ncol = 1, pt.size = 0, split.by = "KO", combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(int.SPEM.cells, features = c("Tff3", "Cd44", "Barx1", "Tacstd2", "Sox4", "Trp63"), ncol = 1, pt.size = 0, split.by = "KO", combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(int.SPEM.cells, features = c("Lgr5", "Bmi1", "Dmbt1", "Sox2", "Tnfrsf19"), ncol = 1, pt.size = 0, split.by = "KO", combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

DimPlot(int.SPEM.cells, label = TRUE)

DefaultAssay(int.SPEM.cells) <- "integrated"
library(dplyr)
meta.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(subset(int.SPEM.cells, downsample=100), features = top20$gene, angle = 90, size=4, lines.width = 1) + scale_fill_gradientn(colors = c("blue","white","red"))+ theme(axis.text.y = element_text(size = 10))

saveRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/Final2.rds")

#TrajectoryAnalysis
library(patchwork)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
DefaultAssay(int.SPEM.cells) <- "RNA"
readRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/Final2.rds")
DimPlot(int.SPEM.cells, reduction = "umap", label = TRUE)
levels(int.SPEM.cells)

#int.SPEM.cells$celltype <- Idents(int.SPEM.cells)
#int.SPEM <- int.SPEM.cells[, int.SPEM.cells$celltype %in% c("Low-SPEM", "Pit-SPEM", "Chief-SPEM", "Super-SPEM", "Stem-SPEM", "Neck-SPEM",
#                                                           "Neck-SPEM2", "Low-SPEM2")]

DefaultAssay(int.SPEM.cells) <- "RNA"
int.SPEM.cells$celltype.KO <- paste(Idents(int.SPEM.cells), int.SPEM.cells$KO, sep = "_")
int.SPEM.cells$celltype <- Idents(int.SPEM.cells)
Idents(int.SPEM.cells) <- "celltype.KO"

DimPlot(int.SPEM.cells, reduction = "umap", split.by = "KO")
int.SPEM.cells$celltype <- Idents(int.SPEM.cells)
SPEM.AIG <- int.SPEM.cells[, int.SPEM.cells$celltype.KO %in% c("Low-SPEM_AIG", "Pit-SPEM_AIG", "Chief-SPEM_AIG", "Super-SPEM_AIG", 
                                                               "Stem-SPEM_AIG", "Neck-SPEM_AIG", "Neck-SPEM2_AIG", "Low-SPEM2_AIG")]
SPEM.HDT <- int.SPEM.cells[,  int.SPEM.cells$celltype.KO %in% c("Low-SPEM_HDT", "Pit-SPEM_HDT", "Chief-SPEM_HDT", "Super-SPEM_HDT", 
                                                                "Stem-SPEM_HDT", "Neck-SPEM_HDT","Neck-SPEM2_HDT", "Low-SPEM2_HDT")]
int.SPEM <- int.SPEM.cells[, int.SPEM.cells$celltype.KO %in% c("Low-SPEM_AIG", "Pit-SPEM_AIG", "Chief-SPEM_AIG", "Super-SPEM_AIG", 
                                                               "Stem-SPEM_AIG", "Neck-SPEM_AIG", "Neck-SPEM2_AIG", "Low-SPEM2_AIG", 
                                                               "Low-SPEM_HDT", "Pit-SPEM_HDT", "Chief-SPEM_HDT", "Super-SPEM_HDT", 
                                                               "Stem-SPEM_HDT", "Neck-SPEM_HDT","Neck-SPEM2_HDT", "Low-SPEM2_HDT")]
#AIG
SPEM.AIG.cds <- as.cell_data_set(SPEM.AIG)
## Calculate size factors using built-in function in monocle3
SPEM.AIG.cds <- estimate_size_factors(SPEM.AIG.cds)
## Add gene names into CDS
SPEM.AIG.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(SPEM.AIG[["RNA"]])

SPEM.AIG.cds <- preprocess_cds(SPEM.AIG.cds, preprocess_method = "PCA")
SPEM.AIG.cds <- reduce_dimension(SPEM.AIG.cds, reduction_method = "UMAP") #also preprocess_cds with PCA by default
SPEM.AIG.cds <- cluster_cells(cd = SPEM.AIG.cds, reduction_method = "UMAP")
SPEM.AIG.cds <- learn_graph(SPEM.AIG.cds, use_partition = TRUE)
colData(SPEM.AIG.cds)
write.csv(colData(SPEM.AIG.cds), file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/PseudotimesummaryAIG.csv")

#HDT
SPEM.HDT.cds <- as.cell_data_set(SPEM.HDT)
## Calculate size factors using built-in function in monocle3
SPEM.HDT.cds <- estimate_size_factors(SPEM.HDT.cds)
## Add gene names into CDS
SPEM.HDT.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(SPEM.HDT[["RNA"]])

SPEM.HDT.cds <- preprocess_cds(SPEM.HDT.cds, preprocess_method = "PCA")
SPEM.HDT.cds <- reduce_dimension(SPEM.HDT.cds, reduction_method = "UMAP") #also preprocess_cds with PCA by default
SPEM.HDT.cds <- cluster_cells(cd = SPEM.HDT.cds, reduction_method = "UMAP")
SPEM.HDT.cds <- learn_graph(SPEM.HDT.cds, use_partition = TRUE)
colData(SPEM.HDT.cds)
write.csv(colData(SPEM.HDT.cds), file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/PseudotimesummaryHDT.csv")

#BOTH
int.SPEM.cds <- as.cell_data_set(int.SPEM)
## Calculate size factors using built-in function in monocle3
int.SPEM.cds <- estimate_size_factors(int.SPEM.cds)
## Add gene names into CDS
int.SPEM.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(int.SPEM[["RNA"]])

int.SPEM.cds <- preprocess_cds(int.SPEM.cds, preprocess_method = "PCA")
int.SPEM.cds <- reduce_dimension(int.SPEM.cds, reduction_method = "UMAP") #also preprocess_cds with PCA by default
int.SPEM.cds <- cluster_cells(cd = int.SPEM.cds, reduction_method = "UMAP")
int.SPEM.cds <- learn_graph(int.SPEM.cds, use_partition = TRUE)
colData(int.SPEM.cds)
write.csv(colData(int.SPEM.cds), file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/CSV/PseudotimesummaryBOTH.csv")

plot_cells(
  cds = SPEM.AIG.cds,
  color_cells_by = "celltype", label_cell_groups = FALSE,
  label_leaves = TRUE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE)

plot_cells(
  cds = SPEM.HDT.cds,
  color_cells_by = "celltype", label_cell_groups = FALSE,
  label_leaves = TRUE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE)

plot_cells(
  cds = int.SPEM.cds,
  color_cells_by = "celltype", label_cell_groups = FALSE,
  label_leaves = TRUE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE)

SPEM.AIG.cds <- order_cells(cds = SPEM.AIG.cds, reduction = "UMAP")
SPEM.HDT.cds <- order_cells(cds = SPEM.HDT.cds, reduction = "UMAP")
int.SPEM.cds <- order_cells(cds = int.SPEM.cds, reduction = "UMAP")

plot_cells(
  cds = SPEM.AIG.cds,
  color_cells_by = "celltype", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = SPEM.AIG.cds,
  color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = SPEM.AIG.cds,
  color_cells_by = "pseudotime",
  genes = c("Gkn3", "Tff2", "Chil4", "Spp1", "Stmn1", "Birc5", "Mki67", "Dmbt1", "Bhlha15"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

#HFD Genes
plot_cells(
  cds = SPEM.AIG.cds,
  color_cells_by = "pseudotime",
  genes = c("Chil4", "Bpifb1", "Myh1", "Pdia2", "Gif", "Tnnc2", "Ltf"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

plot_cells(
  cds = SPEM.HDT.cds,
  color_cells_by = "celltype", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = SPEM.HDT.cds,
  color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = SPEM.HDT.cds,
  color_cells_by = "pseudotime",
  genes = c("Gkn3", "Tff2", "Chil4", "Spp1", "Stmn1", "Birc5", "Mki67", "Dmbt1", "Bhlha15"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

#HFD Genes
plot_cells(
  cds = SPEM.HDT.cds,
  color_cells_by = "pseudotime",
  genes = c("Chil4", "Bpifb1", "Myh1", "Pdia2", "Gif", "Tnnc2", "Ltf"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

plot_cells(
  cds = int.SPEM.cds,
  color_cells_by = "celltype", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = int.SPEM.cds,
  color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  group_label_size = 5, cell_size = 1,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = int.SPEM.cds,
  color_cells_by = "pseudotime",
  genes = c("Gkn3", "Tff2", "Chil4", "Spp1", "Stmn1", "Birc5", "Mki67", "Dmbt1", "Bhlha15"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

#HFD Genes
plot_cells(
  cds = int.SPEM.cds,
  color_cells_by = "pseudotime",
  genes = c("Chil4", "Bpifb1", "Myh1", "Pdia2", "Gif", "Tnnc2", "Ltf"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE)

#finds the top20 markers in each pseudotime cluster
marker_test_res <- top_markers(SPEM.AIG.cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(20, pseudo_R2)

head(top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id)))
plot_genes_by_group(SPEM.AIG.cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

SPEM.AIG.cds <- AddMetaData(
  object = SPEM.AIG.cds,
  metadata = SPEM.AIG.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "SPEM.AIG"
)

marker_test_res <- top_markers(SPEM.HDT.cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(20, pseudo_R2)

SPEM.HDT.cds <- AddMetaData(
  object = SPEM.HDT.cds,
  metadata = SPEM.HDT.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "SPEM.HDT.cds"
)

head(top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id)))
plot_genes_by_group(SPEM.HDT.cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

marker_test_res <- top_markers(int.SPEM.cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(20, pseudo_R2)

head(top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id)))
plot_genes_by_group(int.SPEM.cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

int.SPEM.cds <- AddMetaData(
  object = int.SPEM.cds,
  metadata = int.SPEM.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "int.SPEM.cds"
)

FeaturePlot(int.SPEM.cells, c("SPEM.AIG.cds", "SPEM.HDT.cds", "int.SPEM.cds"), label = TRUE, pt.size = 0.1) & scale_color_viridis_c()
saveRDS(int.SPEM.cells, file = "/Users/dipaololabmacbook/Desktop/scRNA21/DataIntegration3a/SPEM/afterPseudotime.rds")
FeaturePlot(int.SPEM.cells, c("SPEM.AIG.cds", "SPEM.HDT.cds", "int.SPEM.cds"), label = TRUE, pt.size = 0.1)

