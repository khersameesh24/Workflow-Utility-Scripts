#include library seurat and dplyr
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
#load datasets
hspc.data.control <- Read10X(data.dir = "/home/RNA_Seq/Sc_Seq/scRNA_seq_data/seq_data/Control/10x_data/outs/filtered_feature_bc_matrix/")
hspc.control <- CreateSeuratObject(counts = hspc.data.control , project = "Control", min.features = 200)
hspc.control

hspc.data.knockout <- Read10X(data.dir = "/home/RNA_Seq/Sc_Seq/scRNA_seq_data/seq_data/Knockout/10x_data/outs/filtered_feature_bc_matrix/")
hspc.knockout <- CreateSeuratObject(counts = hspc.data.knockout , project = "Knockout", min.features = 200)
hspc.knockout

#The percentage of reads that map to the mitochondrial genome(all genes starting with mt-)
hspc.control[["percent.mt"]] <- PercentageFeatureSet(hspc.control, pattern = "^mt-")
hspc.knockout[["percent.mt"]] <- PercentageFeatureSet(hspc.knockout, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(hspc.control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
boxplot.control <-boxplot(hspc.control$percent.mt, xlab="Mitochondrial genes", ylab="Percentage")
feature.control <- boxplot(hspc.control$nFeature_RNA, xlab="#Features", ylab="Percentage")

VlnPlot(hspc.knockout, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
boxplot.knockout <-boxplot(hspc.knockout$percent.mt, xlab="Mitochondrial genes", ylab="Percentage")
feature.knockout <- boxplot(hspc.knockout$nFeature_RNA, xlab="#Features", ylab="Percentage")

# FeatureScatter Plot
plot1 <- FeatureScatter(hspc.control, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hspc.control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(hspc.knockout, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hspc.knockout, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Removing Outliers
hspc.control <- subset(hspc.control, subset = percent.mt > 0.5 & percent.mt < 10.58)
hspc.control

hspc.knockout <- subset(hspc.knockout, subset = percent.mt > 0.581 & percent.mt < 10.84)
hspc.knockout

#Normalizing the data
hspc.control <- NormalizeData(hspc.control)
hspc.knockout <- NormalizeData(hspc.knockout)

#Identification of highly variable features 
hspc.control <- FindVariableFeatures(hspc.control)
hspc.knockout <- FindVariableFeatures(hspc.knockout)

# Identify the 30 most highly variable genes
top30.control <- head(VariableFeatures(hspc.control), 30)
top30.knockout <- head(VariableFeatures(hspc.knockout), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hspc.control)
plot2 <- LabelPoints(plot = plot1, points = top30.control[1:10], repel = TRUE)

plot1 <- VariableFeaturePlot(hspc.knockout)
plot2 <- LabelPoints(plot = plot1, points = top30.knockout[1:10], repel = TRUE)
CombinePlots(plots = list(plot1, plot2))  #Dimension issue

#Scaling the data (linear transformation) 
all.genes.control <- rownames(hspc.control)
hspc.control <- ScaleData(hspc.control, features = all.genes.control) 

all.genes.knockout <- rownames(hspc.knockout)
hspc.knockout <- ScaleData(hspc.knockout, features = all.genes.knockout) 

#Perform linear dimensional reduction 
hspc.control <- RunPCA(hspc.control, features = VariableFeatures(object = hspc.control)) 
hspc.knockout <- RunPCA(hspc.knockout, features = VariableFeatures(object = hspc.knockout))

#print(hspc.control[["pca"]], dims = 1:5, nfeatures = 5)
#print(hspc.knockout[["pca"]], dims = 1:5, nfeatures = 5)

# # PCA related plots
VizDimLoadings(hspc.control, dims = 1:2, reduction = "pca")
DimPlot(hspc.control, reduction = "pca")
DimHeatmap(hspc.control, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(hspc.control, dims = 1:15, cells = 500, balanced = TRUE)

VizDimLoadings(hspc.knockout, dims = 1:2, reduction = "pca")
DimPlot(hspc.knockout, reduction = "pca")
DimHeatmap(hspc.knockout, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(hspc.knockout, dims = 1:15, cells = 500, balanced = TRUE)

# #Determine the ‘dimensionality’ of the dataset
hspc.control <- JackStraw(hspc.control, num.replicate = 100)
hspc.control <- ScoreJackStraw(hspc.control, dims = 1:20)

hspc.knockout <- JackStraw(hspc.knockout, num.replicate = 100)
hspc.knockout <- ScoreJackStraw(hspc.knockout, dims = 1:20)

#JackStrawPlot - visualization tool for comparing the distribution of p-values 
JackStrawPlot(hspc.control, dims = 1:15)
ElbowPlot(hspc.control)

JackStrawPlot(hspc.knockout, dims = 1:15)
ElbowPlot(hspc.knockout)

#Cluster Cells -- graph-based clustering approach
hspc.control <- FindNeighbors(hspc.control, dims = 1:15)
hspc.control <- FindClusters(hspc.control, resolution = 0.5)

hspc.knockout <- FindNeighbors(hspc.knockout, dims = 1:15)
hspc.knockout <- FindClusters(hspc.knockout, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP/tSNE)
hspc.control <- RunUMAP(hspc.control, dims = 1:15)
hspc.knockout <- RunUMAP(hspc.knockout, dims = 1:15)
hspc.control <- RunTSNE(hspc.control, dims = 1:15)
hspc.knockout <- RunTSNE(hspc.knockout, dims = 1:15)

#----------------------###CombineObjects####-------------------------------------------#

hspc.combined <- merge(x = hspc.control , y = hspc.knockout)
hspc.combined <- SplitObject(hspc.combined, split.by = "orig.ident")
hspc.combined <- lapply(X = hspc.combined, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
})
hspc.anchors <- FindIntegrationAnchors(object.list = hspc.combined,dims = 1:15)
hspc.combined <- IntegrateData(anchorset = hspc.anchors, dims = 1:15)

DefaultAssay(hspc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
hspc.combined <- ScaleData(hspc.combined, verbose = FALSE)
hspc.combined <- RunPCA(hspc.combined, npcs = 30, verbose = FALSE)
# UMAP and Clustering
hspc.combined <- RunUMAP(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- RunTSNE(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- FindNeighbors(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- FindClusters(hspc.combined, resolution = 0.5)

hspc.combined <- RunUMAP(hspc.combined, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)

# Plot the clusters
DimPlot(data, group.by = "RNA_snn_res.1")


# Visualization UMAP/TSNE
p1 <- DimPlot(hspc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(hspc.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
p1 <- DimPlot(hspc.combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(hspc.combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

DimPlot(hspc.combined, reduction = "umap", split.by = "orig.ident")

#To identify canonical cell type marker genes that are conserved across conditions
DefaultAssay(hspc.combined) <- "RNA"
conserved.markers <- FindConservedMarkers(hspc.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)

#explore marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(hspc.combined, features = c("Cd34",
"Procr",
"Slamf1",
"Cd48",
"Alcam",
"Meis1",
"Erg",
"Trpc6",
"Mpl",
"Hoxb5",
"Fgd5",
"Ctnnal1",
"Serpina3f"), min.cutoff = "q9",split.by = "orig.ident")

#Name the clusters according to the cell types
hspc.combined <- RenameIdents(hspc.combined,'0' = "GMP−1",
'1' = "pMo−1",
'2' = "pNeu",
'3' = "pMo−2",
'4' = "GMP−2",
'5' = "GMP−3",
'6' = "MEP",
'7' ="Ery",
'8' = "Myeloid−Biased HSC",
'9' = "pBa",
'10' = "pDC",
'11' = "Lymphoid−Biased HSC",
'12' = "pT",
'13' = "pB−1",
'14' = "pB−2",
'15' = "pB−3",
'16' = "Platelet−Biased HSC",
'17' = "GMP−4")
DimPlot(hspc.combined, label = TRUE)


#Identify differential expressed genes across conditions
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
c1 <- subset(hspc.combined, idents = "C_1")
Idents(c1) <- "knockout"
avg.c1 <- log1p(AverageExpression(c1, verbose = FALSE)$RNA)
avg.c1$gene <- rownames(avg.c1)


c2 <- subset(hspc.combined, idents = "C_2")
Idents(c2) <- "knockout"
avg.c2 <- log1p(AverageExpression(c2, verbose = FALSE)$RNA)
avg.c2$gene <- rownames(avg.c2)

genes.to.label = c("Itgb3", "Trp53", "Klf5", "Cops5", "Gata6")
p1 <- ggplot(avg.c1, aes(CTRL, KO)) + geom_point() + ggtitle("C_1")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.c2, aes(CTRL, KO)) + geom_point() + ggtitle("C_2")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)

hspc.combined$orig.ident <- paste(Idents(hspc.combined), hspc.combined$orig.ident, sep = "_")
hspc.combined$orig.ident <- Idents(hspc.combined)
Idents(hspc.combined) <- "orig.ident"
differential gene expression between cluster 0 of control and Knockout
diff.exp.genes <- FindMarkers(hspc.combined, ident.1 = "Knockout_0", ident.2 = "Control_0", verbose = FALSE)
head(diff.exp.genes, n = 15)



