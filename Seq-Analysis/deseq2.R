# Analyze RNA-seq data with DeSeq2
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfil
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)

setwd("~/jnc/Heatmaps/results/DESEQ2_WT-LT_vs_KO-LT")

# read the counts table
countsTable <- read.table("WT_LT_vs_KO_LT_readcounts.tsv", sep = "\t", header = T)
countsTable$Gene_ID <- paste(countsTable$Gene_ID, 1:nrow(countsTable), sep = "_")
rownames(countsTable) <- countsTable$Gene_ID
countsTable <- countsTable[, -1]
head(countsTable)

# specify the conditions
coldata <- read.table("conditions.tsv")
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# convert the counts table to counts-matrix
cts <- as.matrix(countsTable)
head(cts)

# construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)

# specify explicitly the comparisons to be made
dds$condition <- factor(dds$condition, levels = c("Knockout", "Wildtype"))

# differential expression analysis
dds <- DESeq(dds)

# extracts results table with log2 fold changes, p-values and adjusted-p-values
res <- results(dds)
# coefficient/contrast could be specified to build a result table for
# res <- results(dds, contrast=c("condition","Knockout","Wildtype"))

# Information about which variables and tests were used can be found
# by calling the function mcols on the results object.
mcols(res)$description

# Summarize basic tallies using the summary fn()
summary(res)

# shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
resLFC <- lfcShrink(dds, coef = "condition_Wildtype_vs_Knockout", type = "apeglm")

# order the results table by the smallest p value
resOrdered <- res[order(res$pvalue),]

# how many adjusted p-values are less than 0.1
sum(res$padj < 0.1, na.rm = TRUE)

# plotMA shows the log2 fold changes attributable to a given variable 
# over the mean of normalized counts for all the samples in the DESeqDataSet
plotMA(res, ylim = c(-2, 2))
plotMA(resLFC, ylim = c(-2, 2))

# write data to a csv file
write.csv(as.data.frame(resOrdered), file = "WT_vs_KO_LT_results.csv")

# filter out entries based on a cutoff
system("sed 's/,/\t/g' WT_vs_KO_LT_results.csv | sed 's/\"//g' | awk '($3>2||$3<-2) && $7<0.05' > significants_WT_vs_KO_LT_results.tsv")

# read from the filtered results
resSig <- read.table("significants_WT_vs_KO_LT_results.tsv", sep = "\t", header = F)
# specify column names
colnames(resSig) <- c("id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# assay function is used to extract the matrix of normalized values
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

norm <- as.data.frame(assay(vsd))
write.table(norm, "normalized_WT_vs_KO_LT.tsv", sep = "\t")
norm$id <- rownames(norm)

geneFilt <- merge(norm, resSig, by = "id")
gene_DE <- geneFilt[, - (8:13)]

rownames(gene_DE) <- gene_DE$id
gene_DE <- gene_DE[, -1]
colnames(gene_DE) <- c("WT-1", "WT-2", "WT-3", "ArjKO-1", "ArjKO-2", "ArjKO-3")

gene <- as.matrix(gene_DE)
gene <- apply(gene_DE, 2, as.numeric)
rownames(gene) <- rownames(gene_DE)

annotation <- data.frame(label = c(rep("Wildtype", 3), rep("Knockout", 3)))
rownames(annotation) <- colnames(as.character(gene))

# specify the color pallete
colors <- colorRampPalette(brewer.pal(9, "RdYlBu"))(25)

# scale the matrix
gene_scaled = t(scale(t(gene)))
# save heatmap as a jpg
jpeg(filename = "WT_vs_ArjKO_LT_heatmap_most_significant_genes.jpg",
     width = 464,
     height = 360,
     pointsize = 12,
     quality = 100,
     bg = "white",
     res = NA)
# generate heatmap
heatmap(gene_scaled,
        labRow = F,
        col = colors,
        cexCol = 0.9,
        keep.dendro = T)
dev.off()

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
# heatmap of this distance matrix gives us an overview over
# similarities and dissimilarities between samples
jpeg(filename = "WT_vs_ArjKO_LT_sample_to_sample_distances.jpg",
     width = 464,
     height = 360,
     pointsize = 12,
     quality = 100,
     bg = "white",
     res = NA)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
# Principal component plot of the samples
# shows the samples in the 2D plane spanned by their first two principal components.
# useful for visualizing the overall effect of experimental covariates and batch effects
jpeg(filename = "WT_vs_ArjKO_LT_PCA_plot.jpg",
     width = 464,
     height = 360,
     pointsize = 12,
     quality = 100,
     bg = "white",
     res = NA)
# plotPCA(vsd, intgroup = c("condition", "type"))
# customize the PCA plot using the ggplot function
pcaData <- plotPCA(vsd, intgroup = c("condition", "type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()
