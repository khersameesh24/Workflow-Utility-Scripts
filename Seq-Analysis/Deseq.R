#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")

#loading DESeq
#.libPaths("/opt/R-3.4.2/library")
library("DESeq")
library(pheatmap)
library(gplots)


#read the counts file.
countsTable <- read.table("WT_LT_vs_KO_LT_readcounts.filtered.tsv",sep="\t",header=T)
countsTable$Gene_ID = paste(countsTable$Gene_ID,1:nrow(countsTable),sep ="_")

##specify that the first column is the gene name.
rownames(countsTable) <- countsTable$Gene_ID
countsTable <- countsTable[,-1]
head(countsTable)

##specify the conditions
conditions <- factor(c(rep("Normal",3),rep("Knockout",3)))

##create the main DESeq object
countDataSet <- newCountDataSet(countsTable,conditions)

##adjust for the difference in reads mapped to in each sample
countDataSet <- estimateSizeFactors(countDataSet)
sizeFactors(countDataSet)
head(counts(countDataSet))

##to plot heatmap
normalized = counts(countDataSet,normalized=TRUE)
norm <- cbind("id"=rownames(normalized),normalized)

countDataSet <-estimateDispersions(countDataSet)

##plot dispersion
pdf("disp_DESeq.pdf")
plotDispEsts(countDataSet)
dev.off()

##all differential expression values
DEVal=nbinomTest(countDataSet,"Normal","Knockout")

##plot the log2 fold changes (cancer v normal) against the mean normalised counts
pdf("plotMA_DESeq.pdf")
plotMA(DEVal,col = ifelse(DEVal$pval<=0.05, "gray32", "red"), linecol = "red")
dev.off()

##plot histogram of p values
pdf("hist_pVal_DESeq.pdf")
hist(DEVal$pval, breaks=100, col="skyblue", border="slateblue")
dev.off()

##plot Volcano
pdf("VolcanoPlot.pdf")
plot(DEVal$log2FoldChange, log10(DEVal$pval), cex = 0.5, col="blue")
dev.off()

##saving data
write.csv(DEVal,file="DE_LT.filtered.csv")

#-------------------------------------------------------------------------------------------------------------
## Filter from All_DE_miRNA_normal_vs_glio_ouput.csv based on Pval <= 0.05 and log2foldchange >= 2 and <= -2 using awk:
system(awk '($6>2 || $6<-2) && $7<0.05' All_DE_RNA_normal_vs_glio_output.csv >  Significantly_DE_RNA_normal_vs_glio_output.csv)

significant <- read.table("sig_DE_LT.filtered.tsv",sep="\t",header=F)
colnames(significant) <- c("id","baseMean","baseMeanA","baseMeanB","foldChange","log2FoldChange","pval","padj")
#-------------------------------------------------------------------------------------------------------------
  
geneFilt=merge(norm,significant,by="id")
gene_DE <- geneFilt[,-(8:14)]

rownames(gene_DE) <- gene_DE$id
gene_DE <- gene_DE[,-1]
colnames(gene_DE) <- c("WT1","WT3","WT4","KO2","KO3","KO4")


gene <- as.matrix(gene_DE)
gene <- apply(gene_DE,2,as.numeric)

dist <- dist(t(gene))
plot(hclust(dist(gene[1:10])))
##plot heatmap
pdf("HeatMap.pdf")
heatmap(t(gene))
dev.off()


