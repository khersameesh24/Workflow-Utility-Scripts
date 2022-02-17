library(gplots)
library(pheatmap)
gc()
set.seed(42)
imagefile = "heatmap_WT-LT_vs_KO-LT.pdf"
dat1<- read.table("heatmap.csv",sep=",",header=TRUE)
dat1_filtered <- as.matrix(data.frame(dat1$WT_LT,dat1$KO_LT))
colnames(dat1_filtered) <- c("WT_LT","KO_LT")
rownames(dat1_filtered) <- dat1$Geneids	
pdf(imagefile,width=5, height=10)
pheatmap(dat1_filtered ,scale = "column", cluster_row=FALSE , cellwidth =40, cellheight = 10, color = colorRampPalette(c("green","red"))(50),margins=c(3,25),fontsize = 5)
dev.off()

