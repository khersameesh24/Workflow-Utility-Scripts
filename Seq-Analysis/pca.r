data = read.csv (file, sep = "" , header = "")
pca_res = prcomp(data, center=T ,scale =T)
plot(pca_res,type="I")
pca_df = as.data.frame(pca_res$x)
ggplot(pca_df, aes(x=PC1,y=PC2,col =Class)+geom_point(alpha=0.5))

