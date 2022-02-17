library(monocle3)
library(BiocGenerics)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
cds <- as.cell_data_set(hspc.combined)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
p1 <- plot_cells(cds, show_trajectory_graph = T)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = T)
#pdf("Overlap_trajectory")
wrap_plots(p1, p2)
#dev.off()

cds <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(hspc.combined)
cds <- learn_graph(cds)
#pdf("Overlap_trajcetroy_labelled")
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)
#dev.off()

max.avp <- which.max(unlist(FetchData(hspc.combined, vars = 'ident')))
max.avp <- colnames(hspc.combined)[max.avp]
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = max.avp)
#pdf("Overlap_trajectory_reduced_labelled")
plot_cells(cds, color_cells_by = "partition", label_cell_groups = T, label_leaves = T,label_branch_points = T)
#dev.off()

hspc.combined.sub <- as.Seurat(cds)
hspc.combined.sub <- as.Seurat(cds, assay = "integrated")
FeaturePlot(hspc.combined.sub, "monocle3_pseudotime")

cds_ctrl_sub <- as.Seurat(cds_ctrl)
> cds_ctrl_sub <- as.Seurat(cds_ctrl, assay = "integrated")
> FeaturePlot(cds_ctrl_sub, "monocle3_pseudotime")


