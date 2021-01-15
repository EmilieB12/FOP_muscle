#Load object
#subset myogenic cells
subset = subset(ips.integrated, idents = c("3", "10", "11"))
DimPlot(object = subset, reduction = "umap", label = TRUE)


library(ggplot2)
library(patchwork)

plot <- DimPlot(subset, reduction = "umap")
select.cells <- CellSelector(plot = plot)

Idents(subset, cells = select.cells) <- "doublets"
levels(subset)

#Subset out all the doublets
subset2 = subset(subset, idents = c("3", "10", "11"))
DimPlot(object = subset2, reduction = "umap", label = TRUE)
subset2 <- FindNeighbors(object = subset2, reduction = "pca", dims = 1:15)
subset2 <- FindClusters(subset2, resolution = 0.4)

##make new UMAP Plot from subset
subset2 <- RunUMAP(object = subset2, dims = 1:15)

DimPlot(object = subset2, reduction = "umap")
DimPlot(object = subset2, reduction = "umap", group.by = "ips")

saveRDS(subset2, file = "subset2_ips.rds")


##########Identify differential expressed genes across conditions############
DefaultAssay(subset2) <- "RNA"
subset2 <- ScaleData(object = subset2, features = rownames(subset2))

subset2.markers <- FindAllMarkers(subset2, only.pos = TRUE, min.pct = 0.25, logfc.threshmuscle = 0.25)
subset2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- subset2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(object = subset2, features = top5$gene)

write.csv(subset2.markers, file=paste("cluster_markers_subset_ips.csv", sep=""))


subset2$celltype.ips <- paste(Idents(subset2), subset2$ips, sep = "_")
subset2$celltype <- Idents(subset2)
Idents(subset2) <- "celltype.ips"

DEips.cluster1 <- FindMarkers(subset2, ident.1 = "1_ctrl", ident.2 = "1_fop", verbose = FALSE)
head(DEips.cluster1, n = 15)
write.csv(DEips.cluster1, file=paste("DEips.cluster1E.csv", sep=""))





#or for primary vs ips
subset = subset(muscle, idents = c("2", "7", "10"))
DimPlot(object = subset, reduction = "umap", legend = TRUE, label = TRUE)


plot <- DimPlot(subset, reduction = "umap")
select.cells <- CellSelector(plot = plot)

Idents(subset, cells = select.cells) <- "non-myogenic cells"
levels(subset)

#Subset out all the doublets
subset2 = subset(subset, idents = c("2", "7", "10"))
DimPlot(object = subset2, reduction = "umap", legend = TRUE, label = TRUE)

