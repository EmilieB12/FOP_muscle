library(ggplot2)
library(dplyr)
library(cowplot)
library(Matrix)
library(Seurat)

########Merged Seurat v3######
ctrl.data <- Read10X(data.dir = "1323/filtered_gene_bc_matrices/GRCh38/")
fop.data <- Read10X(data.dir = "F3/filtered_gene_bc_matrices/GRCh38/")

# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "ips_ctrl", min.cells = 5)
ctrl$ips <- "ctrl"
fop <- CreateSeuratObject(counts = fop.data, project = "ips_fop", min.cells = 5)
fop$ips <- "fop"

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
ctrl[["percent.mt"]] <- PercentageFeatureSet(object = ctrl, pattern = "^MT-")
fop[["percent.mt"]] <- PercentageFeatureSet(object = fop, pattern = "^MT-")

ips <- merge(x = ctrl, y = fop)
VlnPlot(ips, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "ips")
plot1 <- FeatureScatter(ips, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ips, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

ips <- subset(ips, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
ips.list <- SplitObject(ips, split.by = "ips")

for (i in 1:length(ips.list)) {
  ips.list[[i]] <- NormalizeData(ips.list[[i]], verbose = FALSE)
  ips.list[[i]] <- FindVariableFeatures(ips.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}

ips.anchors <- FindIntegrationAnchors(object.list = ips.list, dims = 1:30)
ips.integrated <- IntegrateData(anchorset = ips.anchors, dims = 1:30)

DefaultAssay(object = ips.integrated) <- "integrated"
##########Scaling the data#####
all.genes <- rownames(x = ips.integrated)
ips.integrated <- ScaleData(object = ips.integrated, features = all.genes, vars.to.regress =  "percent.mt")
####Regress cell cycle### 
##Alternate Workflow, regressing out differences between g2m and s
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ips.integrated <- CellCycleScoring(ips.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ips.integrated$CC.Difference <- ips.integrated$S.Score - ips.integrated$G2M.Score
ips.integrated <- ScaleData(ips.integrated, vars.to.regress = "CC.Difference", features = rownames(ips.integrated))

# Run the standard workflow for visualization and clustering
ips.integrated <- RunPCA(ips.integrated, npcs = 30, verbose = FALSE)
ips.integrated <- RunUMAP(ips.integrated, reduction = "pca", dims = 1:30)
ips.integrated <- FindNeighbors(object = ips.integrated, reduction = "pca", dims = 1:30)
ips.integrated <- FindClusters(object = ips.integrated, resolution = 0.5)

DimPlot(object = ips.integrated, reduction = "umap")

saveRDS(ips.integrated, file = "ips.integrated_v3.rds")



##########Identify differential expressed genes across conditions############
DefaultAssay(ips.integrated) <- "RNA"
ips.integrated <- ScaleData(object = ips.integrated, features = rownames(ips.integrated))

ips.integrated.markers <- FindAllMarkers(ips.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshmuscle = 0.25)
ips.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- ips.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(object = ips.integrated, features = top5$gene)

write.csv(ips.integrated.markers, file=paste("cluster_markers_ips_06.csv", sep=""))




