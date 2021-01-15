library (ggplot2)
library(cowplot)
library(Matrix)
library(Seurat)
library(dplyr)

setwd("/Volumes/Genetic/1 Seq/SingleCell/02_iPS-SatelliteDiffn")
FOP3hg19.data <- Read10X(data.dir = "F3-hg19/outs/filtered_gene_bc_matrices/hg19/")
WT1323hg19.data <- Read10X(data.dir = "1323-hg19/outs/filtered_gene_bc_matrices/hg19/")
setwd("/Volumes/SSD/scRNAseq/Pom10x")
vastus.data <- Read10X(data.dir = "56-vastus/hg19/")

# Set up sample integration
vastus <- CreateSeuratObject(counts = vastus.data, project = "1 Primary SC", min.cells = 7)
vastus$muscle <- "1 Primary SC"
wt <- CreateSeuratObject(counts = WT1323hg19.data, project = "2 WT_iPS", min.cells = 7)
wt$muscle <- "2 WT_iPS"
fop <- CreateSeuratObject(counts = FOP3hg19.data, project = "3 FOP_iPS", min.cells = 7)
fop$muscle <- "3 FOP_iPS"

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
vastus[["percent.mt"]] <- PercentageFeatureSet(object = vastus, pattern = "^MT-")
wt[["percent.mt"]] <- PercentageFeatureSet(object = wt, pattern = "^MT-")
fop[["percent.mt"]] <- PercentageFeatureSet(object = fop, pattern = "^MT-")


muscle <- merge(x = vastus, y = list(wt, fop))
VlnPlot(muscle, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "muscle")
plot1 <- FeatureScatter(muscle, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(muscle, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


# Filter unwanted cells
muscle <- subset(muscle, subset =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
muscle.list <- SplitObject(muscle, split.by = "muscle")

# Normalize object, find variable feature
for (i in 1:length(muscle.list)) {
  muscle.list[[i]] <- NormalizeData(muscle.list[[i]], verbose = FALSE)
  muscle.list[[i]] <- FindVariableFeatures(muscle.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

# Integration of all the muscle cell datasets
muscle.anchors <- FindIntegrationAnchors(object.list = list(vastus, wt, fop), dims = 1:30)
muscle.integrated <- IntegrateData(anchorset = muscle.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(object = muscle.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
#all.genes <- rownames(muscle.integrated)
#muscle.combined <- ScaleData(muscle.integrated, features = all.genes)
muscle.combined <- ScaleData(muscle.integrated, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)

##regressing out differences between g2m and s
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
muscle.combined <- CellCycleScoring(muscle.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
muscle.combined$CC.Difference <- muscle.combined$S.Score - muscle.combined$G2M.Score
muscle.combined <- ScaleData(muscle.combined, vars.to.regress = "CC.Difference", features = rownames(muscle.combined))

muscle.combined <- ScaleData(muscle.integrated, vars.to.regress = c("percent.mt"), verbose = FALSE)
muscle.combined <- ScaleData(muscle.integrated, verbose = FALSE)
# Run the standard workflow for visualization and clustering
muscle.combined <- RunPCA(muscle.combined, npcs = 30, verbose = FALSE)
muscle.combined <- RunUMAP(muscle.combined, reduction = "pca", dims = 1:30)
muscle.combined <- FindNeighbors(object = muscle.combined, reduction = "pca", dims = 1:30)
muscle.combined <- FindClusters(object = muscle.combined, resolution = 0.5)


DimPlot(object = muscle.combined, reduction = "umap")
DimPlot(object = muscle.combined, reduction = "umap", group.by = "muscle")

DefaultAssay(muscle.combined) <- "RNA"
FeaturePlot(object = muscle.combined, features = c("PAX7", "MYF5", "MYOD1", "MYOG"), min.cutoff = "q9")
FeaturePlot(object = muscle.combined, features = c("MYOG", "MEF2C", "ASPN", "MAP2", "NEUROD4", "ASCL1"), min.cutoff = "q9")





##########Identify conserved cell type markers##########
# find markers for every cluster compared to all remaining cells, report only the positive ones
muscle.combined.markers <- FindAllMarkers(object = muscle.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
muscle.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

top5 <- muscle.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(object = muscle.combined, features = top5$gene)
write.csv(muscle.combined.markers, file=paste("primary_ips combined_v3_markers.csv", sep=""))


