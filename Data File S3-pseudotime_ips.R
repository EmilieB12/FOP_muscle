library(stats4)
library(splines)
library(VGAM)
library(parallel)
library(irlba)
library(Matrix)
library(DDRTree)
library(BiocGenerics)
library(Biobase)
library(ggplot2)
library(Seurat)
library(monocle)

#Load object
#Rename seurat object
seuratX <- subset2

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(seuratX@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seuratX@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.915


#Run ordering algorithm
var_genes <-seuratX[["RNA"]]@var.features
ordering_genes <- var_genes

HSMM <- setOrderingFilter(HSMM, ordering_genes)
print(dim(exprs(HSMM)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM <- reduceDimension(HSMM,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=2,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

# First decide what you want to color your cells by
print(head(pData(HSMM)))

## order cells change colors and theta to match your plot
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, 
                     color_by = "seurat_clusters",
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 1)

plot_cell_trajectory(HSMM, color_by = "muscle",
                     show_branch_points = FALSE,
                     show_tree = TRUE, 
                     cell_size = 1.5) + scale_color_manual(values = c("dodgerblue1", "black", "deeppink"))



#Plots
plot_cell_trajectory(HSMM, markers = c("PAX7", "MYL1", "COL1A1", "DCN"), use_color_gradient = T, cell_size = 1, show_branch_points = FALSE)
plot_cell_trajectory(HSMM, color_by = "Pseudotime", show_branch_points = FALSE)

marker.genes <- c("PAX7", "MYF5", "MYOD1", "MYOG")
plot_genes_branched_pseudotime(HSMM[marker.genes,],
                               branch_point = 1, panel_order = c("PAX7", "MYF5", "MYOD1", "MYOG"),
                               color_by = "seurat_clusters",
                               cell_size = 0.5,
                               ncol = 2)



#Save object
save(HSMM, file = "Pseudo_myosubset_primaryVSips_may.rds")

#order the cells to use the desired "root state" of the trajectory
#Here setting the root state to state 2 since I want to compare branch 1 to branch 2 in later analysis
#you must already have called orderCells() once to use this argument.
HSMM <- orderCells(HSMM, root_state = 3)

#Analyzing branches in single-cell trajectories
#Modules of genes that co-vary across pseudotime
BEAM_res2 <- BEAM(HSMM, branch_point = 1, cores = 1)
BEAM_res2 <- BEAM_res2[order(BEAM_res2$qval),]
BEAM_res2 <- BEAM_res2[,c("gene_short_name", "pval", "qval")]
write.csv(BEAM_res2, 
          file=paste("BEAMres2_markers_primaryVSips_may.csv", sep=""))

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("PAX7", "MYF5", "MYOD1", "MYOG", "MYL1", "ACVR2A", "COL1A1", "COL1A2", "COL3A1", "DCN", "OGN", "BGN", "IGF2", "DLK1", "CAV1", "SPRY1", "HEY1", "HES1", "BMP4", "ID1", "ID3", "CHODL", "APOE", "MMP2", "TGFB2")))

plot_genes_branched_heatmap(HSMM[to_be_tested,],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)


                                         
                                           
                                           
