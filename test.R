rm(list=ls())
cat('\f')
library('ggplot2')
library("shiny")
library("shinydashboard")
library("Seurat")
library('monocle')
library('kohonen')
library("viridis") 
library('switchde')
source('base/share.R')


dd   <- as.matrix(read.csv( "data/lung_exprs_data.csv",row.names = 1))
idx1 <- (rowSums(dd) > 0)
idx2 <- (apply(dd,1,sd) > 0)
dd   <- dd[idx1&idx2, ]
## log data
dd <- log2(dd+1)
gdata_expr <- dd



gdata_phenotype <- read.csv( "data/lung_phenotype_data.csv",row.names = 1)
gene_set        <- readGMT("data/ss.gmt")

### dimision reduction
dd <- CreateSeuratObject(gdata_expr,meta.data =gdata_phenotype)
dd <- NormalizeData(object = dd)
dd <- FindVariableFeatures(dd,selection.method = "vst")
dd <- ScaleData(object = dd)
dd <- RunPCA(object = dd)
dd <- RunTSNE(object = dd)
dd <- RunUMAP(object = dd,dims = 1:5)
dd <- FindNeighbors(dd)
dd <- FindClusters(object = dd)
dd <- FindClusters(object = dd,resolution = 1)
dd <- FindClusters(object = dd,resolution = 1.2)
dd <- FindClusters(object = dd,resolution = 0.6)
dd <- FindClusters(object = dd,resolution = 0.4)

DimPlot(object = dd,reduction = "umap",label = TRUE)

sdata <- dd


############################################# Trajectory analysis ##########################







