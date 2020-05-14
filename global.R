r_encoding <- "UTF-8"

r_path <- '.'
## pkgs used
#r_pkgs <- c('ggplot2',"shiny", "pheatmap", "shinydashboard","NormqPCR",'limma','reshape')
#sapply(r_pkgs, require, character.only = TRUE)
#options(repos = BiocInstaller::biocinstallRepos())
library('ggplot2')
library("shiny")
#library("pheatmap")
library("shinydashboard")
library("Seurat")
library('monocle')
library('kohonen')
library("viridis") #may as use viridis to make the gene expression heatmap prettier honestly
#library('scales')
#library('matrixStats')
#library('qusage')
library('switchde')
options(shiny.maxRequestSize = 100*1024^2) 