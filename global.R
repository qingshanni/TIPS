r_encoding <- "UTF-8"
r_path <- '.'

library('ggplot2')
library("shiny")
library("shinydashboard")
library("remotes")
library("fs")
library('Seurat', lib.loc = path_home_r("R/R-3.6.3/old_library/"))
library('monocle')
library('kohonen')
library("viridis") 
library('switchde')
options(shiny.maxRequestSize = 100*1024^2) 