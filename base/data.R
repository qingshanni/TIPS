
##############   data   ###############################

gdata_expr_raw <- reactive({
    inFile <- input$param_file_exprData
    req(inFile)
    
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, 
                 expr = {
                      dd <- as.matrix(read.csv( inFile$datapath,row.names = 1))
                  })
  dd
  
})


gdata_expr <- reactive({
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 
                 if(input$param_file_inputSource){
                   dd <- as.matrix(read.csv( "data/lung_exprs_data.csv",row.names = 1))
                   dd <- log2(dd+1)
                 }else{
                   inFile <- input$param_file_exprData
                   req(inFile)
                   dd <- gdata_expr_raw()
                   
                   ####### filter 
                   idx1 <- T
                   idx2 <- T
                   if(input$param_file_filter_sum0){
                     idx1 <- (rowSums(dd) > 0)
                   }
                   if(input$param_file_filter_sd0){
                     idx2 <- (apply(dd,1,sd) > 0)
                   }            
                   dd <- dd[idx1&idx2, ]
                   
                   ###  log data
                   if(input$param_file_exprDatalog){
                     dd <- log2(dd+1)
                   }
                 }
                 
               })
  dd
  
})


gdata_phenotype <- reactive({
  if(input$param_file_inputSource){
    dd <- read.csv( "data/lung_phenotype_data.csv",row.names = 1)
  }else{
    inFile <- input$param_file_phenoData
    req(inFile)
    dd <- read.csv( inFile$datapath,row.names = 1)
  }
  dd
})


gene_set <- reactive({
  if(input$param_file_inputSource){
    dd <- readGMT("data/ss.gmt")
  }else{
    tp <- input$parm_file_GeneSets_type
    
    if(  tp == "in"){
      
      fn <- file.path('data',paste(input$parm_GeneSets_dataset,'.gmt' ,sep=''))
      
    }else{
      inFile <- input$parm_file_GeneSets
      req(inFile)
      
      fn <- inFile$datapath
    }
    
    dd <- readGMT(fn)
  }
  dd
})



##############################  dimreduce  ##################################
####    seurat version   2.3.4

SeuratData <- reactive({
  
  if(getInputValue(input$btn_dimreduce,0) == 0){
    return(NULL)
  }
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 isolate({
                 
                 dd <- CreateSeuratObject(gdata_expr(),meta.data =gdata_phenotype() )

                 if (getInputValue(input$parm_dimreduce_norm,TRUE)){
                   dd <- NormalizeData(object = dd)
                 }
                 
                 xlow    <- input$parm_dimreduce_xlow
                 
                 
                 
                 xhigh <- Inf
                 if(input$parm_dimreduce_xhightype == 'usr'){
                   xhigh   <- input$parm_dimreduce_xhigh
                 }
                 ycutoff <- input$parm_dimreduce_ycutoff

                 dd <- FindVariableGenes(object = dd, x.low.cutoff = xlow, x.high.cutoff = xhigh , y.cutoff = ycutoff)
                 dd <- ScaleData(object = dd)
                 dd <- RunPCA(object = dd)
                 dd <- RunTSNE(object = dd, genes.use = dd@var.genes, perplexity = input$parm_dimreduce_tSNE_ppl)
                 dd <- RunUMAP(object = dd, genes.use = dd@var.genes, reduction.use = "umap", n_neighbors = input$parm_dimreduce_umap_kNeighbors)
                 dd <- FindClusters(object = dd, reduction.type = "umap", dims.use = 1:2, print.output = F)
                 dd <- FindClusters(object = dd, reduction.type = "umap", dims.use = 1:2, print.output = F, resolution = 1)
                 dd <- FindClusters(object = dd, reduction.type = "umap", dims.use = 1:2, print.output = F, resolution = 1.2)
                 dd <- FindClusters(object = dd, reduction.type = "umap", dims.use = 1:2, print.output = F, resolution = 0.6)
                 dd <- FindClusters(object = dd, reduction.type = "umap", dims.use = 1:2, print.output = F, resolution = 0.4)
                 
                 })
               })
  dd
})

 
####   monocle version 2.6.4

scData <- reactive({
  if(getInputValue(input$btn_trajectory,0) == 0){
    return(NULL)
  }

  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 isolate({
                 
                 if(input$parm_traj_NegBinom){
                   fdata           <-  data.frame( gene_short_name = rownames(gdata_expr()))
                   rownames(fdata) <- rownames(gdata_expr())
                   fd <- new('AnnotatedDataFrame', data = fdata) 
                   pd <- new('AnnotatedDataFrame', data = SeuratData()@meta.data) 
                   sc.data  <- newCellDataSet(gdata_expr(),phenoData = pd, featureData = fd, 
                                              lowerDetectionLimit = 0.1)
                   sc.data <- estimateSizeFactors(sc.data)
                   sc.data <- estimateDispersions(sc.data)
                 }else{  
                   sc.data <- importCDS(SeuratData(), import_all = T)
                 }  
                   
                   sc.data <- detectGenes(sc.data, min_expr = 0.1)
                   sc.data <- setOrderingFilter(sc.data, SeuratData()@var.genes)
                   if(input$parm_traj_norm){
                     sc.data <- reduceDimension(sc.data, reduction_method = "DDRTree")
                   }else{
                     sc.data <- reduceDimension(sc.data, reduction_method = "DDRTree", norm_method = "none")
                   }
                   sc.data <- orderCells(sc.data, reverse = input$parm_traj_reverse)                   
                })
               })
  sc.data
})

# ##########################################################
# SeuratData <- reactive({
#   
#   #  input$btn_trajectory
#   
#   withProgress(message = 'Calculation in progress',
#                detail = 'This may take a while...', value = 0, 
#                expr = {
#                  #                 isolate({
#                  
#                  dd <- CreateSeuratObject(counts = gdata_expr(),meta.data =gdata_phenotype() )
#                  dd <- NormalizeData(object = dd)
#                  dd <- FindVariableFeatures(object = dd)
#                  dd <- ScaleData(object = dd)
#                  dd <- RunPCA(object = dd)
#                  dd <- RunTSNE(object = dd)
#                  dd <- RunUMAP(object = dd,features = dd@assays[[1]]@var.features)
# 
#                  #                 })
#                })
#   dd
# })
# 
# 
# ##########################################################
# 
# 
# scData <- reactive({
#   
# #  input$btn_trajectory
#   
#   withProgress(message = 'Calculation in progress',
#                detail = 'This may take a while...', value = 0, 
#                expr = {
# #                 isolate({
#                    
#                    fdata           <-  data.frame( gene_short_name = rownames(gdata_expr()))
#                    rownames(fdata) <- rownames(gdata_expr())
#                    fd <- new('AnnotatedDataFrame', data = fdata) 
#                    pd <- new('AnnotatedDataFrame', data = gdata_phenotype()) 
#                    sc.data  <- newCellDataSet(gdata_expr(),phenoData = pd, featureData = fd, 
#                                               lowerDetectionLimit = 0.1,
#                                               expressionFamily = negbinomial.size())
#                    sc.data <- estimateSizeFactors(sc.data)
#                    sc.data <- estimateDispersions(sc.data)
#                    
#                    sc.data <- detectGenes(sc.data, min_expr = 0.1)
#                    sc.data <- reduceDimension(sc.data, reduction_method = "DDRTree")
#                    sc.data <- orderCells(sc.data)
# #                 })
#                })
#   sc.data
# })


PseudotimeListAll <- reactive({
  
  if(getInputValue(input$btn_geneset,0) == 0){
    return(NULL)
  }
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0,
               expr = {
                 isolate({
                   #### filter by gene number cutoff
                   gene_num_cutoff <- input$parm_GeneSets_gcutoff  ###
                   gene_set.new    <- lapply(gene_set(), function(genes){intersect(genes, rownames(scData())) })
                   gene_set.new.n  <- lapply(gene_set.new,length)
                   gene_set.new    <- gene_set.new[gene_set.new.n >= gene_num_cutoff]
                   
                   
                   ##### 
                   gene_set.idx <- NULL   
                   v.Pseudotime <- list()
                   k <- 0
                   for(i in 1:length(gene_set.new)){
                     tryCatch({
                        
                        cat('Handle: ', i, '\n')
                        genes      <- gene_set.new[[i]]
                        X_mono <- setOrderingFilter(scData(), genes)
                        X_mono <- reduceDimension(X_mono, reduction_method = "DDRTree", norm_method = "none")
                        X_mono <- orderCells(X_mono)
                        v1 <- pData(phenoData(X_mono))$Pseudotime  
                        k  <- k + 1
                        gene_set.idx <- c(gene_set.idx,i)
                        v.Pseudotime[[k]] <- v1
                        incProgress(i/length(gene_set.new))
                       },error=function(e){})
                   }
                   gene_set.new <- gene_set.new[gene_set.idx]
                   gv           <- pData(phenoData(scData()))$Pseudotime 
                   Correlation  <- sapply(v.Pseudotime,function(v1){ cor(gv,v1)})
                   vsort <- sort(abs(Correlation),decreasing=T, index.return=T)
                   
                   gene_set.new <- gene_set.new[vsort$ix]
                   v.Pseudotime <- v.Pseudotime[vsort$ix]
                   Correlation  <- vsort$x
                   
                   SetSize     <- sapply(gene_set.new,length)
                   log10pval   <- sapply(v.Pseudotime,function(v1){ -log10(cor.test(gv,v1)$p.value) })
                   
                   corData     <- data.frame( Pathway     = names(gene_set.new), 
                                              Correlation = Correlation, 
                                              SetSize     = SetSize, 
                                              log10pval   = log10pval)
                   rownames(corData) <- names(gene_set.new)
                 })
               })
  
  list(gene_set = gene_set.new, v.Pseudotime = v.Pseudotime,corData = corData )
})


PseudotimeListAllDis <- reactive({
  
  if(getInputValue(input$btn_geneset,0) == 0){
    return(NULL)
  }
  
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0,
               expr = {
                 all <- PseudotimeListAll()
                 
                 gtime <-  pData(phenoData(scData()))$Pseudotime
                 ptime <-  all$v.Pseudotime
                 
                 n  <- input$parm_GeneSets_nbin
                 ghist     <- hist(gtime, breaks=getBreaks(gtime,n), plot=F)
                 Correlation  <- sapply(ptime,function(p){ 
                   phist <- hist(p, breaks=getBreaks(p,n), plot=F)
                   cor(ghist$density, phist$density)
                 })
                 log10pval  <- sapply(ptime,function(p){ 
                   phist <- hist(p, breaks=getBreaks(p,n), plot=F)
                   -log10(cor.test(ghist$density, phist$density)$p.value)
                 })
                 
                 corData     <- data.frame( Pathway     = all$corData$Pathway, 
                                            SetSize     = all$corData$SetSize, 
                                            Correlation = all$corData$Correlation,
                                            log10pval   = all$corData$log10pval,
                                            Correlation.dis = Correlation,
                                            log10pval.dis   =log10pval,
                                            stringsAsFactors = F
                 )
                 rownames(corData) <- rownames(all$corData)
                 
               })
  list(gene_set = all$gene_set, v.Pseudotime = all$v.Pseudotime,corData = corData )
})



# PseudotimeList <- reactive({
#  
#   if(getInputValue(input$btn_geneset,0) == 0){
#     return(NULL)
#   }
#   withProgress(message = 'Calculation in progress',
#                detail = 'This may take a while...', value = 0, 
#                expr = {
#                 isolate({
#                     gene_num_topn   <- min(input$parm_GeneSets_topn,length(PseudotimeListAll()$gene_set))
#                     
#                     lst <- list(gene_set     = PseudotimeListAll()$gene_set[1:gene_num_topn], 
#                                 v            = PseudotimeListAll()$v.Pseudotime[1:gene_num_topn],
#                                 corData      = PseudotimeListAll()$corData[1:gene_num_topn,] )
#                  })
#                })
#   lst
# })







# GenesetPathCor<- reactive({
#   
#   vs          <- PseudotimeList()$v
#   gene_set    <- PseudotimeList()$gene_set
#   v           <- pData(phenoData(scData()))$Pseudotime
#   Correlation <- PseudotimeList()$Correlation
#   SetSize     <- sapply(gene_set,length)
#   log10pval   <- sapply(vs,function(v1){ -log10(cor.test(v,v1)$p.value) })
#   
#   data.frame( Pathway=names(gene_set), Correlation = abs(Correlation), SetSize=SetSize, log10pval = log10pval)
#   
# })





####  correlation of gene expression and pseudotime
PtimeCorExpr <- reactive({
  
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 
                 expr  <-  exprs(scData())
                 ptime <- pData(phenoData(scData()))$Pseudotime  
                 
                 vcor  <- apply(expr,1,function(v){ cor(ptime,v) })
                 
                 
                 ### sort #
                 vsort <- sort(abs(vcor),decreasing=T, index.return=T)
                 
               })
  list(vcor = vcor, vsort = vsort)
})


exprSelectMX <- reactive({
  
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 
                 n      <- input$parm_GeneCor_genenum
                 expr   <-  exprs(scData())
                 vsort  <-  PtimeCorExpr()$vsort
                 vcor   <-  PtimeCorExpr()$vcor
                 
                 expr_sel <- expr[vsort$ix[1:n],]
                 vcor_sel <- vcor[vsort$ix[1:n]]
                 
                 vssort   <- sort(vcor_sel,decreasing=T, index.return=T)
                 expr_sel <- expr_sel[vssort$ix,]
                 
               })
  expr_sel
})



data.som <- reactive({
  # if(getInputValue(input$btn_som,0) == 0){
  #   return(NULL)
  # }
  
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 
                 expr  <-  exprs(scData())
                 ptime <- pData(phenoData(scData()))$Pseudotime  
                 gene_set <- PseudotimeListAll()$gene_set
                 if( input$parm_som_type == "top"){
                       n            <- min(input$parm_som_pathnum, length(gene_set))
                       gene_set.new <- gene_set[1:n]
                   
                 }else{
                   if (is.null(input$parm_som_pathselect)) {
                     stop("Please select geneset!")
                   }
                   gene_set.new <- gene_set[input$parm_som_pathselect]
                   
                 }
                 expr.sum <- sapply(gene_set.new,function(genes){
                 colSums(expr[genes,])
                 }) 
                 mx <- cbind(ptime,expr.sum) 
                 df <- scale(data.frame(mx))
                 colnames(df) <- c("Pseudotime", paste("Score", names(gene_set.new),sep = "_" ))
                 dd.som <- som(scale(df), grid = somgrid(input$parm_som_gridx, input$parm_som_gridy, "rectangular"))
               })
  dd.som
})



## overall switchgene table
data.switch.org <-  reactive({
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 mx     <- gdata_expr()
                 
                 #### filter genes by ratio of nonzero expression
                 nn <- apply(mx,1,function(v){ sum(v>0)})  ## number of nonzero expression
                 mx <- mx[nn > input$parm_switchHist_nonzero_threshold * dim(mx)[2] / 100,]
                 
                 ptime  <- pData(phenoData(scData()))$Pseudotime
                 mx_sde <- switchde(mx, ptime,lower_threshold=input$parm_switchHist_lower_threshold) #overall switchgene table, exportable as one output
                 rownames(mx_sde) <- mx_sde$gene
               })
  mx_sde
})


data.switch <-  reactive({
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 sde  <- data.switch.org()
                 gene <- sde$gene
                 idx  <- (sde$qval <= input$parm_switchHist_qcutoff) & (abs(sde$k) >=input$parm_switchHist_kcutoff) & (sde$mu0 >input$parm_switchHist_ucutoff)
                 sde  <- sde[idx, ]   ### filter by q value
                 rownames(sde)<-gene[idx]
               })
  sde
})






## pathway switchgene
pathway.switch <-  reactive({
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, 
               expr = {
                 
                 sde      <- data.switch()
                 gset     <- gene_set()
                 
                 #save(sde,gset,file='tt.RData')
                 #### filter by gene number cutoff
                 gene_num_cutoff <- input$parm_switchScater_ncutoff  ###
                 gene_set.new    <- lapply(gset, function(genes){intersect(genes, rownames(sde)) })
                 gene_set.new.n  <- lapply(gene_set.new,length)
                 gene_set.new    <- gene_set.new[gene_set.new.n >= gene_num_cutoff]
                 
                 pathway_sde <- lapply(gene_set.new,function(genes){ 
                   s <- sde[genes,] 
                   rownames(s) <- s$gene
                   s
                 })
                 names(pathway_sde) <- names(gene_set.new)
                 
               })
  save(sde,gset,pathway_sde,file="tt.RData")
  
  pathway_sde
})


#####  




