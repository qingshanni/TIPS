####### main:  data #####

output$showExprsTable <- renderTable({
  n <- min(nrow(gdata_expr()),5)
  m <- min(ncol(gdata_expr()),5)
  gdata_expr()[1:n,1:m]
}, bordered = T, rownames = T)

output$showExprsSummary <- renderUI({
  div( h3('Expression data'),
       p("Rows: ", nrow(gdata_expr()), " Columns: ", ncol(gdata_expr())),
       p("The first 5 rows and columns are shown below"),
       tableOutput("showExprsTable")
  )
})

output$showPhenotypeTable <- renderTable({
  n <- min(nrow(gdata_phenotype()),5)
  m <- min(ncol(gdata_phenotype()),5)
  gdata_phenotype()[1:n,1:m]
}, bordered = T, rownames = T)

output$showPhenotypeSummary <- renderUI({
  div( hr(), 
       h3('Phenotype data'),
       p("Rows: ", nrow(gdata_phenotype()), " Columns: ", ncol(gdata_phenotype())),
       p("The first 5 rows and columns are shown below"),
       tableOutput("showPhenotypeTable")
  )
})

output$showGeneSetSummary <- renderUI({
  div( hr(), 
       h3('GeneSet data'),
       p("Set Numbers: ", length(gene_set())),
       p("The first 3 geneset are shown below"),
       p(strong(names(gene_set())[1]),'(', length(gene_set()[[1]])  ,')'),p(paste(gene_set()[[1]], collapse = ',')),
       p(strong(names(gene_set())[2]),'(', length(gene_set()[[2]])  ,')'),p(paste(gene_set()[[2]], collapse = ',')),
       p(strong(names(gene_set())[3]),'(', length(gene_set()[[3]])  ,')'),p(paste(gene_set()[[3]], collapse = ','))
  )
})





output$ui_gdata <- renderUI({
  
  fluidRow( 
        
        column(width = 3, wellPanel( checkboxInput("param_file_inputSource", "Load demo data", FALSE)),
                          conditionalPanel(condition = "input.param_file_inputSource == 0",
                                   wellPanel(fileInput('param_file_exprData', label = h5('Choose sequencing data(.csv)'),
                                                       accept = c('text/csv')),
                                             fileInput('param_file_phenoData', label = h5('Choose phenotype data(.csv)'),
                                                       accept = c('text/csv')),
                                             radioButtons("parm_file_GeneSets_type", "Pathway type:",
                                                          c("Build in" = "in",
                                                            "User input" = "usr"),inline=T,selected ="usr"),
                                             conditionalPanel(
                                               condition = "input.parm_file_GeneSets_type == 'in'",
                                               selectInput("parm_GeneSets_dataset", "Select dataset:",
                                                           c("microRNA msigdb"="EnrichedmicroRNA_msigdb",
                                                             "Hallmarks msigdb" = "Hallmarks_msigdb", 
                                                             "TFs msigdb" = "TFs_msigdb",
                                                             "Biocarta"  ="Pathways_Biocarta",
                                                             "KEGG"  ="Pathways_KEGG",                                                        
                                                             "Reactome"  ="Pathways_Reactome"))),
                                             
                                             conditionalPanel(
                                               condition = "input.parm_file_GeneSets_type == 'usr'",
                                               fileInput('parm_file_GeneSets', label = h5('Choose pathway genes data(.gmt)'),
                                                         accept = c('text/gmt'))
                                             ))
                                   ,
                                   wellPanel( checkboxInput("param_file_exprDatalog", "Log sequencing data", TRUE),
                                              checkboxInput("param_file_filter_sum0", "Filter genes with all zero value", TRUE),
                                              checkboxInput("param_file_filter_sd0", "Filter genes with zero variance ", TRUE)
                                   )
                                   )
               ),
        column(width = 9,  box( uiOutput("showExprsSummary") , 
                                uiOutput("showPhenotypeSummary"), 
                                uiOutput("showGeneSetSummary") ,                
                                width = NULL,  collapsible = F,title = "Data", status = "primary", solidHeader = TRUE)) 
        ) 
})


############################################# Dimensional reduction analysis ##########################
output$dimreducePlot <- renderPlot({
  if(is.null(SeuratData())) {
    return(NULL)
  } 
  
  if(input$parm_dimreduce_bycolor %in% colnames(SeuratData()@meta.data) ){
    p <- DimPlot(object = SeuratData(), reduction = input$parm_dimreduce_md, group.by = input$parm_dimreduce_bycolor,pt.size = input$parm_dimreduce_ptsize)
  }else{
    p <- DimPlot(object = SeuratData(), reduction = input$parm_dimreduce_md,pt.size = input$parm_dimreduce_ptsize)
  }
  p
},height=function(){input$parm_dimreduce_figh},width = function(){input$parm_dimreduce_figw})

output$parm_dimreduce_dl_SeuratData <- downloadHandler(
  filename = function() { 
    'data_seurat.rds'
  },
  content  = function(file) {
    saveRDS(SeuratData(),file = file )
  },
  contentType =  ".rds"
)

output$parm_dimreduce_dl_plot <- downloadHandler(
  filename = function() { 
    'dimreduce_plot.pdf'
  },
  content  = function(file) {
    pdf(file,width = input$parm_dimreduce_figw/75, height=input$parm_dimreduce_figh/75)
    if(input$parm_dimreduce_bycolor %in% colnames(SeuratData()@meta.data) ){
      p <- DimPlot(object = SeuratData(), reduction = input$parm_dimreduce_md, group.by = input$parm_dimreduce_bycolor,pt.size = input$parm_dimreduce_ptsize)
    }else{
      p <- DimPlot(object = SeuratData(), reduction = input$parm_dimreduce_md,pt.size = input$parm_dimreduce_ptsize)
    }
    p
    dev.off()
  },
  contentType =  NA
)



output$ui_parm_dimreduce_bycolor <- renderUI({
    selectizeInput(inputId = "parm_dimreduce_bycolor", label = "Color by:", 
               choices = {if(!is.null(SeuratData())) c('[NONE]',colnames(SeuratData()@meta.data))}, 
               selected = '[NONE]', multiple = F)
})



output$ui_dimreduce <- renderUI({

  fluidRow( 
    
    column(width = 3,  wellPanel(        checkboxInput("parm_dimreduce_norm", "Normalize the data", TRUE),
                                         
                                         numericInput("parm_dimreduce_xlow", "x.low.cutoff:", 0.1),
                                         
                                         radioButtons("parm_dimreduce_xhightype", "x.high.cutoff:",
                                                      c("Not set" = "Inf",
                                                        "Input by user" = "usr"),inline=T),
                                         conditionalPanel(
                                           condition = "input.parm_dimreduce_xhightype == 'usr'",
                                           numericInput("parm_dimreduce_xhigh", NULL, 100000000)),
                                         
                                         numericInput("parm_dimreduce_ycutoff", "y.cutoff:", 1),
                                         
                                         numericInput("parm_dimreduce_tSNE_ppl", "tSNE perplexity:", 20, min = 1,step = 1),
                                         numericInput("parm_dimreduce_umap_kNeighbors", "kNeighbors:", 30, min = 1,step = 1),
                                         selectizeInput(inputId = "parm_dimreduce_md", label = "dimensional reduction method:", 
                                                        choices = c("UMAP" = 'umap',"tSNE" = 'tsne', "PCA" = 'pca'), 
                                                        selected = "umap", multiple = F),

                                         actionButton("btn_dimreduce", "Compute")
                                 ),
                       wellPanel(       sliderInput("parm_dimreduce_figw", "Figure width:", min = 100, max = 1500, value = 600),
                                        sliderInput("parm_dimreduce_figh", "Figure  Hight:", min = 100, max = 1000, value = 400),
                                        uiOutput('ui_parm_dimreduce_bycolor'),
                                        numericInput("parm_dimreduce_ptsize", "Point Size:", 1, min = 1,step = 1)
                                 )
            ),
    column(width = 9,  box(    plotOutput("dimreducePlot",width = "100%", height = "100%"),                         
                               width = NULL, height=NULL, collapsible = TRUE,title = " ", status = "primary", solidHeader = TRUE), 
                       downloadButton('parm_dimreduce_dl_SeuratData', 'Download data (rds)'),
                       downloadButton('parm_dimreduce_dl_plot', 'Download figure(pdf)'))
  ) 
})


############################################# Trajectory analysis ##########################
output$cellTrajectoryPlot1 <- renderPlot({
  if(is.null(scData())) {
    return(NULL)
  } 
 
  if(input$parm_traj_pd %in% colnames(pData(scData())) ){
    p <- plot_cell_trajectory(scData(), color_by = input$parm_traj_pd, cell_size = input$parm_traj_ptsize)
    
  }else{
    p <- plot_cell_trajectory(scData(),cell_size = input$parm_traj_ptsize)
  }
  p
},height=function(){input$parm_traj_figh},width = function(){input$parm_traj_figw})

output$cellTrajectory_dl_monoData <- downloadHandler(
  filename = function() { 
    'data_monocle.rds'
  },
  content  = function(file) {
    saveRDS(scData(),file = file )
  },
  contentType =  ".rds"
)


output$parm_cellTrajectory_dl_plot <- downloadHandler(
  filename = function() { 
    'cellTrajectory_plot.pdf'
  },
  content  = function(file) {

    if(input$parm_traj_pd %in% colnames(pData(scData())) ){
      p <- plot_cell_trajectory(scData(), color_by = input$parm_traj_pd, cell_size = input$parm_traj_ptsize)
      
    }else{
      p <- plot_cell_trajectory(scData(), cell_size = input$parm_traj_ptsize)
    }
    p
    ggsave(file, device = "pdf",width= input$parm_traj_figw/75, height=input$parm_traj_figh/75, units="in")
  },
  contentType =  NA
)

output$ui_parm_traj_pd <- renderUI({
  selectizeInput(inputId = "parm_traj_pd", label = "Color by:", 
                 #  choices = c('Pseudotime',colnames(gdata_phenotype())), 
                  choices = {if(!is.null(scData())) colnames(pData(scData()))}, 
                 selected = 'Pseudotime', multiple = F)
  
})


output$ui_trajectory <- renderUI({
  
   
  
  fluidRow(    
    
    column(width = 3,  wellPanel(        
                                         checkboxInput("parm_traj_norm", "Normalize the data", TRUE),
                                         checkboxInput("parm_traj_NegBinom", "NegBinom (recommended for droplet-based and/or very sparse data)", TRUE),
                                         checkboxInput("parm_traj_reverse", "Reverse the beginning and end points of the learned biological process", FALSE),
                                         actionButton("btn_trajectory", "Compute")
    ),
                    wellPanel(       sliderInput("parm_traj_figw", "Figure width:", min = 100, max = 1500, value = 600),
                                     sliderInput("parm_traj_figh", "Figure  Hight:", min = 100, max = 1000, value = 400),
                                     numericInput("parm_traj_ptsize", "Point Size:", 1, min = 1,step = 1),
                                     uiOutput('ui_parm_traj_pd')
                              )),
    column(width = 9,  box(    plotOutput("cellTrajectoryPlot1",width = "100%", height = "100%"),                         
                               width = NULL, height=NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                               downloadButton('cellTrajectory_dl_monoData', 'Download data (rds)'),
                               downloadButton('parm_cellTrajectory_dl_plot', 'Download figure(pdf)')) 
  ) 
})

###########################################   geneset  Pathway correlations ################################
####  compute Pseudotime using genes in each pathway，
####  Then compute correlation between the each pathway Pseudotime with the global Pseudotime


output$HeatmapPlot <- renderPlot({
    
  if(is.null(PseudotimeListAll())) {
    return(NULL)
  } 
  
  if( input$parm_GeneSets_seltype == "top"){
    
    gene_num_topn   <- min(input$parm_GeneSets_topn,length(PseudotimeListAllDis()$gene_set))
    df              <- PseudotimeListAllDis()$corData[1:gene_num_topn,] 
    
  }else{
    sels            <- input$parm_GeneSets_pathselect
    df              <- PseudotimeListAllDis()$corData[sels,] 
  }
  p1 <- ggplot(df,aes(x = Correlation, y = Pathway, size = SetSize, color = log10pval))+geom_point() 
  p2 <- ggplot(df,aes(x = Correlation.dis, y = Pathway, size = SetSize, color = log10pval.dis))+geom_point()
  p3 <- ggplot(df,aes(x = Correlation.order, y = Pathway, size = SetSize, color = log10pval.order))+geom_point()
  
  
  p4 <- ggplot(df,aes(x = Correlation, y = Correlation.dis)) + theme_bw() + geom_point()
  p5 <- ggplot(df,aes(x = Correlation, y = Correlation.order)) + theme_bw() + geom_point()
  p6 <- ggplot(df,aes(x = Correlation.order, y = Correlation.dis)) + theme_bw() + geom_point()
  
  plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2, ncol=3,labels = c("A", "B","C","D","E","F"))
  
},height=function(){input$parm_GeneSets_figh},width = function(){input$parm_GeneSets_figw})


output$Heatmap_dl_Data <- downloadHandler(
  filename = function() { 
    'data.csv'
  },
  content  = function(file) {
    write.csv(PseudotimeListAllDis()$corData, file = file, row.names = F, quote =F)
  },
  contentType =  "text/csv"
)


output$parm_GeneSets_dl_plot <- downloadHandler(
  filename = function() { 
    'GeneSets_plot.pdf'
  },
  content  = function(file) {
    if( input$parm_GeneSets_seltype == "top"){
      
      gene_num_topn   <- min(input$parm_GeneSets_topn,length(PseudotimeListAllDis()$gene_set))
      df              <- PseudotimeListAllDis()$corData[1:gene_num_topn,] 
      
    }else{
      sels            <- input$parm_GeneSets_pathselect
      df              <- PseudotimeListAllDis()$corData[sels,] 
    }
    p1 <- ggplot(df,aes(x = Correlation, y = Pathway, size = SetSize, color = log10pval))+geom_point() 
    p2 <- ggplot(df,aes(x = Correlation.dis, y = Pathway, size = SetSize, color = log10pval.dis))+geom_point()
    p3 <- ggplot(df,aes(x = Correlation.order, y = Pathway, size = SetSize, color = log10pval.order))+geom_point()
    
    
    p4 <- ggplot(df,aes(x = Correlation, y = Correlation.dis)) + theme_bw() + geom_point()
    p5 <- ggplot(df,aes(x = Correlation, y = Correlation.order)) + theme_bw() + geom_point()
    p6 <- ggplot(df,aes(x = Correlation.order, y = Correlation.dis)) + theme_bw() + geom_point()
    
    p123 <- plot_grid(p1,p2,p3,nrow = 1, ncol=3,labels = c("A", "B","C"))
    p456 <- plot_grid(p4,p5,p6,nrow = 1, ncol=3,labels = c("D", "E","F"))
    
    p<- plot_grid(p123,p456, nrow = 2, ncol=1, labels = c("", ""))
    
    ggsave(file, device = "pdf",width= input$parm_GeneSets_figw/75, height=input$parm_GeneSets_figh/75, units="in")
  },
  contentType =  NA
)


output$ui_geneset_param <- renderUI({
  
  selectInput("parm_GeneSets_pathselect", "Select geneset:", 
              {
                if(is.null(PseudotimeListAll())){
                  pathlist <- NULL
                }else{
                  pathlist <- names(PseudotimeListAll()$gene_set)
                }
                pathlist
              },multiple=T)
})

output$ui_geneset <- renderUI({
  
  fluidRow( 
    column(width = 3, wellPanel(       
      numericInput("parm_GeneSets_gcutoff", "Gene number cutoff:", 20, min = 20),
      actionButton("btn_geneset", "Compute")),
      wellPanel(       
        radioButtons("parm_GeneSets_seltype", "Selection type:",
                     c("Top N" = "top",
                       "User selected" = "usr"),inline=T),
        conditionalPanel(
          condition = "input.parm_GeneSets_seltype == 'top'",
          numericInput("parm_GeneSets_topn", "Top N gene sets (Correlation to pseudo time):", 10, min = 1,step = 1)),
        
        conditionalPanel(
          condition = "input.parm_GeneSets_seltype == 'usr'",
          uiOutput("ui_geneset_param")), 
          numericInput("parm_GeneSets_nbin", "Number of bins for compute distribution correlation:", 10, min = 1,step = 1),
          sliderInput("parm_GeneSets_figw", "Figure width:", min = 100, max = 1500, value = 800),
          sliderInput("parm_GeneSets_figh", "Figure  Hight:", min = 100, max = 1000, value = 600)
      )),
    column(width = 9,  box(    plotOutput("HeatmapPlot",width = "100%", height = "100%"),
                               width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
           downloadButton('Heatmap_dl_Data', 'Download data (csv)'),
           downloadButton('parm_GeneSets_dl_plot', 'Download figure (pdf)') )
  ) 
})



###########################################   Gene correlations  ################################

#### the correlation between Pseudotime and Gene expression

output$GeneCorPlot <- renderPlot({
  

  gene  <- input$parm_GeneCorPlot_gene
  expr  <-  exprs(scData())[gene,]
  ptime <- pData(phenoData(scData()))$Pseudotime
  
  df <- data.frame(ptime = ptime, expr = expr)
  
  ggplot(df,aes(x = ptime, y = expr))+geom_point()+ 
       xlab( "Pseudotime" ) + ylab( "Gene expression" )
  
},height=function(){input$parm_GeneCorPlot_figh},width = function(){input$parm_GeneCorPlot_figw})


output$GeneCorHeatmap <- renderPlot({
  
  mx <- exprSelectMX()
  x  <- colnames(mx) 
  y  <- rownames(mx) 
  data <- expand.grid(Cell=x, Gene=y)
  data$Expression <- as.vector(t(mx))
  ggplot(data, aes(x=Cell, y=Gene, fill= Expression)) + geom_tile(colour = "grey50") +
    scale_fill_gradient(low="white", high="red") + theme(axis.text.x = element_blank()) 
  
},height=function(){input$parm_GeneCorHeatmap_figh},width = function(){input$parm_GeneCorHeatmap_figw})




output$parm_GeneCorPlot_dl_plot <- downloadHandler(
  filename = function() { 
    'GeneCorPlot_plot.pdf'
  },
  content  = function(file) {
    gene  <- input$parm_GeneCorPlot_gene
    expr  <-  exprs(scData())[gene,]
    ptime <- pData(phenoData(scData()))$Pseudotime
    
    df <- data.frame(ptime = ptime, expr = expr)
    
    ggplot(df,aes(x = ptime, y = expr))+geom_point()+ xlab( "Pseudotime" ) + ylab( "Gene expression" )
    ggsave(file, device = "pdf",width= input$parm_GeneCorPlot_figw/75, height=input$parm_GeneCorPlot_figh/75, units="in")
  },
  contentType =  NA
)


output$parm_GeneCorHeatmap_dl_plot <- downloadHandler(
  filename = function() { 
    'GeneCorHeatmap_plot.pdf'
  },
  content  = function(file) {
    mx <- exprSelectMX()
    x  <- colnames(mx) 
    y  <- rownames(mx) 
    data <- expand.grid(Cell=x, Gene=y)
    data$Expression <- as.vector(t(mx))
    ggplot(data, aes(x=Cell, y=Gene, fill= Expression)) + geom_tile(colour = "grey50") +
      scale_fill_gradient(low="white", high="red") + theme(axis.text.x = element_blank()) 
   ggsave(file, device = "pdf",width= input$parm_GeneCorHeatmap_figw/75, height=input$parm_GeneCorHeatmap_figh/75, units="in")

  },
  contentType =  NA
)


output$ui_gene_cor <- renderUI({
  
  genenames <- rownames(scData())
  n <- length(genenames)
  
  
  div(fluidRow( 
    column(width = 3, wellPanel(div('Analysis the correlation between Pseudotime and Gene expression')),
                      wellPanel(    numericInput("parm_GeneCor_genenum", "Gene number cutoff:", 50, min = 5,max=n, step = 1)
                                ),
                      wellPanel(       sliderInput("parm_GeneCorHeatmap_figw", "Figure width:", min = 100, max = 1500, value = 600),
                                       sliderInput("parm_GeneCorHeatmap_figh", "Figure  Hight:", min = 100, max = 1000, value = 400)
           )),
    column(width = 9,  box(plotOutput("GeneCorHeatmap",width = "100%", height = "100%"),
                               width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                       downloadButton('parm_GeneCorHeatmap_dl_plot', 'Download figure(pdf)')) 
  ) ,hr(),p(),
  fluidRow( 
    column(width = 3,  wellPanel( selectizeInput(inputId = "parm_GeneCorPlot_gene", label = "Select gene:", 
                                                    choices = genenames, selected = genenames[1], multiple = F)),
                      wellPanel(  sliderInput("parm_GeneCorPlot_figw", "Figure width:", min = 100, max = 1500, value = 600),
                                  sliderInput("parm_GeneCorPlot_figh", "Figure  Hight:", min = 100, max = 1000, value = 400))
    ),
    column(width = 9,  box(plotOutput("GeneCorPlot",width = "100%", height = "100%"),                          
                           width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                       downloadButton('parm_GeneCorPlot_dl_plot', 'Download figure(pdf)')) 
  ) )
  
})


#############################  SOM ##################################
##  https://clarkdatalabs.github.io/soms/SOM_NBA




output$somPlot <- renderPlot({
  
  if(is.null( pData(scData()))){ return(NULL) }

  phall <- pData(scData()) 
  som   <- data.som()

  path <- input$parm_som_pathway
  ph   <- input$parm_som_ph
  codes <- som$codes
  if( is.list(codes))
    codes <- codes[[1]]
  end
  
  op <- par(mfrow = c(2, 2))
  plot(som, main = "TIPS SOM") 
  plot(som, type = "property", property = codes[,paste('Score_',path,sep='')], main=path) #allows emphasis of the change in a given pathway
  if(ph == '[NONE]'){
    plot(som, type = "mapping", pchs = 20, main = "Density")
  }else{
    ft         <- as.factor(phall[,ph])
    levels(ft) <- hue_pal()(length(levels(ft)))  #so this replaces N factors with N default ggplot colors
    
    plot(som, type = "mapping", pchs = 20, main = ph, col = as.character(ft)) #this now maps the metadata onto the points
  }
  plot(som, type = "count", main = "Node Counts", palette.name = viridis) #heatmap of the number of cells in each cluster
  
  par(op) 
  
},height=function(){input$parm_som_figh},width = function(){input$parm_som_figw})


output$parm_som_dl_plot <- downloadHandler(
  filename = function() { 
    'GeneCorHeatmap_plot.pdf'
  },
  content  = function(file) {
    pdf(file,width = input$parm_som_figw/75, height=input$parm_som_figh/75)
    
    phall <- pData(scData()) 
    som   <- data.som()
    
    path <- input$parm_som_pathway
    ph   <- input$parm_som_ph
    
    codes <- som$codes
    if( is.list(codes))
      codes <- codes[[1]]
    end
    
    op <- par(mfrow = c(2, 2))
    plot(som, main = "TIPS SOM") 
    plot(som, type = "property", property = codes[,paste('Score_',path,sep='')], main=path) #allows emphasis of the change in a given pathway
    if(ph == '[NONE]'){
      plot(som, type = "mapping", pchs = 20, main = "Density")
    }else{
      ft         <- as.factor(phall[,ph])
      levels(ft) <- hue_pal()(length(levels(ft)))  #so this replaces N factors with N default ggplot colors
      
      plot(som, type = "mapping", pchs = 20, main = ph, col = as.character(ft)) #this now maps the metadata onto the points
    }
    plot(som, type = "count", main = "Node Counts", palette.name = viridis) #heatmap of the number of cells in each cluster
    
    par(op) 
    
    dev.off()
  },
  contentType =  NA
)


output$ui_som <- renderUI({
  

  div(fluidRow( 
    column(width = 3,  wellPanel(   numericInput("parm_som_gridy", "Row numbers:", 4, min = 2,step = 1),
                                    numericInput("parm_som_gridx", "Column numbers:", 6, min = 2,step = 1),
                                    radioButtons("parm_som_type", "Selection type:",
                                                 c("Top N" = "top",
                                                   "User selected" = "usr"),inline=T),
                                    conditionalPanel(
                                      condition = "input.parm_som_type == 'top'",
                                      numericInput("parm_som_pathnum", NULL, min(10,length(PseudotimeListAll()$gene_set)), min = 2,max=length(PseudotimeListAll()$gene_set),step = 1)),
                                    
                                    conditionalPanel(
                                      condition = "input.parm_som_type == 'usr'",
                                      selectInput("parm_som_pathselect", NULL,
                                                  names(PseudotimeListAll()$gene_set), multiple=T)),
                                    selectInput("parm_som_pathway", "Select a gene set:",
                                                names(PseudotimeListAll()$gene_set),selected =names(PseudotimeListAll()$gene_set)[1], multiple=F),
                                    selectInput("parm_som_ph", "Select ph:",
                                                c('[NONE]',colnames(pData(scData()))),selected ="[NONE]", multiple=F)
#                                    ,actionButton("btn_som", "Compute")
                                    ),
                      wellPanel(  sliderInput("parm_som_figw", "Figure width:", min = 100, max = 1500, value = 600),
                                  sliderInput("parm_som_figh", "Figure  Hight:", min = 100, max = 1000, value = 400)) 
           ),
    column(width = 9,  box(plotOutput("somPlot",width = "100%", height = "100%"), 
                           width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                           downloadButton('parm_som_dl_plot', 'Download figure(pdf)')) 
  ))
  
})


#############################  switch  ##################################
##  https://github.com/kieranrcampbell/switchde


#########   Hist

output$switchPlotHist <- renderPlot({
  
  if(is.null( data.switch())){ return(NULL) }
  
  mx_sde <- data.switch()

  h <- hist(mx_sde$t0, plot = FALSE)
  plot(h, xlab = "Time", ylab = "Frequency",
       main = "", col = "white")
  
  
},height=function(){input$parm_switchHist_figh},width = function(){input$parm_switchHist_figw})


output$switchHist_dl_plot <- downloadHandler(
  filename = function() { 
    'switchHist_plot.pdf'
  },
  content  = function(file) {
    pdf(file,width = input$parm_switchHist_figw/75, height=input$parm_switchHist_figh/75)
    
    mx_sde <- data.switch()
    h <- hist(mx_sde$t0, plot = FALSE)
    plot(h, xlab = "Time", ylab = "Frequency", main = "", col = "white")
    
    dev.off()
  },
  contentType =  NA
)

output$switchHist_dl_Data <- downloadHandler(
  filename = function() { 
    'data.rds'
  },
  content  = function(file) {
    saveRDS(data.switch(),file = file )
  },
  contentType =  "rds"
)



output$ui_switchHist <- renderUI({

  div(fluidRow( 
    column(width = 3, wellPanel(numericInput("parm_switchHist_lower_threshold", "lower_threshold:", 0.01),
                                numericInput("parm_switchHist_nonzero_threshold","Threshold ratio of nonzero expression (%):", 10, min=0, max=80),
                                numericInput("parm_switchHist_qcutoff", "q value cutoff:", 0.05),
                                numericInput("parm_switchHist_kcutoff", "k cutoff:", 0.5,min = 0),
                                numericInput("parm_switchHist_ucutoff", "u cutoff:", 0,min = 0)), 
                      wellPanel(  sliderInput("parm_switchHist_figw", "Figure width:", min = 100, max = 1500, value = 600),
                       sliderInput("parm_switchHist_figh", "Figure  Hight:", min = 100, max = 1000, value = 400)) 
    ),
    column(width = 9,  box(plotOutput("switchPlotHist",width = "100%", height = "100%"), 
                           width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                       downloadButton('switchHist_dl_Data', 'Download data (.rds)'),
                       downloadButton('switchHist_dl_plot', 'Download figure (pdf)')) 
  ))
  
})



#########   Scater

output$switchPlotScater <- renderPlot({
  
  if(is.null( data.switch())){ return(NULL) }
  
  sdes <- pathway.switch()
  
  N    <- sapply(sdes, function(s) dim(s)[1])    #gene number
  t0   <- sapply(sdes, function(s) median(s$t0))    #median timepoint for switch
  k    <- sapply(sdes, function(s) median(abs(s$k))) #median switch intensity
  
  df <- data.frame(t0=t0, Pathway = names(sdes), Intensity=k, N=N)
  ggplot(df,aes(x = t0, y = Pathway, size = N, color = Intensity))+geom_point()
  
  
},height=function(){input$parm_switchScater_figh},width = function(){input$parm_switchScater_figw})





output$switchScater_dl_Data <- downloadHandler(
  filename = function() { 
    'data.csv'
  },
  content  = function(file) {
    sdes <- pathway.switch()
    
    N    <- sapply(sdes, function(s) dim(s)[1])    #gene number
    t0   <- sapply(sdes, function(s) median(s$t0))    #median timepoint for switch
    k    <- sapply(sdes, function(s) median(abs(s$k))) #median switch intensity
    df   <- data.frame(Pathway = names(sdes),Median_t0=t0,  Median_Intensity=k, Gene_Number=N,stringsAsFactors = F)
    write.csv(df, file = file, row.names = F, quote =F)
  },
  contentType =  "text/csv"
)

output$switchScater_dl_plot <- downloadHandler(
  filename = function() { 
    'switchScater_plot.pdf'
  },
  content  = function(file) {
    sdes <- pathway.switch()
    
    N    <- sapply(sdes, function(s) dim(s)[1])    #gene number
    t0   <- sapply(sdes, function(s) median(s$t0))    #median timepoint for switch
    k    <- sapply(sdes, function(s) median(abs(s$k))) #median switch intensity
    
    df <- data.frame(t0=t0, Pathway = names(sdes), Intensity=k, N=N)
    ggplot(df,aes(x = t0, y = Pathway, size = N, color = Intensity))+geom_point()
    ggsave(file, device = "pdf",width= input$parm_GeneCorPlot_figw/75, height=input$parm_GeneCorPlot_figh/75, units="in")
  },
  contentType =  NA
)




output$ui_switchScater <- renderUI({
  
  
  div(fluidRow( 
    column(width = 3,  wellPanel(  numericInput("parm_switchScater_ncutoff", "gene number cutoff:", 20)), 
                       wellPanel(  sliderInput("parm_switchScater_figw", "Figure width:", min = 100, max = 1500, value = 600),
                                   sliderInput("parm_switchScater_figh", "Figure  Hight:", min = 100, max = 1000, value = 400)) 
    ),
    column(width = 9,  box(plotOutput("switchPlotScater",width = "100%", height = "100%"), 
                           width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                      downloadButton('switchScater_dl_Data', 'Download data (.csv)'),
                      downloadButton('switchScater_dl_plot', 'Download figure (pdf)')) 
  ))
  
})


#########   Line

output$switchPlotLine <- renderPlot({
  
  if(is.null( data.switch())){ return(NULL) }
  if(is.null( input$parm_switchLine_genes)){ return(NULL) }
  
  ptime  <- pData(phenoData(scData()))$Pseudotime
 
  sde1  <- pathway.switch()[[input$parm_switchLine_pathway]]
  genes <- input$parm_switchLine_genes
  mx.new<- gdata_expr()[genes,,drop=F]
  dm    <- dim(mx.new)
  v     <- as.vector(t(mx.new))
  dd    <- data.frame(Expression = v, Pseudotime = rep(ptime,dm[1]),gene = rep(genes,each= dm[2] ))
  min_max <- range(ptime)
  data.line <- NULL
  
  for(gene in genes){
    df <- sigmoid(ptime,params = extract_pars(sde1, gene))
    data.line <- rbind(data.line,data.frame(df,gene = rep(gene, dim(df)[1])))
  }
  
  p1 <- ggplot() + theme_bw() +
    geom_point(data=dd, aes(x = Pseudotime, y = Expression,color = gene,fill=gene), shape = 21,size= input$parm_switchLine_pz) +  
    geom_line(data=data.line, aes(x = x, y = y,color = gene),size= input$parm_switchLine_linew)
  
  
  p2 <- ggplot(sde1, aes(x=t0, y=..density..)) +
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) + geom_density() 
  
  
  
  if(input$parm_switchLine_groupby != '[NONE]'){
    idx <- order(pData(scData())[,input$parm_switchLine_groupby])
    mx.new <- mx.new[,idx,drop=F]
  }
  
  x <- colnames(mx.new) 
  y <- rownames(mx.new) 
  data <- expand.grid(Cell=x, Gene=y)
  data$Expression <- as.vector(t(mx.new))
  p3 <- ggplot(data, aes(x=Cell, y=Gene, fill= Expression)) + geom_tile(colour = "grey50") +
    scale_fill_gradient(low="white", high="red") + theme(axis.text.x = element_blank())
  
#  p3 <- pheatmap(log2(mx.new + 1), cluster_rows = F, cluster_cols = F, show_colnames = F)
  
  
  p12 <- plot_grid(p1,p2,nrow = 1, ncol=2,labels = c("A", "B"))
  plot_grid(p12,p3, nrow = 2, ncol=1, labels = c("", "C"))
  
},height=function(){input$parm_switchLine_figh},width = function(){input$parm_switchLine_figw})


output$switchPlotLine_dl_plot <- downloadHandler(
  filename = function() { 
    'switch_plot.pdf'
  },
  content  = function(file) {
    if(is.null( data.switch())){ return(NULL) }
    if(is.null( input$parm_switchLine_genes)){ return(NULL) }
    
    ptime  <- pData(phenoData(scData()))$Pseudotime
    
    sde1  <- pathway.switch()[[input$parm_switchLine_pathway]]
    genes <- input$parm_switchLine_genes
    mx.new<- gdata_expr()[genes,,drop=F]
    dm    <- dim(mx.new)
    v     <- as.vector(t(mx.new))
    dd    <- data.frame(Expression = v, Pseudotime = rep(ptime,dm[1]),gene = rep(genes,each= dm[2] ))
    min_max <- range(ptime)
    data.line <- NULL
    
    for(gene in genes){
      df <- sigmoid(ptime,params = extract_pars(sde1, gene))
      data.line <- rbind(data.line,data.frame(df,gene = rep(gene, dim(df)[1])))
    }
    
    p1 <- ggplot() + theme_bw() +
      geom_point(data=dd, aes(x = Pseudotime, y = Expression,color = gene,fill=gene), shape = 21,size= input$parm_switchLine_pz) +  
      geom_line(data=data.line, aes(x = x, y = y,color = gene),size= input$parm_switchLine_linew)
    
    
    p2 <- ggplot(sde1, aes(x=t0, y=..density..)) +
      geom_histogram(fill="cornsilk", colour="grey60", size=.2) + geom_density() 
    
    
    
    if(input$parm_switchLine_groupby != '[NONE]'){
      idx <- order(pData(scData())[,input$parm_switchLine_groupby])
      mx.new <- mx.new[,idx,drop=F]
    }
    
    x <- colnames(mx.new) 
    y <- rownames(mx.new) 
    data <- expand.grid(Cell=x, Gene=y)
    data$Expression <- as.vector(t(mx.new))
    p3 <- ggplot(data, aes(x=Cell, y=Gene, fill= Expression)) + geom_tile(colour = "grey50") +
      scale_fill_gradient(low="white", high="red") + theme(axis.text.x = element_blank())
    
    #  p3 <- pheatmap(log2(mx.new + 1), cluster_rows = F, cluster_cols = F, show_colnames = F)
    
    
    p12 <- plot_grid(p1,p2,nrow = 1, ncol=2,labels = c("A", "B"))
    plot_grid(p12,p3, nrow = 2, ncol=1, labels = c("", "C"))
    ggsave(file, device = "pdf",width= input$parm_GeneCorPlot_figw/75, height=input$parm_GeneCorPlot_figh/75, units="in")
  },
  contentType =  NA
)



output$ui_switchLine_param <- renderUI({
  sde <- pathway.switch()[[input$parm_switchLine_pathway]]
  genes <- rownames(sde)
  
  
  idx   <- order(sde$qval, decreasing = F)
  gname <- paste(genes,'(q:',format(sde$qval,digits=3),')',sep='')
  names(genes) <- gname
  
  genes <- genes[idx]
  
  selectInput("parm_switchLine_genes", "Select genes:", genes, multiple=T)
})






output$ui_switchLine <- renderUI({
  
  
  div(fluidRow( 
    column(width = 3,  
                       wellPanel(  selectInput("parm_switchLine_pathway", "Select pathway:", names(pathway.switch()),multiple=F),
                                   uiOutput("ui_switchLine_param"),
                                   
                                   selectizeInput(inputId = "parm_switchLine_groupby", label = "Grouped by:", 
                                                  choices = c('[NONE]',colnames(pData(scData()))), 
                                                  selected = NULL, multiple = F)
                                ),  
           
                       wellPanel(  sliderInput("parm_switchLine_pz", "Point size:", min = 1, max = 10, value = 1),
                                   sliderInput("parm_switchLine_linew", "Line width:", min = 0.1, max = 2, value = 0.5),
                                   sliderInput("parm_switchLine_figw", "Figure width:", min = 100, max = 1500, value = 800),
                                   sliderInput("parm_switchLine_figh", "Figure  Hight:", min = 100, max = 1000, value = 600)) 
    ),
    column(width = 9,  box(plotOutput("switchPlotLine",width = "100%", height = "100%"), 
                           width = NULL,  collapsible = F,title = " ", status = "primary", solidHeader = TRUE),
                      downloadButton('switchPlotLine_dl_plot', 'Download figure (pdf)') ) 
  ))
  
})





