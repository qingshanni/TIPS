
output$ui_intro <- renderUI({
  div(
    h2("Introduction"),
    wellPanel(
    p("Recent advances in bioinformatics analyses have led to the development of novel tools enabling
      the capture and trajectory mapping of single-cell RNA sequencing (scRNAseq) data. However,
      there is a lack of methods to assess the contributions of biological pathways and transcription
      factors to an overall developmental trajectory mapped from scRNAseq data. In this manuscript,
      we present a simplified approach for trajectory inference of pathway significance (TIPS) that
      leverages existing knowledgebases of functional pathways and other gene lists to provide further
      mechanistic insights into a biological process. TIPS identifies key pathways which contribute to
      a process of interest, as well as the individual genes that best reflect these changes. TIPS also
      provides insight into the relative timing of pathway changes, as well as a suite of visualizations
      to enable simplified data interpretation of scRNAseq libraries generated using a wide range of
      techniques. TIPS may serve as a useful tool to help biologists perform deeper
      functional analyses and visualization of their single-cell data. A schematic overview of TIPS is 
      showed in the following figure.")),
    br(),br(),br(),br(),
    img(src='flowchart.jpg', style="display: block; margin-left: auto; margin-right: auto")
  )
})
