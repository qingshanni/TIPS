output$ui_intro_text <- renderText({
  paste('  Real-time quantitative polymerase-chain-reaction (qPCR) is a standard technique', 
        'in most laboratories used for various applications in basic research. Analysis of qPCR', 
        'data is a crucial part of the entire experiment.\n ',
        '  For the convenience of users, we developed qPCRshiny, a software for qPCR data analysis.',
        'The workflow of qPCR data analysis is as following',
    sep = ' '
  )
  
  
})

output$ui_intro <- renderUI({
  div(
    h2("Introduction"),
    wellPanel(
    verbatimTextOutput("ui_intro_text"),
    img(src='workflow.png', style="display: block; margin-left: auto; margin-right: auto")
    )
  )
})
