library(shinydashboard)

sidebar <-   dashboardSidebar(
  sidebarMenu(
    menuItem("Introduction", tabName = "tab_intro", icon = icon("eye"),selected=TRUE),
    menuItem("Data upload",   tabName = "tab_gdata", icon = icon("dashboard")),
    menuItem("Data quality",   tabName = "tab_filter", icon = icon("dashboard")),
    menuItem("Dimensional reduction",   tabName = "tab_dimreduce", icon = icon("dashboard")),
    menuItem("Trajectory analysis",   tabName = "tab_trajectory", icon = icon("dashboard")),
    menuItem("Pathway correlation analysis",   tabName = "tab_geneset", icon = icon("dashboard")),
    menuItem("Gene correlation analysis",   tabName = "tab_gene_cor", icon = icon("dashboard")),
    menuItem("SOM analysis",   tabName = "tab_som", icon = icon("dashboard")),
    menuItem("Switch analysis",   icon = icon("dashboard"),
             menuSubItem("Switch time distribution",  tabName = "tab_switchHist", icon = icon("dashboard")),
             menuSubItem("Switch t0",  tabName = "tab_switchScater", icon = icon("dashboard")),
             menuSubItem("Switch genes",  tabName = "tab_switchLine", icon = icon("dashboard"))
             )

  ))  

boardpage <-  dashboardBody( 
  tabItems(
#    tabItem(tabName = "tab_intro", includeMarkdown("intro.md")),    
    tabItem(tabName = "tab_intro", uiOutput("ui_intro")),   
    tabItem(tabName = "tab_gdata", uiOutput("ui_gdata")),
    tabItem(tabName = "tab_filter", uiOutput("ui_filter")),
    tabItem(tabName = "tab_dimreduce", uiOutput("ui_dimreduce")),    
    tabItem(tabName = "tab_trajectory", uiOutput("ui_trajectory")),
    tabItem(tabName = "tab_geneset", uiOutput("ui_geneset")),
    tabItem(tabName = "tab_gene_cor", uiOutput("ui_gene_cor")),
    tabItem(tabName = "tab_som", uiOutput("ui_som")),
    tabItem(tabName = "tab_switchHist", uiOutput("ui_switchHist")),
    tabItem(tabName = "tab_switchScater", uiOutput("ui_switchScater")),
    tabItem(tabName = "tab_switchLine", uiOutput("ui_switchLine"))
    
  )
)


dashboardPage(
  dashboardHeader(title = 'TIPS'),
  sidebar,
  boardpage,
  title = 'TIPS analysis'
)

