

library(shiny)
library(rclipboard)
library(shinydashboard)

library(ggplot2)
# library(data.table)
# library(ggrepel)
# library(viridis)
# library(ggnewscale)
# library(RColorBrewer)
# library(grid)
# library(gridExtra) 
library(dplyr)
# library(ggrastr)
# library(ggpubr)



library(Seurat)
library(Matrix)
# library(DT)
# library(patchwork)


# library(ggsankey) #devtools::install_github("davidsjoberg/ggsankey")
# library(highcharter)

# library("BiocParallel")
# register(MulticoreParam(6))

# col_vector <- scCustFx::ColorTheme()$col_vector
col_vector = c("#7FC97F", "#38170B", "#BEAED4", "#BF1B0B", "#FFC465", "#386CB0", 
               "#66ADE5", "#F0027F", "#252A52", "#BF5B17", "#999999", "#666666", 
               "#E69F00", "#1B9E77", "#56B4E9", "#D95F02", "#009E73", "#7570B3", 
               "#F0E442", "#E7298A", "#0072B2", "#66A61E", "#D55E00", "#E6AB02",
               "#CC79A7", "#A6761D", "#e6194b", "#666666", "#3cb44b", "#A6CEE3", 
               "#ffe119", "#1F78B4", "#4363d8", "#B2DF8A", "#f58231", "#33A02C",
               "#911eb4", "#FB9A99", "#46f0f0", "#E31A1C", "#FDBF6F", "#FF7F00",
               "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", 
               "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", 
               "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", 
               "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
               "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
               "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
               "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
               "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
               "#CCEBC5", "#FFED6F")

# ggtheme = scCustFx:::theme_simple


# library(scCustFx)
# library(CellMembrane)

library(biomaRt)




pathi = getwd()


# Define UI for application that draws a histogram
ui <- dashboardPage(skin="red",
                    dashboardHeader(title = "MoSpaTA v0.1A"),
                    #https://rstudio.github.io/shinydashboard/appearance.html#icons
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Introduction", tabName = "MainDash", icon = icon("dashboard"), selected = T),
                        menuItem("Load Data", tabName = "LoadData", icon = icon("save")),
                        menuItem("Gene Expression", tabName = "GeneExpr", icon = icon("dna")),
                        menuItem("MoDSTA (scRNA-seq) MetaData", tabName = "MoDSTAmeta", icon = icon("wrench")),
                        menuItem("STseqMT1 (Spatial) MetaData", tabName = "STseqMT1meta", icon = icon("wrench")),
                        menuItem("SDA Component Browser", tabName = "SDABrowser", icon = icon("dna")),
                        
                        # menuItem("Save Out", tabName = "SaveOut", icon = icon("save")),
                        menuItem("@eisamahyari", icon = icon("heart"), 
                                 href = "https://eisascience.github.io")
                      )
                    ),
                    
                    dashboardBody(
                      # useShinyjs(),
                      tags$head(
                        tags$style(HTML("
                                        .content-wrapper {
                                        background-color: black !important;
                                        }
                                        .main-sidebar {
                                        background-color: black !important;
                                        }
                                        .multicol .shiny-options-group{
                                        -webkit-column-count: 5; /* Chrome, Safari, Opera */
                                        -moz-column-count: 5;    /* Firefox */
                                        column-count: 5;
                                        -moz-column-fill: balanced;
                                        -column-fill: balanced;
                                        }
                                        .checkbox{
                                        margin-top: 0px !important;
                                        -webkit-margin-after: 0px !important; 
                                        }
                                        "))),
                      tabItems(
                        
                        # Tab Items ------------
                        
                        ## Main tab ------------
                        tabItem(tabName = "MainDash",
                                h2("Main Dashboard"),
                                fluidRow(
                                  # valueBoxOutput("InfoBox_Main", width = 6),
                                  
                                  # box(,
                                  #     width = 5, background = "olive"
                                  #     
                                  # )
                                )),
                        
                        ## Load Data tab ------------
                        tabItem(tabName = "LoadData",
                                h2("Load In Data"),
                                
                                fluidRow(
                                  
                                  valueBoxOutput("InfoBox_Load", width = 6),
                                  
                                  box(actionButton("MoDSTA", "Mouse Deveopmental scRNASeq Testis Atlas"),
                                      actionButton("STseqMouse1", "Mouse STseq Testis 1"),
                                      actionButton("SDAres", "SDA on MoDSTA"),
                                      actionButton("BiomaRtLoad", "Load BiomaRt Data"),
                                      actionButton("AllInputs", "** Load All **"),
                                      width = 10, background = "olive")
                                  
                                  )),
                        
                        ## Gene Expression tab ------------
                        tabItem(tabName = "GeneExpr",
                                h2("Gene Expression"),
                                fluidRow(
                                  box(
                                    title = "Inputs", status = "warning", solidHeader = TRUE,
                                    "Multiple formatting of gene sets accepted", 
                                    br(), "List can be seperated by comma e.g. from ", 
                                    br(), "   or spaces e.g. from Excel", 
                                    br(), "Also, single or double quotes or not",
                                    #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                                    textInput("GeneSet", "A set of genes", "'Prm1', 'Tnp1'"),
                                    width = 10
                                  ),
                                  
                                  # box(
                                  #   title = "MoDSTA UMAP Unintegrated", status = "primary", solidHeader = TRUE,
                                  #   collapsible = TRUE,
                                  #   plotOutput("GeneExpr_MoDSTA_UMAP"),
                                  #   width = 5
                                  # ),
                                  
                                  
                                  # box(
                                  #   title = "Celltypes MoDSTA", status = "primary", solidHeader = TRUE,
                                  #   collapsible = TRUE,
                                  #   plotOutput("CellType_MoDSTA_UMAPintg"),
                                  #   width = 10
                                  # ),
                                  
                                  box(
                                    title = "Expression on The Mouse Developmental Single-cell Testis Atlas (MoDSTA) UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_MoDSTA_UMAPintg"),
                                    width = 10
                                  ),
                                  
                                  box(
                                    title = "Expression on the Stereo-seq Mouse Testis Sample 1 UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_STseqMT1_UMAP"),
                                    width = 10
                                  ),
                                  box(
                                    title = "Spatial Expression on the Stereo-seq Mouse Testis Sample 1", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_STseqMT1_Spatial"),
                                    width = 10
                                  )
                                  
                                  
                                  
                                )),
                        
                        
                        ## MoDSTA Meta tab ------------
                        tabItem(tabName = "MoDSTAmeta",
                                h2("Meta Data Associated to The Mouse Developmental Single-cell Testis Atlas"),
                                fluidRow(
                                  
                                  
                                  box(
                                    title = "Celltypes MoDSTA", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("CellType_MoDSTA_UMAPintg"),
                                    width = 10
                                  )
                                  
                                )),
                        ## STseqMT1 Meta tab ------------
                        tabItem(tabName = "STseqMT1meta",
                                h2("Meta Data Associated to Spatial (Stereo-seq) of Adult Mouse Testis Sample 1"),
                                fluidRow(
                                  
                                  
                                  box(
                                    title = "Unuspervised Clustering on Spatial (Stereo-seq) of Adult Mouse Testis Sample 1 UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("UnSupClusters_STseqMT1_UMAP"),
                                    width = 10
                                  )
                                  
                                )),
                        
                        
                        
                        ## SDABrowser tab ------------
                        tabItem(tabName = "SDABrowser",
                                h2("SDA Component Analysis"),
                                fluidRow(
                                  
                                  valueBoxOutput("InfoBox_SDABrowser", width = 6),
                                  
                                  box(
                                    title = "SDA Run and Component Selector", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    uiOutput("Select.SDArun"),
                                    
                                    width = 5
                                  ),
                                  box(
                                    title = "Projection", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                            
                                    actionButton("Apply2MoDSTA", "Apply to MoDSTA"),
                                    actionButton("Apply2STseqMT1", "Apply to TseqMT1"),
                                    
                                    
                                    width = 5
                                  ),
                                  box(
                                    title = "Navigation", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    uiOutput("Select.SDAcomponentN"),
                                    
                                    actionButton("PreviousComp", "Prev"),
                                    actionButton("NextComp", "Next"),
                                    
                                    
                                    width = 5
                                  ),
                                  box(
                                    title = "Score Projection on MoDSTA UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAloadings_coord"),
                                    plotOutput("SDAscore_STseqMT1_Spatial"),
                                    plotOutput("SDAscore_MoDSTA_UMAPintg"),
                                    plotOutput("SDAscore_STseqMT1_UMAP"),
                                    
                                    
                                    
                                    width = 10
                                  )
                                  
                                ))
                        
                        
                        
                      ) #end of tabItems
                    ) #end of body
) #end UI

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ## environment defaults
  envv=reactiveValues(y=NULL)
  
  source("app_OE_load.R",local = TRUE)
  source("app_OE_SDA.R",local = TRUE)
  source("app_InfoBox.R",local = TRUE)
  source("app_Fxs.R",local = TRUE)
  source("app_Figs_UMAP.R",local = TRUE)
  

  
  ### SDA local folder
  output$Select.SDArun <-
    renderUI(expr = selectInput(inputId = 'sda.run',
                                label = 'SDA Run Name',
                                choices = names(envv$SDARedDataLS$loadings) ))
  
  
  # Reactive expression to track the selected run name
  selected_SDArun <- reactive({
    input$sda.run
  })

  ### SDA local folder
  # output$Select.SDAcomponentN <-
  #   renderUI(expr = selectInput(inputId = 'sda.comp.N',
  #                               label = 'SDA Component Number',
  #                               choices = paste0("Comp ", 1: nrow(envv$SDARedDataLS$loadings[[output$Select.SDArun]]$loadings)) ))

  output$Select.SDAcomponentN <- renderUI({
    if (!is.null(selected_SDArun()) && !is.null(envv$SDARedDataLS$loadings[[selected_SDArun()]])) {
      n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
      selectInput(inputId = 'sda.comp.N',
                  label = 'SDA Component Number',
                  choices = paste0("Comp_", 1:n_comps))
    } else {
      selectInput(inputId = 'sda.comp.N',
                  label = 'SDA Component Number',
                  choices = character(0))
    }
  })
 
  # Reactive expression to track the selected run name
  selected_compN <- reactive({
    input$sda.comp.N
  })
  
  # if(length(selected_compN)>1) {
  #   envv$CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
  #   envv$CompNname = paste0("sda.", selected_SDArun, ".V", CompN)
  # }
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
