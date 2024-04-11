

library(shiny)
library(rclipboard)
library(shinydashboard)

library(ggplot2)
library(data.table)
library(ggrepel)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(grid)
library(gridExtra) 
require(dplyr)
library(ggrastr)
library(ggpubr)



library(Seurat)
library(Matrix)
library(DT)
library(patchwork)


# library(ggsankey) #devtools::install_github("davidsjoberg/ggsankey")
# library(highcharter)

library("BiocParallel")
register(MulticoreParam(6))

col_vector <- scCustFx::ColorTheme()$col_vector

ggtheme = scCustFx:::theme_simple


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
                                      actionButton("SDAres", "SDA on MoSTA"),
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
                                  
                                  
                                  box(
                                    title = "Celltypes MoDSTA", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("CellType_MoDSTA_UMAPintg"),
                                    width = 10
                                  ),
                                  
                                  box(
                                    title = "Expr MoDSTA UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_MoDSTA_UMAPintg"),
                                    width = 10
                                  ),
                                  
                                  box(
                                    title = "Expr STSeq MT1 UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_STseqMT1_UMAP"),
                                    width = 10
                                  ),
                                  box(
                                    title = "Expr STSeq MT1 Spatial", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_STseqMT1_Spatial"),
                                    width = 10
                                  ),
                                  
                                  box(
                                    title = "UnsupClus STSeq MT1 UMAP", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("UnSupClusters_STseqMT1_UMAP"),
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
  source("app_InfoBox.R",local = TRUE)
  source("app_Fxs.R",local = TRUE)
  source("app_Figs_UMAP.R",local = TRUE)
  

  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
