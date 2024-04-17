

library(shiny)
library(rclipboard)
library(shinydashboard)

library(knitr)

library(ggplot2)
# library(data.table)
library(ggrepel)
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


VersionID = "0.3.2A"

pathi = getwd()


# Define UI for application that draws a histogram
ui <- dashboardPage(skin="yellow",
                    dashboardHeader(title = paste0("MoSpaTA v", VersionID)),
                    #https://rstudio.github.io/shinydashboard/appearance.html#icons
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Introduction", tabName = "MainDash", icon = icon("dashboard"), selected = T),
                        menuItem("Load Data", tabName = "LoadData", icon = icon("save")),
                        menuItem("Gene Expression", tabName = "GeneExpr", icon = icon("dna")),
                        menuItem("MoDSTA (scRNA-seq) MetaData", tabName = "MoDSTAmeta", icon = icon("wrench")),
                        menuItem("STseqMT1 (Spatial) MetaData", tabName = "STseqMT1meta", icon = icon("wrench")),
                        menuItem("SDA Component Browser", tabName = "SDABrowser", icon = icon("dna")),
                        menuItem("SDA Gene Search", tabName = "SDAgeneSearch", icon = icon("dna")),
                        
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
                                        background-color: grey2 !important;
                                        }
                                        markdown-content h1, .markdown-content p {
                                          font-family: Arial, sans-serif;
                                        }
                                        .main-sidebar {
                                        background-color: grey2 !important;
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
                                        .markdown-content { 
                                            width: 100%; 
                                            max-width: none; 
                                          }
                                          .markdown-content * { 
                                            max-width: 100%;
                                          }
                                        "))),
                      tabItems(
                        
                        # Tab Items ------------
                        
                        ## Main tab ------------
                        tabItem(tabName = "MainDash",
                                h2("Main Dashboard"),
                                fluidRow(
                               # box( uiOutput('MainPageHTML'),
                               #      width = 10)
                                  #(Version 0.3A)
                                  column(12,
                                         tags$h1(paste0("MoSpaTA Version ", VersionID )),
                                         tags$p("Welcome to the Mouse Spatial Omic Testis Atlas (MoSpaTA), an interactive tool for exploring spatial and single-cell transcriptomic data of mouse testis."),
                                         tags$p("This early alpha release of MoSpaTA aims to foster discussion and interest in developing a unified spatial omic atlas of the testis."),
                                         tags$p(tags$strong("Conrad Lab 2024")),
                                         tags$p(tags$em("Eisa Mahyari (@eisamahyari)")),
                                         tags$a(href="https://github.com/eisascience/MoSpaTA", "See: github.com/eisascience/MoSpaTA"),
                                         tags$h2("Underlying Data"),
                                         tags$p("MoSpaTA integrates the following key resources:"),
                                         tags$ul(
                                           tags$li(tags$strong("STseq MT1:"), " This dataset represents a spatial transcriptomic atlas of a healthy adult wild-type (WT) mouse. The cell-level data are derived from algorithms that define cell boundaries using a DAPI image, covering a total of approximately 35,000 cells to optimize boundary selection and enhance the performance of this interactive atlas."),
                                           tags$li(tags$strong("MoDSTA (Mouse Developmental Single-cell Testis Atlas):"), " Comprising 177,891 testis cells from 20 publicly available 10X Genomics samples, this dataset includes cells from 2 embryonic, 4 postnatal, and 14 WT adult samples. For performance optimization, MoDSTA has been randomly downsampled to 30,000 cells. Future versions will include integrated analyses combining STseq MT1 and MoDSTA datasets."),
                                           tags$li(tags$strong("Machine Learning on MoDSTA:"), " Using SDA (soft clustering), our machine-learning models learn transcriptional patterns within high-dimensional transcriptomic spaces. Several SDA models have been trained on subsets of MoDSTA to enable specialized analysis.")
                                         ),
                                         tags$h2("User Manual"),
                                         tags$p("Navigate MoSpaTA using the left-side navigation bar to access various features."),
                                         tags$h3("Load Data Tab"),
                                         tags$p("Utilize the five buttons available for data management, with the 'Load All' button recommended for general use. If data loading issues arise, try loading each dataset individually for troubleshooting."),
                                         tags$h3("Gene Expression Tab"),
                                         tags$p("Explore gene expression by typing the names of single or multiple genes. View gene expression maps derived from UMAP analyses of MoDSTA and STseq MT1 datasets, along with their spatial representations in the testis."),
                                         tags$h3("MoDSTA MetaData Tab"),
                                         tags$p("Access visual representations and metadata related to cell type labels from the MoDSTA dataset."),
                                         tags$h3("STseq MT1 MetaData Tab"),
                                         tags$p("Explore metadata and related visual content for the STseq MT1 dataset."),
                                         tags$h3("SDA Component Browser"),
                                         tags$p("After loading the data (see Load Data tab), you can explore various SDA runs. There are several SDA runs that were performed on MoDSTA; some are runs on subsets of MoDSTA, others are simply different seeds of the combo objects. The aim is to have machine-learning train on every part of this data to find signatures of potential interest. These signatures (also called components) are effectively multi-geneic weighted signatures that fine-tune to a signature; providing surgical dissection of the signature. "),
                                         tags$p("First, an SDA run needs to be selected, e.g., sda_FullCombo."),
                                         tags$p("Second, click on the two Apply Buttons, to project the SDA run in to each dataset. This produces a score, for each cell of each data for each component in each selected run. A run may have 20-100 components, which were selected a priori heuristically. "),
                                         tags$p("Third, navigate the components either using the drop down or the next/prev buttons. "),
                                         tags$p("The genes and their weights are shown in the first figure. The top weighted genes are shown. The weights have direction, so the positive genes score towards the positive direction and the negative loaded genes, score towards the negative direction. The key point is the grouping of genes in each direction, as a 'soft' cluster of transcriptional correlation between those identified genes. In other words, the component has found a transcriptional signature between them. Next, How this scores cells in MoDSTA and STseqMT1 is visualized."),
                                         tags$h3("SDA Runs Table:"),
                                         tags$table(class = "table table-striped",
                                                    tags$thead(
                                                      tags$tr(
                                                        tags$th("ID"),
                                                        tags$th("Description"),
                                                        tags$th("Name")
                                                      )
                                                    ),
                                                    tags$tbody(
                                                      tags$tr(tags$td("569507"), tags$td("Other Minor Clusters of the Testis"), tags$td("Minor Clusters")),
                                                      tags$tr(tags$td("569509"), tags$td("Leydig Cell Clusters"), tags$td("Leydig Clusters")),
                                                      tags$tr(tags$td("569516"), tags$td("Germ Cell Clusters"), tags$td("GC Clusters")),
                                                      tags$tr(tags$td("569520"), tags$td("Sertoli Cell Clusters"), tags$td("SC Clusters")),
                                                      tags$tr(tags$td("572405"), tags$td("WT Adults only Combo"), tags$td("WT8Week Full Combo")),
                                                      tags$tr(tags$td("573391"), tags$td("Juviniles and Prenatal only Combo"), tags$td("JuvPreNatal Combo")),
                                                      tags$tr(tags$td("574419"), tags$td("Full Combo"), tags$td("Full Combo DS-seed1")),
                                                      tags$tr(tags$td("574426"), tags$td("Full Combo"), tags$td("Full Combo DS-seed2")),
                                                      tags$tr(tags$td("574430"), tags$td("Full Combo"), tags$td("Full Combo DS-seed3")),
                                                      tags$tr(tags$td("WTCombo"), tags$td("Adults Only"), tags$td("Adults Combo DS-seed1")),
                                                      tags$tr(tags$td("WTCombo2"), tags$td("Adults Only"), tags$td("Adults Combo DS-seed2")),
                                                      tags$tr(tags$td("568747"), tags$td("Full Combo"), tags$td("Full Combo DS-seed4"))
                                                    )),
                                         tags$h3(""),
                                         tags$p(".")
                                  )
                                )
                                ),
                        
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
                                      width = 10, background = "black")
                                  
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
                                  
                                  box(
                                    title = "Threshold Expression, 0 to 100%", status = "primary", solidHeader = TRUE,
                                    #   collapsible = TRUE,
                                    sliderInput("ExprThresh", "% of Expression: (default 0, no cutting)",
                                                min = 0, max = 1, value = 0, step = 0.1),
                                    width = 5
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
                                    title = "Spatial Expression on the Stereo-seq Mouse Testis Sample 1", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExpr_STseqMT1_Spatial"),
                                    width = 10
                                  ), 
                                  
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
                                    plotOutput("DatasetName_MoDSTA_UMAPintg"),
                                    plotOutput("DatasetGrouping_MoDSTA_UMAPintg"),
                                    plotOutput("UnsupClusters_MoDSTA_UMAPintg"),
                                    plotOutput("Breed_MoDSTA_UMAPintg"),
                                    
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
                                    plotOutput("UnSupClusters_SpatialFacet_STseqMT1_UMAP"),
                                    
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
                                    title = "Threshold SDA Score, 0 to 100%", status = "primary", solidHeader = TRUE,
                                    #   collapsible = TRUE,
                                    sliderInput("SDAThresh", "% of SDA Score: (default 0, no cutting)",
                                                min = 0, max = 1, value = 0, step = 0.1),
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
                                  
                                )),
                        
                        ## SDA gene Search  tab ------------
                        tabItem(tabName = "SDAgeneSearch",
                                h2("Search SDA components top loaded genes for your gene"),
                                fluidRow(
                                  
                                  box(
                                    title = "Single Gene Input", status = "warning", solidHeader = TRUE,
                                    textInput("Gene2SearchSDA", "A single gene", "Prm1"),
                                    width = 5
                                  ),
                                  box(
                                    title = "Positive Hits", status = "info", solidHeader = TRUE,
                                    textOutput("displayPosHits"),  # UI element to display positive hits
                                    width = 6
                                  ),
                                  box(
                                    title = "Negative Hits", status = "info", solidHeader = TRUE,
                                    textOutput("displayNegHits"),  # UI element to display negative hits
                                    width = 6
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
  source("app_Figs_SDA.R",local = TRUE)
  
  # output$MainPageHTML <- renderUI({
  #   # includeHTML('MainPage.html')
  #   # print(getwd())
  #   # tags$iframe(seamless="seamless", src= "./MainPage.html", width=800, height=800)
  #   # tags$iframe(seamless = "seamless", src = "./inst/app/MainPage.html", width = "100%", height = "800")
  #   
  #   HTML(markdown::markdownToHTML(knit('MainPage.md', quiet = TRUE)))
  #   
  # })
  # 
  
  ### SDA local folder
  output$Select.SDArun <-
    renderUI(expr = selectInput(inputId = 'sda.run',
                                label = 'SDA Run Name',
                                choices = names(envv$SDARedDataLS$loadings) ))
  
  # observe({
  #   print(input$sda.run)
  # })
  
  # Reactive expression to track the selected run name
  selected_SDArun <- reactive({
    # req(input$sda.run)
    # print(input$sda.run)
    input$sda.run
  })
 

  ### SDA local folder
  # output$Select.SDAcomponentN <-
  #   renderUI(expr = selectInput(inputId = 'sda.comp.N',
  #                               label = 'SDA Component Number',
  #                               choices = paste0("Comp ", 1: nrow(envv$SDARedDataLS$loadings[[output$Select.SDArun]]$loadings)) ))

  output$Select.SDAcomponentN <- renderUI({
    
    if (!is.null(selected_SDArun())) {
      n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
      selectInput(inputId = 'sda.comp.N',
                  label = 'SDA Component Number',
                  choices = paste0("Comp_", 1:n_comps))
    } else {
      selectInput(inputId = 'sda.comp.N',
                  label = 'SDA Component Number',
                  choices = character(0))
    }
    # req(selected_SDArun())  
    # n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
    # selectInput(inputId = 'sda.comp.N',
    #             label = 'SDA Component Number',
    #             choices = paste0("Comp_", 1:n_comps))
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
