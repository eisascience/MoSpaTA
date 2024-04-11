



## Gene Expr MoDSTA Unintegrated UMAP ------

output$GeneExpr_MoDSTA_UMAP <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  if(length(grep(",", GeneSet)) == 0){
    
    if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
      GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
    } else {
      GeneSet <- unlist(strsplit(GeneSet, " "))
    }
    
    
  } else {
    GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
    
  }
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MoDSTA)]
  
  print(GeneSet)
  
  if(is.null(envv$MoDSTA) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load MoDSTA dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      FeaturePlot(envv$MoDSTA , reduction = "umap", features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      envv$MoDSTA$SummedExpr = colSums(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
      
      FeaturePlot(envv$MoDSTA , reduction = "umap", features = "SummedExpr",
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }

    
    
  }
  
})





## Gene Expr MoDSTA Integrated UMAP ------

output$GeneExpr_MoDSTA_UMAPintg <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  if(length(grep(",", GeneSet)) == 0){
    
    if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
      GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
    } else {
      GeneSet <- unlist(strsplit(GeneSet, " "))
    }
    
    
  } else {
    GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
    
  }
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MoDSTA)]
  
  print(GeneSet)
  
  if(is.null(envv$MoDSTA) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load MoDSTA dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      FeaturePlot(envv$MoDSTA , reduction = "umapSup.harmony", features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      envv$MoDSTA$SummedExpr = colSums(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
      
      FeaturePlot(envv$MoDSTA , reduction = "umapSup.harmony", features = "SummedExpr",
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    
  }
  
})


## Gene Expr STseq MT1 UMAP ------

output$GeneExpr_STseqMT1_UMAP <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  if(length(grep(",", GeneSet)) == 0){
    
    if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
      GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
    } else {
      GeneSet <- unlist(strsplit(GeneSet, " "))
    }
    
    
  } else {
    GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
    
  }
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
  
  print(GeneSet)
  
  if(is.null(envv$MSTseqCells1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq MT1 dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      FeaturePlot(envv$MSTseqCells1 , reduction = "umap", features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      envv$MSTseqCells1$SummedExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
      
      FeaturePlot(envv$MSTseqCells1 , reduction = "umap", features = "SummedExpr",
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    
  }
  
})


output$GeneExpr_STseqMT1_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  if(length(grep(",", GeneSet)) == 0){
    
    if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
      GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
    } else {
      GeneSet <- unlist(strsplit(GeneSet, " "))
    }
    
    
  } else {
    GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
    
  }
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
  
  print(GeneSet)
  
  if(is.null(envv$MSTseqCells1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq MT1 dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      SpatialFeaturePlot(envv$MSTseqCells1 , features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      envv$MSTseqCells1$SummedExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
      
      SpatialFeaturePlot(envv$MSTseqCells1 , features = "SummedExpr",
                  # raster = T, raster.dpi = c(800, 800),
                   max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    
  }
  
})




## Celltype on Integrated UMAP ------

output$CellType_MoDSTA_UMAPintg  <- renderPlot({
  
  DimPlot(envv$MoDSTA , reduction = "umapSup.harmony", 
          group.by = "Pheno1",
          repel = T,
          cols = col_vector, label = T,
          raster = T, raster.dpi = c(800, 800))  +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "right",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("CellTypes")
  
  
  
})



## UnsupClusters on STseq MT1 UMAP ------

output$UnSupClusters_STseqMT1_UMAP  <- renderPlot({
  
  DimPlot(envv$MSTseqCells1 , reduction = "umap", 
          group.by = "SCT_snn_res.0.1",
          repel = T,
          cols = col_vector, label = T#,raster = T, raster.dpi = c(800, 800)
          )  +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "right",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("Transcriptional Groups")
  
  
  
})




## SDA Score on Integratred UMAP

output$SDAscore_MoDSTA_UMAPintg <- renderPlot({
  
  
  
  # envv$CompN = as.numeric( strsplit(input$sda.comp.N, "_")[[1]][2])
  # envv$CompNname = paste0("sda.", selected_SDArun, ".V", envv$CompN)
  # CompNname = envv$CompNname
  # 
  # print(CompNname)
  

  selected_SDArun = input$sda.run
  selected_compN = input$sda.comp.N
  

  
  
  if(is.null(envv$MoDSTA) | is.null(envv$commands$SDAproj_MoDSTA)){
    plot(x=0, y=0, main="Load MoDSTA dataset")
  } else {
    
    # # Split the component number from its label
    # compN <- as.numeric(strsplit(selected_compN(), " ")[[1]][2])
    # 
    # # Construct the dynamic variable name
    # compNname <- paste0("sda.", selected_SDArun(), ".V", compN)
    
    CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
    CompNname = paste0("sda.", selected_SDArun, ".V", CompN)
    
    FeaturePlot(envv$MoDSTA , reduction = "umapSup.harmony", features = CompNname, order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
      theme_classic(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      ) &
      ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
    
    
    
    
  }
  
})




## SDA Score on STseqMT1 UMAP

output$SDAscore_STseqMT1_UMAP <- renderPlot({
  
  

  
  selected_SDArun = input$sda.run
  selected_compN = input$sda.comp.N
  
  
  
  
  if(is.null(envv$MoDSTA) | is.null(envv$commands$SDAproj_MSTseqCells1)){
    plot(x=0, y=0, main="Load MSTseq1 Cells dataset")
  } else {
    
    # # Split the component number from its label
    # compN <- as.numeric(strsplit(selected_compN(), " ")[[1]][2])
    # 
    # # Construct the dynamic variable name
    # compNname <- paste0("sda.", selected_SDArun(), ".V", compN)
    
    CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
    CompNname = paste0("sda.", selected_SDArun, ".V", CompN)
    
    FeaturePlot(envv$MSTseqCells1 , reduction = "umap", features = CompNname, order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
      theme_classic(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      )   &
      ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
    
    
    
    
  }
  
})





## SDA Score on STseqMT1 Spatial

output$SDAscore_STseqMT1_Spatial <- renderPlot({
  
  
  
  
  selected_SDArun = input$sda.run
  selected_compN = input$sda.comp.N
  
  
  
  
  if(is.null(envv$MoDSTA) | is.null(envv$commands$SDAproj_MSTseqCells1)){
    plot(x=0, y=0, main="Load MSTseq1 Cells dataset")
  } else {
    
    # # Split the component number from its label
    # compN <- as.numeric(strsplit(selected_compN(), " ")[[1]][2])
    # 
    # # Construct the dynamic variable name
    # compNname <- paste0("sda.", selected_SDArun(), ".V", compN)
    
    CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
    CompNname = paste0("sda.", selected_SDArun, ".V", CompN)
    
    SpatialFeaturePlot(envv$MSTseqCells1, features = CompNname, 
                       max.cutoff = 'q99', min.cutoff = 'q01')  +
      theme_classic(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      )    &
      ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
    
    
    
    
  }
  
})


## SDA Score on STseqMT1 Spatial

output$SDAloadings_coord <- renderPlot({
  
  
  
  
  selected_SDArun = input$sda.run
  selected_compN = input$sda.comp.N
  
  
  
  
  if(is.null(envv$MoDSTA) | is.null(envv$commands$SDAproj_MSTseqCells1)){
    plot(x=0, y=0, main="Load MSTseq1 Cells dataset")
  } else {
    
    # # Split the component number from its label
    # compN <- as.numeric(strsplit(selected_compN(), " ")[[1]][2])
    # 
    # # Construct the dynamic variable name
    # compNname <- paste0("sda.", selected_SDArun(), ".V", compN)
    
    CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
    CompNname = paste0("sda.", selected_SDArun, ".V", CompN)
    
    
    
    # RedN = unlist(lapply(strsplit(SDAcomps[CompN], "\\."), function(x){gsub("V", "", x[2])}))  
    
    
    gg = plot_loadings_coordinates(SDARedDataLS=envv$SDARedDataLS,
                                   mart= envv$mart, genes = envv$genes, 
                                   highlight_genes =NULL, 
                                   reduction = selected_SDArun, 
                                   dimN = CompN, 
                                   #invertWeights = ifelse(xi %in% c(1), T, F),
                                   data_set = "mmusculus_gene_ensembl", #mmulatta_gene_ensembl
                                   includeUnMapped = T, 
                                   geneLocPath="./data/GeneCoordinants_FullCombo_March202024.rds"
    ); gg
    

    
  }
  
})


