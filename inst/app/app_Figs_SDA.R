
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
                                   # mart= envv$mart, 
                                   genes = envv$genes, 
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
    
    
    if(input$SDAThresh %in% c(0, 1) ){
      
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
    
    } else {
      
      
      MaxVal = max(envv$MSTseqCells1@meta.data[, CompNname])
      
      
      envv$MSTseqCells1$SDAScoreBL = ifelse( envv$MSTseqCells1@meta.data[, CompNname]  > MaxVal * input$SDAThresh,
                                              "Above", "Below")
      
      SpatialDimPlot(envv$MSTseqCells1 , 
                     group.by = "SDAScoreBL",
      ) +facet_wrap(~SDAScoreBL, ncol=4) +
        theme_classic(base_size = 14) +
        guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
        theme(legend.position = "none",
              axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        ) + ggtitle(paste0("Threshold cut above: ", round(MaxVal * input$SDAThresh,3)))
      
      
      
      
    }
    
    
  }
  
})









  