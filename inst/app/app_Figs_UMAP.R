



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
    
    if(input$ExprThresh %in% c(0, 1) ){
      
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
      
    } else {
      # when threshold is not 0 or 1
      
      if(length(GeneSet)==1) {
        envv$MSTseqCells1$SummedExpr = GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ]
        
      } else {
        envv$MSTseqCells1$SummedExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      envv$MSTseqCells1$SummedExprBL = ifelse(envv$MSTseqCells1$SummedExpr>max(envv$MSTseqCells1$SummedExpr) * input$ExprThresh,
                                              "Above", "Below")
      
      SpatialDimPlot(envv$MSTseqCells1 , 
                     group.by = "SummedExprBL",
      ) +facet_wrap(~SummedExprBL, ncol=4) +
        theme_classic(base_size = 14) +
        guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
        theme(legend.position = "none",
              axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        ) + ggtitle(paste0("Threshold cut above: ", round(max(envv$MSTseqCells1$SummedExpr) * input$ExprThresh,3)))
      

      
      
      
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

## DatasetName on Integrated UMAP ------

output$DatasetName_MoDSTA_UMAPintg  <- renderPlot({
  
  DimPlot(envv$MoDSTA , reduction = "umapSup.harmony", 
          group.by = "DatasetName",
          repel = T,
          cols = col_vector, label = F,
          raster = T, raster.dpi = c(800, 800))  +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "right",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("DatasetName")
  
  
  
})

## DatasetGrouping on Integrated UMAP ------

output$DatasetGrouping_MoDSTA_UMAPintg  <- renderPlot({
  
  DimPlot(envv$MoDSTA , reduction = "umapSup.harmony", 
          group.by = "SubjectId",
          repel = T,
          cols = col_vector, label = F,
          raster = T, raster.dpi = c(800, 800))  +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "right",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("DatasetGrouping")
  
  
  
})

## Breed on Integrated UMAP ------

output$Breed_MoDSTA_UMAPintg  <- renderPlot({
  
  DimPlot(envv$MoDSTA , reduction = "umapSup.harmony", 
          group.by = "Breed",
          repel = T,
          cols = col_vector, label = F,
          raster = T, raster.dpi = c(800, 800))  +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "right",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("Breed")
  
  
  
})

## UnsupClusters on Integrated UMAP ------

output$UnsupClusters_MoDSTA_UMAPintg  <- renderPlot({
  
  DimPlot(envv$MoDSTA , reduction = "umapSup.harmony", 
          group.by = "RNA_snn_res.0.6",
          repel = T,
          cols = col_vector, label = F,
          raster = T, raster.dpi = c(800, 800))  +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=2)) +
    theme(legend.position = "right",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("RNA_snn_res.0.6")
  
  
  
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


## UnsupClusters Faceted of spatial on STseq MT1 UMAP ------

output$UnSupClusters_SpatialFacet_STseqMT1_UMAP  <- renderPlot({
  
  SpatialDimPlot(envv$MSTseqCells1 , 
          group.by = "SCT_snn_res.0.1",
  ) +facet_wrap(~SCT_snn_res.0.1, ncol=4) +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle("Transcriptional Groups")
  
  
  
})



