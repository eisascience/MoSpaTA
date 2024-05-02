

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