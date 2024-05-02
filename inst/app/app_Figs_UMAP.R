



## Gene Expr MoDSTA Unintegrated UMAP ------
# 
# output$GeneExpr_MoDSTA_UMAP <- renderPlot({
#   
#   GeneSet <- input$GeneSet
#   
#   
#   if(length(grep(",", GeneSet)) == 0){
#     
#     if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
#       GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
#     } else {
#       GeneSet <- unlist(strsplit(GeneSet, " "))
#     }
#     
#     
#   } else {
#     GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
#     
#   }
#   
#   GeneSet = GeneSet[GeneSet %in% rownames(envv$MoDSTA)]
#   
#   print(GeneSet)
#   
#   if(is.null(envv$MoDSTA) | length(GeneSet) < 1){
#     plot(x=0, y=0, main="Load MoDSTA dataset")
#   } else {
#     
#     if(length(GeneSet)==1) {
#       
#       FeaturePlot(envv$MoDSTA , reduction = "umap", features = GeneSet[1], 
#                   # raster = T, raster.dpi = c(800, 800),
#                   order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
#         theme_classic(base_size = 14) +
#         theme(axis.line = element_blank(),
#               axis.text.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank() #,plot.title = element_blank()
#         )  &
#         ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
#       
#     } else{
#       
#       if(input$GeneExprOpr == "sum"){
#         envv$MoDSTA$MultiGeneExpr = colSums(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#         
#       } else {
#         envv$MoDSTA$MultiGeneExpr = colMeans(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#         
#       }
#       
#       
#       FeaturePlot(envv$MoDSTA , reduction = "umap", features = "MultiGeneExpr",
#                   # raster = T, raster.dpi = c(800, 800),
#                   order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
#         theme_classic(base_size = 14) +
#         theme(axis.line = element_blank(),
#               axis.text.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank() #,plot.title = element_blank()
#         ) + ggtitle(ifelse(input$GeneExprOpr == "sum", "Summed Expression", "Mean of Expression Set"))  &
#         ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
#       
#     }
# 
#     
#     
#   }
#   
# })
# 

## Gene Expr MoDSTA Integrated UMAP ------

output$GeneExpr_MoDSTA_UMAPintg <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  
  if(input$ExprThresh %in% c(0, 1) ){

    Render_GeneExpr_MoDSTA_UMAPintg(envv, input, GeneSet)

  } else {

    Threshold_GeneExpr_MoDSTA_UMAPintg(envv, input, GeneSet)


  }
  
})


## Gene Expr STseq MT1 UMAP ------

output$GeneExpr_STseqMT1_UMAP <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  gg1 = Render_GeneExpr_STseqMT1_UMAP(envv, input, GeneSet)
  
  gg1
 
  
})



## Gene Expr SlideSeqV1 WT1 UMAP ------

output$GeneExpr_SlideSeqV1MT1_UMAP <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  gg1 = Render_GeneExpr_SlideSeqV1MT1_UMAP(envv, input, GeneSet)
  
  gg1
  
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



