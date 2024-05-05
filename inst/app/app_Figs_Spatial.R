

output$GeneExpr_SlideSeqV1MT1_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  

  
  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_SlideSeqV1MT1_Spatial(envv, input, GeneSet)
    
  } else {
    
    # gg1
    
    Threshold_GeneExpr_SlideSeqV1MT1_Spatial(envv, input, GeneSet)
      
    }
    
 
  
})


output$GeneExpr_SlideSeqV1MT2_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  
  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_SlideSeqV1MT2_Spatial(envv, input, GeneSet)
    
  } else {
    
    # gg1
    
    Threshold_GeneExpr_SlideSeqV1MT2_Spatial(envv, input, GeneSet)
    
  }
  
  
  
})



output$GeneExpr_SlideSeqV1MT3_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  
  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_SlideSeqV1MT3_Spatial(envv, input, GeneSet)
    
  } else {
    
    # gg1
    
    Threshold_GeneExpr_SlideSeqV1MT3_Spatial(envv, input, GeneSet)
    
  }
  
  
  
})



output$GeneExpr_STseqMT1_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  

  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_STseqMT1_Spatial(envv, input, GeneSet)
    
  } else {
    
    
    Threshold_GeneExpr_STseqMT1_Spatial(envv, input, GeneSet)
    
   
    
  }
  
  
})
