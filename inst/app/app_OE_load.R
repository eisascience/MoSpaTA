observeEvent(input$MoDSTA, {
  
  envv$InfoBox_Load = "Loaded Mouse scRNA-Seq Testis Atlas (MoSTA)"
  
  envv = Load_MoDSTA(envv)
  

})

observeEvent(input$STseqMouse1, {
  
  envv$InfoBox_Load = "Loaded Stereo-seq Spatial omics mouse testis 1"
  
  envv = Load_STseqMT1(envv)
  
})

observeEvent(input$SDAres, {
  
  envv$InfoBox_Load = "Loaded SDA results trained on Mouse scRNA-Seq Testis Atlas (MoSTA)"
  
  envv = Load_SDA(envv)
  
})

observeEvent(input$AllInputs, {
  
  envv$InfoBox_Load = "Loaded all inputs, MoSTA, Spatial, and SDA"
  
  envv = Load_MoDSTA(envv)
  envv = Load_STseqMT1(envv)
  envv = Load_SDA(envv)
  
})