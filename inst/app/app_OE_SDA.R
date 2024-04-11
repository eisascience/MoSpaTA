observeEvent(input$Apply2MoDSTA, {
  
  envv$InfoBox_SDABrowser = "Projected SDA results trained on Mouse scRNA-Seq Testis Atlas (MoSTA)"
  
  selected_SDArun = input$sda.run
  # selected_compN = input$sda.comp.N
  
  if(!is.null(selected_SDArun)){
    print(selected_SDArun)
    
    envv$MoDSTA = ImputeSDA2SerV2(SerObj = envv$MoDSTA,
                                  sda_loadings = envv$SDARedDataLS$loadings[[selected_SDArun]]$loadings,
                                  # keepComps = unique(c(1, 2, CompN)),
                                  sdaObjID = selected_SDArun, plot=F, MakeReduc = F)
    print("Projection complete")
    envv$SDAcomps = colnames(envv$MoDSTA@meta.data)[grep("sda.", colnames(envv$MoDSTA@meta.data))]
    envv$SDAcomps = naturalsort::naturalsort(envv$SDAcomps)
    envv$commands$SDAproj_MoDSTA = T
    # envv$CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
    # envv$CompNname = paste0("sda.", selected_SDArun, ".V", envv$CompN)
    
  } else{ 
    print("selected_SDArun is NULL")
    }
  
  
  
})

observeEvent(input$Apply2STseqMT1, {
  
  envv$InfoBox_SDABrowser = "Projected SDA results trained on Mouse scRNA-Seq Testis Atlas (MoSTA)"
  
  selected_SDArun = input$sda.run
  
  if(!is.null(selected_SDArun)){
    
  envv$MSTseqCells1 = ImputeSDA2SerV2(SerObj = envv$MSTseqCells1 ,
                                      sda_loadings = envv$SDARedDataLS$loadings[[selected_SDArun]]$loadings,
                                      # keepComps = unique(c(1, 2, CompN)),
                                      sdaObjID = selected_SDArun, plot=F, MakeReduc = F, assay = "SCT")
  print("Projection complete")
  envv$commands$SDAproj_MSTseqCells1 = T
  
  # envv$SDAcomps = colnames(envv$MSTseqCells1@meta.data)[grep("sda.", colnames(envv$MSTseqCells1@meta.data))]
  # envv$SDAcomps = naturalsort::naturalsort(envv$SDAcomps)
  } else{ 
    print("selected_SDArun is NULL")
  }
  
})

observeEvent(input$NextComp, {
  selected_compN = input$sda.comp.N
  CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
 
  n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
  
  
  if (CompN < n_comps) {
    new_comp <- CompN + 1
    updateSelectInput(session, "sda.comp.N", selected = paste0("Comp_", new_comp))
  }
  
})


observeEvent(input$PreviousComp, {
  selected_compN = input$sda.comp.N
  CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
  
  n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
  
  
  if (CompN > 1) {
    new_comp <- CompN - 1
    updateSelectInput(session, "sda.comp.N", selected = paste0("Comp_", new_comp))
  }
  
})
