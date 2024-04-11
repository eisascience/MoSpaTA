
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