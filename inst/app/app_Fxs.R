
# Envv update fx ---------


#' This function loads MoSTA
#' 
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Load_MoDSTA <- function(envv, input, session){
  
  print(getwd())
  
  if(is.null(envv$MoDSTA)) {
    
    envv$MoDSTA = readRDS("./data/MTDCombo_SerRedObj30K.rds")
  
  
  # envv$MoDSTA[["umapSup.harmony"]] = readRDS("../MouseTestis/MouseTestisAtlas2024/data/MTDCombo_HarmonyUnsupcellScores.rds")
  envv$MoDSTA$Pheno0 = unlist(lapply(strsplit(envv$MoDSTA$Pheno1, "_"), function(x){x[1]}))
  
  envv$MoDSTA$Breed_Pheno0 = paste0(gsub("C57BL/", "C57BL", envv$MoDSTA$Breed), "_", envv$MoDSTA$Pheno0)
  
  envv$MoDSTA$Condition = gsub("Embryonic D", "E", gsub("Postnatal ", "", gsub("WT 8W", "WT", envv$MoDSTA$SubjectId))) 
  
  envv$MoDSTA$Condition2 = gsub("D7", "Postnatal", gsub("D2", "Postnatal", gsub("E18.5", "Prenatal", envv$MoDSTA$Condition)))
  
  
  envv$MoDSTA$Breed_Cond2 = naturalsort::naturalfactor(paste0(envv$MoDSTA$Condition2, "_", gsub("C57BL/", "C57BL", envv$MoDSTA$Breed)))
  } else {
    print("MoDSTA already leaded ")
  }
  
  return(envv)
  
}


#' This function loads STseq Spatial omics mouse testis 1
#' 
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Load_STseqMT1 <- function(envv){
  

  if(is.null(envv$MSTseqCells1)) {
    
  envv$MSTseqCells1 = readRDS( "./data/B02622D6_CellsAdjusted.ser.proc_SerRedObj.rds")
  
  } else {
    print("STseq MT1 already leaded ")
  }
  
  return(envv)
  
}

#' This function loads SDA trained on MoSTA
#' 
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Load_SDA <- function(envv, input, session){
  
  if(is.null(envv$SDARedDataLS)) {
    
    envv$SDARedDataLS = readRDS("./data/SDARedDataLS_Apr82024.rds")
    names(envv$SDARedDataLS$loadings) = paste0("sda_", names(envv$SDARedDataLS$loadings))
    
  } else {
    print("SDA results already leaded ")
  }
  return(envv)
}