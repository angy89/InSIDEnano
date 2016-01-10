
clique_3_search = function(ADJ,type="NDD"){
  mcl = list()
  
  if(type=="NDD"){
    if(DEBUGGING)
      cat("clique_3_search function type--> NDD \n")
    x_ADJ = colnames(ADJ)[colnames(ADJ) %in% nano]
    y_ADJ = colnames(ADJ)[colnames(ADJ) %in% drugs]
    z_ADJ = colnames(ADJ)[colnames(ADJ) %in% disease]
  }
  if(type=="NDC"){
    if(DEBUGGING)
      cat("clique_3_search function type--> NDC \n")    
    x_ADJ = colnames(ADJ)[colnames(ADJ) %in% nano]
    y_ADJ = colnames(ADJ)[colnames(ADJ) %in% drugs]
    z_ADJ = colnames(ADJ)[colnames(ADJ) %in% chemical]
  }
  if(type=="DCD"){
    if(DEBUGGING)
      cat("clique_3_search function type--> DCD \n")
    x_ADJ = colnames(ADJ)[colnames(ADJ) %in% drugs]
    y_ADJ = colnames(ADJ)[colnames(ADJ) %in% chemical]
    z_ADJ = colnames(ADJ)[colnames(ADJ) %in% disease]
  }
  
  
  i = 1
  for(n in x_ADJ){
  n_drugs = y_ADJ[ADJ[n,y_ADJ]!=0]
    for(dr in n_drugs){
      n_dis = z_ADJ[ADJ[dr,z_ADJ]!=0]
      for(di in n_dis){
        
        and_sum = (ADJ[n,dr]!=0) & (ADJ[dr,di]!=0) & (ADJ[di,n]!=0)
        if(and_sum){
          cliq = c(n,dr,di)
          names(cliq) = c(n,dr,di)
          mcl[[i]]=cliq
          i = i + 1
        }
      }
    }
  }
  
  if(DEBUGGING)
    cat("Nro cliques: ",length(mcl),"\n")
  return(mcl)
}

NDDC_clique_search = function(ADJ,nano,drugs,disease,chemical){
  mcl = list()
  nano_ADJ = colnames(ADJ)[colnames(ADJ) %in% nano]
  drug_ADJ = colnames(ADJ)[colnames(ADJ) %in% drugs]
  dis_ADJ = colnames(ADJ)[colnames(ADJ) %in% disease]
  chem_ADJ = colnames(ADJ)[colnames(ADJ) %in% chemical]
  
  if(DEBUGGING){
    cat("NDDC_clique_search function!\n")
    cat("ADJ--> ",dim(ADJ),"\n")
    
    cat("nano_ADJ--> ",length(nano_ADJ),"\n")
    cat("drug_ADJ--> ",length(drug_ADJ),"\n")
    cat("dis_ADJ--> ",length(dis_ADJ),"\n")
    cat("chem_ADJ--> ",length(chem_ADJ),"\n")
    
    
  }
  
  i = 1
  for(n in nano_ADJ){
    n_drugs =     drug_ADJ[ADJ[n,drug_ADJ]!=0]
    for(dr in n_drugs){
      n_dis = dis_ADJ[ADJ[dr,dis_ADJ]!=0]
      for(di in n_dis){
        n_chem = chem_ADJ[ADJ[di,chem_ADJ]!=0]
        for(c in n_chem){
          and_sum = (ADJ[n,dr]!=0) & (ADJ[n,di]!=0) & (ADJ[n,c]!=0) & (ADJ[dr,di]!=0) & (ADJ[dr,c]!=0) & (ADJ[di,c]!=0)
          
          if(and_sum){
            cliq = c(n,dr,di,c)
            names(cliq) = c(n,dr,di,c)
            mcl[[i]]=cliq
            i = i +1
          }
        }
        
      }
    }
  }
  if(DEBUGGING)
    cat("NDDC_clique_search function Number of cliques --> ",length(mcl),"\n")
  return(mcl)
}