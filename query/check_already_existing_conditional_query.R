check_already_existing_conditional_query = function(input,output,LOG_CONDITIONAL){
  if(DEBUGGING){
    cat("--->START check_already_existing_conditional_query function ---> \n")
  }
  n_log = dim(LOG_CONDITIONAL)[1]
  library(plyr)
  
  nano_query = paste(c(input$nano_input,""),collapse="_")
  drug_query = input$drug_input
  chemical_query = input$chemical_input
  disease_query = paste(c(input$disease_input,""),collapse="_")
  th_p = input$th_slider2
  query_th = input$percentuale_somma
  nElem_cliques = input$nroCliques
  CLIQUE_TYPE = input$clique_type
  
  if(DEBUGGING){
    cat(input$disease_input[1]," <--> ",input$disease_input[2],"\n")
    cat(nchar(input$disease_input[1])," <--> ",nchar(input$disease_input[2]),"\n")
    
    cat("class(input$disease_input) ",class(input$disease_input),"\n")
    cat("into function Inputs--> ", nano_query,drug_query,chemical_query,disease_query,th_p,query_th,nElem_cliques,CLIQUE_TYPE,"\n")
  }
  
row_to_find <- data.frame(Nano=nano_query,Drugs = "",Disease=disease_query,Chemical="",Th=th_p,n1=query_th,n2 = nElem_cliques,type=CLIQUE_TYPE)
  #row_to_find <- data.frame(Nano="MWCNT",Drugs = "",Disease="Asthma",Chemical="",Th=80,n1=2,n2 = 2,type="ALL")
  xx = match_df(LOG_CONDITIONAL[,1:8], row_to_find)
  
  if(DEBUGGING){
    cat("dim(xx) ",dim(xx),"\n")
  }
  
  if(dim(xx)[1]>0){
    row_idx = as.numeric(rownames(xx))
    file_name = paste(APP_PATH,"Log_folder/",as.character(LOG_CONDITIONAL[row_idx,9]),sep="")
    cat("---> END check_already_existing_conditional_query function ---> \n")
    
    return(file_name)
  }else{
    cat("---> END check_already_existing_conditional_query function ---> \n")
    
    return(NULL)
  }
}