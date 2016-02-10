conditional_query = function(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,LOG_CONDITIONAL,items_list){
  output$info2_1 <- renderUI({HTML(info_text)}) 
  CLIQUE_TYPE = input$clique_type
  if(DEBUGGING)
    cat("CLIQUE_TYPE", CLIQUE_TYPE, "\n")
  query_th = input$percentuale_somma
  if(DEBUGGING)
    cat("query_th", query_th, "\n")
  nElem_cliques = input$nroCliques
  if(DEBUGGING)
    cat("nElem_cliques", nElem_cliques, "\n")
  
  if(DEBUGGING){
    cat("CONDITIONAL QUERY THslider2 --> ",input$th_slider2,"\n")
  }
  
  validate( need(input$percentuale_somma != "", "Please select the treshold"))
  validate( need(input$nroCliques != "", "Please select the number of nodes"))
  
  CQN = conditional_query_nodes(input,output)
  query_nodes=CQN$query_nodes
  disease_list=CQN$disease_list
  selected_nodes=CQN$selected_nodes
  
  if(DEBUGGING)
    cat("query_nodes:",query_nodes,"\n")
  
  withProgress(message = 'Progress...', min = 1,max = 11, {
    incProgress(1, detail = "Evaluating input list...")
    free_query_UI_node_of_interest_output(input,output,disease_list) #in free_query_UI.R
    
    incProgress(1, detail = "Thresholding...")
    THS = find_thresholds(W_ADJ,th_p = input$th_slider2/100 ) #in query_utilities.R
    
    incProgress(1, detail = "Removing edges under threshold...")
    W_ADJ = apply_thresholds(W_ADJ,THS) #in query_utilities.R
    
    nano_qn_e = CQN$nano_query
    drug_qn_e = CQN$drug_query
    dis_qn_e = CQN$disease_query
    chem_qn_e = CQN$chemical_query
    type_qn = CQN$type_qn
    combination_vect = CQN$combination_vect
    selected_nodes = CQN$selected_nodes
    
    
    if(DEBUGGING){
      message("In conditional_query::: nano_qn_e:",nano_qn_e,"\n")
      message("In conditional_query::: drug_qn_e:",drug_qn_e,"\n")
      message("In conditional_query::: dis_qn_e:",dis_qn_e,"\n")
      message("In conditional_query::: chem_qn_e:",chem_qn_e,"\n")
      message("In conditional_query::: combination_vect:",combination_vect,"\n")
    }
    
    
    CQI = conditional_query_items(input,output,nano_qn_e,drug_qn_e,
                                       dis_qn_e,chem_qn_e,type_qn,combination_vect,
                                       query_nodes,W_ADJ,query_th,selected_nodes)
    
    W_ADJ=CQI$W_ADJ
    combinations=CQI$combinations
    toREM=CQI$toREM
    Col_Sum_list=CQI$Col_Sum_list
    
    if(dim(W_ADJ)[1] == 0){
      if(DEBUGGING)
        cat("ARGHHH --> No items to display! Try to decreases the threshold.", max(unlist(Col_Sum_list)),"\n") 
    }
    
    incProgress(1, detail = "Creating graph...")

    CQSC = conditional_query_subgraph_creation(input,output,W_ADJ,toREM,info_text,node_type)
    graph_gw=CQSC$graph_gw
    graph_s=CQSC$graph_s
    ADJ_S = CQSC$ADJ_S
    info_text=CQSC$info_text  
    
    incProgress(1, detail = "Searching for cliques")
    
    estimated_tyme = estimate_query_tyme(input,output,graph_s)
      
    if(CLIQUE_TYPE == "ALL"){
      mcl_NDD = clique_3_search(ADJ_S,type="NDD")
      mcl_NDC = clique_3_search(ADJ_S,type="NDC")
      mcl_DCD = clique_3_search(ADJ_S,type="DCD")
      mcl_NDDC = NDDC_clique_search(ADJ_S,nano,drugs,disease,chemical)
      mcl = c(mcl_NDD,mcl_NDC,mcl_DCD,mcl_NDDC)
      good_cliques = mcl
    }#both cliques of length 3 and 4
    if(CLIQUE_TYPE == "NDCD"){ #clique da 4 oggetti
      mcl = NDDC_clique_search(ADJ = ADJ_S,nano,drugs,disease,chemical)
      good_cliques = mcl
    }# end IF clique di 4
    
    if(CLIQUE_TYPE == "NDD"){ #clique Nano-Drug-Disease
      mcl = clique_3_search(ADJ = ADJ_S,type="NDD")
      good_cliques = mcl
    }# end IF Nano-Drug-Disease
    if(CLIQUE_TYPE == "NDC"){ #clique Nano-Drug-Chemical
      mcl = clique_3_search(ADJ = ADJ_S,type="NDC")
      good_cliques = mcl
    }# end IF Nano-Drug-Chemical
    if(CLIQUE_TYPE == "DCD"){ #clique Drug-Chemical-Disease
      mcl = clique_3_search(ADJ = ADJ_S,type="DCD")
      good_cliques = mcl
    }# end IF Drug-Chemical-Disease
    
    incProgress(1, detail = "Evaluating cliques")
    
    is_good2 = unlist(lapply(X = good_cliques,FUN = function(obj){
      length(which(selected_nodes %in% names(obj)))
    }))
    
    idx2 = which(is_good2>=nElem_cliques)
    good_cliques = good_cliques[idx2]
    if(DEBUGGING)
      cat("Nro good cliques-->: ",length(good_cliques),"\n")
    
    if(length(good_cliques) == 0){
      info_text = "No results! \n"
    }
    
    validate(need(length(good_cliques)>0, paste("No clique with this threshold")))
    info_text = paste(info_text, "Number of cliques:",length(good_cliques),"<br/>")
    incProgress(1, detail = "Clustering cliques")
    
    cliques_groups=conditional_clustering_cliques(input,output,good_cliques,ADJ_S)
    info_text = paste(info_text, "Number of clique groups:",length(cliques_groups),"<br/>")
    plot_list_pattern(input,output,cliques_groups)
    
    if(DEBUGGING)
      cat("Nro of cliques groups: ",length(cliques_groups),"\n")
    
    BQL = build_query_list(input,output,cliques_groups,nano,drugs,chemical,disease,good_cliques)
    MM_list = BQL$MM_list
    MList=BQL$MList
    
    incProgress(1, detail = "Building tables")
    building_table_and_proxy(input,output,MM_list)
  
    incProgress(1, detail = "Saving Results")
    dim(LOG_CONDITIONAL)[1] -> log_counter
    log_counter + 1 -> log_counter
    
    save(CLIQUE_TYPE,query_th,nElem_cliques,query_nodes,disease_list,selected_nodes,
         info_text,W_ADJ,MList,MM_list,graph_gw,g,g_geni2,graph_s,ADJ_S,chemMat,
         good_cliques,join10,estimated_tyme,cliques_groups,
         file=paste(LOCAL_PATH,"Log_folder/",log_counter,".RData",sep=""))
    
    if(DEBUGGING){
      nano_query = input$nano_input
      drug_query = input$drug_input
      chemical_query = input$chemical_input
      disease_query = input$disease_input
      th_p = input$th_slider2
      query_th = input$percentuale_somma
      nElem_cliques = input$nroCliques
      CLIQUE_TYPE = input$clique_type
      
      cat("Inputs--> ", nano_query,drug_query,chemical_query,disease_query,th_p,query_th,nElem_cliques,CLIQUE_TYPE,"\n")
    }
    
    LOG_CONDITIONAL = as.matrix(LOG_CONDITIONAL)
    LOG_CONDITIONAL = rbind(LOG_CONDITIONAL,c(paste(c(input$nano_input,""),collapse="_"),
                                              "",
                                              paste(c(input$disease_input,""),collapse="_"),
                                              "",
                                              input$th_slider2,input$percentuale_somma,input$nroCliques,input$clique_type,paste(log_counter,".RData",sep="")))
    cat("New row to add to the log\n")
    cat("Row--> ",c(paste(c(input$nano_input,""),collapse="_"),
                    "",
                    paste(c(input$disease_input,""),collapse="_"),
                    "",
                    input$th_slider2,input$percentuale_somma,input$nroCliques,input$clique_type,paste(log_counter,".RData",sep="")),"\n")
    
    LOG_CONDITIONAL = as.data.frame(LOG_CONDITIONAL)
    
    save(LOG_CONDITIONAL,file=paste(LOCAL_PATH,"LOG.RData",sep=""))
    
    incProgress(1, detail = "Preparing Output")
    
    clique_graph_cq_plot(input,output,MList,MM_list,proxy,graph_s)#in conditional_query_output.R
    #create_histrogram_condition(input,output,MList)#in pattern_UI.R; genera le combinazioni di items per l'istogramma      
    enrich_clique(input,output,MList,MM_list,proxy,graph_s,g,items_list,"CONDITIONAL")
    genes_data_table_output(input,output,MList,MM_list,proxy,graph_s,g,g_geni2,items_list,"CONDITIONAL") #in conditional_query_output.R
    barplot_pattern_conditional_query(input,output,MList,graph_gw) #in conditional_query_output.R
    #save(ADJ_S,chemMat,good_cliques,join10,file = "/home/aserra/InsideNano/www/immagine_per_debugging.RData")
    plot_force_based_subnetwork_query_resutls(input,output,ADJ_S,chemMat,good_cliques,join10) #in qury_outputs.R
    plot_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10)#in qury_outputs.R
    plot_gene_subnetwork(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
    plot_gene_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
    
    })
  
}