free_query = function(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,items_list){
  output$infoFree <- renderUI({HTML(info_free_query_text)}) 
  
  withProgress(message = 'Progress...', min = 1,max = 10, {
  
    validate(need(input$disease != "", "Please select a Disease"))
    
    SNC = select_node_query(input,output,disease_list,selected_nodes)
    selected_nodes = SNC$selected_nodes
    disease_list = SNC$disease_list
    
    incProgress(1, detail = "Evaluating input list...")
    free_query_UI_node_of_interest_output(input,output,disease_list)
    
    incProgress(1, detail = "Thresholding...")
    THS = find_thresholds(W_ADJ,th_p ) #in query_utilities.R
    
    incProgress(1, detail = "Removing edges under threshold...")
    W_ADJ = apply_thresholds(W_ADJ,THS) #in query_utilities.R
    
    incProgress(1, detail = "Creating graph...")
    graph_gw = creating_graph(W_ADJ,node_type) #in query_utilities.R
    
    incProgress(1, detail = "Creating subgraph of selected nodes")
    SSNFQ = subgraph_selected_nodes_free_query(graph_gw,selected_nodes) #in query_utilities.R
    graph_s=SSNFQ$graph_s
    ADJ_S=SSNFQ$ADJ_S
  
    info_free_query_text = paste(info_free_query_text, "Number of elements in the matrix:",nrow(ADJ_S),"<br/>")
  
    #save(ADJ_S=ADJ_S,graph_s=graph_s,file="./asthma_query.RData")
    
    incProgress(1, detail = "Searching for cliques")
    #mcl = cliques_search(graph_s,min=4,max = 4) #in query_utilities.R
    mcl = NDDC_clique_search(ADJ=ADJ_S,nano,drugs,disease,chemical)
  
    if(length(mcl)==0){
      info_free_query_text = "No results!"
    }else{
      info_free_query_text = paste(info_free_query_text, "Number of cliques:",length(mcl),"<br/>")
    }
  
    validate(need(length(mcl)>0, paste("No clique with this threshold")))

    incProgress(1, detail = "Evaluating cliques")
    #good_cliques = cliques_evaluation(mcl) #in query_utilities.R
    good_cliques = mcl
    
    incProgress(1, detail = "Clustering cliques")
    
    cliques_groups = cliques_clustering(good_cliques,ADJ_S) #in query_utilities.R
    info_free_query_text = paste(info_free_query_text, "Number of clique groups:",length(cliques_groups),"<br/>")
  
    output_cliques_pattern(input,output,cliques_groups) #in patterns_UI.R
    MList = group_cliques_list(cliques_groups,good_cliques)
    MM_list = clickable_cliques_list(MList,cliques_groups)
    
    incProgress(1, detail = "Building tables")
    proxy = building_table_and_proxy(input,output,MM_list)
    
    incProgress(1, detail = "Preparing Output")
    
    #save_free_query_results(input,output,selected_nodes,disease_list,THS,W_ADJ,graph_gw,SSNFQ,graph_s,ADJ_S,mcl,good_cliques,cliques_groups,MList,MM_list,"free_query_results")
    if(DEBUGGING){
      cat("Plotting free query ggplot totale\n")
    }
    free_query_UI_ggplot_totale(input,output,MList,graph_gw) # in free_query_UI.R
    #plot_clique_graph(input,output,MM_list,graph_s,proxy) #in query_outputs.R
    if(DEBUGGING){
      cat("Enrich cliques\n")
    }
    plot_clique_graph(input,output,MM_list,graph_s,proxy)
      
    enrich_clique(input,output,MList,MM_list,proxy,graph_s,g,items_list,"FREE")
    
    #create_histrogram_condition(input,output,MList)#in pattern_UI.R; genera le combinazioni di items per l'istogramma
    if(DEBUGGING){
      cat("Plotting barplot for cliques\n")
    }
    barplot_patter_output(input,output,MList,graph_gw) #in query_outputs.R
    if(DEBUGGING){
      cat("Preparing genes data table\n")
    }
    genes_data_table_output(input,output,MList,MM_list,proxy,graph_s,g,g_geni2,items_list,"FREE") #in conditional_query_output.R
    if(DEBUGGING){
      cat("Preparing force based subnetwork\n")
    }
    plot_force_based_subnetwork_query_resutls(input,output,ADJ_S,chemMat,good_cliques,join10) #in qury_outputs.R
    if(DEBUGGING){
      cat("Preparing force based subnetwork statistics\n")
    }
    plot_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10)
    if(DEBUGGING){
      cat("Preparing force based geene subnetwork\n")
    }
    plot_gene_subnetwork(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)
    if(DEBUGGING){
      cat("Preparing force based gene subnetwork statistics\n")
    }
    plot_gene_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)
    #save.image("Free_query_debugging.RData")
    
  })#end progress bar
  
  if(DEBUGGING){
    cat("End Free Query Function!\n")
  }
}