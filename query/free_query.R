free_query = function(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,items_list){
  withProgress(message = 'Progress...', min = 1,max = 9, {
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
    save(ADJ_S=ADJ_S,graph_s=graph_s,file="./asthma_query.RData")
    
    incProgress(1, detail = "Searching for cliques")
    #mcl = cliques_search(graph_s,min=4,max = 4) #in query_utilities.R
    mcl = NDDC_clique_search(ADJ=ADJ_S,nano,drugs,disease,chemical)


    incProgress(1, detail = "Evaluating cliques")
    #good_cliques = cliques_evaluation(mcl) #in query_utilities.R
    good_cliques = mcl
    
    incProgress(1, detail = "Clustering cliques")
    
    cliques_groups = cliques_clustering(good_cliques,ADJ_S) #in query_utilities.R
    output_cliques_pattern(input,output,cliques_groups) #in patterns_UI.R
    MList = group_cliques_list(cliques_groups,good_cliques)
    MM_list = clickable_cliques_list(MList,cliques_groups)
    
    incProgress(1, detail = "Building tables")
    proxy = building_table_and_proxy(input,output,MM_list)
    
  })#end progress bar
  
  save_free_query_results(input,output,selected_nodes,disease_list,THS,W_ADJ,graph_gw,SSNFQ,graph_s,ADJ_S,mcl,good_cliques,cliques_groups,MList,MM_list,"free_query_results")
  free_query_UI_ggplot_totale(input,output,MList,graph_gw) # in free_query_UI.R
  #plot_clique_graph(input,output,MM_list,graph_s,proxy) #in query_outputs.R
  enrich_clique(input,output,MList,MM_list,proxy,graph_s,items_list)
  #create_histrogram_condition(input,output,MList)#in pattern_UI.R; genera le combinazioni di items per l'istogramma
  barplot_patter_output(input,output,MList,graph_gw) #in query_outputs.R
  genes_data_table_output(input,output,MList,MM_list,proxy,graph_s,g,g_geni2,items_list) #in conditional_query_output.R
  plot_force_based_subnetwork_query_resutls(input,output,ADJ_S,chemMat,good_cliques,join10) #in qury_outputs.R
  plot_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10)
  plot_gene_subnetwork(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)
  plot_gene_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)
}