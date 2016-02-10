load_conditional_query = function(input,output,RData_file,items_list,g){
  
  
  withProgress(message = 'Progress...', min = 1,max = 6, {
    incProgress(1, detail = "Loading file...")
    
    load(file = RData_file)
    
    output$info2_1 <- renderUI({HTML(info_text)}) 
    free_query_UI_node_of_interest_output(input,output,disease_list) #in free_query_UI.R
    
    output$extimatedTime = renderUI({
      HTML(paste("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size<br/> <strong>Estimated Iteration: ",estimated_tyme,"<strong/><br/>"))
    })
    
    output$NetworkPattern <- renderUI({
      cliques_LL = list()
      
      for(i in 1:length(cliques_groups)){
        cliques_LL[[paste0("Type",i)]] = paste0("M",i)
      }
      selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
    })
    
    output$clique_data_table = DT::renderDataTable({
      type = input$NetworkPattern
      validate(
        need(input$NetworkPattern != "", "Please select a pattern type")
      )
      type = as.integer(gsub(pattern = "M",x =type,replacement = ""))  
      
      DT::datatable(data =  MM_list[[type]],
                    options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE,scrollY = "400px", scrollCollapse = TRUE,paging=FALSE),
                    escape=FALSE,
                    selection = "single")
      
    })
    
    incProgress(1, detail = "Preparing Output...")
    
    proxy = dataTableProxy("clique_data_table")
    incProgress(1, detail = "Preparing Barplot output...")
    
    barplot_pattern_conditional_query(input,output,MList,graph_gw) #in conditional_query_output.R
    clique_graph_cq_plot(input,output,MList,MM_list,proxy,graph_s)#in conditional_query_output.R
    incProgress(1, detail = "Preparing Genes Data table Output...")
    
    genes_data_table_output(input,output,MList,MM_list,proxy,graph_s,g,g_geni2,items_list,"CONDITIONAL") #in conditional_query_output.R
    incProgress(1, detail = "Enrich Cliques...")
    
    enrich_clique(input,output,MList,MM_list,proxy,graph_s,g,items_list,"CONDITIONAL")
    
    incProgress(1, detail = "Plot subnetwork and statistics...")
    
    plot_force_based_subnetwork_query_resutls(input,output,ADJ_S,chemMat,good_cliques,join10) #in qury_outputs.R
    plot_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10)#in qury_outputs.R
    plot_gene_subnetwork(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
    plot_gene_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
    
  })
  
}