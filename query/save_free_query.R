#object to save
#selected_nodes
#disease_list
#THS
#W_ADJ
#graph_gw
#SSNFQ
save_free_query_results = function(input,output,selected_nodes,disease_list,THS,W_ADJ,graph_gw,SSNFQ,graph_s,ADJ_S,mcl,good_cliques,cliques_groups,MList,MM_list,file_name){
  output$downloadData <- downloadHandler(
    filename = function() { paste(file_name, '.RData', sep='') },
    content = function(file_name) {
      save(selected_nodes,disease_list,THS,W_ADJ,graph_gw,SSNFQ,graph_s,ADJ_S,mcl,good_cliques,cliques_groups,MList,MM_list, file=file_name)
    }
  )
}