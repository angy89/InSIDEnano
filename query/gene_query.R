gene_query = function(input,output,disease_list,selected_items_,W_ADJ,th_p = input$th_slider3/100,node_type,chemMat,join10,g,g_geni2,gene_input){

   # selected_genes = input$gene_input
  selected_genes = genes_input[1]  
  if(DEBUGGING){
      cat("class(Selected genes): ",class(selected_genes),"\n")
      cat("Selected genes: ",selected_genes,"\n")
    }
    
    if("ALL" %in% selected_genes){
      genes_to_be_queried = genes_input
    }else{
      genes_to_be_queried = selected_genes 
    }
    
    selected_items_ = igraph::neighbors(graph = g,v = genes_to_be_queried)
    selected_items_ = names(selected_items_)
    selected_items__idx = which(selected_items_ %in% rownames(W_ADJ))
    selected_items_ = selected_items_[selected_items__idx]
    
    if(DEBUGGING){
      cat("selected item based on gene:", selected_items_,"\n")
    }
  
    gene_nano = selected_items_[which(selected_items_ %in% nano)]
    gene_drug = selected_items_[which(selected_items_ %in% drugs)]
    gene_chem = selected_items_[which(selected_items_ %in% chemical)]
    gene_dis = selected_items_[which(selected_items_ %in% disease)]
    
    gene_query_UI_part2(input,output,gene_nano,gene_drug,gene_chem,gene_dis)
      
  
}