conditional_table_proxy = function(input,output,MM_list){
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
  proxy = dataTableProxy("clique_data_table")
  return(proxy)
}


plot_list_pattern = function(input,output,cliques_groups){
  output$NetworkPattern <- renderUI({
    cliques_LL = list()
    
    for(i in 1:length(cliques_groups)){
      cliques_LL[[paste0("Type",i)]] = paste0("M",i)
    }
    selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
  })
}

clique_graph_cq_plot = function(input,output,MList,MM_list,proxy,graph_s){
  output$xx = renderPlot({
    type = input$NetworkPattern
    clique_type = input$clique_type
    
    validate(
      need(input$NetworkPattern != "", "Please select a pattern type")
    )
    if(DEBUGGING)
      cat("Il tipo selezionato ??: ",type,"\n")
    type = as.integer(gsub(pattern = "M",x =type,replacement = ""))
    if(DEBUGGING)
      cat("Il numero selezionato ??: ",type,"\n")
    if(DEBUGGING)
      cat("class input$clique_data_table: ",class(input$clique_data_table))
    s = input$clique_data_table_rows_selected
    if(DEBUGGING)
      cat("s vale: ",s,"\n")
    
    if(clique_type == "NDCD"){
      g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "NDD"){
      g_= internal_render_plotNDD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "NDC"){
      g_= internal_render_plotNDC(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "DCD"){
      g_= internal_render_plotDCD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "ALL"){
      if(("" %in% MList[[type]][1,]) == FALSE){
        g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
      }else{
        col_idx = which(MList[[type]][1,] %in% "")
        if(col_idx == 4){
          g_= internal_render_plotNDD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
        }
        if(col_idx == 1){
          g_= internal_render_plotNDC(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
        }
        if(col_idx == 2){
          g_= internal_render_plotDCD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
        }
      }
    }
    if(DEBUGGING)
      cat("g_class: ",class(g_),"\n")
    
    if(vcount(g_) == 4){
      #plot(c(-1,0,1,0),c(0,1,0,-1))
      my_layout = matrix(c(-1,0,0,1,1,0,0,-1),4,2,byrow = TRUE)
    }
    if(vcount(g_)==3){
      #plot(c(-1,0,1),c(0,1,0))
      my_layout = matrix(c(-1,0,0,0.5,1,0),3,2,byrow = TRUE)
    }
    
    plot(g_,vertex.color = V(g_)$color,
         vertex.size = 50,edge.width = abs(E(g_)$weight)+2,
         vertex.label.color = "black",layout = my_layout)
    legend(x = "bottom",legend = c("Positive Correlation","Negative Correlation"),fill = c("red","darkgreen"))
    
  })
}

genes_data_table_output = function(input,output,MList,MM_list,proxy,graph_s,g,g_geni2){
  output$genes_data_table = DT::renderDataTable({
    type = input$NetworkPattern
    clique_type = input$clique_type
    
    validate(
      need(input$NetworkPattern != "", "Please select a pattern type")
    )
    if(DEBUGGING)
      cat("Il tipo selezionato ??: ",type,"\n")
    type = as.integer(gsub(pattern = "M",x =type,replacement = ""))
    if(DEBUGGING)
      cat("Il numero selezionato ??: ",type,"\n")
    if(DEBUGGING)
      cat("class input$clique_data_table: ",class(input$clique_data_table))
    s = input$clique_data_table_rows_selected
    if(DEBUGGING)
      cat("s vale: ",s,"\n")
    
    if(clique_type == "NDCD"){
      g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "NDD"){
      g_= internal_render_plotNDD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "NDC"){
      g_= internal_render_plotNDC(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "DCD"){
      g_= internal_render_plotDCD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(clique_type == "ALL"){
      if(("" %in% MList[[type]][1,]) == FALSE){
        g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
      }else{
        col_idx = which(MList[[type]][1,] %in% "")
        if(col_idx == 4){
          g_= internal_render_plotNDD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
        }
        if(col_idx == 1){
          g_= internal_render_plotNDC(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
        }
        if(col_idx == 2){
          g_= internal_render_plotDCD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
        }
      }
    }
    if(DEBUGGING){
      cat("g_class: ",class(g_),"\n")
      cat("graph names ",V(g_)$name,"\n")
      
    }
    
    
    validate(
      need(length(s)!=0 , "Please select a clique")
    )          
    
    idx_g = which(V(g)$node_type == "gene")
    GMAT = igraph::get.adjacency(g)
    genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,V(g_)$name]))
    GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]
    idx_gg = which(GG_genes %in% V(g_geni2)$name)
    geni_toPlot = igraph::induced.subgraph(g_geni2,GG_genes[idx_gg])
    
    GENE_INFO = cbind(V(geni_toPlot)$name,V(geni_toPlot)$group)
    colnames(GENE_INFO) = c("Gene Symbol","KEGG Pathway")
    
    for(row_i in 1:dim(GENE_INFO)[1]){
      for(col_j in 1:dim(GENE_INFO)[2]){
        GENE_INFO[row_i,col_j] = paste('<a target="_blank" href=\"https://www.google.com/?q=',GENE_INFO[row_i,col_j],'">',GENE_INFO[row_i,col_j],'</a>',sep="")
      }
    }
    
    DT::datatable(data =  GENE_INFO,
                  options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE,scrollY = "400px", scrollCollapse = TRUE,paging=FALSE),
                  escape=FALSE,
                  selection = "single")
    
  })
}


barplot_pattern_conditional_query = function(input,output,MList,graph_gw){
  output$ggplot = renderPlot({
    type = input$NetworkPattern #type of clique
    triple_type = input$plotTripel
    
    validate(
      need(input$NetworkPattern != "", "Please select a pattern type")
    )
    
    clique_type = input$clique_type
    
    if(clique_type == "ALL"){
      if(DEBUGGING)
        cat("Il tipo selezionato: ",type,"\n")
      i = as.integer(gsub(pattern = "M",x =type,replacement = ""))
      if(DEBUGGING)
        cat("Il numero selezionato: ",i,"\n")
      
      Mi = MList[[i]]
      Mi = as.data.frame(Mi)
      colnames(Mi)=c("Disease","Nano","Drug","Chemical")
      
      if(triple_type==1){
        columns_ = c("Disease","Nano","Drug")
        classes_elem = c("disease","nano","drugs","chemical")
      }
      if(triple_type==2){
        columns_ = c("Disease","Nano","Chemical")
        classes_elem = c("disease","nano","chemical")
        
      }
      if(triple_type==3){
        columns_ = c("Disease","Nano")
        classes_elem = c("disease","nano")
        
      }
      if(triple_type==4){
        columns_ = c("Disease","Drug")
        classes_elem = c("disease","drugs")
        
      }
      if(triple_type==5){
        columns_ = c("Disease","Chemical")
        classes_elem = c("disease","chemical")
        
      }
      if(triple_type==6){
        columns_ = c("Chemical","Nano")
        classes_elem = c("chemical","nano")
        
      }
      if(triple_type==7){
        columns_ = c("Chemical","Drug")
        classes_elem = c("chemical","drugs")
        
      }
      if(triple_type==8){
        columns_ = c("Nano","Drug")
        classes_elem = c("nano","drugs")
      }
      
      col_idx = which(as.matrix(Mi[1,]) %in% "")
      if(length(col_idx)>0){
        if(col_idx == 4){
          validate(
            need(("Chemical" %in% columns_) == FALSE, "Chemical not selected")
          ) 
        }
        if(col_idx == 3){
          validate(
            need(("Drug" %in% columns_) == FALSE, "Drug not selected")
          ) 
        }
        if(col_idx == 2){
          validate(
            need(("Nano" %in% columns_) == FALSE, "Nano not selected")
          ) 
        }
        if(col_idx == 1){
          validate(
            need(("Disease" %in% columns_) == FALSE, "Disease not selected")
          ) 
        }
      }
      
      
      Mi_count = count(Mi, vars=columns_)
      nano_id = unique(Mi_count$Nano)
      drug_id = unique(Mi_count$Drug)
      disease_id = unique(Mi_count$Disease)
      
      Mi_count = Mi_count[order(Mi_count$freq,decreasing = FALSE),]
      nElem_c = length(columns_)
      
      d_id = input$InterestingNodes_items
      d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
      column_index = which(classes_elem %in% d_node_type)
      
      validate(
        need(input$percentuale != "", "Please select a percentage")
      )
      
      if(length(columns_)>2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        dim(Mi_count)[1] -> up_b
        low_b = up_b - nElem
        Mi_count2 = Mi_count$freq[low_b:up_b]
        names(Mi_count2) =  paste(Mi_count[low_b:up_b,indexes_c[1]],Mi_count[low_b:up_b,indexes_c[2]],sep="-")
        
        
        #               Mi_count2 = Mi_count$freq[1:nElem]
        #               names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep=" ")
        
      }
      if(length(columns_) == 2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        dim(Mi_count)[1] -> up_b
        low_b = up_b - nElem
        Mi_count2 = Mi_count$freq[low_b:up_b]
        names(Mi_count2) =  paste(Mi_count[low_b:up_b,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        
        
        #               Mi_count2 = Mi_count$freq[1:nElem]
        #               #names(Mi_count2) =  paste(Mi_count[1:nElem,columns_[2:nElem_c]],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        #               names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        
      }
      
      #index = which(Mi_count[1:nElem,column_index] %in% d_id)
      index = which(Mi_count[low_b:up_b,column_index] %in% d_id)
      
      
      validate(
        need(length(index)>0, paste("No clique for", d_id ,"with this threshold"))
      )
      
      mar.default = c(5,4,4,2) + 0.1
      par(mar = mar.default+ c(0,15,0,0))
      barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
              names.arg=names(Mi_count2)[index], cex.names=1, las=1)
      
      
    }# end if node_type == ALL
    
    
    if(clique_type == "NDCD"){
      if(DEBUGGING)
        cat("Il tipo selezionato: ",type,"\n")
      i = as.integer(gsub(pattern = "M",x =type,replacement = ""))
      if(DEBUGGING)
        cat("Il numero selezionato: ",i,"\n")
      
      Mi = MList[[i]]
      Mi = as.data.frame(Mi)
      colnames(Mi)=c("Disease","Nano","Drug","Chemical")
      
      if(triple_type==1){
        columns_ = c("Disease","Nano","Drug")
        classes_elem = c("disease","nano","drugs","chemical")
      }
      if(triple_type==2){
        columns_ = c("Disease","Nano","Chemical")
        classes_elem = c("disease","nano","chemical")
        
      }
      if(triple_type==3){
        columns_ = c("Disease","Nano")
        classes_elem = c("disease","nano")
        
      }
      if(triple_type==4){
        columns_ = c("Disease","Drug")
        classes_elem = c("disease","drugs")
        
      }
      if(triple_type==5){
        columns_ = c("Disease","Chemical")
        classes_elem = c("disease","chemical")
        
      }
      if(triple_type==6){
        columns_ = c("Chemical","Nano")
        classes_elem = c("chemical","nano")
        
      }
      if(triple_type==7){
        columns_ = c("Chemical","Drug")
        classes_elem = c("chemical","drugs")
        
      }
      if(triple_type==8){
        columns_ = c("Nano","Drug")
        classes_elem = c("nano","drugs")
      }
      
      Mi_count = count(Mi, vars=columns_)
      nano_id = unique(Mi_count$Nano)
      drug_id = unique(Mi_count$Drug)
      disease_id = unique(Mi_count$Disease)
      
      Mi_count = Mi_count[order(Mi_count$freq),]
      nElem_c = length(columns_)
      
      d_id = input$InterestingNodes_items
      d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
      column_index = which(classes_elem %in% d_node_type)
      
      validate(
        need(input$percentuale != "", "Please select a percentage")
      )
      
      if(length(columns_)>2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep=" ")
        
      }
      if(length(columns_) == 2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        #names(Mi_count2) =  paste(Mi_count[1:nElem,columns_[2:nElem_c]],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        
      }
      
      index = which(Mi_count[1:nElem,column_index] %in% d_id)
      
      
      validate(
        need(length(index)>0, paste("No clique for", d_id ,"with this threshold"))
      )
      
      mar.default = c(5,4,4,2) + 0.1
      par(mar = mar.default+ c(0,15,0,0))
      barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
              names.arg=names(Mi_count2)[index], cex.names=1, las=1)
      
      
    }# end if node_type == NDCD
    if(clique_type == "NDD"){ #Nano Drug Disease
      if(DEBUGGING)
        cat("Il tipo selezionato: ",type,"\n")
      i = as.integer(gsub(pattern = "M",x =type,replacement = ""))
      if(DEBUGGING)
        cat("Il numero selezionato: ",i,"\n")
      
      Mi = MList[[i]]
      Mi = as.data.frame(Mi)
      colnames(Mi)=c("Disease","Nano","Drug")
      
      if(triple_type==1){
        columns_ = c("Disease","Nano","Drug")
        classes_elem = c("disease","nano","drugs","chemical")
      }
      if(triple_type==2){
        columns_ = c("Disease","Nano","Chemical")
        classes_elem = c("disease","nano","chemical")
        
      }
      if(triple_type==3){
        columns_ = c("Disease","Nano")
        classes_elem = c("disease","nano")
        
      }
      if(triple_type==4){
        columns_ = c("Disease","Drug")
        classes_elem = c("disease","drugs")
        
      }
      if(triple_type==5){
        columns_ = c("Disease","Chemical")
        classes_elem = c("disease","chemical")
        
      }
      if(triple_type==6){
        columns_ = c("Chemical","Nano")
        classes_elem = c("chemical","nano")
        
      }
      if(triple_type==7){
        columns_ = c("Chemical","Drug")
        classes_elem = c("chemical","drugs")
        
      }
      if(triple_type==8){
        columns_ = c("Nano","Drug")
        classes_elem = c("nano","drugs")
      }
      validate(
        need(("Chemical" %in% columns_) == FALSE, "Chemical not selected")
      )
      
      Mi_count = count(Mi, vars=columns_)
      nano_id = unique(Mi_count$Nano)
      drug_id = unique(Mi_count$Drug)
      disease_id = unique(Mi_count$Disease)
      
      Mi_count = Mi_count[order(Mi_count$freq),]
      nElem_c = length(columns_)
      
      d_id = input$InterestingNodes_items
      d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
      column_index = which(classes_elem %in% d_node_type)
      
      validate(
        need(input$percentuale != "", "Please select a percentage")
      )
      
      if(length(columns_)>2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep=" ")
        
      }
      if(length(columns_) == 2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        
      }
      
      index = which(Mi_count[1:nElem,column_index] %in% d_id)
      
      validate(
        need(length(index)>0, paste("No clique for", d_id ,"with this threshold"))
      )
      
      mar.default = c(5,4,4,2) + 0.1
      par(mar = mar.default+ c(0,15,0,0))
      barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
              names.arg=names(Mi_count2)[index], cex.names=1, las=1)
      
      
    }# end if node_type == NDD
    
    
    if(clique_type == "NDC"){ #Nano Drug Chemical
      if(DEBUGGING)
        cat("Il tipo selezionato: ",type,"\n")
      i = as.integer(gsub(pattern = "M",x =type,replacement = ""))
      if(DEBUGGING)
        cat("Il numero selezionato: ",i,"\n")
      
      Mi = MList[[i]]
      Mi = as.data.frame(Mi)
      colnames(Mi)=c("Chemical","Nano","Drug")
      
      if(triple_type==1){
        columns_ = c("Disease","Nano","Drug")
        classes_elem = c("disease","nano","drugs","chemical")
      }
      if(triple_type==2){
        columns_ = c("Disease","Nano","Chemical")
        classes_elem = c("disease","nano","chemical")
        
      }
      if(triple_type==3){
        columns_ = c("Disease","Nano")
        classes_elem = c("disease","nano")
        
      }
      if(triple_type==4){
        columns_ = c("Disease","Drug")
        classes_elem = c("disease","drugs")
        
      }
      if(triple_type==5){
        columns_ = c("Disease","Chemical")
        classes_elem = c("disease","chemical")
        
      }
      if(triple_type==6){
        columns_ = c("Chemical","Nano")
        classes_elem = c("chemical","nano")
        
      }
      if(triple_type==7){
        columns_ = c("Chemical","Drug")
        classes_elem = c("chemical","drugs")
        
      }
      if(triple_type==8){
        columns_ = c("Nano","Drug")
        classes_elem = c("nano","drugs")
      }
      validate(
        need(("Drug" %in% columns_) == FALSE, "Chemical not selected")
      )
      
      Mi_count = count(Mi, vars=columns_)
      nano_id = unique(Mi_count$Nano)
      drug_id = unique(Mi_count$Drug)
      disease_id = unique(Mi_count$Disease)
      
      Mi_count = Mi_count[order(Mi_count$freq),]
      nElem_c = length(columns_)
      
      d_id = input$InterestingNodes_items
      d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
      column_index = which(classes_elem %in% d_node_type)
      
      validate(
        need(input$percentuale != "", "Please select a percentage")
      )
      
      if(length(columns_)>2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep=" ")
        
      }
      if(length(columns_) == 2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        
      }
      
      index = which(Mi_count[1:nElem,column_index] %in% d_id)
      
      validate(
        need(length(index)>0, paste("No clique for", d_id ,"with this threshold"))
      )
      
      mar.default = c(5,4,4,2) + 0.1
      par(mar = mar.default+ c(0,15,0,0))
      barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
              names.arg=names(Mi_count2)[index], cex.names=1, las=1)
      
      
    }# end if node_type == NDC
    
    if(clique_type == "DCD"){ #Drug Chemical Disease
      if(DEBUGGING)
        cat("Il tipo selezionato: ",type,"\n")
      i = as.integer(gsub(pattern = "M",x =type,replacement = ""))
      if(DEBUGGING)
        cat("Il numero selezionato: ",i,"\n")
      
      Mi = MList[[i]]
      Mi = as.data.frame(Mi)
      colnames(Mi)=c("Disease","Chemical","Drug")
      
      if(triple_type==1){
        columns_ = c("Disease","Nano","Drug")
        classes_elem = c("disease","nano","drugs","chemical")
      }
      if(triple_type==2){
        columns_ = c("Disease","Nano","Chemical")
        classes_elem = c("disease","nano","chemical")
        
      }
      if(triple_type==3){
        columns_ = c("Disease","Nano")
        classes_elem = c("disease","nano")
        
      }
      if(triple_type==4){
        columns_ = c("Disease","Drug")
        classes_elem = c("disease","drugs")
        
      }
      if(triple_type==5){
        columns_ = c("Disease","Chemical")
        classes_elem = c("disease","chemical")
        
      }
      if(triple_type==6){
        columns_ = c("Chemical","Nano")
        classes_elem = c("chemical","nano")
        
      }
      if(triple_type==7){
        columns_ = c("Chemical","Drug")
        classes_elem = c("chemical","drugs")
        
      }
      if(triple_type==8){
        columns_ = c("Nano","Drug")
        classes_elem = c("nano","drugs")
      }
      validate(
        need(("Nano" %in% columns_) == FALSE, "Chemical not selected")
      )
      
      Mi_count = count(Mi, vars=columns_)
      nano_id = unique(Mi_count$Nano)
      drug_id = unique(Mi_count$Drug)
      disease_id = unique(Mi_count$Disease)
      
      Mi_count = Mi_count[order(Mi_count$freq),]
      nElem_c = length(columns_)
      
      d_id = input$InterestingNodes_items
      d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
      column_index = which(classes_elem %in% d_node_type)
      
      validate(
        need(input$percentuale != "", "Please select a percentage")
      )
      
      if(length(columns_)>2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep=" ")
        
      }
      if(length(columns_) == 2){
        indexes_c = 1:(dim(Mi_count)[2]-1)
        indexes_c = indexes_c[-column_index]
        nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
        Mi_count2 = Mi_count$freq[1:nElem]
        names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
        
      }
      
      index = which(Mi_count[1:nElem,column_index] %in% d_id)
      
      validate(
        need(length(index)>0, paste("No clique for", d_id ,"with this threshold"))
      )
      
      mar.default = c(5,4,4,2) + 0.1
      par(mar = mar.default+ c(0,15,0,0))
      barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
              names.arg=names(Mi_count2)[index], cex.names=1, las=1)
      
      
    }# end if node_type == DCD
    
    
  })
}


# output$Subnetwork_plot = renderForceNetwork({ #in go2
#   if(DEBUGGING)
#     cat("Subnetwork plot function\n")
#   ADJ_toPlot = ADJ_S
#   nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#   chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#   drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#   disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#   
#   chemMat_index = which(chemMat[,1] %in% chemical_s)
#   chemMat = chemMat[chemMat_index,]
#   
#   d = input$InterestingNodes
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   validate(
#     need(input$InterestingNodes != "", "SELECT A NODE!")
#   )
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   neig = lapply(X = good_cliques,FUN = function(i_cliques){
#     v_names = names(i_cliques)
#     if(sum(v_names %in% d)>0){
#       return(v_names)
#     }
#   })
#   
#   vds = unique(unlist(neig))
#   
#   edges_type = input$sub_Edges
#   if(DEBUGGING)
#     cat("Edges_type ",edges_type,"\n")
#   type = input$sub_checkGroup
#   
#   validate(
#     need(input$sub_checkGroup != "", "Please select an object group")
#   )
#   if(DEBUGGING)
#     cat("Tipi selezionati",length(type),"\n",type,"\n")
#   ATC_code = input$sub_ATCGroup
#   if(ATC_code == "ALL"){
#     ATC_code = ATC_letter_vector
#   }
#   if(DEBUGGING)
#     cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
#   chem_group = input$sub_ChemicalGroup
#   if(chem_group == "ALL"){
#     chem_group = chemical_group_vector
#   }
#   if(DEBUGGING)
#     cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
#   
#   items = c() 
#   items_type = c()
#   items_lengt = c()
#   #items_color = c()
#   for(i in 1:length(type)){
#     if(type[i]=="nano"){
#       items = c(items,nano_s)
#       items_type = c(items_type,"nano")
#       items_lengt = c(items_lengt,length(nano_s))
#       #items_color = c(items_color,"pink")
#     }
#     if(type[i]=="drugs"){
#       validate(
#         need(input$sub_ATCGroup != "", "Please select an ATC group")
#       )
#       already_selected = c()
#       for(j in 1:length(ATC_code)){
#         ATC_lev1 = substr(join10$code,1,1)
#         ATC_index = which(ATC_lev1 %in% ATC_code[j])
#         new_drugs = unique(join10$name[ATC_index])
#         new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#         index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#         already_selected = new_drugs[index_new_drugs]
#         
#         toRem = which(new_drugs[index_new_drugs] %in% items)
#         if(length(toRem)>0){
#           index_new_drugs = index_new_drugs[-toRem]
#         }
#         
#         items = c(items,new_drugs[index_new_drugs])
#         items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#         items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#       }
#     }
#     if(type[i]=="chemical"){
#       validate(
#         need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#       )
#       
#       already_selected = c()
#       for(j in 1:length(chem_group)){
#         chem_index = which(chemMat[,2] %in% chem_group[j])
#         new_chem = unique(chemMat[chem_index,1])
#         index_new_chem = which((new_chem %in% already_selected)==FALSE)
#         already_selected = new_chem[index_new_chem]
#         
#         toRem = which(new_chem[index_new_chem] %in% items)
#         if(length(toRem)>0){
#           index_new_chem = index_new_chem[-toRem]
#         }
#         
#         items = c(items,new_chem[index_new_chem])
#         items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#         items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#       }
#       
#     }
#     if(type[i]=="disease"){
#       good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#       items = c(items,good_disease)
#       items_type = c(items_type,"disease")
#       items_lengt = c(items_lengt,length(good_disease))
#     }
#   }
#   if(DEBUGGING){
#     cat("item_type: ",items_type,"\n")
#     cat("items_lengt: ",items_lengt,"\n")
#   }
#   if(edges_type == "P"){
#     ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
#   }
#   if(edges_type == "N"){
#     ADJ_toPlot[which(ADJ_toPlot>0)] = 0
#   }
#   
#   if(DEBUGGING)
#     cat("Graph creation... \n")
#   gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot[items,items],mode = "undirected",weighted = TRUE)
#   V(gto_plot)$type = rep(items_type,items_lengt)
#   gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
#   gto_plot = igraph::simplify(gto_plot)
#   
#   if(DEBUGGING)
#     cat("Data frame creation... \n")
#   data_frame = from_igraph_to_data_frame(gto_plot)
#   
#   MyClickScript <- 
#     '      d3.select(this).select("circle").transition()
#   .duration(750)
#   .attr("r", 30)'
#   
#   if(DEBUGGING){
#     cat("Force Networks plot... \n")
#     cat("Repulserad ",input$sub_repulserad,"\n")
#     cat("Edges length: ",input$sub_length,"\n")
#   }
#   forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
#                Source = "source", Target = "target",
#                Value = "value", NodeID = "name",Nodesize="size",
#                Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
#                legend = TRUE,
#                clickAction = MyClickScript,
#                charge = -input$sub_repulserad,
#                linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,
#                                         "}")))
# })
# 
# 
# output$dig_dist = renderPlot({
#   if(DEBUGGING)
#     cat("Subnetwork plot function\n")
#   ADJ_toPlot = ADJ_S
#   nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#   chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#   drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#   disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#   
#   chemMat_index = which(chemMat[,1] %in% chemical_s)
#   chemMat = chemMat[chemMat_index,]
#   
#   d = input$InterestingNodes
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   validate(
#     need(input$InterestingNodes != "", "SELECT A NODE!")
#   )
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   neig = lapply(X = good_cliques,FUN = function(i_cliques){
#     v_names = names(i_cliques)
#     if(sum(v_names %in% d)>0){
#       return(v_names)
#     }
#   })
#   
#   vds = unique(unlist(neig))
#   
#   edges_type = input$sub_Edges
#   if(DEBUGGING)
#     cat("Edges_type ",edges_type,"\n")
#   type = input$sub_checkGroup
#   
#   validate(
#     need(input$sub_checkGroup != "", "Please select an object group")
#   )
#   if(DEBUGGING)
#     cat("Tipi selezionati",length(type),"\n",type,"\n")
#   ATC_code = input$sub_ATCGroup
#   if(ATC_code == "ALL"){
#     ATC_code = ATC_letter_vector
#   }
#   if(DEBUGGING)
#     cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
#   chem_group = input$sub_ChemicalGroup
#   if(chem_group == "ALL"){
#     chem_group = chemical_group_vector
#   }
#   if(DEBUGGING)
#     cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
#   
#   items = c() 
#   items_type = c()
#   items_lengt = c()
#   #items_color = c()
#   for(i in 1:length(type)){
#     if(type[i]=="nano"){
#       items = c(items,nano_s)
#       items_type = c(items_type,"nano")
#       items_lengt = c(items_lengt,length(nano_s))
#       #items_color = c(items_color,"pink")
#     }
#     if(type[i]=="drugs"){
#       validate(
#         need(input$sub_ATCGroup != "", "Please select an ATC group")
#       )
#       already_selected = c()
#       for(j in 1:length(ATC_code)){
#         ATC_lev1 = substr(join10$code,1,1)
#         ATC_index = which(ATC_lev1 %in% ATC_code[j])
#         new_drugs = unique(join10$name[ATC_index])
#         new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#         index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#         already_selected = new_drugs[index_new_drugs]
#         
#         toRem = which(new_drugs[index_new_drugs] %in% items)
#         if(length(toRem)>0){
#           index_new_drugs = index_new_drugs[-toRem]
#         }
#         
#         items = c(items,new_drugs[index_new_drugs])
#         items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#         items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#       }
#     }
#     if(type[i]=="chemical"){
#       validate(
#         need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#       )
#       
#       already_selected = c()
#       for(j in 1:length(chem_group)){
#         chem_index = which(chemMat[,2] %in% chem_group[j])
#         new_chem = unique(chemMat[chem_index,1])
#         index_new_chem = which((new_chem %in% already_selected)==FALSE)
#         already_selected = new_chem[index_new_chem]
#         
#         toRem = which(new_chem[index_new_chem] %in% items)
#         if(length(toRem)>0){
#           index_new_chem = index_new_chem[-toRem]
#         }
#         
#         items = c(items,new_chem[index_new_chem])
#         items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#         items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#       }
#       
#     }
#     if(type[i]=="disease"){
#       good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#       items = c(items,good_disease)
#       items_type = c(items_type,"disease")
#       items_lengt = c(items_lengt,length(good_disease))
#     }
#   }
#   if(DEBUGGING){
#     cat("item_type: ",items_type,"\n")
#     cat("items_lengt: ",items_lengt,"\n")
#   }
#   if(edges_type == "P"){
#     ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
#   }
#   if(edges_type == "N"){
#     ADJ_toPlot[which(ADJ_toPlot>0)] = 0
#   }
#   
#   if(DEBUGGING)
#     cat("Graph creation... \n")
#   gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot[items,items],mode = "undirected",weighted = TRUE)
#   V(gto_plot)$type = rep(items_type,items_lengt)
#   gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
#   gto_plot = igraph::simplify(gto_plot)
#   
#   degree_dist = range01(igraph::degree(gto_plot))
#   closeness_dist = range01(igraph::closeness(gto_plot,weights = abs(E(gto_plot)$weight)))
#   betweenness_dist = range01(igraph::betweenness(gto_plot,weights = abs(E(gto_plot)$weight)))
#   eccentricity_dist = range01(igraph::eccentricity(gto_plot))
#   
#   y_m = max(degree_dist,closeness_dist,betweenness_dist,eccentricity_dist)
#   
#   names(degree_dist) = V(gto_plot)$name
#   
#   plot(degree_dist, type="l", col="red", ylim = c(0,y_m), ylab = "",xlab="",lwd= 4)
#   lines(closeness_dist, col="green",lwd= 4)
#   lines(betweenness_dist, col="yellow",lwd= 4)
#   lines(eccentricity_dist, col = "purple",lwd= 4)
#   legend(x = "top",legend = c("Degree","Closeness","Betwenness","Eccentricity"),fill=c("red","green","yellow","purple"))
# })
# 
# output$Subnetwork_plot_statistic = renderPlot({
#   if(DEBUGGING)
#     cat("STATISTICHE SOTTONETWORK...\n")
#   ADJ_toPlot = ADJ_S
#   nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#   chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#   drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#   disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#   
#   chemMat_index = which(chemMat[,1] %in% chemical_s)
#   chemMat = chemMat[chemMat_index,]
#   
#   d = input$InterestingNodes
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   validate(
#     need(input$InterestingNodes != "", "SELECT A NODE!")
#   )
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   neig = lapply(X = good_cliques,FUN = function(i_cliques){
#     v_names = names(i_cliques)
#     if(sum(v_names %in% d)>0){
#       return(v_names)
#     }
#   })
#   
#   vds = unique(unlist(neig))
#   
#   edges_type = input$sub_Edges
#   if(DEBUGGING)
#     cat("Edges_type ",edges_type,"\n")
#   type = input$sub_checkGroup
#   if(DEBUGGING)
#     cat("type: ",type,"\n")
#   
#   validate(
#     need(input$sub_checkGroup != "", "Please select an object group")
#   )
#   if(DEBUGGING)
#     cat("Tipi selezionati",length(type),"\n",type,"\n")
#   ATC_code = input$sub_ATCGroup
#   if(DEBUGGING)
#     cat("ATC_code ",ATC_code,"\n")
#   
#   if(ATC_code == "ALL"){
#     ATC_code = ATC_letter_vector
#   }
#   if(DEBUGGING)
#     cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
#   
#   chem_group = input$sub_ChemicalGroup
#   if(chem_group == "ALL"){
#     chem_group = chemical_group_vector
#   }
#   if(DEBUGGING)
#     cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
#   
#   items = c() 
#   items_type = c()
#   items_lengt = c()
#   #items_color = c()
#   for(i in 1:length(type)){
#     if(type[i]=="nano"){
#       items = c(items,nano_s)
#       items_type = c(items_type,"nano")
#       items_lengt = c(items_lengt,length(nano_s))
#       #items_color = c(items_color,"pink")
#     }
#     if(type[i]=="drugs"){
#       validate(
#         need(input$sub_ATCGroup != "", "Please select an ATC group")
#       )
#       already_selected = c()
#       for(j in 1:length(ATC_code)){
#         ATC_lev1 = substr(join10$code,1,1)
#         ATC_index = which(ATC_lev1 %in% ATC_code[j])
#         new_drugs = unique(join10$name[ATC_index])
#         new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#         index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#         already_selected = new_drugs[index_new_drugs]
#         
#         toRem = which(new_drugs[index_new_drugs] %in% items)
#         if(length(toRem)>0){
#           index_new_drugs = index_new_drugs[-toRem]
#         }
#         
#         items = c(items,new_drugs[index_new_drugs])
#         items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#         items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#       }
#     }
#     if(type[i]=="chemical"){
#       validate(
#         need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#       )
#       
#       already_selected = c()
#       for(j in 1:length(chem_group)){
#         chem_index = which(chemMat[,2] %in% chem_group[j])
#         new_chem = unique(chemMat[chem_index,1])
#         index_new_chem = which((new_chem %in% already_selected)==FALSE)
#         already_selected = new_chem[index_new_chem]
#         
#         toRem = which(new_chem[index_new_chem] %in% items)
#         if(length(toRem)>0){
#           index_new_chem = index_new_chem[-toRem]
#         }
#         
#         items = c(items,new_chem[index_new_chem])
#         items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#         items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#       }
#       
#     }
#     if(type[i]=="disease"){
#       good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#       items = c(items,good_disease)
#       items_type = c(items_type,"disease")
#       items_lengt = c(items_lengt,length(good_disease))
#     }
#   }
#   
#   slices = items_lengt
#   lbls = items_type
#   pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(items_type))),
#       radius = 1, cex = 0.7)
#   
# })
# 
# 
# output$gene_Subnetwork_plot = renderForceNetwork({
#   if(DEBUGGING)
#     cat("Subnetwork plot\n")
#   ADJ_toPlot = ADJ_S
#   nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#   chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#   drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#   disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#   
#   chemMat_index = which(chemMat[,1] %in% chemical_s)
#   chemMat = chemMat[chemMat_index,]
#   
#   d = input$InterestingNodes
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   validate(
#     need(input$InterestingNodes != "", "SELECT A NODE!")
#   )
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   neig = lapply(X = good_cliques,FUN = function(i_cliques){
#     v_names = names(i_cliques)
#     if(sum(v_names %in% d)>0){
#       return(v_names)
#     }
#   })
#   
#   vds = unique(unlist(neig))
#   
#   edges_type = input$sub_Edges
#   if(DEBUGGING)
#     cat("Edges_type ",edges_type,"\n")
#   type = input$sub_checkGroup
#   
#   validate(
#     need(input$sub_checkGroup != "", "Please select an object group")
#   )
#   
#   ATC_code = input$sub_ATCGroup
#   if(ATC_code == "ALL"){
#     ATC_code = ATC_letter_vector
#   }
#   chem_group = input$sub_ChemicalGroup
#   if(chem_group == "ALL"){
#     chem_group = chemical_group_vector
#   }
#   
#   items = c() 
#   items_type = c()
#   items_lengt = c()
#   #items_color = c()
#   for(i in 1:length(type)){
#     if(type[i]=="nano"){
#       items = c(items,nano_s)
#       items_type = c(items_type,"nano")
#       items_lengt = c(items_lengt,length(nano_s))
#       #items_color = c(items_color,"pink")
#     }
#     if(type[i]=="drugs"){
#       validate(
#         need(input$sub_ATCGroup != "", "Please select an ATC group")
#       )
#       already_selected = c()
#       for(j in 1:length(ATC_code)){
#         ATC_lev1 = substr(join10$code,1,1)
#         ATC_index = which(ATC_lev1 %in% ATC_code[j])
#         new_drugs = unique(join10$name[ATC_index])
#         new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#         index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#         already_selected = new_drugs[index_new_drugs]
#         
#         toRem = which(new_drugs[index_new_drugs] %in% items)
#         if(length(toRem)>0){
#           index_new_drugs = index_new_drugs[-toRem]
#         }
#         
#         items = c(items,new_drugs[index_new_drugs])
#         items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#         items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#       }
#     }
#     if(type[i]=="chemical"){
#       validate(
#         need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#       )
#       
#       already_selected = c()
#       for(j in 1:length(chem_group)){
#         chem_index = which(chemMat[,2] %in% chem_group[j])
#         new_chem = unique(chemMat[chem_index,1])
#         index_new_chem = which((new_chem %in% already_selected)==FALSE)
#         already_selected = new_chem[index_new_chem]
#         
#         toRem = which(new_chem[index_new_chem] %in% items)
#         if(length(toRem)>0){
#           index_new_chem = index_new_chem[-toRem]
#         }
#         items = c(items,new_chem[index_new_chem])
#         items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#         items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#       }
#     }
#     if(type[i]=="disease"){
#       good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#       items = c(items,good_disease)
#       items_type = c(items_type,"disease")
#       items_lengt = c(items_lengt,length(good_disease))
#     }
#   }
#   
#   idx_g = which(V(g)$node_type == "gene")
#   GMAT = igraph::get.adjacency(g)
#   genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,items]))
#   GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]
#   idx_gg = which(GG_genes %in% V(g_geni2)$name)
#   geni_toPlot = igraph::induced.subgraph(g_geni2,GG_genes[idx_gg])
#   geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
#   
#   validate(
#     need(vcount(geni_toPlot)>0, "Empty Network")
#   )
#   
#   #         gto_plot = igraph::induced.subgraph(graph = g,vids = c(V(g)$name[idx_g],items))
#   #         gto_plot = igraph::delete.vertices(gto_plot,items)
#   #         idx_gg = which(V(gto_plot)$name %in% V(g_geni2)$name)
#   #         geni_toPlot = igraph::induced.subgraph(g_geni2,V(gto_plot)$name[idx_gg])
#   #         geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
#   
#   data_frame = get.data.frame(x = geni_toPlot,what = "both")
#   edges = data_frame$edges
#   edges$value = round(abs(edges$weight * 10),digits = 0)
#   colnames(edges) = c("source","target","weight","value")
#   edges$value = edges$value + 0.2
#   vertices = data_frame$vertices
#   vertices$size = igraph::degree(geni_toPlot)
#   colnames(vertices) = c("name","group","type","size")
#   
#   if(dim(edges)[1]>0){
#     for(i in 1:dim(edges)[1]){
#       edges[i,"source"] = which(vertices[,"name"] %in% edges[i,"source"]) - 1
#       edges[i,"target"] = which(vertices[,"name"] %in% edges[i,"target"]) - 1
#     }  
#   }
#   
#   vertices$name = as.factor(vertices$name)
#   vertices$group = as.factor(vertices$group)
#   vertices$size = as.numeric(vertices$size)
#   vertices$type = as.factor(vertices$type)
#   edges$source = as.integer(edges$source)
#   edges$target  = as.integer(edges$target)
#   edges$value = as.integer(edges$value)
#   
#   MyClickScript <- 
#     '      d3.select(this).select("circle").transition()
#   .duration(750)
#   .attr("r", 30)'
#   
#   forceNetwork(Links = edges, Nodes = vertices,
#                Source = "source", Target = "target",
#                Value = "value", NodeID = "name",Nodesize="size",
#                Group = "group",zoom = TRUE,opacity = 0.95,fontSize = 8,
#                legend = TRUE,
#                clickAction = MyClickScript,
#                charge = -input$sub_repulserad,
#                linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,
#                                         "}")))
# })
# 
# output$gene_Subnetwork_plot_statistics = renderPlot({
#   if(DEBUGGING)
#     cat("Subnetwork plot\n")
#   ADJ_toPlot = ADJ_S
#   nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#   chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#   drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#   disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#   
#   chemMat_index = which(chemMat[,1] %in% chemical_s)
#   chemMat = chemMat[chemMat_index,]
#   
#   d = input$InterestingNodes
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   validate(
#     need(input$InterestingNodes != "", "SELECT A NODE!")
#   )
#   if(DEBUGGING)
#     cat("INTERESTING NODES: ",d,"\n")
#   
#   neig = lapply(X = good_cliques,FUN = function(i_cliques){
#     v_names = names(i_cliques)
#     if(sum(v_names %in% d)>0){
#       return(v_names)
#     }
#   })
#   
#   vds = unique(unlist(neig))
#   
#   edges_type = input$sub_Edges
#   if(DEBUGGING)
#     cat("Edges_type ",edges_type,"\n")
#   type = input$sub_checkGroup
#   
#   validate(
#     need(input$sub_checkGroup != "", "Please select an object group")
#   )
#   
#   ATC_code = input$sub_ATCGroup
#   if(ATC_code == "ALL"){
#     ATC_code = ATC_letter_vector
#   }
#   chem_group = input$sub_ChemicalGroup
#   if(chem_group == "ALL"){
#     chem_group = chemical_group_vector
#   }
#   
#   items = c() 
#   items_type = c()
#   items_lengt = c()
#   #items_color = c()
#   for(i in 1:length(type)){
#     if(type[i]=="nano"){
#       items = c(items,nano_s)
#       items_type = c(items_type,"nano")
#       items_lengt = c(items_lengt,length(nano_s))
#       #items_color = c(items_color,"pink")
#     }
#     if(type[i]=="drugs"){
#       validate(
#         need(input$sub_ATCGroup != "", "Please select an ATC group")
#       )
#       already_selected = c()
#       for(j in 1:length(ATC_code)){
#         ATC_lev1 = substr(join10$code,1,1)
#         ATC_index = which(ATC_lev1 %in% ATC_code[j])
#         new_drugs = unique(join10$name[ATC_index])
#         new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#         index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#         already_selected = new_drugs[index_new_drugs]
#         
#         toRem = which(new_drugs[index_new_drugs] %in% items)
#         if(length(toRem)>0){
#           index_new_drugs = index_new_drugs[-toRem]
#         }
#         
#         items = c(items,new_drugs[index_new_drugs])
#         items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#         items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#       }
#     }
#     if(type[i]=="chemical"){
#       validate(
#         need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#       )
#       
#       already_selected = c()
#       for(j in 1:length(chem_group)){
#         chem_index = which(chemMat[,2] %in% chem_group[j])
#         new_chem = unique(chemMat[chem_index,1])
#         index_new_chem = which((new_chem %in% already_selected)==FALSE)
#         already_selected = new_chem[index_new_chem]
#         
#         toRem = which(new_chem[index_new_chem] %in% items)
#         if(length(toRem)>0){
#           index_new_chem = index_new_chem[-toRem]
#         }
#         items = c(items,new_chem[index_new_chem])
#         items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#         items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#       }
#     }
#     if(type[i]=="disease"){
#       good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#       items = c(items,good_disease)
#       items_type = c(items_type,"disease")
#       items_lengt = c(items_lengt,length(good_disease))
#     }
#   }
#   
#   idx_g = which(V(g)$node_type == "gene")
#   GMAT = igraph::get.adjacency(g)
#   genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,items]))
#   GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]
#   idx_gg = which(GG_genes %in% V(g_geni2)$name)
#   geni_toPlot = igraph::induced.subgraph(g_geni2,GG_genes[idx_gg])
#   geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
#   
#   validate(
#     need(vcount(geni_toPlot)>0, "Empty Network")
#   )
#   
#   
#   data_frame = get.data.frame(x = geni_toPlot,what = "both")
#   edges = data_frame$edges
#   edges$value = round(abs(edges$weight * 10),digits = 0)
#   colnames(edges) = c("source","target","weight","value")
#   
#   vertices = data_frame$vertices
#   vertices$size = igraph::degree(geni_toPlot)
#   colnames(vertices) = c("name","group","type","size")
#   
#   for(i in 1:dim(edges)[1]){
#     edges[i,"source"] = which(vertices[,"name"] %in% edges[i,"source"]) - 1
#     edges[i,"target"] = which(vertices[,"name"] %in% edges[i,"target"]) - 1
#   }
#   
#   vertices$name = as.factor(vertices$name)
#   vertices$group = as.factor(vertices$group)
#   vertices$size = as.numeric(vertices$size)
#   vertices$type = as.factor(vertices$type)
#   
#   edges$source = as.integer(edges$source)
#   edges$target  = as.integer(edges$target)
#   edges$value = as.integer(edges$value)
#   
#   
#   
#   
#   slices = table(vertices$group)
#   lbls = names(table(vertices$group))
#   #           df = data.frame(as.vector(slices),lbls)
#   #           colnames(df) = c("slices","KEGG")
#   # p <- ggplot(df, aes(x=1, y=slices, fill=KEGG)) +
#   #   ggtitle("Network Composition") +
#   #   # black border around pie slices
#   #   geom_bar(stat="identity", color='black') +
#   #   # remove black diagonal line from legend
#   #   guides(fill=guide_legend(override.aes=list(colour=NA))) +
#   #   # polar coordinates
#   #   coord_polar(theta='y') +
#   #   # label aesthetics
#   #   theme(axis.ticks=element_blank(),  # the axis ticks
#   #         axis.title=element_blank(),  # the axis labels
#   #         axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
#   #         axis.text.x=element_text(color='black')) +
#   #   # pie slice labels
#   #   scale_y_continuous(
#   #     breaks=cumsum(df$slices) - df$slices/2,
#   #     labels=df$KEGG
#   #   )
#   # 
#   # print(p)
#   
#   
#   #radial.pie(slices,labels=lbls,main="Network statistics")
#   
#   
#   pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(slices)),
#       radius = 1,cex = 0.7,las=1)
#   
# })



# output$drugs_chemical_DT = DT::renderDataTable({
#   type = input$NetworkPattern
#   clique_type = input$clique_type
#   
#   validate(
#     need(input$NetworkPattern != "", "Please select a pattern type")
#   )
#   if(DEBUGGING)
#     cat("Il tipo selezionato ??: ",type,"\n")
#   type = as.integer(gsub(pattern = "M",x =type,replacement = ""))
#   if(DEBUGGING)
#     cat("Il numero selezionato ??: ",type,"\n")
#   if(DEBUGGING)
#     cat("class input$clique_data_table: ",class(input$clique_data_table))
#   s = input$clique_data_table_rows_selected
#   if(DEBUGGING)
#     cat("s vale: ",s,"\n")
#   
#   if(clique_type == "NDCD"){
#     g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#   }
#   if(clique_type == "NDD"){
#     g_= internal_render_plotNDD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#   }
#   if(clique_type == "NDC"){
#     g_= internal_render_plotNDC(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#   }
#   if(clique_type == "DCD"){
#     g_= internal_render_plotDCD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#   }
#   if(clique_type == "ALL"){
#     if(("" %in% MList[[type]][1,]) == FALSE){
#       g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#     }else{
#       col_idx = which(MList[[type]][1,] %in% "")
#       if(col_idx == 4){
#         g_= internal_render_plotNDD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#       }
#       if(col_idx == 1){
#         g_= internal_render_plotNDC(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#       }
#       if(col_idx == 2){
#         g_= internal_render_plotDCD(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
#       }
#     }
#   }
#   if(DEBUGGING){
#     cat("g_class: ",class(g_),"\n")
#     cat("graph names ",V(g_)$name,"\n")
#     
#   }
#   
#   
#   if(DEBUGGING){
#     cat("Name ---> ",V(g_)$name,"\n")
#   }
#   
#   validate(
#     need(length(s)!=0 , "Please select a clique")
#   )   
#   
#   if(length(V(g_)$name) == 3){
#     cat("3 VERTICI --> \n")
#     
#     cli_drug = V(g_)$name[which(V(g_)$name %in% drugs)]  
#     idx_drug = which(join10$name == cli_drug)
#     idx_drug_letter = which(ATC_letter_vector == join10[idx_drug,3])
#     ATC_l1 = names(ATC_choice_list)[idx_drug_letter + 1]
#     drug_row = c(cli_drug,ATC_l1)
#     INFO = rbind(drug_row,c("",""))
#     
#   }
#   
#   if(length(V(g_)$name) == 4){
#     cat("4 VERTICI ---> \n")
#     
#     cli_drug = V(g_)$name[which(V(g_)$name %in% drugs)]  
#     idx_drug = which(join10$name == cli_drug)
#     idx_drug_letter = which(ATC_letter_vector == join10[idx_drug,3])
#     ATC_l1 = names(ATC_choice_list)[idx_drug_letter + 1]
#     drug_row = c(cli_drug,ATC_l1)
#     
#     
#     cli_chem = V(g_)$name[which(V(g_)$name %in% chemical)]  
#     idx_chem = which(chemMat[,1] == cli_chem) 
#     chem_row = chemMat[idx_chem,]
#     
#     INFO = rbind(drug_row,chem_row)
#     
#   }
#   
#   
#   colnames(INFO) = c("Phenotypic Entity","Class")
#   rownames(INFO) = NULL
#   DT::datatable(data =  INFO,
#                 options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE,scrollY = "400px", scrollCollapse = TRUE,paging=FALSE),
#                 escape=FALSE,
#                 selection = "single")
#   
# })   