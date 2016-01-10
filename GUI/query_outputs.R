barplot_patter_output = function(input,output,MList,graph_gw){
  output$ggplot = renderPlot({
    type = input$NetworkPattern #type of clique
    triple_type = input$plotTripel
    validate(
      need(input$percentuale != "", "Please select a percentage")
    )
    validate(
      need(input$NetworkPattern != "", "Please select a pattern type")
    )
    
    i = as.integer(gsub(pattern = "M",x =type,replacement = ""))
    
    if(DEBUGGING){
      cat("Il tipo selezionato: ",type,"\n")
      cat("Il numero selezionato: ",i,"\n")
    }
    
    Mi = MList[[i]]
    Mi = as.data.frame(Mi)
    colnames(Mi)=c("Disease","Nano","Drug","Chemical")
    
    if(triple_type==1){
      columns_ = c("Disease","Nano","Drug")
    }
    if(triple_type==2){
      columns_ = c("Disease","Nano","Chemical")
    }
    if(triple_type==3){
      columns_ = c("Disease","Nano")
    }
    if(triple_type==4){
      columns_ = c("Disease","Drug")
    }
    if(triple_type==5){
      columns_ = c("Disease","Chemical")
    }
    
    Mi_count = count(Mi, vars=columns_)
    nano_id = unique(Mi_count$Nano)
    drug_id = unique(Mi_count$Drug)
    disease_id = unique(Mi_count$Disease)
    
    Mi_count = Mi_count[order(Mi_count$freq,decreasing = FALSE),]
    
    if(DEBUGGING)
      cat("Mi_count FREQUENCIES ",Mi_count$freq,"\n")
    
    nElem_c = length(columns_)
    
    d_id = input$InterestingNodes_items
    d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
    column_index = which(c("disease","nano","drugs","chemical") %in% d_node_type)
    
    
    if(DEBUGGING){
      cat("INPUT PERCENTUALE CAMBIATA: ",input$percentuale,"\n")
    }
    
    if(length(columns_)>2){
      indexes_c = 1:(dim(Mi_count)[2]-1)
      indexes_c = indexes_c[-column_index]
      nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
      dim(Mi_count)[1] -> up_b
      low_b = up_b - nElem
      Mi_count2 = Mi_count$freq[low_b:up_b]
      names(Mi_count2) =  paste(Mi_count[low_b:up_b,indexes_c[1]],Mi_count[low_b:up_b,indexes_c[2]],sep="-")
      #Mi_count2 = Mi_count$freq[1:nElem]
      #names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep="-")
    }
    if(length(columns_) == 2){
      indexes_c = 1:(dim(Mi_count)[2]-1)
      indexes_c = indexes_c[-column_index]
      nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
      dim(Mi_count)[1] -> up_b
      low_b = up_b - nElem
      Mi_count2 = Mi_count$freq[low_b:up_b]
      names(Mi_count2) =  paste(Mi_count[low_b:up_b,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
      
      #Mi_count2 = Mi_count$freq[1:nElem]
      #names(Mi_count2) =  paste(Mi_count[1:nElem,columns_[2:nElem_c]],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
    }
    
    
    #index = which(Mi_count[1:nElem,column_index] %in% d_id)
    index = which(Mi_count[low_b:up_b,column_index] %in% d_id)
    
    
    validate(
      need(length(index)>0, "No cliques for the selected node!")
    )
    
    #par(mfrow = c(1,length(disease_id)))
    #for(d_id in disease_id){
    mar.default = c(5,4,4,2) + 0.1
    par(mar = mar.default+ c(0,15,0,0))
    barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
            names.arg=names(Mi_count2)[index], cex.names=1, las=1)
    #}
  })
}

plot_clique_graph = function(input,output,MM_list,graph_s,proxy){
  output$xx = renderPlot({
    s = input$clique_data_table_rows_selected
    type = input$NetworkPattern
    
    validate(
      need(input$NetworkPattern != "", "Please select a pattern type")
    )
    
    type = as.integer(gsub(pattern = "M",x =type,replacement = ""))
    
    
    if(DEBUGGING){
      cat("Il tipo selezionato : ",type,"\n")
      cat("Il numero selezionato : ",type,"\n")
      cat("class input$clique_data_table: ",class(input$clique_data_table))
    }
    
    #s = input$clique_data_table_rows_selected
    if(DEBUGGING){
      cat("s = input$clique_data_table_rows_selected: ",input$clique_data_table_rows_selected,"\n")
    }
    g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    
    if(DEBUGGING){
      cat("s vale: ",s,"\n")
      cat("g_class: ",class(g_),"\n")
    }
    
    plot(g_,vertex.color = V(g_)$color,
         vertex.size = 50,edge.width = abs(E(g_)$weight)+2,
         vertex.label.color = "black")
    legend(x = "bottom",legend = c("Positive Correlation","Negative Correlation"),fill = c("red","darkgreen"))
    
  })
}

plot_force_based_subnetwork_query_resutls = function(input,output,ADJ_S,chemMat,good_cliques,join10){
  output$Subnetwork_plot = renderForceNetwork({ 
    if(DEBUGGING)
      cat("Subnetwork plot function\n")
    ADJ_toPlot = ADJ_S
    nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
    chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
    drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
    disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
    
    chemMat_index = which(chemMat[,1] %in% chemical_s)
    chemMat = chemMat[chemMat_index,]
    
    d = input$InterestingNodes
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    validate(
      need(input$InterestingNodes != "", "SELECT A NODE!")
    )
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    neig = lapply(X = good_cliques,FUN = function(i_cliques){
      v_names = names(i_cliques)
      if(sum(v_names %in% d)>0){
        return(v_names)
      }
    })
    
    vds = unique(unlist(neig))
    
    edges_type = input$sub_Edges
    if(DEBUGGING)
      cat("Edges_type ",edges_type,"\n")
    type = input$sub_checkGroup
    
    validate(
      need(input$sub_checkGroup != "", "Please select an object group")
    )
    if(DEBUGGING)
      cat("Tipi selezionati",length(type),"\n",type,"\n")
    ATC_code = input$sub_ATCGroup
    if(ATC_code == "ALL"){
      ATC_code = ATC_letter_vector
    }
    if(DEBUGGING)
      cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
    chem_group = input$sub_ChemicalGroup
    if(chem_group == "ALL"){
      chem_group = chemical_group_vector
    }
    if(DEBUGGING)
      cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
    
    items = c() 
    items_type = c()
    items_lengt = c()
    #items_color = c()
    for(i in 1:length(type)){
      if(type[i]=="nano"){
        items = c(items,nano_s)
        items_type = c(items_type,"nano")
        items_lengt = c(items_lengt,length(nano_s))
        #items_color = c(items_color,"pink")
      }
      if(type[i]=="drugs"){
        validate(
          need(input$sub_ATCGroup != "", "Please select an ATC group")
        )
        already_selected = c()
        for(j in 1:length(ATC_code)){
          ATC_lev1 = substr(join10$code,1,1)
          ATC_index = which(ATC_lev1 %in% ATC_code[j])
          new_drugs = unique(join10$name[ATC_index])
          new_drugs = new_drugs[which(new_drugs %in% drug_s)]
          index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
          already_selected = new_drugs[index_new_drugs]
          
          toRem = which(new_drugs[index_new_drugs] %in% items)
          if(length(toRem)>0){
            index_new_drugs = index_new_drugs[-toRem]
          }
          
          items = c(items,new_drugs[index_new_drugs])
          items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
          items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
        }
      }
      if(type[i]=="chemical"){
        validate(
          need(input$sub_ChemicalGroup != "", "Please select a chemical group")
        )
        
        already_selected = c()
        for(j in 1:length(chem_group)){
          chem_index = which(chemMat[,2] %in% chem_group[j])
          new_chem = unique(chemMat[chem_index,1])
          index_new_chem = which((new_chem %in% already_selected)==FALSE)
          already_selected = new_chem[index_new_chem]
          
          toRem = which(new_chem[index_new_chem] %in% items)
          if(length(toRem)>0){
            index_new_chem = index_new_chem[-toRem]
          }
          
          items = c(items,new_chem[index_new_chem])
          items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
          items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
        }
        
      }
      if(type[i]=="disease"){
        good_disease = disease[disease %in% colnames(ADJ_toPlot)]
        items = c(items,good_disease)
        items_type = c(items_type,"disease")
        items_lengt = c(items_lengt,length(good_disease))
      }
    }
    if(DEBUGGING){
      cat("item_type: ",items_type,"\n")
      cat("items_lengt: ",items_lengt,"\n")
    }
    if(edges_type == "P"){
      ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
    }
    if(edges_type == "N"){
      ADJ_toPlot[which(ADJ_toPlot>0)] = 0
    }
    
    if(DEBUGGING)
      cat("Graph creation... \n")
    
    gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot[items,items],mode = "undirected",weighted = TRUE)
    V(gto_plot)$type = rep(items_type,items_lengt)
    gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
    gto_plot = igraph::simplify(gto_plot)
    
    if(DEBUGGING)
      cat("Data frame creation... \n")
    
    data_frame = from_igraph_to_data_frame(gto_plot)
    
    MyClickScript <- 
      '      d3.select(this).select("circle").transition()
  .duration(750)
  .attr("r", 30)'
    
    if(DEBUGGING){
      cat("Force Networks plot... \n")
      cat("Repulserad ",input$sub_repulserad,"\n")
      cat("Edges length: ",input$sub_length,"\n")
    }
    forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
                 Source = "source", Target = "target",
                 Value = "value", NodeID = "name",Nodesize="size",
                 Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
                 legend = TRUE,
                 clickAction = MyClickScript,
                 charge = -input$sub_repulserad,
                 linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,
                                          "}")))
  })
}

plot_subnetwork_statistics = function(input,output,ADJ_S,chemMat,good_cliques,join10){
  output$Subnetwork_plot_statistic = renderPlot({
    if(DEBUGGING)
      cat("STATISTICHE SOTTONETWORK...\n")
    ADJ_toPlot = ADJ_S
    nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
    chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
    drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
    disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
    
    chemMat_index = which(chemMat[,1] %in% chemical_s)
    chemMat = chemMat[chemMat_index,]
    
    d = input$InterestingNodes
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    validate(
      need(input$InterestingNodes != "", "SELECT A NODE!")
    )
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    neig = lapply(X = good_cliques,FUN = function(i_cliques){
      v_names = names(i_cliques)
      if(sum(v_names %in% d)>0){
        return(v_names)
      }
    })
    
    vds = unique(unlist(neig))
    
    edges_type = input$sub_Edges
    if(DEBUGGING)
      cat("Edges_type ",edges_type,"\n")
    type = input$sub_checkGroup
    if(DEBUGGING)
      cat("type: ",type,"\n")
    
    validate(
      need(input$sub_checkGroup != "", "Please select an object group")
    )
    if(DEBUGGING)
      cat("Tipi selezionati",length(type),"\n",type,"\n")
    ATC_code = input$sub_ATCGroup
    if(DEBUGGING)
      cat("ATC_code ",ATC_code,"\n")
    
    if(ATC_code == "ALL"){
      ATC_code = ATC_letter_vector
    }
    if(DEBUGGING)
      cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
    
    chem_group = input$sub_ChemicalGroup
    if(chem_group == "ALL"){
      chem_group = chemical_group_vector
    }
    if(DEBUGGING)
      cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
    
    items = c() 
    items_type = c()
    items_lengt = c()
    #items_color = c()
    for(i in 1:length(type)){
      if(type[i]=="nano"){
        items = c(items,nano_s)
        items_type = c(items_type,"nano")
        items_lengt = c(items_lengt,length(nano_s))
        #items_color = c(items_color,"pink")
      }
      if(type[i]=="drugs"){
        validate(
          need(input$sub_ATCGroup != "", "Please select an ATC group")
        )
        already_selected = c()
        for(j in 1:length(ATC_code)){
          ATC_lev1 = substr(join10$code,1,1)
          ATC_index = which(ATC_lev1 %in% ATC_code[j])
          new_drugs = unique(join10$name[ATC_index])
          new_drugs = new_drugs[which(new_drugs %in% drug_s)]
          index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
          already_selected = new_drugs[index_new_drugs]
          
          toRem = which(new_drugs[index_new_drugs] %in% items)
          if(length(toRem)>0){
            index_new_drugs = index_new_drugs[-toRem]
          }
          
          items = c(items,new_drugs[index_new_drugs])
          items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
          items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
        }
      }
      if(type[i]=="chemical"){
        validate(
          need(input$sub_ChemicalGroup != "", "Please select a chemical group")
        )
        
        already_selected = c()
        for(j in 1:length(chem_group)){
          chem_index = which(chemMat[,2] %in% chem_group[j])
          new_chem = unique(chemMat[chem_index,1])
          index_new_chem = which((new_chem %in% already_selected)==FALSE)
          already_selected = new_chem[index_new_chem]
          
          toRem = which(new_chem[index_new_chem] %in% items)
          if(length(toRem)>0){
            index_new_chem = index_new_chem[-toRem]
          }
          
          items = c(items,new_chem[index_new_chem])
          items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
          items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
        }
        
      }
      if(type[i]=="disease"){
        good_disease = disease[disease %in% colnames(ADJ_toPlot)]
        items = c(items,good_disease)
        items_type = c(items_type,"disease")
        items_lengt = c(items_lengt,length(good_disease))
      }
    }
    
    slices = items_lengt
    lbls = items_type
    pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(items_type))),
        radius = 1, cex = 0.7)
    
  })
}

plot_gene_subnetwork = function(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2){
  output$gene_Subnetwork_plot = renderForceNetwork({
    if(DEBUGGING)
      cat("Subnetwork plot\n")
    ADJ_toPlot = ADJ_S
    nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
    chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
    drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
    disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
    
    chemMat_index = which(chemMat[,1] %in% chemical_s)
    chemMat = chemMat[chemMat_index,]
    
    d = input$InterestingNodes
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    validate(
      need(input$InterestingNodes != "", "SELECT A NODE!")
    )
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    neig = lapply(X = good_cliques,FUN = function(i_cliques){
      v_names = names(i_cliques)
      if(sum(v_names %in% d)>0){
        return(v_names)
      }
    })
    
    vds = unique(unlist(neig))
    
    edges_type = input$sub_Edges
    if(DEBUGGING)
      cat("Edges_type ",edges_type,"\n")
    type = input$sub_checkGroup
    
    validate(
      need(input$sub_checkGroup != "", "Please select an object group")
    )
    
    ATC_code = input$sub_ATCGroup
    if(ATC_code == "ALL"){
      ATC_code = ATC_letter_vector
    }
    chem_group = input$sub_ChemicalGroup
    if(chem_group == "ALL"){
      chem_group = chemical_group_vector
    }
    
    items = c() 
    items_type = c()
    items_lengt = c()
    #items_color = c()
    for(i in 1:length(type)){
      if(type[i]=="nano"){
        items = c(items,nano_s)
        items_type = c(items_type,"nano")
        items_lengt = c(items_lengt,length(nano_s))
        #items_color = c(items_color,"pink")
      }
      if(type[i]=="drugs"){
        validate(
          need(input$sub_ATCGroup != "", "Please select an ATC group")
        )
        already_selected = c()
        for(j in 1:length(ATC_code)){
          ATC_lev1 = substr(join10$code,1,1)
          ATC_index = which(ATC_lev1 %in% ATC_code[j])
          new_drugs = unique(join10$name[ATC_index])
          new_drugs = new_drugs[which(new_drugs %in% drug_s)]
          index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
          already_selected = new_drugs[index_new_drugs]
          
          toRem = which(new_drugs[index_new_drugs] %in% items)
          if(length(toRem)>0){
            index_new_drugs = index_new_drugs[-toRem]
          }
          
          items = c(items,new_drugs[index_new_drugs])
          items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
          items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
        }
      }
      if(type[i]=="chemical"){
        validate(
          need(input$sub_ChemicalGroup != "", "Please select a chemical group")
        )
        
        already_selected = c()
        for(j in 1:length(chem_group)){
          chem_index = which(chemMat[,2] %in% chem_group[j])
          new_chem = unique(chemMat[chem_index,1])
          index_new_chem = which((new_chem %in% already_selected)==FALSE)
          already_selected = new_chem[index_new_chem]
          
          toRem = which(new_chem[index_new_chem] %in% items)
          if(length(toRem)>0){
            index_new_chem = index_new_chem[-toRem]
          }
          items = c(items,new_chem[index_new_chem])
          items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
          items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
        }
      }
      if(type[i]=="disease"){
        good_disease = disease[disease %in% colnames(ADJ_toPlot)]
        items = c(items,good_disease)
        items_type = c(items_type,"disease")
        items_lengt = c(items_lengt,length(good_disease))
      }
    }
    
    idx_g = which(V(g)$node_type == "gene")
    GMAT = igraph::get.adjacency(g)
    genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,items]))
    GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]
    idx_gg = which(GG_genes %in% V(g_geni2)$name)
    geni_toPlot = igraph::induced.subgraph(g_geni2,GG_genes[idx_gg])
    geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
    
    validate(
      need(vcount(geni_toPlot)>0, "Empty Network")
    )
    
    #         gto_plot = igraph::induced.subgraph(graph = g,vids = c(V(g)$name[idx_g],items))
    #         gto_plot = igraph::delete.vertices(gto_plot,items)
    #         idx_gg = which(V(gto_plot)$name %in% V(g_geni2)$name)
    #         geni_toPlot = igraph::induced.subgraph(g_geni2,V(gto_plot)$name[idx_gg])
    #         geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
    
    data_frame = get.data.frame(x = geni_toPlot,what = "both")
    edges = data_frame$edges
    edges$value = round(abs(edges$weight * 10),digits = 0)
    colnames(edges) = c("source","target","weight","value")
    edges$value = edges$value + 0.2
    vertices = data_frame$vertices
    vertices$size = igraph::degree(geni_toPlot)
    colnames(vertices) = c("name","group","type","size")
    
    if(dim(edges)[1]>0){
      for(i in 1:dim(edges)[1]){
        edges[i,"source"] = which(vertices[,"name"] %in% edges[i,"source"]) - 1
        edges[i,"target"] = which(vertices[,"name"] %in% edges[i,"target"]) - 1
      }  
    }
    
    vertices$name = as.factor(vertices$name)
    vertices$group = as.factor(vertices$group)
    vertices$size = as.numeric(vertices$size)
    vertices$type = as.factor(vertices$type)
    edges$source = as.integer(edges$source)
    edges$target  = as.integer(edges$target)
    edges$value = as.integer(edges$value)
    
    MyClickScript <- 
      '      d3.select(this).select("circle").transition()
    .duration(750)
    .attr("r", 30)'
    
    forceNetwork(Links = edges, Nodes = vertices,
                 Source = "source", Target = "target",
                 Value = "value", NodeID = "name",Nodesize="size",
                 Group = "group",zoom = TRUE,opacity = 0.95,fontSize = 8,
                 legend = TRUE,
                 clickAction = MyClickScript,
                 charge = -input$sub_repulserad,
                 linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,
                                          "}")))
  })      
}

plot_gene_subnetwork_statistics = function(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2){
  output$gene_Subnetwork_plot_statistics = renderPlot({
    if(DEBUGGING)
      cat("Subnetwork plot\n")
    ADJ_toPlot = ADJ_S
    nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
    chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
    drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
    disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
    
    chemMat_index = which(chemMat[,1] %in% chemical_s)
    chemMat = chemMat[chemMat_index,]
    
    d = input$InterestingNodes
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    validate(
      need(input$InterestingNodes != "", "SELECT A NODE!")
    )
    if(DEBUGGING)
      cat("INTERESTING NODES: ",d,"\n")
    
    neig = lapply(X = good_cliques,FUN = function(i_cliques){
      v_names = names(i_cliques)
      if(sum(v_names %in% d)>0){
        return(v_names)
      }
    })
    
    vds = unique(unlist(neig))
    
    edges_type = input$sub_Edges
    if(DEBUGGING)
      cat("Edges_type ",edges_type,"\n")
    type = input$sub_checkGroup
    
    validate(
      need(input$sub_checkGroup != "", "Please select an object group")
    )
    
    ATC_code = input$sub_ATCGroup
    if(ATC_code == "ALL"){
      ATC_code = ATC_letter_vector
    }
    chem_group = input$sub_ChemicalGroup
    if(chem_group == "ALL"){
      chem_group = chemical_group_vector
    }
    
    items = c() 
    items_type = c()
    items_lengt = c()
    #items_color = c()
    for(i in 1:length(type)){
      if(type[i]=="nano"){
        items = c(items,nano_s)
        items_type = c(items_type,"nano")
        items_lengt = c(items_lengt,length(nano_s))
        #items_color = c(items_color,"pink")
      }
      if(type[i]=="drugs"){
        validate(
          need(input$sub_ATCGroup != "", "Please select an ATC group")
        )
        already_selected = c()
        for(j in 1:length(ATC_code)){
          ATC_lev1 = substr(join10$code,1,1)
          ATC_index = which(ATC_lev1 %in% ATC_code[j])
          new_drugs = unique(join10$name[ATC_index])
          new_drugs = new_drugs[which(new_drugs %in% drug_s)]
          index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
          already_selected = new_drugs[index_new_drugs]
          
          toRem = which(new_drugs[index_new_drugs] %in% items)
          if(length(toRem)>0){
            index_new_drugs = index_new_drugs[-toRem]
          }
          
          items = c(items,new_drugs[index_new_drugs])
          items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
          items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
        }
      }
      if(type[i]=="chemical"){
        validate(
          need(input$sub_ChemicalGroup != "", "Please select a chemical group")
        )
        
        already_selected = c()
        for(j in 1:length(chem_group)){
          chem_index = which(chemMat[,2] %in% chem_group[j])
          new_chem = unique(chemMat[chem_index,1])
          index_new_chem = which((new_chem %in% already_selected)==FALSE)
          already_selected = new_chem[index_new_chem]
          
          toRem = which(new_chem[index_new_chem] %in% items)
          if(length(toRem)>0){
            index_new_chem = index_new_chem[-toRem]
          }
          items = c(items,new_chem[index_new_chem])
          items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
          items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
        }
      }
      if(type[i]=="disease"){
        good_disease = disease[disease %in% colnames(ADJ_toPlot)]
        items = c(items,good_disease)
        items_type = c(items_type,"disease")
        items_lengt = c(items_lengt,length(good_disease))
      }
    }
    
    idx_g = which(V(g)$node_type == "gene")
    GMAT = igraph::get.adjacency(g)
    genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,items]))
    GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]
    idx_gg = which(GG_genes %in% V(g_geni2)$name)
    geni_toPlot = igraph::induced.subgraph(g_geni2,GG_genes[idx_gg])
    geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
    
    validate(
      need(vcount(geni_toPlot)>0, "Empty Network")
    )
    
    
    data_frame = get.data.frame(x = geni_toPlot,what = "both")
    edges = data_frame$edges
    edges$value = round(abs(edges$weight * 10),digits = 0)
    colnames(edges) = c("source","target","weight","value")
    
    vertices = data_frame$vertices
    vertices$size = igraph::degree(geni_toPlot)
    colnames(vertices) = c("name","group","type","size")
    
    for(i in 1:dim(edges)[1]){
      edges[i,"source"] = which(vertices[,"name"] %in% edges[i,"source"]) - 1
      edges[i,"target"] = which(vertices[,"name"] %in% edges[i,"target"]) - 1
    }
    
    vertices$name = as.factor(vertices$name)
    vertices$group = as.factor(vertices$group)
    vertices$size = as.numeric(vertices$size)
    vertices$type = as.factor(vertices$type)
    
    edges$source = as.integer(edges$source)
    edges$target  = as.integer(edges$target)
    edges$value = as.integer(edges$value)
    
    slices = table(vertices$group)
    lbls = names(table(vertices$group))  
    pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(slices)),
        radius = 1,cex = 0.7,las=1)
    
  })
}



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
