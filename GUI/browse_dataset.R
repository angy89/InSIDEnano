plot_kegg_path_network = function(input,output,KEGG_ADJ){
  output$KEGG_Network = renderForceNetwork({
    
    rownames(KEGG_ADJ) = colnames(KEGG_ADJ) = gsub(x = rownames(KEGG_ADJ),pattern = "KEGG_",replacement = "")
    
    th = input$KEGG_Similarity
    th = quantile(KEGG_ADJ,th/100)
    KEGG_ADJ[KEGG_ADJ<th] = 0
    
    toRem =which(colSums(KEGG_ADJ)==0)
    if(length(toRem)>0){
      KEGG_ADJ = KEGG_ADJ[-toRem,-toRem]
    }
    
    KEGG_path_graph = graph.adjacency(adjmatrix = KEGG_ADJ,mode = "undirected",weighted = TRUE)
    V(KEGG_path_graph)$size = igraph::degree(KEGG_path_graph)
    V(KEGG_path_graph)$group = "KEGG Pathways"
    
    KEGG_path_data_frame = igraph::get.data.frame(KEGG_path_graph,what = "both")
    
    vertices = KEGG_path_data_frame$vertices
    edges = KEGG_path_data_frame$edges
    
    edges$from = match(edges$from,vertices$name) - 1
    edges$to = match(edges$to,vertices$name) - 1
    
    vertices$name = as.factor(vertices$name)
    vertices$group = as.factor(vertices$group)
    vertices$size = as.numeric(vertices$size)
    edges$from = as.integer(edges$from)
    edges$to  = as.integer(edges$to)

    forceNetwork(Links = edges, Nodes = vertices,
                 Source = "from", Target = "to",
                 Value = "weight", NodeID = "name",Nodesize="size",
                 zoom = TRUE,opacity = 0.85,fontSize = 10,Group = "group",
                 legend = TRUE, 
                # clickAction = JS(Shiny.onInputChange("selected_KEGG", d.name)),
                 charge = -input$KEGG_repulserad,
                 linkDistance = JS(paste0("function(d){return d.value*",input$KEGG_length,"}")))
  })
}

plot_gene_network = function(input,output,g,g_geni2){
  output$geneNetwork = renderForceNetwork({
    
    groups_path = input$Patway_g
    validate(
      need(input$Patway_g != "", "Please select a Pathway")
    )
    
    if(DEBUGGING)
      cat("Number of groups ",length(groups_path),"\n")
    
   
    #ADJ_g = igraph::as_adjacency_matrix(g_geni2,attr = "weight",type = "both")
    
    good_index = which(V(g_geni2)$group %in% groups_path)  
    geni_toPlot = igraph::induced.subgraph(graph = g_geni2,vids = V(g_geni2)$name[good_index])
    geni_toPlot = igraph::delete.vertices(graph = geni_toPlot,v = which(igraph::degree(geni_toPlot)<1))
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
    
    MyClickScript <- 'd3.select(this).select("circle").transition().duration(750).attr("r", 30)'
    
    MyClickScript2 <- 'd3.select(this).select(dLink.target).select("circle").transition().duration(750).attr("r",30),
    d3.select(dLink.target).select("circle").transition().duration(750).attr("r",30)'
    
    forceNetwork(Links = edges, Nodes = vertices,
                 Source = "source", Target = "target",
                 Value = "value", NodeID = "name",Nodesize="size",
                 zoom = TRUE,opacity = 0.85,fontSize = 10,Group = "group",
                 legend = TRUE, height = input$gene_width,width =input$gene_height, 
                 clickAction = MyClickScript,charge = -input$gene_repulserad,
                 linkDistance = JS(paste0("function(d){return d.value*",input$gene_length,"}")))
  })
  
  
}



plot_item_network_pie = function(input,output,W2_ADJ){
  output$pie_chart <- renderPlot({
    ADJ_toPlot = W2_ADJ
    edges_type = input$Edges
    
    type = input$checkGroup
    ATC_code = input$ATCGroup
    
    if(ATC_code == "ALL"){
      ATC_code = ATC_letter_vector
    }
    
    chem_group = input$ChemicalGroup
    
    if(chem_group == "ALL"){
      chem_group = chemical_group_vector
    }
    if(DEBUGGING)
      cat(type," ",length(type)," ",ATC_code," ",chem_group,"\n")
    
    validate(
      need(input$checkGroup != "", "Please select an object group")
    )
    
    items = c()
    items_type = c()
    items_lengt = c()
    #items_color = c()
    for(i in 1:length(type)){
      if(type[i]=="nano"){
        items = c(items,nano)
        items_type = c(items_type,"nano")
        items_lengt = c(items_lengt,length(nano))
      }
      if(type[i]=="drugs"){
        #         validate(
        #           need(input$ATCGroup != "", "Please select an ATC group")
        #         )
        already_selected = c()
        for(j in 1:length(ATC_code)){
          ATC_lev1 = substr(join10$code,1,1)
          ATC_index = which(ATC_lev1 %in% ATC_code[j])
          new_drugs = unique(join10$name[ATC_index])
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
        #         validate(
        #           need(input$ChemicalGroup != "", "Please select a chemical group")
        #         )
        
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
        validate(
          need(input$input_dis != "", "Please select a disease.")
        )
        
        good_disease_to_plot = input$input_dis
        if(good_disease_to_plot =="ALL"){
          good_disease_to_plot = good_disease
        }
        if(DEBUGGING)
          cat("Good disease: ",good_disease_to_plot,"\n")
        
        #good_disease = disease[disease %in% colnames(ADJ_toPlot)]
        items = c(items,good_disease_to_plot)
        items_type = c(items_type,"disease")
        items_lengt = c(items_lengt,length(good_disease_to_plot))
      }
    }
    
    th = input$slider1
    
    pos_index = which(ADJ_toPlot>0)
    neg_index = which(ADJ_toPlot<0)
    
    th_p_val = quantile(ADJ_toPlot[pos_index],th/100)
    th_n_val = quantile(ADJ_toPlot[neg_index],1-(th/100))
    
    ADJ_toPlot[which(ADJ_toPlot>0 & ADJ_toPlot<th_p_val)] = 0
    ADJ_toPlot[which(ADJ_toPlot<0 & ADJ_toPlot>th_n_val)] = 0
    
    if(edges_type == "P"){
      ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
    }
    if(edges_type == "N"){
      ADJ_toPlot[which(ADJ_toPlot>0)] = 0
    }
    if(DEBUGGING)
      cat("Edges_type: ",edges_type,"\n")
    
    gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot[items,items],mode = "undirected",weighted = TRUE)
    V(gto_plot)$type = rep(items_type,items_lengt)
    gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
    gto_plot = igraph::simplify(gto_plot)
    
    slices = table(V(gto_plot)$type)
    lbls = names(table(V(gto_plot)$type))
    #slices = items_lengt
    #  lbls = items_type
    pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(items_type))))
    
  })
}

plot_item_network = function(input,output,W2_ADJ){
  output$itemNetwork = renderForceNetwork({
    ADJ_toPlot = W2_ADJ
    
    edges_type = input$Edges
    
    type = input$checkGroup
    ATC_code = input$ATCGroup
    
    if(ATC_code == "ALL"){
      ATC_code = ATC_letter_vector
    }
    
    chem_group = input$ChemicalGroup
    if(chem_group == "ALL"){
      chem_group = chemical_group_vector
    }
    
    if(DEBUGGING)
      cat(type," ",length(type)," ",ATC_code," ",chem_group,"\n")
    
    validate(
      need(input$checkGroup != "", "Please select an object group")
      #need(length(input$checkGroup)>1, "Please select a couple of object groups")
      
    )
    
    items = c()
    items_type = c()
    items_lengt = c()
    #items_color = c()
    for(i in 1:length(type)){
      if(type[i]=="nano"){
        items = c(items,nano)
        items_type = c(items_type,"nano")
        items_lengt = c(items_lengt,length(nano))
      }
      if(type[i]=="drugs"){
        validate(
          need(input$ATCGroup != "", "Please select an ATC group")
        )
        already_selected = c()
        for(j in 1:length(ATC_code)){
          ATC_lev1 = substr(join10$code,1,1)
          ATC_index = which(ATC_lev1 %in% ATC_code[j])
          new_drugs = unique(join10$name[ATC_index])
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
          need(input$ChemicalGroup != "", "Please select a chemical group")
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
        validate(
          need(input$input_dis != "", "Please select a disease.")
        )
        
        good_disease_to_plot = input$input_dis
        if("ALL" %in% good_disease_to_plot){
          good_disease_to_plot = disease[disease %in% colnames(ADJ_toPlot)]
        }
        if(DEBUGGING)
          cat("Good disease: ",good_disease_to_plot,"\n")
        
        items = c(items,good_disease_to_plot)
        items_type = c(items_type,"disease")
        items_lengt = c(items_lengt,length(good_disease_to_plot))
      }
    }
    
    th = input$slider1
    
    pos_index = which(ADJ_toPlot>0)
    neg_index = which(ADJ_toPlot<0)
    
    th_p_val = quantile(ADJ_toPlot[pos_index],th/100)
    th_n_val = quantile(ADJ_toPlot[neg_index],1-(th/100))
    
    ADJ_toPlot[which(ADJ_toPlot>0 & ADJ_toPlot<th_p_val)] = 0
    ADJ_toPlot[which(ADJ_toPlot<0 & ADJ_toPlot>th_n_val)] = 0
    
    if(edges_type == "P"){
      ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
    }
    if(edges_type == "N"){
      ADJ_toPlot[which(ADJ_toPlot>0)] = 0
    }
    if(DEBUGGING)
      cat("Edges_type: ",edges_type,"\n")
    
    gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot[items,items],mode = "undirected",weighted = TRUE)
    V(gto_plot)$type = rep(items_type,items_lengt)
    gto_plot = igraph::delete.vertices(gto_plot,which(igraph::degree(gto_plot)<1))
    gto_plot = igraph::simplify(gto_plot)
    
    data_frame = from_igraph_to_data_frame(gto_plot)
    MyClickScript <- 
      '      d3.select(this).select("circle").transition()
    .duration(750)
    .attr("r", 30)'
    
    MyClickScript5<- 
      '            d3.select(this).select("circle").attr("opacity", 1);
    '
    
    MyClickScript3 <- 'd3.select(this).node().__data__;
    node.style("opacity", function (o) {
    return neighboring(d, o) | neighboring(o, d) ? 1 : 0.1;
    });'
    
    MyClickScript6 <- 'd3.fisheye.circular()
    .radius(120)'
    
    MyClickScript4 <- 'alert("You clicked " + d.name + " which is in row " + (d.index + 1) + " degree " + d.weight + 
    " type " + d.type + 
    " of your original R data frame");'
    if(DEBUGGING)
      cat("Render force network...\n")
    
    forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
                 Source = "source", Target = "target",
                 Value = "value", NodeID = "name",Nodesize="size",
                 Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
                 #radiusCalculation = JS("d.nodesize + 6"),
                 legend = TRUE,# height = my_screen_h,width =my_screen_w/2, 
                 clickAction = MyClickScript3,charge = -input$repulserad,
                 linkDistance = JS(paste0("function(d){return d.value*",input$length,"}")))
    
  })
}