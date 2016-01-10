select_node_query = function(input,output,disease_list,selected_nodes){
  for(i in input$disease){
    disease_list[[i]] = i
    selected_nodes = c(selected_nodes,i)
  }
  if(DEBUGGING){
    cat("selected_nodes", selected_nodes, "\n")  
    cat("disease_list", length(disease_list), "\n")  
  }
  return(list(disease_list = disease_list,selected_nodes=selected_nodes))
}

apply_thresholds = function(W_ADJ,THS){
  if(DEBUGGING) cat("Apply threshold function \n")
  W2 = W_ADJ
  W_ADJ[nano,disease][which(W_ADJ[nano,disease]>0 & W_ADJ[nano,disease]<THS$th_ndis_p)] = 0 #nano disease
  W_ADJ[nano,disease][which(W_ADJ[nano,disease]<0 & W_ADJ[nano,disease]>THS$th_ndis_n)] = 0
  
  W_ADJ[nano,drugs][which(W_ADJ[nano,drugs]>0 & W_ADJ[nano,drugs]<THS$th_nd_p)] = 0 #nano drugs
  W_ADJ[nano,drugs][which(W_ADJ[nano,drugs]<0 & W_ADJ[nano,drugs]>THS$th_nd_n)] = 0
  
  W_ADJ[disease,drugs][which(W_ADJ[disease,drugs]>0 & W_ADJ[disease,drugs]<THS$th_dd_p)] = 0 #disease drugs
  W_ADJ[disease,drugs][which(W_ADJ[disease,drugs]<0 & W_ADJ[disease,drugs]>THS$th_dd_n)] = 0
  
  W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]>0 & W_ADJ[nano,chemical]<THS$th_nc_p)] = 0 #nano chemical
  W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]<0 & W_ADJ[nano,chemical]>THS$th_nc_n)] = 0
  
  W_ADJ[drugs,chemical][which(W_ADJ[drugs,chemical]>0 & W_ADJ[drugs,chemical]<THS$th_dd_p)] = 0 #drugs chemical
  W_ADJ[drugs,chemical][which(W_ADJ[drugs,chemical]<0 & W_ADJ[drugs,chemical]>THS$th_dd_n)] = 0
  
  W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]>0 & W_ADJ[nano,chemical]<THS$th_drc_p)] = 0 #nano chemical
  W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]<0 & W_ADJ[nano,chemical]>THS$th_drc_n)] = 0
  
  W_ADJ[disease,chemical][which(W_ADJ[disease,chemical]>0 & W_ADJ[disease,chemical]<THS$th_dc_p)] = 0 #disease chemical
  W_ADJ[disease,chemical][which(W_ADJ[disease,chemical]<0 & W_ADJ[disease,chemical]>THS$th_dc_n)] = 0
  
  W_ADJ[nano,nano] = 0
  W_ADJ[drugs,drugs] = 0
  W_ADJ[disease,disease] = 0
  W_ADJ[chemical,chemical] = 0
  
  W_ADJ[selected_nodes,selected_nodes] = W2[selected_nodes,selected_nodes]
  return(W_ADJ)
}

find_thresholds = function(W_ADJ,th_p){
  if(DEBUGGING) cat("Find threshold function \n")
  
  diag(W_ADJ) = 0
  th_n = 1-th_p
  
  if(DEBUGGING){
    cat("dim(W_ADJ)", dim(W_ADJ), "\n")
    cat("th ", th_p, "\n")
  }
  
  W = W_ADJ[nano,disease]
  th_ndis_p = quantile(W[which(W>0)],th_p) #threshold per nano disease positiva
  th_ndis_n = quantile(W[which(W<0)],th_n) #threshold per nano disease negativa
  W = W_ADJ[nano,drugs]
  th_nd_p = quantile(W[which(W>0)],th_p)#threshold per nano drugs positiva
  th_nd_n = quantile(W[which(W<0)],th_n)#threshold per nano drugs negativa
  W = W_ADJ[disease,drugs]
  th_dd_p = quantile(W[which(W>0)],th_p)#threshold per disease drugs positiva
  th_dd_n = quantile(W[which(W<0)],th_n)#threshold per disease drugs negativa
  W = W_ADJ[nano,chemical]
  th_nc_p = quantile(W[which(W>0)],th_p)#threshold per nano chemical positiva
  th_nc_n = quantile(W[which(W<0)],th_n)#threshold per nano chemical negativa
  W = W_ADJ[drugs,chemical]
  th_drc_p = quantile(W[which(W>0)],th_p)#threshold per drugs chemical positiva
  th_drc_n = quantile(W[which(W<0)],th_n)#threshold per drugs chemical negativa
  W = W_ADJ[disease,chemical]
  th_dc_p = quantile(W[which(W>0)],th_p)#threshold per disease chemical positiva
  th_dc_n = quantile(W[which(W<0)],th_n)#threshold per disease chemical negativa
  
  THS = list(th_ndis_p=th_ndis_p,
             th_ndis_n=th_ndis_n,
             th_nd_p=th_nd_p,
             th_nd_n=th_nd_n,
             th_dd_p=th_dd_p,
             th_dd_n=th_dd_n,
             th_nc_p=th_nc_p,
             th_nc_n=th_nc_n,
             th_drc_p=th_drc_p,
             th_drc_n=th_drc_n,
             th_dc_p=th_dc_p,
             th_dc_n=th_dc_n)
  
  return(THS)
}

creating_graph = function(W_ADJ,node_type){

  if(DEBUGGING)
    cat("creating graph \n")
  
  graph_gw = graph.adjacency(W_ADJ,mode="undirected",weighted=TRUE)
  
  if(DEBUGGING)
    cat(class(graph_gw)," ", ecount(graph_gw), " ", vcount(graph_gw), "\n")
  
  V(graph_gw)$type = node_type
  graph_gw = igraph::delete.vertices(graph_gw,which(igraph::degree(graph_gw)<1))
  graph_gw = igraph::simplify(graph_gw)
  return(graph_gw)
}

subgraph_selected_nodes_free_query = function(graph_gw,selected_nodes){
  if(DEBUGGING){
    cat("selected_nodes", selected_nodes, "\n")
    cat(class(graph_gw)," ", ecount(graph_gw), " ", vcount(graph_gw), "\n")
  }
  nn = igraph::neighborhood(graph = graph_gw,order = 1,nodes = selected_nodes,mode = "all")
  
  if(DEBUGGING)
    cat("length(nn) ",length(nn), "\n")
  
  nn_v = c()
  for(i in 1:length(nn)){
    nn_v=c(nn_v,names(nn[[i]]))
  }
  nn_vv = unique(nn_v)
  
  graph_s = igraph::induced.subgraph(graph = graph_gw,vids = c(selected_nodes,nn_vv))
  
  if(DEBUGGING){
    cat("GRAPH_S ", vcount(graph_s), ecount(graph_s))
    cat("Graph_gw and Graph_s classes", class(graph_gw),class(graph_s),"\n")
  }
  
  index_nano = which(V(graph_s)$name %in% nano)
  index_drug = which(V(graph_s)$name %in% drugs)
  index_chem = which(V(graph_s)$name %in% chemical)
  index_dis = which(V(graph_s)$name %in% disease)
  
  tipes = rep("nano",length(V(graph_s)$name))
  tipes[index_drug] = "drugs"
  tipes[index_chem] = "chemical"
  tipes[index_dis] = "disease"
  
  col = rep("pink",length(V(graph_s)$name))
  col[index_drug] = "yellow"
  col[index_chem] = "violet"
  col[index_dis] = "skyblue"
  
  V(graph_s)$type = tipes
  V(graph_s)$color = col
  graph_s = igraph::set.edge.attribute(graph_s,name = "color",index = E(graph_s)[which(sign(E(graph_s)$weight)<0)],value = "darkgreen")
  graph_s = igraph::set.edge.attribute(graph_s,name = "color",index = E(graph_s)[which(sign(E(graph_s)$weight)>0)],value = "red")
  
  ADJ_S = get.adjacency(graph_s,attr = "weight")
  ADJ_S = as.matrix(ADJ_S)
  
  return(list(graph_s=graph_s,ADJ_S=ADJ_S))
}

cliques_search = function(graph_s,min=4,max = 4){
  nano_ADJ = V(graph_s)$name[V(graph_s)$name %in% nano]
  drug_ADJ = V(graph_s)$name[V(graph_s)$name %in% drugs]
  dis_ADJ = V(graph_s)$name[V(graph_s)$name %in% disease]
  chem_ADJ = V(graph_s)$name[V(graph_s)$name %in% chemical]
  
  if(DEBUGGING){
    cat("NDDC_clique_search function!\n")
    cat("NELEM--> ",vcount(graph_s),"\n")
    
    cat("nano_ADJ--> ",length(nano_ADJ),"\n")
    cat("drug_ADJ--> ",length(drug_ADJ),"\n")
    cat("dis_ADJ--> ",length(dis_ADJ),"\n")
    cat("chem_ADJ--> ",length(chem_ADJ),"\n")
    
    
  }
  mcl = cliques(graph=graph_s, min=min, max=max)
  
  if(DEBUGGING)
    cat("Nro cliques: ",length(mcl),"\n")
  
  return(mcl)
}

cliques_evaluation = function(mcl){
  is_good = lapply(X = mcl,FUN = function(obj){
    vertici = names(obj)
    ni = di = ddis = dc = 0
    
    ni = sum(vertici %in% nano)
    if(ni == 1){
      di = sum(vertici %in% drugs)
      if(di == 1){
        ddis = sum(vertici %in% disease)
        if(ddis == 1){
          dc = sum(vertici %in% chemical)
          if(dc==1){
            return(TRUE)
          }
        }
      }
    }
    
    return(FALSE)
    
  })
  is_good = unlist(is_good)
  sum(is_good)
  if(DEBUGGING){
    cat("Number of good cliques: ",sum(is_good),"\n")
  }
  idx = which(is_good==TRUE)
  good_cliques = mcl[idx]
  return(good_cliques)
}

cliques_clustering = function(good_cliques,ADJ_S){
  a = lapply(X = good_cliques,FUN = function(obj){
    vertices = names(obj)
    v_nano = vertices[which(vertices %in% nano)]
    v_drug = vertices[which(vertices %in% drugs)]
    v_chem = vertices[which(vertices %in% chemical)]
    v_dis = vertices[which(vertices %in% disease)]
    
    intersect(v_chem,v_drug) -> ii
    if(length(ii)>0){
      which(v_chem %in% v_drug) -> index_ii
      v_chem = v_chem[-index_ii]
    }
    
    row = sign(c(ADJ_S[v_nano,v_dis],
                 ADJ_S[v_nano,v_chem],
                 ADJ_S[v_nano,v_drug],
                 ADJ_S[v_drug,v_dis],
                 ADJ_S[v_drug,v_chem],
                 ADJ_S[v_dis,v_chem]))
  })
  
  MAT = do.call(rbind, a)
  uniqueMAT = unique(MAT)
  
  cliques_groups = lapply(1:nrow(uniqueMAT),function(i){
    unlist(lapply(1:nrow(MAT),function(j){
      if(sum(uniqueMAT[i,]!=MAT[j,])==0){
        j
      }
    }))
  })
  
  if(DEBUGGING)
    cat("Nro of cliques groups: ",length(cliques_groups),"\n")
  
  return(cliques_groups)
}

group_cliques_list = function(cliques_groups,good_cliques){
  
  MList = lapply(cliques_groups,FUN=function(obj){
    idx = obj
    good_cliques_i = good_cliques[idx]
    vertices_list = lapply(good_cliques_i,FUN = names)
    ord_vertices = lapply(vertices_list,FUN = function(vertices){
      v_nano = vertices[which(vertices %in% nano)]
      v_drug = vertices[which(vertices %in% drugs)]
      v_chem = vertices[which(vertices %in% chemical)]
      v_dis = vertices[which(vertices %in% disease)]
      
      intersect(v_chem,v_drug) -> ii
      if(length(ii)>0){
        which(v_chem %in% v_drug) -> index_ii
        v_chem = v_chem[-index_ii]
      }
      
      c(v_dis,v_nano,v_drug,v_chem)
    })
    M = do.call(rbind,ord_vertices)
  })
  
  return(MList)
}

clickable_cliques_list = function(MList,cliques_groups){
  MM_list = list()
  nType = length(cliques_groups)
  
  for(i in 1:nType){
    Mi = MList[[i]]
    
    Mi_link = apply(X = Mi,MARGIN = 1,FUN = function(row_i){
      #xx = paste('<a target="_blank",href=\"https://www.google.com/?q=',row_i,' class=','"btn btn-primary"','>',row_i,'</a>',sep="")
      xx = paste('<a target="_blank" href=\"https://www.google.com/?q=',row_i,'">',row_i,'</a>',sep="")
      
      new_row = c(xx[1],xx[2],xx[3],xx[4])
    })
    Mi = t(Mi_link)
    rownames(Mi) = 1:dim(Mi)[1]
    Mi = as.data.frame(Mi)
    colnames(Mi)=c("Disease","Nano","Drug","Chemical")
    MM_list[[i]] = Mi
    #M_output_list[[i]] = DT::renderDataTable(Mi,options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE),escape = FALSE)
  }
  
  return(MM_list)
}

conditional_query_nodes = function(input,output){
  xx = paste(input$nano_input,input$drug_input, input$chemical_input, input$disease_input,sep="")
  if(DEBUGGING){
    cat("concatenazione: ",xx,"\n")
    cat("length(xx): ",length(xx),"\n")
  }
  if(length(xx)==0){
    output$info2_1 <- renderUI({
      HTML("Please insert at least one object for the query!")
    }) 
    validate(need(length(xx)>0, "Please insert at least one object for the query!"))
  }
  
  nano_query = input$nano_input
  if(length(nano_query) !=0 ){
    if(nano_query=="ALL"){
      nano_query = nano
    }
  }
  
  drug_query = input$drug_input
  if(length(drug_query) !=0 ){
    if(drug_query=="ALL"){
      drug_query = drugs
    }
    if(drug_query=="A" || drug_query=="C" || drug_query=="D" || 
       drug_query=="G" || drug_query=="H" || drug_query=="J" || 
       drug_query=="L" || drug_query=="M" || drug_query=="N" ||
       drug_query=="P" || drug_query=="R" || drug_query=="S"){
      index_D= which(join10$ATC_lev1 ==drug_query)
      drug_query = unique(join10[index_D,1])
    }
  }
  
  chemical_query = input$chemical_input
  if(length(chemical_query) != 0){
    if(chemical_query=="ALL"){
      chemical_query = chemical
    }
    if(chemical_query %in% names(table(chemMat[,2]))){
      index_C = which(chemMat[,2] == chemical_query)
      chemical_query = unique(chemMat[index_C,1])
    }
  }
  
  disease_query = input$disease_input
  if(length(disease_query)!=0){
    if(disease_query=="ALL"){
      disease_query = disease
    }
  }
  query_nodes = c(nano_query,drug_query,chemical_query,disease_query)
  for(i in query_nodes){
    disease_list[[i]] = i
    selected_nodes = c(selected_nodes,i)
  }
  
  return(list(query_nodes=query_nodes,disease_list=disease_list,selected_nodes=selected_nodes))
}

conditional_query_items = function(input,output,query_nodes,W_ADJ,query_th){
  
  nano_qn = which(query_nodes %in% nano)
  drug_qn = which(query_nodes %in% drugs)
  dis_qn = which(query_nodes %in% disease)
  chem_qn = which(query_nodes %in% chemical)
  
  type_qn = rep("nano",length(query_nodes))
  type_qn[drug_qn] = "drugs"
  type_qn[dis_qn] = "disease"
  type_qn[chem_qn] = "chemical"
  
  nano_qn_e = query_nodes[which(query_nodes %in% nano)]
  drug_qn_e = query_nodes[which(query_nodes %in% drugs)]
  dis_qn_e = query_nodes[which(query_nodes %in% disease)]
  chem_qn_e = query_nodes[which(query_nodes %in% chemical)]
  
  if(DEBUGGING){
    cat("nano_qn_e:",nano_qn_e,"\n")
    cat("drug_qn_e:",drug_qn_e,"\n")
    cat("dis_qn_e:",dis_qn_e,"\n")
    cat("chem_qn_e:",chem_qn_e,"\n")
    cat("is.null(drug_qn_e ",is.null(drug_qn_e),"\n")
    cat("(length(nano_qn_e) == 0)", (length(nano_qn_e) == 0), "\n")
  }
  #solo chemical
  if((length(nano_qn_e) == 0) & (length(drug_qn_e)==0) & (length(dis_qn_e)==0)){
    if(DEBUGGING)
      cat("Only chemical is not null\n")
    expand.grid(chem_qn_e) -> combinations
  }
  #solo disease
  else if((length(nano_qn_e) == 0) & (length(drug_qn_e)==0) & (length(chem_qn_e)==0)){
    if(DEBUGGING)
      cat("Only disease is not null\n")
    expand.grid(dis_qn_e) -> combinations
  }
  
  #solo nano
  else if((length(dis_qn_e) == 0) & (length(drug_qn_e)==0) & (length(chem_qn_e)==0)){
    if(DEBUGGING)
      cat("Only nano is not null\n")
    
    expand.grid(nano_qn_e) -> combinations
  }
  
  #solo drug
  else if((length(dis_qn_e) == 0) & (length(nano_qn_e)==0) & (length(chem_qn_e)==0)){
    if(DEBUGGING)
      cat("Only drug is not null\n")
    
    expand.grid(drug_qn_e) -> combinations
  }
  
  #solo nano chemical
  else if((length(nano_qn_e)==0) & (length(chem_qn_e)==0)){
    if(DEBUGGING)
      cat("Nano and chemical are null\n")
    
    expand.grid(dis_qn_e,drug_qn_e) -> combinations
  }
  #solo nano drug
  else if((length(nano_qn_e)==0) & (length(drug_qn_e)==0)){
    if(DEBUGGING)
      cat("drug and chemical are null\n")
    
    expand.grid(chem_qn_e,dis_qn_e) -> combinations
  }
  #solo nano disease
  else if((length(nano_qn_e)==0) & (length(dis_qn_e)==0)){
    if(DEBUGGING)
      cat("nano and disease are null\n")
    
    expand.grid(chem_qn_e,drug_qn_e) -> combinations
  }
  
  #solo chem drug
  else if((length(chem_qn_e)==0) & (length(drug_qn_e)==0)){
    if(DEBUGGING)
      cat("drug and chemical are null\n")
    
    expand.grid(nano_qn_e,dis_qn_e) -> combinations
  }
  #solo chem dis
  else if((length(chem_qn_e)==0) & (length(dis_qn_e)==0)){
    if(DEBUGGING)
      cat("dis and chemical are null\n")
    
    expand.grid(nano_qn_e,drug_qn_e) -> combinations
  }
  
  #solo drug dis
  else if((length(drug_qn_e)==0) & (length(dis_qn_e)==0)){
    if(DEBUGGING)
      cat("drug and dis are null\n")
    
    expand.grid(nano_qn_e,chem_qn_e) -> combinations
  }
  
  else if(length(nano_qn_e)==0){
    if(DEBUGGING)
      cat("Nano is null\n")
    
    expand.grid(drug_qn_e,dis_qn_e,chem_qn_e) -> combinations
  }
  else if(length(drug_qn_e)==0){
    if(DEBUGGING)
      cat("drug is null\n")
    
    expand.grid(nano_qn_e,dis_qn_e,chem_qn_e) -> combinations
  }
  else if(length(dis_qn_e)==0){
    if(DEBUGGING)
      cat("dis is null\n")
    
    expand.grid(nano_qn_e,drug_qn_e,chem_qn_e) -> combinations
  }
  else if(length(chem_qn_e)==0){
    if(DEBUGGING)
      cat("chem is null\n")
    
    expand.grid(nano_qn_e,drug_qn_e,dis_qn_e) -> combinations
  }
  
  if(DEBUGGING){
    cat("class(combinations) ",class(combinations),"\n")
    cat("colnames(combinations) ",colnames(combinations),"\n")
  }
  combinations = as.matrix(combinations)
  
  Col_Sum_list = list()
  for(index_qn in 1:dim(combinations)[1]){
    variables_qn = combinations[index_qn,]
    MATRICE = W_ADJ[variables_qn,]
    MATRICE[which(MATRICE!=0)] = 1
    Col_Sum_list[[index_qn]] = colSums(MATRICE)
  }
  
  
  TF_qn = lapply(X = Col_Sum_list,FUN = function(list_qn){
    list_qn < query_th
  })
  
  if(length(TF_qn)>1){
    and_qn = TF_qn[[1]]
    for(index_qn in 2:length(TF_qn)){
      and_qn = and_qn & TF_qn[[index_qn]]
    }
  }else{
    and_qn = unlist(TF_qn)
  }
  
  
  and_qn[selected_nodes] = FALSE
  toREM = which(and_qn == TRUE)
  if(length(toREM)>0){
    W_ADJ = W_ADJ[-toREM,-toREM]
  }
  
  if(DEBUGGING)
    cat("dim(W_ADJ)[1]--> ",dim(W_ADJ)[1],"\n")
  
  if(DEBUGGING)
    cat("Max th--> ",max(unlist(Col_Sum_list)),"\n")
  if(DEBUGGING)
    cat("dim(W_ADJ)[1]--> ",dim(W_ADJ)[1],"\n")
  
  
  return(list(W_ADJ=W_ADJ,toREM=toREM,Col_Sum_list=Col_Sum_list))
  
}

conditional_query_subgraph_creation = function(input,output,W_ADJ,toREM,info_text,graph_gw,node_type){
  
  graph_gw = graph.adjacency(W_ADJ,mode="undirected",weighted=TRUE)
  V(graph_gw)$type = node_type[-toREM]
  graph_gw = igraph::delete.vertices(graph_gw,which(igraph::degree(graph_gw)<1))
  graph_gw = igraph::simplify(graph_gw)
  graph_s = graph_gw
  
  info_text = paste(info_text, "Nodes in the network:", vcount(graph_s),"<br/>")
  info_text = paste(info_text, "Edges in the network:", ecount(graph_s),"<br/>")
  
  
  prova2 <- reactive({
    paste("Nodes in the network:", vcount(graph_s))
  })
  
  index_nano = which(V(graph_s)$name %in% nano)
  index_drug = which(V(graph_s)$name %in% drugs)
  index_chem = which(V(graph_s)$name %in% chemical)
  index_dis = which(V(graph_s)$name %in% disease)
  
  tipes = rep("nano",length(V(graph_s)$name))
  tipes[index_drug] = "drugs"
  tipes[index_chem] = "chemical"
  tipes[index_dis] = "disease"
  
  info_text = paste(info_text, "Nanomaterials:", length(index_nano), "<br/>Drugs:", length(index_drug), 
                    "<br/>Chemical:", length(index_chem),"<br/>Disease:", length(index_dis),"<br/>")
  
  
  col = rep("pink",length(V(graph_s)$name))
  col[index_drug] = "yellow"
  col[index_chem] = "violet"
  col[index_dis] = "skyblue"
  
  V(graph_s)$type = tipes
  V(graph_s)$color = col
  graph_s = igraph::set.edge.attribute(graph_s,name = "color",index = E(graph_s)[which(sign(E(graph_s)$weight)<0)],value = "darkgreen")
  graph_s = igraph::set.edge.attribute(graph_s,name = "color",index = E(graph_s)[which(sign(E(graph_s)$weight)>0)],value = "red")
  
  ADJ_S = get.adjacency(graph_s,attr = "weight")
  ADJ_S = as.matrix(ADJ_S)
  
  return(list(graph_s=graph_s,ADJ_S = ADJ_S,info_text=info_text))
  
}