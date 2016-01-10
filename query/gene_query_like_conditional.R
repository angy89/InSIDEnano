gene_like_conditional_query = function(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2){
  output$info2_1 <- renderUI({HTML(info_text)}) 
  CLIQUE_TYPE = input$gene_query_clique_type
  query_th = input$gene_query_percentuale_somma
  nElem_cliques = input$gene_query_nroCliques

  validate( need(input$gene_query_percentuale_somma != "", "Please select the treshold"))
  validate( need(input$gene_query_nroCliques != "", "Please select the number of nodes"))
  
  CQN = conditional_query_nodes(input,output,nano_query = input$gene_query_nano_input,
                                drug_query = input$gene_query_drug_input,
                                chemical_query = input$gene_query_chemical_input,
                                disease_query = input$gene_query_disease_input)
  query_nodes=CQN$query_nodes
  disease_list=CQN$disease_list
  selected_nodes=CQN$selected_nodes
  
  if(DEBUGGING)
    cat("query_nodes:",query_nodes,"\n")
  
  withProgress(message = 'Progress...', min = 1,max = 9, {
    incProgress(1, detail = "Evaluating input list...")
    free_query_UI_node_of_interest_output(input,output,disease_list) #in free_query_UI.R
    
    incProgress(1, detail = "Thresholding...")
    THS = find_thresholds(W_ADJ,th_p = input$th_slider3/100 ) #in query_utilities.R
    
    incProgress(1, detail = "Removing edges under threshold...")
    W_ADJ = apply_thresholds(W_ADJ,THS) #in query_utilities.R
    
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
    
    if(dim(W_ADJ)[1] < 3){
      info_text = "No results! \n"
      output$info2_1 <- renderUI({
        HTML(info_text)
      }) 
    }
    
    validate(
      need(dim(W_ADJ)[1] > 2, "No items to display! Try to decreases the threshold.")
    )
    if(DEBUGGING)
      cat("Max th--> ",max(unlist(Col_Sum_list)),"\n")
    if(DEBUGGING)
      cat("dim(W_ADJ)[1]--> ",dim(W_ADJ)[1],"\n")
    
    
    if(dim(W_ADJ)[1] == 0){
      if(DEBUGGING)
        cat("ARGHHH --> No items to display! Try to decreases the threshold.", max(unlist(Col_Sum_list)),"\n") 
    }
    
    incProgress(1, detail = "Creating graph...")
    
    
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
    
    incProgress(1, detail = "Searching for cliques")
    
    if(CLIQUE_TYPE == "ALL"){
      
      n_nodi = vcount(graph_s)
      estimated_tyme = (n_nodi^3 * 9) + (n_nodi^4 * 16)
      if(DEBUGGING)
        cat("estimated_tyme: ",estimated_tyme,"\n")
      
      output$extimatedTime = renderUI({
        HTML(paste("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size<br/> <strong>Estimated Iteration: ",estimated_tyme,"<strong/><br/>"))
      })
      
      
      mcl = cliques(graph=graph_s, min=3, max=4)
      
      
      is_good_NDCD = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% nano)
        if(ni == 1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% disease)
            if(ddis == 1){
              dc = sum(vertici %in% chemical)
              if(dc == 1){
                return(TRUE)
              }
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good_NDCD = unlist(is_good_NDCD)
      
      is_good_NDD = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% nano)
        if(ni == 1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% disease)
            if(ddis == 1){
              return(TRUE)
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good_NDD = unlist(is_good_NDD)
      
      is_good_NDC = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% nano)
        if(ni ==1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% chemical)
            if(ddis > 0){
              return(TRUE)
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good_NDC = unlist(is_good_NDC)
      
      is_good_DCD = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% chemical)
        if(ni == 1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% disease)
            if(ddis == 1){
              return(TRUE)
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good_DCD = unlist(is_good_DCD)
      
      is_good  = is_good_NDCD | is_good_NDD | is_good_NDC | is_good_DCD
      
      idx = which(is_good==TRUE)
      good_cliques = mcl[idx]
      
      is_good2 = unlist(lapply(X = good_cliques,FUN = function(obj){
        length(which(selected_nodes %in% names(obj)))
      }))
      
      idx2 = which(is_good2>=nElem_cliques)
      good_cliques = good_cliques[idx2]
      if(DEBUGGING)
        cat("Nro cliques: ",length(good_cliques),"\n")
      
      if(length(good_cliques) == 0){
        info_text = "No results! \n"
      }
      
      validate(
        need(length(good_cliques)>0, paste("No clique with this threshold"))
      )
      
      info_text = paste(info_text, "Number of cliques:",length(good_cliques),"<br/>")
      
      incProgress(1, detail = "Clustering cliques")
      
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
        if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)>0)){
          row = sign(c(ADJ_S[v_nano,v_dis],
                       ADJ_S[v_nano,v_chem],
                       ADJ_S[v_nano,v_drug],
                       ADJ_S[v_drug,v_dis],
                       ADJ_S[v_drug,v_chem],
                       ADJ_S[v_dis,v_chem]))
        }
        if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)==0) & (length(v_dis)>0)){
          row = sign(c(ADJ_S[v_nano,v_dis],
                       0,
                       ADJ_S[v_nano,v_drug],
                       ADJ_S[v_drug,v_dis],
                       0,
                       0
          ))
        }
        if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)==0)){
          row = sign(c(0,
                       ADJ_S[v_nano,v_chem],
                       ADJ_S[v_nano,v_drug],
                       0,
                       ADJ_S[v_drug,v_chem],
                       0
          ))
        }
        if((length(v_nano) > 0) & (length(v_drug)==0) & (length(v_chem)>0) & (length(v_dis)==0)){
          row =  row = sign(c(ADJ_S[v_nano,v_dis],
                              ADJ_S[v_nano,v_chem],
                              0,
                              0,
                              0,
                              ADJ_S[v_dis,v_chem]
          ))
        }
        if((length(v_nano) == 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)>0)){
          row = sign(c(0,
                       0,
                       0,
                       ADJ_S[v_drug,v_dis],
                       ADJ_S[v_drug,v_chem],
                       ADJ_S[v_dis,v_chem]
          ))
        }
        row
        
      })
      if(DEBUGGING)
        cat("length(a): ",length(a),"\n")
      
      MAT = do.call(rbind, a)
      if(DEBUGGING)
        cat("dim MAT ",dim(MAT),"\n")
      
      uniqueMAT = unique(MAT)
      
      
      cliques_groups = lapply(1:nrow(uniqueMAT),function(i){
        unlist(lapply(1:nrow(MAT),function(j){
          if(sum(uniqueMAT[i,]!=MAT[j,])==0){
            j
          }
        }))
      })
      
      info_text = paste(info_text, "Number of clique groups:",length(cliques_groups),"<br/>")
      
      
      output$NetworkPattern <- renderUI({
        cliques_LL = list()
        
        for(i in 1:length(cliques_groups)){
          cliques_LL[[paste0("Type",i)]] = paste0("M",i)
        }
        selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
      })
      if(DEBUGGING)
        cat("Nro of cliques groups: ",length(cliques_groups),"\n")
      
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
          
          if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)>0)){
            row = c(v_dis,v_nano,v_drug,v_chem)
            
          }
          if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)==0) & (length(v_dis)>0)){
            row = c(v_dis,v_nano,v_drug,"")
            
          }
          if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)==0)){
            row = c("",v_nano,v_drug,v_chem)
            
          }
          if((length(v_nano) > 0) & (length(v_drug)==0) & (length(v_chem)>0) & (length(v_dis)==0)){
            row = c(v_dis,v_nano,"",v_chem)
            
          }
          if((length(v_nano) == 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)>0)){
            row = c(v_dis,"",v_drug,v_chem)
            
          }
          row
        })
        M = do.call(rbind,ord_vertices)
      })
      
      
      # M_output_list = list()
      MM_list = list()
      nType = length(cliques_groups)
      
      for(i in 1:nType){
        Mi = MList[[i]]
        
        Mi_link = apply(X = Mi,MARGIN = 1,FUN = function(row_i){
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
      
      incProgress(1, detail = "Building tables")
      
      
    }#both cliques of length 3 and 4
    
    if(CLIQUE_TYPE == "NDCD"){ #clique da 4 oggetti
      mcl = cliques(graph=graph_s, min=4, max=4)
      incProgress(1, detail = "Evaluating cliques")
      is_good = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% nano)
        if(ni == 1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% disease)
            if(ddis == 1){
              dc = sum(vertici %in% chemical)
              if(dc == 1){
                return(TRUE)
              }
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good = unlist(is_good)
      sum(is_good)
      idx = which(is_good==TRUE)
      
      #plot(cl,vertex.color=V(cl)$color)
      
      good_cliques = mcl[idx]
      
      is_good2 = unlist(lapply(X = good_cliques,FUN = function(obj){
        length(which(selected_nodes %in% names(obj)))
      }))
      
      idx2 = which(is_good2>=nElem_cliques)
      good_cliques = good_cliques[idx2]
      if(DEBUGGING)
        cat("Nro cliques: ",length(good_cliques),"\n")
      
      if(length(good_cliques) == 0){
        info_text = "No results! \n"
      }
      
      validate(
        need(length(good_cliques)>0, paste("No clique with this threshold"))
      )
      if(DEBUGGING)
        cat("--> No cliques for the paramenters.\n ")
      
      info_text = paste(info_text, "Number of cliques:",length(good_cliques),"<br/>")
      
      incProgress(1, detail = "Clustering cliques")
      
      
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
      
      #         output$info2_5 <- renderText({
      #           paste("Number of clique groups:",length(cliques_groups))
      #         })
      
      info_text = paste(info_text, "Number of clique groups:",length(cliques_groups),"<br/>")
      
      
      output$NetworkPattern <- renderUI({
        cliques_LL = list()
        
        for(i in 1:length(cliques_groups)){
          cliques_LL[[paste0("Type",i)]] = paste0("M",i)
        }
        selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
      })
      
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
      
      
      # M_output_list = list()
      MM_list = list()
      nType = length(cliques_groups)
      
      for(i in 1:nType){
        Mi = MList[[i]]
        
        Mi_link = apply(X = Mi,MARGIN = 1,FUN = function(row_i){
          #xx = paste('<a href=\"https://www.google.com/#q=',row_i,' class=','"btn btn-primary"','>',row_i,'</a>',sep="")
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
      
      incProgress(1, detail = "Building tables")
    }# end IF clique di 4
    
    if(CLIQUE_TYPE == "NDD"){ #clique Nano-Drug-Disease
      mcl = cliques(graph=graph_s, min=3, max=3)
      incProgress(1, detail = "Evaluating cliques")
      is_good = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% nano)
        if(ni == 1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% disease)
            if(ddis == 1){
              return(TRUE)
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good = unlist(is_good)
      sum(is_good)
      idx = which(is_good==TRUE)
      
      #plot(cl,vertex.color=V(cl)$color)
      
      good_cliques = mcl[idx]
      
      is_good2 = unlist(lapply(X = good_cliques,FUN = function(obj){
        length(which(selected_nodes %in% names(obj)))
      }))
      
      idx2 = which(is_good2>=nElem_cliques)
      good_cliques = good_cliques[idx2]
      if(DEBUGGING)
        cat("Nro cliques: ",length(good_cliques),"\n")
      
      if(length(good_cliques) == 0){info_text = "No results! \n"}
      
      validate(
        need(length(good_cliques)>0, paste("No clique with this threshold"))
      )
      
      info_text = paste(info_text, "Number of cliques:",length(good_cliques),"<br/>")
      incProgress(1, detail = "Clustering cliques")
      
      
      a = lapply(X = good_cliques,FUN = function(obj){
        vertices = names(obj)
        v_nano = vertices[which(vertices %in% nano)]
        v_drug = vertices[which(vertices %in% drugs)]
        v_dis = vertices[which(vertices %in% disease)]
        
        row = sign(c(ADJ_S[v_nano,v_dis],
                     ADJ_S[v_nano,v_drug],
                     ADJ_S[v_drug,v_dis]
        ))
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
      
      info_text = paste(info_text, "Number of clique groups:",length(cliques_groups),"<br/>")
      
      output$NetworkPattern <- renderUI({
        cliques_LL = list()
        
        for(i in 1:length(cliques_groups)){
          cliques_LL[[paste0("Type",i)]] = paste0("M",i)
        }
        selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
      })
      
      
      MList = lapply(cliques_groups,FUN=function(obj){
        idx = obj
        good_cliques_i = good_cliques[idx]
        vertices_list = lapply(good_cliques_i,FUN = names)
        ord_vertices = lapply(vertices_list,FUN = function(vertices){
          v_nano = vertices[which(vertices %in% nano)]
          v_drug = vertices[which(vertices %in% drugs)]
          v_dis = vertices[which(vertices %in% disease)]
          
          c(v_dis,v_nano,v_drug)
        })
        M = do.call(rbind,ord_vertices)
      })
      
      MM_list = list()
      nType = length(cliques_groups)
      
      for(i in 1:nType){
        Mi = MList[[i]]
        
        Mi_link = apply(X = Mi,MARGIN = 1,FUN = function(row_i){
          # xx = paste('<a href=\"https://www.google.com/#q=',row_i,' class=','"btn btn-primary"','>',row_i,'</a>',sep="")
          xx = paste('<a target="_blank" href=\"https://www.google.com/?q=',row_i,'">',row_i,'</a>',sep="")
          
          new_row = c(xx[1],xx[2],xx[3])
        })
        Mi = t(Mi_link)
        rownames(Mi) = 1:dim(Mi)[1]
        Mi = as.data.frame(Mi)
        colnames(Mi)=c("Disease","Nano","Drug")
        MM_list[[i]] = Mi
        #M_output_list[[i]] = DT::renderDataTable(Mi,options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE),escape = FALSE)
      }
      
      incProgress(1, detail = "Building tables")
    }# end IF Nano-Drug-Disease
    
    if(CLIQUE_TYPE == "NDC"){ #clique Nano-Drug-Chemical
      mcl = cliques(graph=graph_s, min=3, max=3)
      incProgress(1, detail = "Evaluating cliques")
      is_good = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% nano)
        if(ni ==1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% chemical)
            if(ddis > 0){
              return(TRUE)
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good = unlist(is_good)
      sum(is_good)
      idx = which(is_good==TRUE)
      good_cliques = mcl[idx]
      
      is_good2 = unlist(lapply(X = good_cliques,FUN = function(obj){
        length(which(selected_nodes %in% names(obj)))
      }))
      
      idx2 = which(is_good2>=nElem_cliques)
      good_cliques = good_cliques[idx2]
      if(DEBUGGING)
        cat("Nro cliques: ",length(good_cliques),"\n")
      
      if(length(good_cliques) == 0){info_text = "No results! \n"}
      
      validate(
        need(length(good_cliques)>0, paste("No clique with this threshold"))
      )
      
      info_text = paste(info_text, "Number of cliques:",length(good_cliques),"<br/>")
      incProgress(1, detail = "Clustering cliques")
      
      
      a = lapply(X = good_cliques,FUN = function(obj){
        vertices = names(obj)
        v_nano = vertices[which(vertices %in% nano)]
        v_drug = vertices[which(vertices %in% drugs)]
        v_chem = vertices[which(vertices %in% chemical)]
        #v_dis = vertices[which(vertices %in% disease)]
        
        intersect(v_chem,v_drug) -> ii
        if(length(ii)>0){
          which(v_chem %in% v_drug) -> index_ii
          v_chem = v_chem[-index_ii]
        }
        
        row = sign(c(
          ADJ_S[v_nano,v_chem],
          ADJ_S[v_nano,v_drug],
          ADJ_S[v_drug,v_chem]
        ))
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
      
      info_text = paste(info_text, "Number of clique groups:",length(cliques_groups),"<br/>")
      
      output$NetworkPattern <- renderUI({
        cliques_LL = list()
        
        for(i in 1:length(cliques_groups)){
          cliques_LL[[paste0("Type",i)]] = paste0("M",i)
        }
        selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
      })
      
      MList = lapply(cliques_groups,FUN=function(obj){
        idx = obj
        good_cliques_i = good_cliques[idx]
        vertices_list = lapply(good_cliques_i,FUN = names)
        ord_vertices = lapply(vertices_list,FUN = function(vertices){
          v_nano = vertices[which(vertices %in% nano)]
          v_drug = vertices[which(vertices %in% drugs)]
          v_chem = vertices[which(vertices %in% chemical)]
          #v_dis = vertices[which(vertices %in% disease)]
          
          intersect(v_chem,v_drug) -> ii
          if(length(ii)>0){
            which(v_chem %in% v_drug) -> index_ii
            v_chem = v_chem[-index_ii]
          }
          
          c(v_chem,v_nano,v_drug)
        })
        M = do.call(rbind,ord_vertices)
      })
      
      MM_list = list()
      nType = length(cliques_groups)
      
      for(i in 1:nType){
        Mi = MList[[i]]
        if(DEBUGGING)
          cat("Mi[1,]", Mi[1,],"\n")
        
        Mi_link = apply(X = Mi,MARGIN = 1,FUN = function(row_i){
          #xx = paste('<a href=\"https://www.google.com/#q=',row_i,' class=','"btn btn-primary"','>',row_i,'</a>',sep="")
          xx = paste('<a target="_blank" href=\"https://www.google.com/?q=',row_i,'">',row_i,'</a>',sep="")
          
          new_row = c(xx[1],xx[2],xx[3])
        })
        Mi = t(Mi_link)
        rownames(Mi) = 1:dim(Mi)[1]
        Mi = as.data.frame(Mi)
        colnames(Mi)=c("Chemical","Nano","Drug")
        MM_list[[i]] = Mi
        #M_output_list[[i]] = DT::renderDataTable(Mi,options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE),escape = FALSE)
      }
      
      incProgress(1, detail = "Building tables")
    }# end IF Nano-Drug-Chemical
    
    if(CLIQUE_TYPE == "DCD"){ #clique Drug-Chemical-Disease
      mcl = cliques(graph=graph_s, min=3, max=3)
      incProgress(1, detail = "Evaluating cliques")
      is_good = lapply(X = mcl,FUN = function(obj){
        vertici = names(obj)
        ni = di = ddis = dc = 0
        
        ni = sum(vertici %in% chemical)
        if(ni == 1){
          di = sum(vertici %in% drugs)
          if(di ==1 ){
            ddis = sum(vertici %in% disease)
            if(ddis == 1){
              return(TRUE)
            }
          }
        }
        
        return(FALSE)
        
      })
      is_good = unlist(is_good)
      sum(is_good)
      idx = which(is_good==TRUE)
      
      good_cliques = mcl[idx]
      
      is_good2 = unlist(lapply(X = good_cliques,FUN = function(obj){
        length(which(selected_nodes %in% names(obj)))
      }))
      
      idx2 = which(is_good2>=nElem_cliques)
      good_cliques = good_cliques[idx2]
      if(DEBUGGING)
        cat("Nro cliques: ",length(good_cliques),"\n")
      
      if(length(good_cliques) == 0){info_text = "No results! \n" }
      
      validate(
        need(length(good_cliques)>0, paste("No clique with this threshold"))
      )
      
      info_text = paste(info_text, "Number of cliques:",length(good_cliques),"<br/>")
      incProgress(1, detail = "Clustering cliques")
      
      
      a = lapply(X = good_cliques,FUN = function(obj){
        vertices = names(obj)
        #v_nano = vertices[which(vertices %in% nano)]
        v_drug = vertices[which(vertices %in% drugs)]
        v_chem = vertices[which(vertices %in% chemical)]
        v_dis = vertices[which(vertices %in% disease)]
        
        intersect(v_chem,v_drug) -> ii
        if(length(ii)>0){
          which(v_chem %in% v_drug) -> index_ii
          v_chem = v_chem[-index_ii]
        }
        
        row = sign(c(
          #ADJ_S[v_nano,v_dis],
          # ADJ_S[v_nano,v_chem],
          #ADJ_S[v_nano,v_drug],
          ADJ_S[v_drug,v_dis],
          ADJ_S[v_drug,v_chem],
          ADJ_S[v_dis,v_chem]
        ))
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
      
      info_text = paste(info_text, "Number of clique groups:",length(cliques_groups),"<br/>")
      
      output$NetworkPattern <- renderUI({
        cliques_LL = list()
        
        for(i in 1:length(cliques_groups)){
          cliques_LL[[paste0("Type",i)]] = paste0("M",i)
        }
        selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
      })
      
      if(DEBUGGING)
        cat("Nro of cliques groups: ",length(cliques_groups),"\n")
      
      MList = lapply(cliques_groups,FUN=function(obj){
        idx = obj
        good_cliques_i = good_cliques[idx]
        vertices_list = lapply(good_cliques_i,FUN = names)
        ord_vertices = lapply(vertices_list,FUN = function(vertices){
          # v_nano = vertices[which(vertices %in% nano)]
          v_drug = vertices[which(vertices %in% drugs)]
          v_chem = vertices[which(vertices %in% chemical)]
          v_dis = vertices[which(vertices %in% disease)]
          
          intersect(v_chem,v_drug) -> ii
          if(length(ii)>0){
            which(v_chem %in% v_drug) -> index_ii
            v_chem = v_chem[-index_ii]
          }
          
          c(v_dis,v_chem,v_drug)
        })
        M = do.call(rbind,ord_vertices)
      })
      
      
      # M_output_list = list()
      MM_list = list()
      nType = length(cliques_groups)
      
      for(i in 1:nType){
        Mi = MList[[i]]
        
        Mi_link = apply(X = Mi,MARGIN = 1,FUN = function(row_i){
          # xx = paste('<a href=\"https://www.google.com/#q=',row_i,' class=','"btn btn-primary"','>',row_i,'</a>',sep="")
          xx = paste('<a target="_blank" href=\"https://www.google.com/?q=',row_i,'">',row_i,'</a>',sep="")
          
          new_row = c(xx[1],xx[2],xx[3])
        })
        Mi = t(Mi_link)
        rownames(Mi) = 1:dim(Mi)[1]
        Mi = as.data.frame(Mi)
        colnames(Mi)=c("Disease","Chemical","Drug")
        MM_list[[i]] = Mi
        #M_output_list[[i]] = DT::renderDataTable(Mi,options = list(order = list(list(1, 'desc')),target = 'row+column',scrollX=TRUE),escape = FALSE)
      }
      
      incProgress(1, detail = "Building tables")
    }# end IF Drug-Chemical-Disease
    
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
    
    
  })
  
  barplot_pattern_conditional_query(input,output,MList,graph_gw) #in conditional_query_output.R
  genes_data_table_output(input,output,MList,MM_list,proxy,graph_s,g,g_geni2) #in conditional_query_output.R
  clique_graph_cq_plot(input,output,MList,MM_list,proxy,graph_s)#in conditional_query_output.R
  plot_force_based_subnetwork_query_resutls(input,output,ADJ_S,chemMat,good_cliques,join10) #in qury_outputs.R
  plot_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10)#in qury_outputs.R
  plot_gene_subnetwork(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
  plot_gene_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
  
}