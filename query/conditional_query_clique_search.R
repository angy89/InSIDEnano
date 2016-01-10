conditional_query_cliques_search = function(input,output,CLIQUE_TYPE,graph_s,ADJ_S,nElem_cliques){
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
  
  return(proxy)

}