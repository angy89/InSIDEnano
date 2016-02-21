pubmed_search = function(input,output,word2){
  output$articleTable<-DT::renderDataTable({
    d1<-input$date1
    d2<-input$date2
    
    cat("Query input: ",word2(),"date1: ",d1,"date2: ",d2,"\n")
    
    word2 = strsplit(x = word2(),split = " ")
    
    validate(need(word2()!="","Please insert query object"))
    validate(need(d1!="","Please insert starting date"))
    validate(need(d2!="","Please insert ending date"))
    
    query_ = ""
    for(i in 1:length(word2[[1]])){
      if(i<length(word2[[1]])){
        query_ = paste(query_,word2[[1]][i]," AND ",sep="")
      }else{
        query_ = paste(query_,word2[[1]][i],"",sep="")
      }
    }
    
    withProgress(message = 'Progress...', min = 1,max = 3, {
      cat("Query input: ",query_,"date1: ",d1,"date2: ",d2,"\n")
      url <- EUtilsQuery(query_, type = "esearch", db = "pubmed")
      lines <- readLines(url, warn = FALSE, encoding = "unknown")
      res <- RISmed:::ParseTags(lines)
      EMPTYCHECK <- length(grep("<eSearchResult><Count>0<\\/Count>", 
                                lines)) != 0
      incProgress(1, detail = "")
      
      #    validate(need(EMPTYCHECK!=TRUE,"No results for the query"))
      
      if (EMPTYCHECK) {
        res$Id <- character(0)
        res$Count <- 0
      }
      res = new("EUtilsSummary", db = "pubmed", count = res$Count, retstart = res$RetStart, 
                retmax = res$RetMax, PMID = res$Id, querytranslation = res$QueryTranslation)
      
      cat("Fetching query...\n")
      fetch <- EUtilsGet(res, type="efetch", db="pubmed")
      incProgress(1, detail = "")
      cat("Articles query...\n")
      
      articles<-data.frame('Abstract'=AbstractText(fetch))
      incProgress(1, detail = "")
      
      validate(need(nrow(articles)!=0,"No results for the query"))
      
      authors = Author(fetch)
      article_title = ArticleTitle(fetch)
      title = (Title(fetch))
      year = YearAccepted(fetch) 
      volume = Volume(fetch)
      
      
      authors_list = list()
      for(i in 1:length(authors)){
        authors_list[[i]] = paste(paste(authors[[i]][,1],authors[[i]][,2]),collapse=",")
      }
      paper_info = data.frame(article_title,title,year,volume,unlist(authors_list))
      
      DT::datatable(paper_info)
            
    })
    
  })
  
  output$wordPlot<-renderPlot({
    d1<-input$date1
    d2<-input$date2
    
    cat("Query input: ",word2(),"date1: ",d1,"date2: ",d2,"\n")
    
    word2 = strsplit(x = word2(),split = " ")
    
    validate(need(word2()!="","Please insert query object"))
    validate(need(d1!="","Please insert starting date"))
    validate(need(d2!="","Please insert ending date"))
    
    query_ = ""
    for(i in 1:length(word2[[1]])){
      if(i<length(word2[[1]])){
        query_ = paste(query_,word2[[1]][i]," AND ",sep="")
      }else{
        query_ = paste(query_,word2[[1]][i],"",sep="")
      }
    }
    
    withProgress(message = 'Progress...', min = 1,max = 3, {
      cat("Query input: ",query_,"date1: ",d1,"date2: ",d2,"\n")
      url <- EUtilsQuery(query_, type = "esearch", db = "pubmed")
      lines <- readLines(url, warn = FALSE, encoding = "unknown")
      res <- RISmed:::ParseTags(lines)
      EMPTYCHECK <- length(grep("<eSearchResult><Count>0<\\/Count>", 
                                lines)) != 0
      incProgress(1, detail = "")
      
      #    validate(need(EMPTYCHECK!=TRUE,"No results for the query"))
      
      if (EMPTYCHECK) {
        res$Id <- character(0)
        res$Count <- 0
      }
      res = new("EUtilsSummary", db = "pubmed", count = res$Count, retstart = res$RetStart, 
                retmax = res$RetMax, PMID = res$Id, querytranslation = res$QueryTranslation)
      
      cat("Fetching query...\n")
      fetch <- EUtilsGet(res, type="efetch", db="pubmed")
      incProgress(1, detail = "")
      cat("Articles query...\n")
      
      articles<-data.frame('Abstract'=AbstractText(fetch))
      abstracts<-as.character(articles$Abstract)
      abstracts<-paste(abstracts, sep="", collapse="") 
      cat("Preparing plot query...\n")
      incProgress(1, detail = "")
      
      validate(need(nrow(articles)!=0,"No results for the query"))
      cat("abstracts : ",abstracts,"\n")
      cat("length(abstracts): ",length(abstracts),"\n")
      wordcloud(abstracts, min.freq=5, max.words=150, colors=brewer.pal(7,"Dark2"))
      
    })
    
  })
}


parse_nano_query_input = function(inserted_nano,nano){
  if(length(inserted_nano)==0){
    return("") 
  }
  if("ALL" %in% inserted_nano){
    inserted_nano = nano 
  }
  inserted_nano = unique(inserted_nano)
  return(inserted_nano)
}

parse_drug_query_input = function(inserted_drug,drugs,ATC_letter_vector){
  if(length(inserted_drug)==0){
    return(NULL) 
  }
  
  drug_parsed = c()
  
  idx = which(inserted_drug %in% ATC_letter_vector)
  categories_inserted = inserted_drug[idx]
  
  if(length(idx)>0){
    remaining_drug = inserted_drug[-idx]
  }
  
  if("ALL" %in% inserted_drug){
    return(drugs)
  }
  if(length(categories_inserted)>0){
    index_D= which(join10$ATC_lev1 %in% categories_inserted)
    drug_parsed = c(drug_parsed,unique(join10[index_D,]$name))
  }
  return(unique(c(drug_parsed,remaining_drug)))
}

parse_chemical_query_input = function(inserted_chemical,chemical,chemical_group_vector){
  chemical_parsed = c()
  
  idx = which(inserted_chemical %in% chemical_group_vector)
  categories_inserted = inserted_chemical[idx]
  
  if(length(idx)>0){
    remaining_chemical = inserted_chemical[-idx]
  }
  
  if("ALL" %in% inserted_chemical){
    return(chemical)
  }
  if(length(categories_inserted)>0){
      index_C = which(chemMat[,2] %in% categories_inserted)
      chemical_query = unique(chemMat[index_C,1])
    
    
    index_D= which(join10$ATC_lev1 %in% categories_inserted)
    drug_parsed = c(drug_parsed,unique(join10[index_D,]$name))
  }
  return(unique(c(drug_parsed,remaining_drug)))
}


find_items_type = function(items,nano,drugs,chemical,disease){
  # "nano"     "drugs"    "chemical" "disease"
  
  items_subtype = c()
  for(i in items){
    if(i %in% nano){
      items_subtype = c(items_subtype,"Nano")
      next
    }
    if(i %in% disease){
      items_subtype = c(items_subtype,"Disease")
      next
    }
    if(i %in% drugs){
      idx = which(join10$name %in% i)
      code = paste(unique(join10[idx,]$ATC_lev1),collapse = ";")
      code = paste("Drug ATC code: ",code,sep="")
      items_subtype = c(items_subtype,code)
      next
    }
    if(i %in% chemical){
      idx = which(chemMat[,1] %in% i)
      if(length(idx)==0){
        code = "Chemical - Unknown subcategory"
      }else{
        code = paste(unique(chemMat[idx,2]),collapse = ";")
      }
      items_subtype = c(items_subtype,code)
      next
    }
  }
  
  names(items_subtype) = items
  
  items_type = rep("nano",length(items))
  items_type[items %in% chemical] = "chemical"
  items_type[items %in% drugs] = "drugs"
  items_type[items %in% disease] = "disease"
  names(items_type) = items
  
  nano_el = items[items %in% nano]
  chemical_el = items[items %in% chemical]
  drugs_el = items[items %in% drugs]
  disease_el = items[items %in% disease]
  
  
  return(list(items_type = items_type,items_subtype=items_subtype,elems = list("nano" = nano_el,
                                                                               "drugs"= drugs_el,
                                                                               "disease" = disease_el,
                                                                               "chemical" = chemical_el)))
  
}


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
  
  ni = paste(input$nano_input,collapse="")
  dri = paste(input$drug_input,collapse="")
  ci = paste(input$chemical_input,collapse="")
  di = paste(input$disease_input,collapse="")
  
  control_i = c(ni!="",dri!="",ci!="",di!="")  
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes. control_i:",control_i)
  }
  
  xx = paste(input$nano_input,input$drug_input, input$chemical_input, input$disease_input,sep="")
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes. Concatenazione: ",xx,"\n")
    message("query_utilities::conditional_query_nodes. length(xx): ",length(xx),"\n")
  }
  if(sum(control_i)<=1){
    output$info2_1 <- renderUI({
      HTML("Please, fill in at least two fields of the form!")
    }) 
    validate(need(sum(control_i)>1, "Please, fill in at least two fields of the form!!"))
  }
  
  nano_query = input$nano_input
  if(length(nano_query) !=0 ){
    if("ALL" %in% nano_query){
      nano_query = nano
    }
  }
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {nano_query = ",nano_query, "}\n")
  }
  
  drug_query = input$drug_input
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {drug_query before checking= ",drug_query, "}\n")
  }
  
  if(length(drug_query) !=0 ){
    if("ALL" %in% drug_query){
      drug_query = drugs
    }
    if(drug_query=="A" || drug_query=="C" || drug_query=="D" || 
       drug_query=="G" || drug_query=="H" || drug_query=="J" || 
       drug_query=="L" || drug_query=="M" || drug_query=="N" ||
       drug_query=="P" || drug_query=="R" || drug_query=="S"){
      index_D= which(join10$ATC_lev1 ==drug_query)
      drug_query = unique(join10[index_D,]$name)
    }
  }
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {drug_query = ",drug_query, "}\n")
  }
  
  chemical_query = input$chemical_input
  if(length(chemical_query) != 0){
    if("ALL" %in% chemical_query){
      chemical_query = chemical
    }
    if(chemical_query %in% names(table(chemMat[,2]))){
      index_C = which(chemMat[,2] == chemical_query)
      chemical_query = unique(chemMat[index_C,1])
    }
  }
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {chemical_query = ",chemical_query, "}\n")
  }
  
  disease_query = input$disease_input
  if(length(disease_query)!=0){
    if("ALL" %in% disease_query){
      disease_query = disease
    }
  }
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {disease_query = ",disease_query, "}\n")
  }
  query_nodes = c(nano_query,drug_query,disease_query,chemical_query)
  type_qn = c(rep("nano",length(nano_query)),
              rep("drugs",length(drug_query)),
              rep("disease",length(disease_query)),
              rep("chemical",length(chemical_query)))
  
  for(i in query_nodes){
    disease_list[[i]] = i
    selected_nodes = c(selected_nodes,i)
  }
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {query_nodes = ",query_nodes, "}\n")
  }
  
  
  combination_vect=c(paste(nano_query,collapse = "")!="",paste(drug_query,collapse = "")!="",
                     paste(disease_query,collapse = "")!="",paste(chemical_query,collapse = "")!="")
  
  names(combination_vect)=c("nano","drug","dis","chem")
  
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes {combination_vect = ",combination_vect, "}\n")
  }
  
  return(list(query_nodes=query_nodes,
              nano_query = nano_query,
              drug_query = drug_query,
              chemical_query = chemical_query,
              disease_query = disease_query,
              combination_vect=combination_vect,
              disease_list=disease_list,selected_nodes=selected_nodes))
}

conditional_query_items = function(input,output,nano_qn_e,drug_qn_e,
                                   dis_qn_e,chem_qn_e,type_qn,combination_vect,
                                   query_nodes,W_ADJ,query_th,selected_nodes){
#   
#   nano_qn_e = c("ZnO","ZnO1")
#   drug_qn_e = c()
#   dis_qn_e=c("Asthma")
#   chem_qn_e = c()
  
  query_elements = list(nano_qn_e = nano_qn_e, drug_qn_e=drug_qn_e,dis_qn_e=dis_qn_e,chem_qn_e=chem_qn_e)

  
#   expand.grid(nano_qn_e,drug_qn_e,dis_qn_e,chem_qn_e)
#   expand.grid(nano_qn_e,dis_qn_e)
#   expand.grid(nano_qn_e)
#   combination_vect = c(TRUE,FALSE,TRUE,FALSE)
  combinations = expand.grid(query_elements[combination_vect])
  combinations = as.matrix(combinations)
  
  
  if(DEBUGGING){
    message("In conditional_query_items::: nano_qn_e:",nano_qn_e,"\n")
    message("In conditional_query_items::: drug_qn_e:",drug_qn_e,"\n")
    message("In conditional_query_items::: dis_qn_e:",dis_qn_e,"\n")
    message("In conditional_query_items::: chem_qn_e:",chem_qn_e,"\n")
    message("In conditional_query_items::: combination_vect:",combination_vect,"\n")
    message("In conditional_query_items::: is.null(drug_qn_e ",is.null(drug_qn_e),"\n")
    message("In conditional_query_items::: (length(nano_qn_e) == 0)", (length(nano_qn_e) == 0), "\n")
    message("In conditional_query_items::: class(combinations) ",class(combinations),"\n")
    message("In conditional_query_items::: colnames(combinations) ",colnames(combinations),"\n")
  }
  
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
  
  
  return(list(W_ADJ=W_ADJ,combinations=combinations,toREM=toREM,Col_Sum_list=Col_Sum_list))
  
}

conditional_query_subgraph_creation = function(input,output,W_ADJ,toREM,info_text,node_type){
  
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
  
  return(list(graph_s=graph_s,ADJ_S = ADJ_S,info_text=info_text,graph_gw=graph_gw))
}

conditional_check_cliques = function(obj,ADJ_S){
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
  
  type_query = c(length(v_nano) * length(v_dis),
                 length(v_nano) * length(v_chem),
                 length(v_nano) * length(v_drug),
                 length(v_drug) * length(v_dis),
                 length(v_drug) * length(v_chem),
                 length(v_dis) * length(v_chem))
  
  type_query[type_query>0] = 1
  
  row = sign(c(ADJ_S[v_nano,v_dis],
               ADJ_S[v_nano,v_chem],
               ADJ_S[v_nano,v_drug],
               ADJ_S[v_drug,v_dis],
               ADJ_S[v_drug,v_chem],
               ADJ_S[v_dis,v_chem]))
  
  row_ = rep(0,length(type_query))
  row_[which(type_query!=0)] = row
  
  row = row_ * type_query
  
#   if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)>0)){
#     row = sign(c(ADJ_S[v_nano,v_dis],
#                  ADJ_S[v_nano,v_chem],
#                  ADJ_S[v_nano,v_drug],
#                  ADJ_S[v_drug,v_dis],
#                  ADJ_S[v_drug,v_chem],
#                  ADJ_S[v_dis,v_chem]))
#   }
#   if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)==0) & (length(v_dis)>0)){
#     row = sign(c(ADJ_S[v_nano,v_dis],
#                  0,
#                  ADJ_S[v_nano,v_drug],
#                  ADJ_S[v_drug,v_dis],
#                  0,
#                  0
#     ))
#   }
#   if((length(v_nano) > 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)==0)){
#     row = sign(c(0,
#                  ADJ_S[v_nano,v_chem],
#                  ADJ_S[v_nano,v_drug],
#                  0,
#                  ADJ_S[v_drug,v_chem],
#                  0
#     ))
#   }
#   if((length(v_nano) > 0) & (length(v_drug)==0) & (length(v_chem)>0) & (length(v_dis)==0)){
#     row =  row = sign(c(ADJ_S[v_nano,v_dis],
#                         ADJ_S[v_nano,v_chem],
#                         0,
#                         0,
#                         0,
#                         ADJ_S[v_dis,v_chem]
#     ))
#   }
#   if((length(v_nano) == 0) & (length(v_drug)>0) & (length(v_chem)>0) & (length(v_dis)>0)){
#     row = sign(c(0,
#                  0,
#                  0,
#                  ADJ_S[v_drug,v_dis],
#                  ADJ_S[v_drug,v_chem],
#                  ADJ_S[v_dis,v_chem]
#     ))
#   }
  return(row)
}

conditional_clustering_cliques = function(input,output,good_cliques,ADJ_S){
  a = lapply(X = good_cliques,FUN = function(obj){
    row = conditional_check_cliques(obj,ADJ_S)
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
  
  return(cliques_groups)
}

build_query_list = function(input,output,cliques_groups,nano,drugs,chemical,disease,good_cliques){
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
  
  return(list(MM_list = MM_list,MList=MList))
  
}

estimate_query_tyme = function(input,output,graph_s){
  n_nodi = vcount(graph_s)
  estimated_tyme = (n_nodi^3 * 9) + (n_nodi^4 * 16)
  if(DEBUGGING)
    cat("estimated_tyme: ",estimated_tyme,"\n")
  
  output$extimatedTime = renderUI({
    HTML(paste("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size<br/> <strong>Estimated Iteration: ",estimated_tyme,"<strong/><br/>"))
  })
  return(estimated_tyme)
}

print_previous_query_table = function(input,output,LOG_CONDITIONAL){
  output$previous_query = DT::renderDataTable({
    toPlot = LOG_CONDITIONAL
    toPlot = apply(X = toPlot,MARGIN = 2,FUN = function(col_){
      gsub(x = col_,pattern = "_",replacement = "")
      
    })
    toPlot = toPlot[-1,]
    col_to_rem = which(colnames(toPlot) %in% "File_name")
    toPlot = toPlot[,-col_to_rem]
    DT::datatable(data =  toPlot,
                #  options = list(target = 'row+column',
                 #                scrollX=TRUE, scrollCollapse = TRUE,paging=FALSE),
                  escape=FALSE,rownames = FALSE,
                  selection = "single")
    
  })
  output$LoadQuery = renderUI({
    actionButton("LoadQuery", label = "Load Query", icon=icon("search",lib = "glyphicon"))
  })
  
}

load_query_from_table = function(input,output,LOG_CONDITIONAL,items_list){
  cat("In Load query from table!\n")
   
  selected_row = input$previous_query_rows_selected
    
  
  RData_file = paste(LOCAL_PATH,"Log_folder/",as.character(LOG_CONDITIONAL[selected_row+1,9]),sep="")
  
  nano_query = LOG_CONDITIONAL[selected_row+1,"Nano"]
  drugs_query = LOG_CONDITIONAL[selected_row+1,"Drugs"]
  disease_query = LOG_CONDITIONAL[selected_row+1,"Disease"]
  chemical_query = LOG_CONDITIONAL[selected_row+1,"Chemical"]
  th_query = LOG_CONDITIONAL[selected_row+1,"Th"]
  n1_query = LOG_CONDITIONAL[selected_row+1,"n1"]
  n2_query = LOG_CONDITIONAL[selected_row+1,"n2"]
  type_query = LOG_CONDITIONAL[selected_row+1,"Type"]
  
  if(nano_query!=""){
    nano_query = strsplit(x = as.character(nano_query),split = "_")
  }
  if(drugs_query!=""){
    drugs_query = strsplit(x = as.character(drugs_query),split = "_")
  }
  if(disease_query!=""){
    disease_query = strsplit(x = as.character(disease_query),split = "_")
  }
  if(chemical_query!=""){
    chemical_query = strsplit(x = as.character(chemical_query),split = "_")
  }
  
  if(DEBUGGING){
    message("In load_query_from_table::nano_query:",nano_query,"\n")
    message("In load_query_from_table::drugs_query:",drugs_query,"\n")
    message("In load_query_from_table::disease_query:",disease_query,"\n")
    message("In load_query_from_table::chemical_query:",chemical_query,"\n")
    message("In load_query_from_table::th_query:",th_query,"\n")
    message("In load_query_from_table::n1_query:",n1_query,"\n")
    message("In load_query_from_table::n2_query:",n2_query,"\n")
    message("In load_query_from_table::type_query:",type_query,"\n")
    
  }
  conditional_query_UI_set_query_values(input,output,nano_query,drugs_query,disease_query,chemical_query,th_query,n1_query,n2_query,type_query)
    
#   load_conditional_query(input=input,output=output,RData_file = RData_file,items_list = items_list,g)  
  
  withProgress(message = 'Progress...', min = 1,max = 6, {
    incProgress(1, detail = "Loading file...")
    
    load(file = RData_file)
    
    output$info2_1 <- renderUI({HTML(info_text)}) 
    
    output$NodesOfInterest <- renderUI({
      selectInput("NodesOfInterest",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
    })
    
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
                    options = list(order = list(list(1, 'desc')),target = 'row+column',
                                   scrollX=TRUE,scrollY = "400px", 
                                   scrollCollapse = TRUE,paging=FALSE),
                    escape=FALSE,
                    selection = "single")
      
    })
    
    incProgress(1, detail = "Preparing Output...")
    
    proxy = dataTableProxy("clique_data_table")
    incProgress(1, detail = "Preparing Barplot output...")
    
   # barplot_pattern_conditional_query(input,output,MList,graph_gw) #in conditional_query_output.R
    barplot_patter_conditional_query_input(input,output,MList)
    barplot_pattern_conditional_query(input,output,MList)
   
    clique_graph_cq_plot(input,output,MList,MM_list,proxy,graph_s)#in conditional_query_output.R
    incProgress(1, detail = "Preparing Genes Data table Output...")
    
    genes_data_table_output(input,output,MList,MM_list,proxy,graph_s,g,g_geni2,items_list,"CONDITIONAL") #in conditional_query_output.R
    incProgress(1, detail = "Enrich Cliques...")
    
    enrich_clique(input,output,MList,MM_list,proxy,graph_s,g,items_list,"CONDITIONAL")
    
    incProgress(1, detail = "Plot subnetwork and statistics...")
    
    plot_force_based_subnetwork_query_resutls(input,output,ADJ_S,chemMat,good_cliques,join10) #in qury_outputs.R
    plot_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10)#in qury_outputs.R
#     if(DEBUGGING){
#       message("Sto per chiamare plot_subnetwork_and_statistics\n")
#     }
#     plot_subnetwork_and_statistics(input,output,good_cliques)
#      
#      if(DEBUGGING){
#        message("Ho chiamato plot_subnetwork_and_statistics\n")
#      }
    plot_gene_subnetwork(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
    plot_gene_subnetwork_statistics(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2)#in qury_outputs.R
    
     
   
  })
    
}