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
      names(Mi_count2) = paste(Mi_count[low_b:up_b,indexes_c[1]],Mi_count[low_b:up_b,indexes_c[2]],sep="-")
      #Mi_count2 = Mi_count$freq[1:nElem]
      #names(Mi_count2) = paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep="-")
    }
    if(length(columns_) == 2){
      indexes_c = 1:(dim(Mi_count)[2]-1)
      indexes_c = indexes_c[-column_index]
      nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
      dim(Mi_count)[1] -> up_b
      low_b = up_b - nElem
      Mi_count2 = Mi_count$freq[low_b:up_b]
      names(Mi_count2) = paste(Mi_count[low_b:up_b,indexes_c],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
      #Mi_count2 = Mi_count$freq[1:nElem]
      #names(Mi_count2) = paste(Mi_count[1:nElem,columns_[2:nElem_c]],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
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

genes_data_table_output = function(input,output,MList,MM_list,proxy,graph_s,g,g_geni2,items_list,query_type){
  output$genes_data_table = DT::renderDataTable({
    type = input$NetworkPattern
    type = as.integer(gsub(pattern = "M",x =type,replacement = ""))
    
    clique_type = input$clique_type
    
    validate(
      need(input$NetworkPattern != "", "Please select a pattern type")
    )
    
    s = input$clique_data_table_rows_selected
    
    if(query_type == "FREE"){
      g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(query_type == "CONDITIONAL"){
      g_ = generate_g_cliques(input,output,graph_s,proxy,MM_list,MList)
    }
    
    validate(need(length(s)!=0 , "Please select a clique"))          
    
    vids = V(g_)$name
  
    il = items_list[vids]
    union_genes = il[[1]]$entrez.genes
    
    for(i in 2:length(il)){
      union_genes = c(union_genes,union(union_genes,il[[i]]$entrez.genes))
    }
    
    ADJ = matrix(0,length(union_genes),length(vids))
    rownames(ADJ) = union_genes
    colnames(ADJ) = names(il)
    
    for(i in vids){
      gi = il[[i]]$entrez.genes
      gw = il[[i]]$edge.weigth
      ADJ[gi,i]=as.numeric(gw)
    }
  
    genes_elem  = rowSums(abs(ADJ))
  
    if(length(vids)==4){
      toRem = which(genes_elem<4)
      if(length(toRem)>0){
        genes_elem = genes_elem[-toRem]
      }
    }
    if(length(vids)==3){
      toRem = which(genes_elem<3)
      if(length(toRem)>0){
        genes_elem = genes_elem[-toRem]
      }
    }
    
    
    ADJ = ADJ[names(genes_elem),]
    
    cat("Evaluating dim(ADJ) :",dim(ADJ),"\n")
    
    validate(
      need(nrow(ADJ)>0, "The intersection of the genes activated by all elements in the clique is empty")
    )
    
    cat("Evaluating dim(ADJ) :",dim(ADJ),"\n")
    
    if(class(ADJ)=="numeric"){
      ADJ = matrix(data = ADJ,nrow = 1,ncol = length(ADJ))
      rownames(ADJ) = names(genes_elem)
      colnames(ADJ) = names(il)
      
    }
    
    x <- org.Hs.egSYMBOL
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    entrez = gsub(x = rownames(ADJ),pattern = "_at",replacement = "")
    xx[entrez] -> MYSYMBOL
    MYSYMBOL = unlist(MYSYMBOL)
    
    GENE_INFO = cbind(MYSYMBOL,ADJ)
    colnames(GENE_INFO)[1] = "Gene Name"
    
    cat("Evaluating dim(GENE_INFO) :",dim(GENE_INFO),"\n")
    
    for(row_i in 1:dim(GENE_INFO)[1]){
      #for(col_j in 1:dim(GENE_INFO)[2]){
      
      #http://www.genecards.org/cgi-bin/carddisp.pl?gene=TFPI2
      GENE_INFO[row_i,1] = paste('<a target="_blank" href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=',GENE_INFO[row_i,1],'">',GENE_INFO[row_i,1],'</a>',sep="")
      
      #GENE_INFO[row_i,1] = paste('<a target="_blank" href=\"https://www.google.com/?q=',GENE_INFO[row_i,1],'">',GENE_INFO[row_i,1],'</a>',sep="")
      #}
      for(i in vids){
        if(i %in% disease){
          if(GENE_INFO[row_i,i]=="1"){
            GENE_INFO[row_i,i] =  '<font color="black"><b>&#9666;&#9656;</b></font>' 
          }
          if(GENE_INFO[row_i,i]=="0"){
            GENE_INFO[row_i,i] =  '<font color="black">-</font>'  
          }
        }else{
          if(GENE_INFO[row_i,i]=="1"){
            GENE_INFO[row_i,i] =  '<font color="red"><b>&#9650;</b></font>' 
          }
          if(GENE_INFO[row_i,i]=="-1"){
            GENE_INFO[row_i,i] =  '<font color="green"><b>&#9660;</b></font>' 
          }
          if(GENE_INFO[row_i,i]=="0"){
            GENE_INFO[row_i,i] =  '<font color="black">-</font>' 
          }
        }
        
      }
    }
    
    DT::datatable(data =  GENE_INFO,
                  options = list(target = 'row+column',scrollX=TRUE,scrollY = "400px", scrollCollapse = TRUE,paging=FALSE),
                  escape=FALSE,rownames = FALSE,
                  selection = "single")
    
  })
}



enrich_clique = function(input,output,MList,MM_list,proxy,graph_s,g,items_list,query_type){  
  output$enriched_clique<- renderNanoCluster({
    
    type = input$NetworkPattern
    #     clique_type = input$clique_type
    #     
    #     validate(
    #       need(input$NetworkPattern != "", "Please select a pattern type")
    #     )
    #     
    type = as.integer(gsub(pattern = "M",x =type,replacement = ""))
    s = input$clique_data_table_rows_selected
    
    validate(
      need(s!="", "Please select a clique")
    )
    
    if(query_type == "FREE"){
      g_= internal_render_plot(MM_list[[type]],gr4=graph_s,s,proxyList = proxy)
    }
    if(query_type == "CONDITIONAL"){
      g_ = generate_g_cliques(input,output,graph_s,proxy,MM_list,MList)
    }
    
    clique_list = list(V(g_)$name)
    #clique_list = list(c("AgNP","salbutamol","11-deoxyprostaglandin E1"))
    message("In enrich_clique: clique_list: ",clique_list,"\n")

    th = input$EnrichTh
   # sets = input$EnrichType
    
    th_l = rep(th/100,length(V(g_)$name))
    message("In enrich_clique: th_l: ",th_l,"\n")
    
    gene_sets_list = list(c1_file,KEGG_file,biocarta_file,reactome_file,
                          c3Mir_file,c3Tft_file,c4_file,c5_file,c6_file,c7_file)
    gene_sets_name = c("Positional Gene Sets","KEGG","Biocarta","Reactome","microRNA targets","Trascription factor targets",
                       "Computational Gene Sets","GO gene sets","Oncogenic Signatures","Immunologic Signatures")
    
    message("In enrich_clique: gene_sets_name: ",gene_sets_name)
    DF = cliques_enrichment(clique_list = clique_list,g = g,g_ = g_,gene_sets_list = gene_sets_list,gene_sets_name = gene_sets_name,th_l = th_l,items_list = items_list)   
    
    edges = DF$edges
    vertices = DF$vertices
    colnames(vertices)[3]="size"
    
    validate(need(nrow(vertices)!=length(V(g_)$name),"No enrichmet for these sets and threshold"))
    
    nanocluster(Links = edges, Nodes = vertices,
                Source = "source", Target = "target",cluster_group = "cliques",
                last_level = 3, groups = names(table(vertices$group)),
                Value = "value", NodeID = "name",Nodesize = "size",
                Group = "group",zoom = TRUE,opacity = 0.95,fontSize = 20,
                legend = TRUE,
                charge = -1000,
                linkDistance = JS(paste0("function(d){return d.value*",10,
                                         "}")))
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
    
    if(vcount(g_) == 4){
      #plot(c(-1,0,1,0),c(0,1,0,-1))
      my_layout = matrix(c(-1,0,0,1,1,0,0,-1),4,2,byrow = TRUE)
    }
    if(vcount(g_)==3){
      #plot(c(-1,0,1),c(0,1,0))
      my_layout = matrix(c(-1,0,0,1,1,0),3,2,byrow = TRUE)
    }
    
      #par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),xpd=TRUE)#, mar = c(0, 0, 0, 0))
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)  
    plot(g_,vertex.color = V(g_)$color,
              vertex.size = 50,edge.width = abs(E(g_)$weight)+2,
              vertex.label.color = "black",layout = my_layout)
    legend("topright", legend = c("Positive Correlation","Negative Correlation"), xpd = TRUE,inset=c(-0.2,0),bty="n",
              fill = c("red","darkgreen"),cex=1)

  })
}

# plot_force_based_subnetwork_query_resutls = function(input,output,ADJ_S,chemMat,good_cliques,join10){
#   output$Subnetwork_plot = renderForceNetwork({ 
#     if(DEBUGGING)
#       cat("Subnetwork plot function\n")
#     ADJ_toPlot = ADJ_S
#     nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#     chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#     drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#     disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#     
#     chemMat_index = which(chemMat[,1] %in% chemical_s)
#     chemMat = chemMat[chemMat_index,]
#     
#     d = input$InterestingNodes
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     validate(
#       need(input$InterestingNodes != "", "SELECT A NODE!")
#     )
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     neig = lapply(X = good_cliques,FUN = function(i_cliques){
#       v_names = names(i_cliques)
#       if(sum(v_names %in% d)>0){
#         return(v_names)
#       }
#     })
#     
#     vds = unique(unlist(neig))
#     
#     edges_type = input$sub_Edges
#     if(DEBUGGING)
#       cat("Edges_type ",edges_type,"\n")
#     type = input$sub_checkGroup
#     
#     validate(
#       need(input$sub_checkGroup != "", "Please select an object group")
#     )
#     if(DEBUGGING)
#       cat("Tipi selezionati",length(type),"\n",type,"\n")
#     ATC_code = input$sub_ATCGroup
#     if(ATC_code == "ALL"){
#       ATC_code = ATC_letter_vector
#     }
#     if(DEBUGGING)
#       cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
#     chem_group = input$sub_ChemicalGroup
#     if(chem_group == "ALL"){
#       chem_group = chemical_group_vector
#     }
#     if(DEBUGGING)
#       cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
#     
#     items = c() 
#     items_type = c()
#     items_lengt = c()
#     #items_color = c()
#     for(i in 1:length(type)){
#       if(type[i]=="nano"){
#         items = c(items,nano_s)
#         items_type = c(items_type,"nano")
#         items_lengt = c(items_lengt,length(nano_s))
#         #items_color = c(items_color,"pink")
#       }
#       if(type[i]=="drugs"){
#         validate(
#           need(input$sub_ATCGroup != "", "Please select an ATC group")
#         )
#         already_selected = c()
#         for(j in 1:length(ATC_code)){
#           ATC_lev1 = substr(join10$code,1,1)
#           ATC_index = which(ATC_lev1 %in% ATC_code[j])
#           new_drugs = unique(join10$name[ATC_index])
#           new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#           index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#           already_selected = new_drugs[index_new_drugs]
#           
#           toRem = which(new_drugs[index_new_drugs] %in% items)
#           if(length(toRem)>0){
#             index_new_drugs = index_new_drugs[-toRem]
#           }
#           
#           items = c(items,new_drugs[index_new_drugs])
#           items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#           items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#         }
#       }
#       if(type[i]=="chemical"){
#         validate(
#           need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#         )
#         
#         already_selected = c()
#         for(j in 1:length(chem_group)){
#           chem_index = which(chemMat[,2] %in% chem_group[j])
#           new_chem = unique(chemMat[chem_index,1])
#           index_new_chem = which((new_chem %in% already_selected)==FALSE)
#           already_selected = new_chem[index_new_chem]
#           
#           toRem = which(new_chem[index_new_chem] %in% items)
#           if(length(toRem)>0){
#             index_new_chem = index_new_chem[-toRem]
#           }
#           
#           items = c(items,new_chem[index_new_chem])
#           items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#           items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#         }
#         
#       }
#       if(type[i]=="disease"){
#         good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#         items = c(items,good_disease)
#         items_type = c(items_type,"disease")
#         items_lengt = c(items_lengt,length(good_disease))
#       }
#     }
#     if(DEBUGGING){
#       cat("item_type: ",items_type,"\n")
#       cat("items_lengt: ",items_lengt,"\n")
#     }
#     if(edges_type == "P"){
#       ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
#     }
#     if(edges_type == "N"){
#       ADJ_toPlot[which(ADJ_toPlot>0)] = 0
#     }
#     
#     if(DEBUGGING)
#       cat("Graph creation... \n")
#     
#     gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot[items,items],mode = "undirected",weighted = TRUE)
#     V(gto_plot)$type = rep(items_type,items_lengt)
#     gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
#     gto_plot = igraph::simplify(gto_plot)
#     
#     if(DEBUGGING)
#       cat("Data frame creation... \n")
#     
#     data_frame = from_igraph_to_data_frame(gto_plot)
#     
#     MyClickScript <- 
#       '      d3.select(this).select("circle").transition()
#   .duration(750)
#   .attr("r", 30)'
#     
#     if(DEBUGGING){
#       cat("Force Networks plot... \n")
#       cat("Repulserad ",input$sub_repulserad,"\n")
#       cat("Edges length: ",input$sub_length,"\n")
#     }
#     forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
#                  Source = "source", Target = "target",
#                  Value = "value", NodeID = "name",Nodesize="size",
#                  Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
#                  legend = TRUE,
#                  clickAction = MyClickScript,
#                  charge = -input$sub_repulserad,
#                  linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,
#                                           "}")))
#   })
# }

plot_force_based_subnetwork_query_resutls = function(input,output,ADJ_S,chemMat,good_cliques,join10){
  output$Subnetwork_plot = renderForceNetwork({ 
    
  nodes = input$NodesOfInterest
  ATC_g  = input$sub_ATCGroup
  chem_g = input$sub_ChemicalGroup
  
  Elems = unique(unlist(good_cliques))
  nElems = length(Elems)
  FIT = find_items_type(Elems,nano,drugs,chemical,disease)
  items_type=FIT$items_type
  items_subtype=FIT$items_subtype
  
  idx_cliques = c()
  for(i in 1:length(good_cliques)){
    nodes_ = FALSE
    drug_ = FALSE
    chem_ = FALSE
    
    clique_elem = names(good_cliques[[i]])
    if(sum(nodes %in% clique_elem)>0){
      nodes_ = TRUE
      
      if("ALL" %in% ATC_g){
        drug_ =TRUE
      }else{
        idx_drug = clique_elem[which(clique_elem %in% drugs)]
        if(length(idx_drug)>0){
          drug_sub_idx = which(names(items_subtype) %in% idx_drug)
          ATC_gg = strsplit(x=gsub(pattern = "Drug ATC code: ",x = items_subtype[drug_sub_idx],replacement = ""),";")
          if(sum(ATC_g %in% ATC_gg)>0){
            drug_ = TRUE
          }
        }else{ # drug non presente nella clique, quindi la prendo lo stesso
          drug_ = TRUE
        }
      }
      if("ALL" %in% chem_g){
        chem_ = TRUE
      }else{
        idx_chem = clique_elem[which(clique_elem %in% chemical)]
        if(length(idx_chem)>0){
          chem_sub_idx = which(names(items_subtype) %in% idx_chem)
          if(sum(chem_g %in% items_subtype[chem_sub_idx])>0){
            chem_ = TRUE
          }
        }else{ # drug non presente nella clique, quindi la prendo lo stesso
          chem_ = TRUE
        }
      }
      if(nodes_ && drug_ && chem_){
        idx_cliques = c(idx_cliques,i)
      }
    }
#     if(nodes_){
#       cat("i: ",i," nodes: ",nodes, " nodes_: ",nodes_," ATC_g: ",ATC_g," drug_: ",drug_," chem_g ",chem_g," chem_: ",chem_,"\n")
#       
#     }
    
  }
  validate(need(length(idx_cliques)!=0,"No one of the nodes of interest is in the cliques!"))
  
  good_cliques = good_cliques[idx_cliques]
  
  if(DEBUGGING){
    message("In plot subnetwork and statistics:: nodes",nodes,"\n")
  }
  
  validate(need(nodes!="","Please select one or more nodes of interest!"))
  
  Elems = unique(unlist(good_cliques))
  nElems = length(Elems)
  
  adj_cliques_mat = matrix(0,nElems,nElems)
  colnames(adj_cliques_mat) = rownames(adj_cliques_mat)=unique(unlist(good_cliques))
  
  for(i in good_cliques){
    adj_cliques_mat[names(i),names(i)] = 1
  }
  
  diag(adj_cliques_mat) = 0
 
  FIT = find_items_type(colnames(adj_cliques_mat),nano,drugs,chemical,disease)
  items_type=FIT$items_type
  elems = FIT$elems
  
  gto_plot = graph.adjacency(adjmatrix = adj_cliques_mat,mode = "undirected",weighted = TRUE)
  V(gto_plot)$type = FIT$items_subtype
  gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
  gto_plot = igraph::simplify(gto_plot)
  
 # plot(gto_plot,vertex.size = igraph::degree(gto_plot)*5,mark.groups=good_cliques[4:7])
  
  
  data_frame = from_igraph_to_data_frame(gto_plot)

  n_clique = list()
  for(i in 1:length(good_cliques)){
    n_clique[[i]] = which(data_frame$vertices[,"name"] %in% names(good_cliques[[i]])) - 1
    
  }
#     forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
#                  Source = "source", Target = "target",
#                  Value = "value", NodeID = "name",Nodesize="size",
#                  Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
#                  legend = TRUE,
#                  charge = -1000,
#                  linkDistance = JS(paste0("function(d){return d.value*",10,"}")))
    
    
  forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
               Source = "source", Target = "target",
               Value = "value", NodeID = "name",Nodesize="size",
               Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
               legend = TRUE,
               clickAction = 'd3.select(this).select("circle").transition().duration(750).attr("r", 30)',
               charge = -input$sub_repulserad,
               linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,"}")))
  })
  
}

plot_subnetwork_statistics = function(input,output,ADJ_S,chemMat,good_cliques,join10){
  output$Subnetwork_plot_statistic = renderPlot({ 
    
    nodes = input$NodesOfInterest
    ATC_g  = input$sub_ATCGroup
    chem_g = input$sub_ChemicalGroup
    
    Elems = unique(unlist(good_cliques))
    nElems = length(Elems)
    FIT = find_items_type(Elems,nano,drugs,chemical,disease)
    items_type=FIT$items_type
    items_subtype=FIT$items_subtype
    
    idx_cliques = c()
    for(i in 1:length(good_cliques)){
      nodes_ = FALSE
      drug_ = FALSE
      chem_ = FALSE
      
      clique_elem = names(good_cliques[[i]])
      if(sum(nodes %in% clique_elem)>0){
        nodes_ = TRUE
        
        if("ALL" %in% ATC_g){
          drug_ =TRUE
        }else{
          idx_drug = clique_elem[which(clique_elem %in% drugs)]
          if(length(idx_drug)>0){
            drug_sub_idx = which(names(items_subtype) %in% idx_drug)
            ATC_gg = strsplit(x=gsub(pattern = "Drug ATC code: ",x = items_subtype[drug_sub_idx],replacement = ""),";")
            if(sum(ATC_g %in% ATC_gg)>0){
              drug_ = TRUE
            }
          }else{ # drug non presente nella clique, quindi la prendo lo stesso
            drug_ = TRUE
          }
        }
        if("ALL" %in% chem_g){
          chem_ = TRUE
        }else{
          idx_chem = clique_elem[which(clique_elem %in% chemical)]
          if(length(idx_chem)>0){
            chem_sub_idx = which(names(items_subtype) %in% idx_chem)
            if(sum(chem_g %in% items_subtype[chem_sub_idx])>0){
              chem_ = TRUE
            }
          }else{ # drug non presente nella clique, quindi la prendo lo stesso
            chem_ = TRUE
          }
        }
        if(nodes_ && drug_ && chem_){
          idx_cliques = c(idx_cliques,i)
        }
      }
#       if(nodes_){
#         cat("i: ",i," nodes: ",nodes, " nodes_: ",nodes_," ATC_g: ",ATC_g," drug_: ",drug_," chem_g ",chem_g," chem_: ",chem_,"\n")
#         
#       }
      
    }
    validate(need(length(idx_cliques)!=0,"No one of the nodes of interest is in the cliques!"))
    
    good_cliques = good_cliques[idx_cliques]
    
    if(DEBUGGING){
      message("In plot subnetwork and statistics:: nodes",nodes,"\n")
    }
    
    validate(need(nodes!="","Please select one or more nodes of interest!"))
    
    Elems = unique(unlist(good_cliques))
    nElems = length(Elems)
    
    adj_cliques_mat = matrix(0,nElems,nElems)
    colnames(adj_cliques_mat) = rownames(adj_cliques_mat)=unique(unlist(good_cliques))
    
    for(i in good_cliques){
      adj_cliques_mat[names(i),names(i)] = 1
    }
    
    diag(adj_cliques_mat) = 0
    
    FIT = find_items_type(colnames(adj_cliques_mat),nano,drugs,chemical,disease)
    items_type=FIT$items_type
    elems = FIT$elems
    
        slices = as.numeric(table(FIT$items_subtype))
      lbls = names(table(FIT$items_subtype))
      data_ = data.frame(slices,lbls,slices/max(slices))
      colnames(data_) = c("Occurrency","Groups","Frequencies")
      
      p =  ggplot(data_, aes(x =Occurrency , y = Groups)) +
        geom_point(aes(size = Occurrency, colour = Occurrency))  
      #theme(axis.text.x = element_text(angle = 90))
      print(p)
    })
  
}

# plot_force_based_subnetwork_query_resutls = function(input,output,ADJ_S,chemMat,good_cliques,join10){
#   output$Subnetwork_plot = renderForceNetwork({ 
#     
#     ADJ_toPlot = ADJ_S
#     nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#     chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#     drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#     disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#     
#     chemMat_index = which(chemMat[,1] %in% chemical_s)
#     chemMat = chemMat[chemMat_index,]
#     
#     d = input$InterestingNodes
#     edges_type = input$sub_Edges
#     type = input$sub_checkGroup
#     ATC_code = input$sub_ATCGroup
#     if(ATC_code == "ALL"){
#       ATC_code = ATC_letter_vector
#     }
#     
#     chem_group = input$sub_ChemicalGroup
#     if(chem_group == "ALL"){
#       chem_group = chemical_group_vector
#     }
#     
#     validate(
#       need(input$InterestingNodes != "", "SELECT A NODE!")
#     )
#     
#     validate(
#       need(input$sub_checkGroup != "", "Please select an object group")
#     )
#     
#     neig = lapply(X = good_cliques,FUN = function(i_cliques){
#       v_names = names(i_cliques)
#       if(sum(v_names %in% d)>0){
#         return(v_names)
#       }
#     })
#     
#     vds = unique(unlist(neig))
#     
#     message("in plot subnetwork vds= ",vds,"\n")
#     
#     FIT = find_items_type(vds,nano,drugs,chemical,disease)
#     items_type=FIT$items_type
#     elems = FIT$elems
#     
#     # "nano"     "drugs"    "chemical" "disease"
#     type_bool = c("nano" %in% type,
#                   "drugs" %in% type,
#                   "chemical" %in% type,
#                   "disease" %in% type)
#     
#     names(type_bool) = c("nano","drugs","chemical","disease")
#     
#     message("in plot subnetwork type_bool= ",type_bool,"\n")
#     
#     
#     class_toMaintain = names(which(type_bool==TRUE))
#     toMaintain = unlist(elems[class_toMaintain])
#     
#     message("in plot subnetwork toMaintain (FIRST)= ",toMaintain,"\n")
#     
#     
#     #Rimuovo per categoria
#     ADJ_toPlot = ADJ_toPlot[toMaintain,toMaintain]
#     
#     
#     #Rimuovo per sottocategoria
#     good_elems = c()
#     
#     if(type_bool["drugs"]){
#       if(length(ATC_code)==1){
#         if(ATC_code!="ALL"){
#         atc_drugs = join10[which(join10$ATC_lev1==ATC_code)]$name
#         atc_drug_q = elems$drugs[elems$drugs %in% atc_drugs]
#         good_elems = c(good_elems,atc_drug_q)
#         }
#       }else{
#           good_elems =c(good_elems,elems$drugs)
#       }
#     }
#     
#     if(type_bool["chemical"]){
#       if(length(chem_group)==1){
#         if(chem_group!="ALL"){
#           chemicals_group = chemMat[chemMat[,2] %in% chem_group,1]
#           atc_drug_q = elems$drugs[elems$drugs %in% atc_drugs]
#           good_elems = c(good_elems,atc_drug_q)
#         }
#       }
#       else{
#         good_elems =c(good_elems,elems$chemical)
#       } 
#     }
#     
#     
#     FIT2 = find_items_type(colnames(ADJ_toPlot),nano,drugs,chemical,disease)
#     items_type2=FIT2$items_type
#     elems2 = FIT2$elems
#     
#     message("in plot subnetwork toMaintain (good_elems)= ",good_elems,"\n")
#     message("in plot subnetwork toMaintain (nano)= ",elems2$nano,"\n")
#     message("in plot subnetwork toMaintain (disease)= ",elems2$disease,"\n")
#   
#     toMaintain2 = c(elems2$nano,elems2$disease,unique(good_elems))
#     
#     message("in plot subnetwork toMaintain (SECOND)= ",toMaintain2,"\n")
#     
#     
#     ADJ_toPlot = ADJ_toPlot[toMaintain2,toMaintain2]
#     
#     if(edges_type == "P"){
#         ADJ_toPlot[which(ADJ_toPlot < 0)] = 0
#     }
#     if(edges_type == "N"){
#         ADJ_toPlot[which(ADJ_toPlot>0)] = 0
#     }
#    
#     FIT3 = find_items_type(colnames(ADJ_toPlot),nano,drugs,chemical,disease)
#     
#     gto_plot = graph.adjacency(adjmatrix = ADJ_toPlot,mode = "undirected",weighted = TRUE)
#     V(gto_plot)$type = FIT3$items_type
#     gto_plot = igraph::delete.vertices(gto_plot,which(degree(gto_plot)<1))
#     gto_plot = igraph::simplify(gto_plot)
#     
#     if(DEBUGGING)
#       cat("Data frame creation... \n")
#     
#     data_frame = from_igraph_to_data_frame(gto_plot)
#     
#     MyClickScript <- 
#       '      d3.select(this).select("circle").transition()
#     .duration(750)
#     .attr("r", 30)'
#     
#     if(DEBUGGING){
#       cat("Force Networks plot... \n")
#       cat("Repulserad ",input$sub_repulserad,"\n")
#       cat("Edges length: ",input$sub_length,"\n")
#     }
#     forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
#                  Source = "source", Target = "target",
#                  Value = "value", NodeID = "name",Nodesize="size",
#                  Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
#                  legend = TRUE,
#                  clickAction = MyClickScript,
#                  charge = -input$sub_repulserad,
#                  linkDistance = JS(paste0("function(d){return d.value*",input$sub_length,
#                                           "}")))
# #     
# #     forceNetwork(Links = data_frame$edges, Nodes = data_frame$vertices,
# #                  Source = "source", Target = "target",
# #                  Value = "value", NodeID = "name",Nodesize="size",
# #                  Group = "type",zoom = TRUE,opacity = 0.95,fontSize = 10,
# #                  legend = TRUE,
# #                  clickAction = MyClickScript,
# #                  charge = -500,
# #                  linkDistance = JS(paste0("function(d){return d.value*",2,
# #                                           "}")))
#   })
# }

# plot_subnetwork_statistics = function(input,output,ADJ_S,chemMat,good_cliques,join10){
#   #output$Subnetwork_plot_statistic = renderPlot({
#   output$Subnetwork_plot_statistic = renderPlotly({
#    if(DEBUGGING)
#       cat("STATISTICHE SOTTONETWORK...\n")
#     ADJ_toPlot = ADJ_S
#     nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#     chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#     drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#     disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#     
#     chemMat_index = which(chemMat[,1] %in% chemical_s)
#     chemMat = chemMat[chemMat_index,]
#     
#     d = input$InterestingNodes
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     validate(
#       need(input$InterestingNodes != "", "SELECT A NODE!")
#     )
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     neig = lapply(X = good_cliques,FUN = function(i_cliques){
#       v_names = names(i_cliques)
#       if(sum(v_names %in% d)>0){
#         return(v_names)
#       }
#     })
#     
#     vds = unique(unlist(neig))
#     
#     edges_type = input$sub_Edges
#     if(DEBUGGING)
#       cat("Edges_type ",edges_type,"\n")
#     type = input$sub_checkGroup
#     if(DEBUGGING)
#       cat("type: ",type,"\n")
#     
#     validate(
#       need(input$sub_checkGroup != "", "Please select an object group")
#     )
#     if(DEBUGGING)
#       cat("Tipi selezionati",length(type),"\n",type,"\n")
#     ATC_code = input$sub_ATCGroup
#     if(DEBUGGING)
#       cat("ATC_code ",ATC_code,"\n")
#     
#     if(ATC_code == "ALL"){
#       ATC_code = ATC_letter_vector
#     }
#     if(DEBUGGING)
#       cat("ATC code",length(ATC_code),"\n",ATC_code,"\n")
#     
#     chem_group = input$sub_ChemicalGroup
#     if(chem_group == "ALL"){
#       chem_group = chemical_group_vector
#     }
#     if(DEBUGGING)
#       cat("Chemical Group",length(chem_group),"\n",chem_group,"\n")
#     
#     items = c() 
#     items_type = c()
#     items_lengt = c()
#     #items_color = c()
#     for(i in 1:length(type)){
#       if(type[i]=="nano"){
#         items = c(items,nano_s)
#         items_type = c(items_type,"nano")
#         items_lengt = c(items_lengt,length(nano_s))
#         #items_color = c(items_color,"pink")
#       }
#       if(type[i]=="drugs"){
#         validate(
#           need(input$sub_ATCGroup != "", "Please select an ATC group")
#         )
#         already_selected = c()
#         for(j in 1:length(ATC_code)){
#           ATC_lev1 = substr(join10$code,1,1)
#           ATC_index = which(ATC_lev1 %in% ATC_code[j])
#           new_drugs = unique(join10$name[ATC_index])
#           new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#           index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#           already_selected = new_drugs[index_new_drugs]
#           
#           toRem = which(new_drugs[index_new_drugs] %in% items)
#           if(length(toRem)>0){
#             index_new_drugs = index_new_drugs[-toRem]
#           }
#           
#           items = c(items,new_drugs[index_new_drugs])
#           items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#           items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#         }
#       }
#       if(type[i]=="chemical"){
#         validate(
#           need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#         )
#         
#         already_selected = c()
#         for(j in 1:length(chem_group)){
#           chem_index = which(chemMat[,2] %in% chem_group[j])
#           new_chem = unique(chemMat[chem_index,1])
#           index_new_chem = which((new_chem %in% already_selected)==FALSE)
#           already_selected = new_chem[index_new_chem]
#           
#           toRem = which(new_chem[index_new_chem] %in% items)
#           if(length(toRem)>0){
#             index_new_chem = index_new_chem[-toRem]
#           }
#           
#           items = c(items,new_chem[index_new_chem])
#           items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#           items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#         }
#         
#       }
#       if(type[i]=="disease"){
#         good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#         items = c(items,good_disease)
#         items_type = c(items_type,"disease")
#         items_lengt = c(items_lengt,length(good_disease))
#       }
#     }
#     
#     slices = items_lengt
#     idx = which(slices!=0)
#     slices = slices[idx]
#     lbls = items_type[idx]
#     
#     x_ = 1:length(lbls)
#     data_ = data.frame(slices,lbls,slices/max(slices))
#     colnames(data_) = c("Occurrency","Groups","Frequencies")
#     
#     plot_ly(data_, x = 1:length(Groups), y = Occurrency, text = Groups,
#             mode = "markers",size = Occurrency/max(Occurrency),
#             color = Frequencies)
#   })
# }

plot_gene_subnetwork = function(input,output,ADJ_S,chemMat,good_cliques,join10,g,g_geni2){
  output$gene_Subnetwork_plot = renderForceNetwork({
#     if(DEBUGGING)
#       cat("Subnetwork plot\n")
#     ADJ_toPlot = ADJ_S
#     nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#     chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#     drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#     disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#     
#     chemMat_index = which(chemMat[,1] %in% chemical_s)
#     chemMat = chemMat[chemMat_index,]
#     
#     d = input$InterestingNodes
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     validate(
#       need(input$InterestingNodes != "", "SELECT A NODE!")
#     )
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     neig = lapply(X = good_cliques,FUN = function(i_cliques){
#       v_names = names(i_cliques)
#       if(sum(v_names %in% d)>0){
#         return(v_names)
#       }
#     })
#     
#     vds = unique(unlist(neig))
#     
#     edges_type = input$sub_Edges
#     if(DEBUGGING)
#       cat("Edges_type ",edges_type,"\n")
#     type = input$sub_checkGroup
#     
#     validate(
#       need(input$sub_checkGroup != "", "Please select an object group")
#     )
#     
#     ATC_code = input$sub_ATCGroup
#     if(ATC_code == "ALL"){
#       ATC_code = ATC_letter_vector
#     }
#     chem_group = input$sub_ChemicalGroup
#     if(chem_group == "ALL"){
#       chem_group = chemical_group_vector
#     }
#     
#     items = c() 
#     items_type = c()
#     items_lengt = c()
#     #items_color = c()
#     for(i in 1:length(type)){
#       if(type[i]=="nano"){
#         items = c(items,nano_s)
#         items_type = c(items_type,"nano")
#         items_lengt = c(items_lengt,length(nano_s))
#         #items_color = c(items_color,"pink")
#       }
#       if(type[i]=="drugs"){
#         validate(
#           need(input$sub_ATCGroup != "", "Please select an ATC group")
#         )
#         already_selected = c()
#         for(j in 1:length(ATC_code)){
#           ATC_lev1 = substr(join10$code,1,1)
#           ATC_index = which(ATC_lev1 %in% ATC_code[j])
#           new_drugs = unique(join10$name[ATC_index])
#           new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#           index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#           already_selected = new_drugs[index_new_drugs]
#           
#           toRem = which(new_drugs[index_new_drugs] %in% items)
#           if(length(toRem)>0){
#             index_new_drugs = index_new_drugs[-toRem]
#           }
#           
#           items = c(items,new_drugs[index_new_drugs])
#           items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#           items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#         }
#       }
#       if(type[i]=="chemical"){
#         validate(
#           need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#         )
#         
#         already_selected = c()
#         for(j in 1:length(chem_group)){
#           chem_index = which(chemMat[,2] %in% chem_group[j])
#           new_chem = unique(chemMat[chem_index,1])
#           index_new_chem = which((new_chem %in% already_selected)==FALSE)
#           already_selected = new_chem[index_new_chem]
#           
#           toRem = which(new_chem[index_new_chem] %in% items)
#           if(length(toRem)>0){
#             index_new_chem = index_new_chem[-toRem]
#           }
#           items = c(items,new_chem[index_new_chem])
#           items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#           items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#         }
#       }
#       if(type[i]=="disease"){
#         good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#         items = c(items,good_disease)
#         items_type = c(items_type,"disease")
#         items_lengt = c(items_lengt,length(good_disease))
#       }
#     }
    
    nodes = input$NodesOfInterest
    ATC_g  = input$sub_ATCGroup
    chem_g = input$sub_ChemicalGroup
    
    Elems = unique(unlist(good_cliques))
    nElems = length(Elems)
    FIT = find_items_type(Elems,nano,drugs,chemical,disease)
    items_type=FIT$items_type
    items_subtype=FIT$items_subtype
    
    idx_cliques = c()
    for(i in 1:length(good_cliques)){
      nodes_ = FALSE
      drug_ = FALSE
      chem_ = FALSE
      
      clique_elem = names(good_cliques[[i]])
      if(sum(nodes %in% clique_elem)>0){
        nodes_ = TRUE
        
        if("ALL" %in% ATC_g){
          drug_ =TRUE
        }else{
          idx_drug = clique_elem[which(clique_elem %in% drugs)]
          if(length(idx_drug)>0){
            drug_sub_idx = which(names(items_subtype) %in% idx_drug)
            ATC_gg = strsplit(x=gsub(pattern = "Drug ATC code: ",x = items_subtype[drug_sub_idx],replacement = ""),";")
            if(sum(ATC_g %in% ATC_gg)>0){
              drug_ = TRUE
            }
          }else{ # drug non presente nella clique, quindi la prendo lo stesso
            drug_ = TRUE
          }
        }
        if("ALL" %in% chem_g){
          chem_ = TRUE
        }else{
          idx_chem = clique_elem[which(clique_elem %in% chemical)]
          if(length(idx_chem)>0){
            chem_sub_idx = which(names(items_subtype) %in% idx_chem)
            if(sum(chem_g %in% items_subtype[chem_sub_idx])>0){
              chem_ = TRUE
            }
          }else{ # drug non presente nella clique, quindi la prendo lo stesso
            chem_ = TRUE
          }
        }
        if(nodes_ && drug_ && chem_){
          idx_cliques = c(idx_cliques,i)
        }
      }
      #     if(nodes_){
      #       cat("i: ",i," nodes: ",nodes, " nodes_: ",nodes_," ATC_g: ",ATC_g," drug_: ",drug_," chem_g ",chem_g," chem_: ",chem_,"\n")
      #       
      #     }
      
    }
    validate(need(length(idx_cliques)!=0,"No one of the nodes of interest is in the cliques!"))
    
    good_cliques = good_cliques[idx_cliques]
    
    if(DEBUGGING){
      message("In plot subnetwork and statistics:: nodes",nodes,"\n")
    }
    
    validate(need(nodes!="","Please select one or more nodes of interest!"))
    
    Elems = unique(unlist(good_cliques))
    nElems = length(Elems)
    
#     adj_cliques_mat = matrix(0,nElems,nElems)
#     colnames(adj_cliques_mat) = rownames(adj_cliques_mat)=unique(unlist(good_cliques))
#     
#     for(i in good_cliques){
#       adj_cliques_mat[names(i),names(i)] = 1
#     }
#     
#     diag(adj_cliques_mat) = 0
#     
#     FIT = find_items_type(colnames(adj_cliques_mat),nano,drugs,chemical,disease)
#     items_type=FIT$items_type
#     elems = FIT$elems
    
    idx_g = which(V(g)$node_type == "gene")
    GMAT = igraph::get.adjacency(g)
    genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,Elems]))
    GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]

    x <- org.Hs.egSYMBOL
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    entrez = gsub(x = GG_genes,pattern = "_at",replacement = "")
    xx[entrez] -> MYSYMBOL

    idx_gg = which(MYSYMBOL %in% V(g_geni2)$name)
    geni_toPlot = igraph::induced.subgraph(g_geni2,unlist(MYSYMBOL[idx_gg]))
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
      '      d3.select(this).select("circle").transition().duration(750)
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
  #output$gene_Subnetwork_plot_statistics = renderPlotly({
    
#     if(DEBUGGING)
#       cat("Subnetwork plot\n")
#     ADJ_toPlot = ADJ_S
#     nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
#     chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
#     drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
#     disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
#     
#     chemMat_index = which(chemMat[,1] %in% chemical_s)
#     chemMat = chemMat[chemMat_index,]
#     
#     d = input$InterestingNodes
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     validate(
#       need(input$InterestingNodes != "", "SELECT A NODE!")
#     )
#     if(DEBUGGING)
#       cat("INTERESTING NODES: ",d,"\n")
#     
#     neig = lapply(X = good_cliques,FUN = function(i_cliques){
#       v_names = names(i_cliques)
#       if(sum(v_names %in% d)>0){
#         return(v_names)
#       }
#     })
#     
#     vds = unique(unlist(neig))
#     
#     edges_type = input$sub_Edges
#     if(DEBUGGING)
#       cat("Edges_type ",edges_type,"\n")
#     type = input$sub_checkGroup
#     
#     validate(
#       need(input$sub_checkGroup != "", "Please select an object group")
#     )
#     
#     ATC_code = input$sub_ATCGroup
#     if(ATC_code == "ALL"){
#       ATC_code = ATC_letter_vector
#     }
#     chem_group = input$sub_ChemicalGroup
#     if(chem_group == "ALL"){
#       chem_group = chemical_group_vector
#     }
#     
#     items = c() 
#     items_type = c()
#     items_lengt = c()
#     #items_color = c()
#     for(i in 1:length(type)){
#       if(type[i]=="nano"){
#         items = c(items,nano_s)
#         items_type = c(items_type,"nano")
#         items_lengt = c(items_lengt,length(nano_s))
#         #items_color = c(items_color,"pink")
#       }
#       if(type[i]=="drugs"){
#         validate(
#           need(input$sub_ATCGroup != "", "Please select an ATC group")
#         )
#         already_selected = c()
#         for(j in 1:length(ATC_code)){
#           ATC_lev1 = substr(join10$code,1,1)
#           ATC_index = which(ATC_lev1 %in% ATC_code[j])
#           new_drugs = unique(join10$name[ATC_index])
#           new_drugs = new_drugs[which(new_drugs %in% drug_s)]
#           index_new_drugs = which((new_drugs %in% already_selected)==FALSE)
#           already_selected = new_drugs[index_new_drugs]
#           
#           toRem = which(new_drugs[index_new_drugs] %in% items)
#           if(length(toRem)>0){
#             index_new_drugs = index_new_drugs[-toRem]
#           }
#           
#           items = c(items,new_drugs[index_new_drugs])
#           items_type = c(items_type,paste("Drugs ATC: ", ATC_code[j],sep=""))
#           items_lengt = c(items_lengt,length(new_drugs[index_new_drugs]))
#         }
#       }
#       if(type[i]=="chemical"){
#         validate(
#           need(input$sub_ChemicalGroup != "", "Please select a chemical group")
#         )
#         
#         already_selected = c()
#         for(j in 1:length(chem_group)){
#           chem_index = which(chemMat[,2] %in% chem_group[j])
#           new_chem = unique(chemMat[chem_index,1])
#           index_new_chem = which((new_chem %in% already_selected)==FALSE)
#           already_selected = new_chem[index_new_chem]
#           
#           toRem = which(new_chem[index_new_chem] %in% items)
#           if(length(toRem)>0){
#             index_new_chem = index_new_chem[-toRem]
#           }
#           items = c(items,new_chem[index_new_chem])
#           items_type = c(items_type,paste("Chemical class: ", chem_group[j],sep=""))
#           items_lengt = c(items_lengt,length(new_chem[index_new_chem]))
#         }
#       }
#       if(type[i]=="disease"){
#         good_disease = disease[disease %in% colnames(ADJ_toPlot)]
#         items = c(items,good_disease)
#         items_type = c(items_type,"disease")
#         items_lengt = c(items_lengt,length(good_disease))
#       }
#     }
    
    nodes = input$NodesOfInterest
    ATC_g  = input$sub_ATCGroup
    chem_g = input$sub_ChemicalGroup
    
    Elems = unique(unlist(good_cliques))
    nElems = length(Elems)
    FIT = find_items_type(Elems,nano,drugs,chemical,disease)
    items_type=FIT$items_type
    items_subtype=FIT$items_subtype
    
    idx_cliques = c()
    for(i in 1:length(good_cliques)){
      nodes_ = FALSE
      drug_ = FALSE
      chem_ = FALSE
      
      clique_elem = names(good_cliques[[i]])
      if(sum(nodes %in% clique_elem)>0){
        nodes_ = TRUE
        
        if("ALL" %in% ATC_g){
          drug_ =TRUE
        }else{
          idx_drug = clique_elem[which(clique_elem %in% drugs)]
          if(length(idx_drug)>0){
            drug_sub_idx = which(names(items_subtype) %in% idx_drug)
            ATC_gg = strsplit(x=gsub(pattern = "Drug ATC code: ",x = items_subtype[drug_sub_idx],replacement = ""),";")
            if(sum(ATC_g %in% ATC_gg)>0){
              drug_ = TRUE
            }
          }else{ # drug non presente nella clique, quindi la prendo lo stesso
            drug_ = TRUE
          }
        }
        if("ALL" %in% chem_g){
          chem_ = TRUE
        }else{
          idx_chem = clique_elem[which(clique_elem %in% chemical)]
          if(length(idx_chem)>0){
            chem_sub_idx = which(names(items_subtype) %in% idx_chem)
            if(sum(chem_g %in% items_subtype[chem_sub_idx])>0){
              chem_ = TRUE
            }
          }else{ # drug non presente nella clique, quindi la prendo lo stesso
            chem_ = TRUE
          }
        }
        if(nodes_ && drug_ && chem_){
          idx_cliques = c(idx_cliques,i)
        }
      }
      #     if(nodes_){
      #       cat("i: ",i," nodes: ",nodes, " nodes_: ",nodes_," ATC_g: ",ATC_g," drug_: ",drug_," chem_g ",chem_g," chem_: ",chem_,"\n")
      #       
      #     }
      
    }
    validate(need(length(idx_cliques)!=0,"No one of the nodes of interest is in the cliques!"))
    
    good_cliques = good_cliques[idx_cliques]
    
    if(DEBUGGING){
      message("In plot subnetwork and statistics:: nodes",nodes,"\n")
    }
    
    validate(need(nodes!="","Please select one or more nodes of interest!"))
    
    Elems = unique(unlist(good_cliques))
    nElems = length(Elems)
    
    idx_g = which(V(g)$node_type == "gene")
    GMAT = igraph::get.adjacency(g)
    genes_items_interaction = rowSums(as.matrix(GMAT[idx_g,Elems]))
    GG_genes = rownames(GMAT)[which(genes_items_interaction>0)]
    
    x <- org.Hs.egSYMBOL
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    entrez = gsub(x = GG_genes,pattern = "_at",replacement = "")
    xx[entrez] -> MYSYMBOL
    
    idx_gg = which(MYSYMBOL %in% V(g_geni2)$name)
    geni_toPlot = igraph::induced.subgraph(g_geni2,unlist(MYSYMBOL[idx_gg]))
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
    
    slices = as.numeric(table(vertices$group))
    lbls = names(table(vertices$group))  
    
#     message("in plotting gene subnetwork statistics --> lbls", lbls,"\n")
# 
#     data_ = data.frame(slices,lbls,slices/max(slices))
#     colnames(data_) = c("Occurrency","Groups","Frequencies")
#     
#     plot_ly(data_, x = 1:length(Groups), y = Occurrency, text = Groups,
#             mode = "markers",size = Occurrency/max(Occurrency),
#             color = Frequencies)
#     
#     slices = as.numeric(table(FIT$items_subtype))
#     lbls = names(table(FIT$items_subtype))
    data_ = data.frame(slices,lbls,slices/max(slices))
    colnames(data_) = c("Occurrency","Groups","Frequencies")
    
    p =  ggplot(data_, aes(x =Occurrency , y = Groups)) +
      geom_point(aes(size = Occurrency, colour = Occurrency))  
    #theme(axis.text.x = element_text(angle = 90))
    print(p)
    
    
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
