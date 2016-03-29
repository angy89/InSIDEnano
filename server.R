shinyServer(function(input, output,session){
  
  check_login(output,FALSE) #plot the login control
  set_observer(input,output,session)  #set the observe event
  
  observeEvent(input$login, {
    validate_login(input,output)
    
    if(input$username != username | input$password != password){
      login_failed(input,output,session)  
    }else{
      login_successful(input,output,session)
      check_login(output,TRUE)
      
      phenotypic_network_UI(input,output)
      gene_network_UI(input,output)
      KEGG_path_network_UI(input,output)
      render_clustering_radial_network(input,output,NANO,CHEMICAL,DISEASE,DRUGS)
      
      render_nano_collapsible_tree(input,output,NANO)
      render_drugs_collapsible_tree(input,output,DRUGS)
      render_disease_collapsible_tree(input,output,DISEASE)
      render_chemical_collapsible_tree(input,output,CHEMICAL)
      
      observe({
        nnode = input$nodeName
        render_nano_gene_topTable(input,output,nnode,toSave,MYSYMBOL)
      })
      
      free_query_UI(input,output)
      items_query_UI(input,output)
      conditionl_query_UI(input,output)
      UI_query(input,output,nano,drugs,chemical,disease)
      
      
      withProgress(message = 'Progress...', min = 1,max = 5, {
        load(paste(APP_PATH,"graph_without_genes_also_intra_classes_edges_with_properties80.RData",sep=""))
        #load(paste(APP_PATH,"node_type.RData",sep=""))
        
        incProgress(1, detail = "Data Loaded 1/5")
        load(paste(APP_PATH,"graph_without_genes_also_intra_classes_edges_network_estimation80_2.RData",sep=""))
        rm(W_ADJ)
        load(paste(APP_PATH,"W_ADJ_clr_all_subset_positive_and_negative.RData",sep="")) #W_ADJ; W2_ADJ
        load(paste(APP_PATH,"network_without_genes_uniform.RData",sep="")) #ADJ01 #ADJ_UNIFORM
        load(paste(APP_PATH,"RankedUniformNetwork.RData",sep="")) #ADJ01_RANK #ADJ_UNIFORM_RANK
        
#         cat("W_ADJ: ",length(which(W_ADJ==0)),"\n")
# #         W_ADJ = aracne(abs(ADJ))
# #         W_ADJ = W_ADJ * sign(ADJ)
#         W_ADJ = ADJ
#         diag(W_ADJ) = 0
        
        
        
        #load(paste(APP_PATH,"W_ADJ.RData",sep=""))
        
        incProgress(1, detail = "Data Loaded 2/5")
        
        load(paste(APP_PATH,"gene_network_KEGG_th99.RData",sep="")) #g_geni2; edges; vertices
        load(paste(APP_PATH,"KEGG_PATH_ADJ.RData",sep="")) #KEGG_ADJ
        incProgress(1, detail = "Data Loaded 3/5")
        
        load(paste(APP_PATH,"big_net_with_chemical_up_down80_2_th_30.RData",sep="")) #g
        incProgress(1, detail = "Data Loaded 4/5")
        
        load(paste(APP_PATH,"items_gene_association_complete_no_genes.RData",sep="")) #items_list
        load(paste(LOCAL_PATH,"LOG.RData",sep=""))#LOG_CONDITIONAL
        
        incProgress(1, detail = "Data Loaded 5/5")
        incProgress(1, detail = "Waiting For input!")
        
      })
      
      cat("Loading LOG file\n")
      
      print_previous_query_table(input,output,LOG_CONDITIONAL)
      proxy_query = dataTableProxy("previous_query")
      
      observeEvent(input$refresh_free,{free_query_UI_refresh(input,output)}) #and observe event refresh_free
      observeEvent(input$Refresh, {conditional_query_UI_refresh(input,output)}) #end Observe Event refresh conditional query
      
      #       gene_network_UI(input,output)
      #       genes_input = gene_query_UI(input,output,g,genes_input)
      
      #gene_network_NANO_DISEASE_query(input,output,g)
      
      toRem = which(chemical %in% drugs)
      chem_idx = which(node_type=="chemical")
      toRem2 = which(colnames(W_ADJ)[chem_idx] %in% intersect(chemical,drugs))
      
      chemical = chemical[-toRem]
      toRem2 = toRem2 + 1292
      W_ADJ = W_ADJ[-toRem2,-toRem2]
      node_type = node_type[-toRem2]
      
      
      path_group = names(table(vertices$group))
      path_group = lapply(path_group,FUN = function(i){i})
      names(path_group) = unlist(path_group)
      
      output$Patway <- renderUI({
        selectInput("Patway_g",label = "Gene Pathway",multiple = TRUE, choices = path_group, selected = input$selected_KEGG)
      })
      
      good_disease = disease[disease %in% colnames(W_ADJ)]
      ii_dis_list = list("All" = "ALL")
      for(ii in good_disease){
        ii_dis_list[[ii]]=ii
      }
      
      output$input_dis <- renderUI({
        selectInput("input_dis",label = "Disease",multiple = TRUE,choices = ii_dis_list,selected = ii_dis_list[[2]])
      })
      if(DEBUGGING)
        cat("waiting for input \n")
      
      observeEvent(input$Go, {
        free_query(input,output,disease_list,selected_nodes,W_ADJ,
                   th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,items_list) #in free_query.R,
      }) #End Listener 1
      
      #ANGELA
      observeEvent(input$Go2, {
        load(paste(LOCAL_PATH,"LOG.RData",sep=""))
        if(DEBUGGING){
          cat("LOG file loaded \n")  
          cat("dim(LOG_CONDITIONAL",dim(LOG_CONDITIONAL),"\n")
        }
        
        log_file_to_load = check_already_existing_conditional_query(input,output,LOG_CONDITIONAL)
        
        if(DEBUGGING){
          cat("Query checked", log_file_to_load, "\n")  
        }
        if(is.null(log_file_to_load)){
          if(DEBUGGING){
            cat("New conditional query \n")  
          }
          conditional_query2(input,output,disease_list,selected_nodes,ADJ,ADJ_UNIFORM,ADJ_UNIFORM_RANK,th_p = input$th_slider2/100,node_type,chemMat,join10,g,g_geni2,LOG_CONDITIONAL,items_list)
          #conditional_query(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider2/100,node_type,chemMat,join10,g,g_geni2,LOG_CONDITIONAL,items_list)
        }else{
          if(DEBUGGING){
            cat("Load conditional query \n")  
          }
          load_conditional_query(session,input,output,log_file_to_load,items_list)
        }  
                
      }) #End Listener 2
      
      observeEvent(input$wordButton,{
        pubmed_search(input,output)#,word2)
        
      })
      
      observeEvent(input$action1_cluster, {
        res = ADJ_matrix(W_ADJ, input,output,nano,drugs,chemical,disease,chemMat,join10)
        ADJ2 = res$ADJ2
        g_clust = res$g_clust
        
        message("class(g_clust)",class(g_clust),"\n")
        
        data_frame = from_igraph_to_data_frame_cluster(g_clust,ADJ2)
        edges = data_frame$edges
        vertices = data_frame$vertices
        
        message("class(g_clust)",class(g_clust),"\n")
        
        output$cluster_output<- renderNanoCluster(
          nanocluster(Links = edges, Nodes = vertices,
                      Source = "source", Target = "target",
                      Value = "value", NodeID = "name",
                      Group = "group",zoom = TRUE,opacity = 0.95,fontSize = 20,
                      legend = TRUE,
                      charge = -input$repulseration_cluster,
                      linkDistance = JS(paste0("function(d){return d.value*",input$edge_length_cluster,
                                               "}")))
        )
        
        message("END Graph plotted!\n")
        
      })
      
      observeEvent(input$Go_couple,{        
        couple_query3(input,output,disease_list,selected_nodes,ADJ,ADJ01_UNIFORM,ADJ_UNIFORM_RANK,
                      th_p = input$th_slider_couple/100,node_type,chemMat,
                      join10,g,g_geni2,
                      items_list)
      })
      
      observeEvent(input$LoadQuery, {
        if(DEBUGGING){
          message("Loading Query...\n")
        }
        
        load_query_from_table(input,output,LOG_CONDITIONAL,items_list)
      })  
      
      observeEvent(input$Refresh_query_table,{
        load(paste(LOCAL_PATH,"LOG.RData",sep=""))
        print_previous_query_table(input,output,LOG_CONDITIONAL)
        proxy_query = dataTableProxy("previous_query")
        
      })

      #       observeEvent(input$Go3, {
      #         if(DEBUGGING) cat("GENE QUERY\n")
      #         gene_query(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider3/100,node_type,chemMat,join10,g,g_geni2,gene_input)
      #         
      #       }) #End Listener 3
      #       
      #       observeEvent(input$Go4,{
      #         gene_like_conditional_query(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider3/100,node_type,chemMat,join10,g,g_geni2)
      #           
      #       })
      
      #Outside Listener
      plot_item_network(input,output,W2_ADJ)
      plot_item_network_pie(input,output,W2_ADJ)
      plot_gene_network(input,output,g,g_geni2)
      plot_kegg_path_network(input,output,KEGG_ADJ)
      
    }
    
  }) #end event login
}) #end server




