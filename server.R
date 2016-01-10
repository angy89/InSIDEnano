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
      free_query_UI(input,output)
      conditionl_query_UI(input,output)
     
      withProgress(message = 'Progress...', min = 1,max = 5, {
        load(paste(APP_PATH,"graph_without_genes_also_intra_classes_edges_with_properties80.RData",sep=""))
        incProgress(1, detail = "Data Loaded 1/4")
        
        load(paste(APP_PATH,"graph_without_genes_also_intra_classes_edges_network_estimation80_2.RData",sep=""))
        incProgress(1, detail = "Data Loaded 2/4")
        
        load(paste(APP_PATH,"gene_network_KEGG_th99.RData",sep=""))
        incProgress(1, detail = "Data Loaded 3/4")
       
        load(paste(APP_PATH,"big_net_with_chemical_up_down80_2.RData",sep=""))
        incProgress(1, detail = "Data Loaded 3/4")

        incProgress(1, detail = "Data Loaded 4/4")
        incProgress(1, detail = "Waiting For input!")
        
      })
      
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
        selectInput("Patway_g",label = "Gene Pathway",multiple = TRUE, choices = path_group, selected = path_group[[1]])
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
                   th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2) #in free_query.R
      }) #End Listener 1
      
      #ANGELA
      observeEvent(input$Go2, {
        load(paste(APP_PATH,"LOG.RData",sep=""))
        
        
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
        conditional_query(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider2/100,node_type,chemMat,join10,g,g_geni2,LOG_CONDITIONAL)
        }else{
          if(DEBUGGING){
            cat("Load conditional query \n")  
          }
          load_conditional_query(input,output,log_file_to_load)
        }  
          
      }) #End Listener 2
      
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
    }
      
  }) #end event login
}) #end server




