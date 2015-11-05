# Define server logic required to draw a histogram
shinyServer(function(input, output,session){
  #initial.ok = 0

  output$checkLOGIN = renderUI({
    HTML("<strong>Please login to proceed with the analisys! <strong/> <br/>")
  })
  
  output$checkLOGIN_free = renderUI({
    HTML("<strong>Please login to proceed with the analisys! <strong/> <br/>")
  })
  
  output$checkLOGIN_gene = renderUI({
    HTML("<strong>Please login to visualize the gene network! <strong/> <br/>")
  })
  
  output$checkLOGIN_entities = renderUI({
    HTML("<strong>Please login to visualize the network of entities! <strong/> <br/>")
  })
  
  output$extimatedTime = renderUI({
    HTML("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size.<br/>")
  })
    
  observe({
    if(input$action1 > 0){
      print('1')
      session$sendCustomMessage("myCallbackHandler", "1")
    }
  })
  observe({
    if(input$action4 > 0){
      print('5')
      session$sendCustomMessage("myCallbackHandler", "5")
    }
  })
  observe({
    if(input$action2 > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler", "2")
    }
  })
  observe({
    if(input$action3 > 0){
      print('3')
      session$sendCustomMessage("myCallbackHandler", "3")
    }
  })
  observe({
    if(input$Demo > 0){
      print('4')
      session$sendCustomMessage("myCallbackHandler", "4")
    }
  })
  
  
  observeEvent(input$login, {
    validate(need(input$username != "", "Please insert login name"))
    validate(need(input$password != "", "Please insert password"))
    
    if(input$username != username | input$password != password){
      if(DEBUGGING)
      cat("Login failed... \n")
      session$sendCustomMessage(type = 'testmessage',
                                message = "Login Failed: Wrong username or password! ")
      
      
    }else{
      session$sendCustomMessage(type = 'testmessage',
                                message = "Login Successfull! ")
      if(DEBUGGING)
      cat("Login successfull... \n")
      
      output$checkLOGIN = renderUI({
        HTML(" ")
      })
      
      output$checkLOGIN_free = renderUI({
        HTML(" ")
      })
      
      output$checkLOGIN_gene = renderUI({
        HTML(" ")
      })
      
      output$checkLOGIN_entities = renderUI({
        HTML(" ")
      })
      
      LOGGED_IN = TRUE
      
      output$items_net_page1 = renderUI({
        sidebarPanel(
          fluidRow(
            column(12, checkboxGroupInput("checkGroup", 
                                          label = "Object Network", 
                                          choices = list("Nanomaterials" = "nano", 
                                                         "Drugs" = "drugs", "Chemical" = "chemical",
                                                         "Disease" = "disease"),selected = "nano"))
          ),
          fluidRow(column(12,selectInput("Edges",
                                         label="Select the type of edges",
                                         choices = list("Positive" = "P",
                                                        "Negative" = "N",
                                                        "Both" = "B"),selected="B"))
          ),
          fluidRow(column(12,uiOutput("input_dis"))),
          
          fluidRow(
            column(12, selectInput("ATCGroup", 
                                   label ="Drugs ATC code", multiple = TRUE,
                                   choices = ATC_choice_list,selected = "A"))
          ),
          
          fluidRow(
            column(12, selectInput("ChemicalGroup", 
                                   label ="Chemical Classes", multiple = TRUE,
                                   choices = chemical_choice_list,
                                   selected = "Amino acids"))
          ),
          fluidRow(
            column(9, sliderInput("slider1", label = " Connection Strength",
                                  min = 1, max = 99, value = 99,step=1))
          ),
          fluidRow(
            column(9, sliderInput("repulserad", label = "Repulseration Strenght",
                                  min = 100, max = 1000, value = 100,step=1))
          ),
          fluidRow(
            column(9, sliderInput("length", label = "Edge length",
                                  min = 2, max = 10, value = 2,step=1))
          )
        ) #end Sidebar
      })
      
      output$items_net_page2 = renderUI({
        
        mainPanel(
          fluidRow(
            column(9,wellPanel(forceNetworkOutput("itemNetwork")))
          ),
          fluidRow(
            column(9,wellPanel(plotOutput("pie_chart")))
          )
        )
      })
      
      
      #Free query UI
      
      output$disease = renderUI({
        selectizeInput(
          'disease', label = "Select the item to start the analysis", choices = c(disease,nano,drugs,chemical), multiple = TRUE,
          options = list(create = TRUE)
        )
      })
      
      output$th_slider = renderUI({
        sliderInput("th_slider", label = "Strenght of Similarity",
                    min = 80, max = 99, value = 99,step=1)
      })
      
      output$plotTripel_total = renderUI({
        selectInput("plotTripel_total",label = "Plot of Association Frequencies",
                    choices = list("Disease-Nano-Drug" = 1,
                                   "Disease-Nano-Chemical"= 2,
                                   "Disease-Nano"=3,
                                   "Disease-Drug"=4,
                                   "Disease-Chemical"=5,
                                   "Chemical-Nano" = 6,
                                   "Chemical-Drug" = 7,
                                   "Nano-Drug" = 8),selected = 1)
      })
      
      
      
      output$Go = renderUI({
        actionButton(inputId = "Go",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
      })
      
      output$ok = renderUI({
        actionButton(inputId = "ok", label = "Stop computation")
      })
      
      output$refresh_free = renderUI({
        actionButton(inputId = "refresh_free", label  = "Refresh",icon("refresh", lib = "glyphicon"))
      })
      
      #Conditional query UI
      
      output$nano_input = renderUI({
        selectizeInput('nano_input', label = "Nanomaterials", choices = c("ALL",nano), multiple = TRUE,
                       options = list(create = TRUE))
      })
      
      output$drug_input = renderUI({
        selectizeInput('drug_input', label = "Drugs", choices = c(unlist(ATC_choice_list),drugs), multiple = TRUE,
                       options = list(create = TRUE))
      })
      
      output$disease_input = renderUI({
        selectizeInput('disease_input', label = "Diseases", choices = c("ALL",disease), multiple = TRUE,
                       options = list(create = TRUE))
      })
      
      output$chemical_input = renderUI({
        selectizeInput('chemical_input', label = "Chemicals", choices = c(unlist(chemical_choice_list),chemical), multiple = TRUE,
                       options = list(create = TRUE))
      })
      
      output$th_slider2 = renderUI({
        sliderInput("th_slider2", label = "Strenght of Similarity", min = 10, max = 99, value = 99,step=1)
      })
      
      output$percentuale_somma = renderUI({
        numericInput(inputId = "percentuale_somma",label = "Number of query items connected to the selected nodes based on the threshold",
                     value = 2,min = 1,max = 10,step=1)
      })
      
      output$nroCliques = renderUI({
        numericInput(inputId = "nroCliques",label = "Number of query items that must be selected in the cliques",
                     value = 2,min = 1,max = 10,step=1)
      })
      
      output$Refresh = renderUI({
        actionButton(inputId = "Refresh",label  = "Refresh",icon("refresh", lib = "glyphicon"))
      })
      
      output$Go2 = renderUI({
        actionButton(inputId = "Go2",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
      })
      
      withProgress(message = 'Progress...', min = 1,max = 5, {
        load(paste(APP_PATH,"graph_without_genes_also_intra_classes_edges_with_properties80.RData",sep=""))
        #load("./graph_without_genes_also_intra_classes_edges_with_properties80.RData") #prima CLR
        incProgress(1, detail = "Data Loaded 1/4")
        
        load(paste(APP_PATH,"graph_without_genes_also_intra_classes_edges_network_estimation80_2.RData",sep=""))
        incProgress(1, detail = "Data Loaded 2/4")
        
        load(paste(APP_PATH,"gene_network_KEGG.RData",sep=""))
        incProgress(1, detail = "Data Loaded 3/4")
       
        load(paste(APP_PATH,"big_net_with_chemical_up_down80_2.RData",sep=""))
        incProgress(1, detail = "Data Loaded 3/4")

        incProgress(1, detail = "Data Loaded 4/4")
        incProgress(1, detail = "Waiting For input!")
        
      })
      
   
      output$clique_type = renderUI({
        selectizeInput("clique_type", label = "Select the type of clique you want to analyze", 
                       choices = list("All" = "ALL","Nano-Drug-Chemical-Disease" = "NDCD",
                                      "Nano-Drug-Disease" = "NDD",
                                      "Nano-Drug-Chemical" = "NDC",
                                      "Drug-Chemical-Disease" = "DCD"),
                       selected="ALL")})
      
      observeEvent(input$refresh_free,{
        output$disease = renderUI({
          selectizeInput(
            'disease', label = "Select the item to start the analysis", choices = c(disease,nano,drugs,chemical), multiple = TRUE,
            options = list(create = TRUE)
          )
        })
        
        output$th_slider = renderUI({
          sliderInput("th_slider", label = "Strenght of Similarity",
                      min = 80, max = 99, value = 99,step=1)
        })
        
        output$NodesOfInterest_totale <- renderUI({
          selectInput("InterestingNodes_totale",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
        })
        
        output$plotTripel_total = renderUI({
          selectInput("plotTripel_total",label = "Plot of Association Frequencies",
                      choices = list("Disease-Nano-Drug" = 1,
                                     "Disease-Nano-Chemical"= 2,
                                     "Disease-Nano"=3,
                                     "Disease-Drug"=4,
                                     "Disease-Chemical"=5,
                                     "Chemical-Nano" = 6,
                                     "Chemical-Drug" = 7,
                                     "Nano-Drug" = 8),selected = 1)
        })
        
      }) #and observe event refresh_free
      
      observeEvent(input$Refresh, {
        output$nano_input = renderUI({
          selectizeInput('nano_input', label = "Nanomaterials", choices = c("ALL",nano), multiple = TRUE,
                         options = list(create = TRUE))
        })
        
        output$extimatedTime = renderUI({
          HTML(paste("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size<br/> 
                     <strong> Estimated Iteration:"))
        })
        
        output$drug_input = renderUI({
          selectizeInput('drug_input', label = "Drugs", choices = c(unlist(ATC_choice_list),drugs), multiple = TRUE,
                         options = list(create = TRUE))
        })
        
        output$disease_input = renderUI({
          selectizeInput('disease_input', label = "Diseases", choices = c("ALL",disease), multiple = TRUE,
                         options = list(create = TRUE))
        })
        
        output$chemical_input = renderUI({
          selectizeInput('chemical_input', label = "Chemicals", choices = c(unlist(chemical_choice_list),chemical), multiple = TRUE,
                         options = list(create = TRUE))
        })
        
        output$th_slider2 = renderUI({
          sliderInput("th_slider2", label = "Strenght of Similarity", min = 10, max = 99, value = 99,step=1)
        })
        
        output$percentuale_somma = renderUI({
          numericInput(inputId = "percentuale_somma",label = "Number of query items connected to the selected nodes based on the threshold",
                       value = 2,min = 1,max = 10,step=1)
        })
        
        output$nroCliques = renderUI({
          numericInput(inputId = "nroCliques",label = "Number of query items that must be selected in the cliques",
                       value = 2,min = 1,max = 10,step=1)
        })
        
        output$Refresh = renderUI({
          actionButton(inputId = "Refresh",label  = "Refresh",icon("refresh", lib = "glyphicon"))
        })
        
        output$Go2 = renderUI({
          actionButton(inputId = "Go2",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
        })
        
        output$info2_1 <- renderUI({
          HTML(" ")
        }) 
        
        output$clique_type = renderUI({
          selectizeInput("clique_type", label = "Cliques type", 
                         choices = list("All" = "ALL","Nano-Drug-Chemical-Disease" = "NDCD",
                                        "Nano-Drug-Disease" = "NDD",
                                        "Nano-Drug-Chemical" = "NDC",
                                        "Drug-Chemical-Disease" = "DCD"),
                         selected="ALL")})
      }) #end Observe Event refresh conditional query
      
      output$gene_net_page1 = renderUI({
        sidebarPanel(
          fluidRow(column(12,uiOutput("Patway"))),
          fluidRow(
            column(9, sliderInput("gene_repulserad", label = "Repulseration Strenght",
                                  min = 100, max = 1000, value = 100,step=1))
          ),
          fluidRow(
            column(9, sliderInput("gene_length", label = "Edge length",
                                  min = 2, max = 10, value = 2,step=1))
          )
        )
      })
      
      output$gene_net_page2 = renderUI({
        mainPanel(
          fluidRow(column(12, wellPanel(forceNetworkOutput("geneNetwork"))))   
        )
      })
 
      
      output$gene_net_page3 = renderUI({
        sidebarPanel(
          fluidRow(
            column(9, 
                   selectizeInput(
                     'gene_nano_query', label = "Select the nanomaterial to start the analysis", choices = c(nano), multiple = FALSE,selected="MWCNT",
                     options = list(create = TRUE)
                   )       
            )
          ),
          fluidRow(
            column(9, 
                   selectizeInput(
                     'gene_disease_query', label = "Select the disease to start the analysis", choices = c(disease), multiple = FALSE,selected="Asthma",
                     options = list(create = TRUE)
                   )       
            )
          ),
          fluidRow(
            column(9, sliderInput("th_gene_query", label = "Percentage of drugs",
                                  min = 1, max = 100, value = 99,step=1))
          )
        )
      })
      
      output$gene_net_page4 = renderUI({
        mainPanel(
          fluidRow(column(12, wellPanel(forceNetworkOutput("geneNetwork_query"))))   
        )
      })
      
      output$geneNetwork_query = renderForceNetwork({
        
        nano_query = input$gene_nano_query
        disease_query = input$gene_disease_query
        drug_perc = input$th_gene_query
        
        validate(
          need(input$gene_nano_query != "", "Please select a nano")
        )
        
        validate(
          need(input$gene_disease_query != "", "Please select a disease")
        )
        
        if(DEBUGGING)
          cat("Query from ",input$gene_nano_query ,input$gene_disease_query,"to: ","\n")
        
      
        geni_toPlot = subgraph_nano_disease = function(g,from_nano = nano_query,to_disease=disease_query,drug_perc = drug_perc/100)
          
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
        
        MyClickScript <- 
          '      d3.select(this).select("circle").transition()
        .duration(750)
        .attr("r", 30)'
        
        MyClickScript2 <- 
          'd3.select(this).select(dLink.target).select("circle").transition().duration(750).attr("r",30),
        d3.select(dLink.target).select("circle").transition().duration(750).attr("r",30)
        '
        
        forceNetwork(Links = edges, Nodes = vertices,
                     Source = "source", Target = "target",
                     Value = "value", NodeID = "name",Nodesize="size",
                     zoom = TRUE,opacity = 0.85,fontSize = 10,Group = "group",
                     legend = TRUE, height = input$gene_width,width =input$gene_height, 
                     clickAction = MyClickScript,charge = -input$gene_repulserad,
                     linkDistance = JS(paste0("function(d){return d.value*",input$gene_length,"}")))
      })
      
      
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
        
        
        
        withProgress(message = 'Progress...', min = 1,max = 9, {
          validate(
            need(input$disease != "", "Please select a Disease")
          )
          
          for(i in input$disease){
            disease_list[[i]] = i
            selected_nodes = c(selected_nodes,i)
          }
          if(DEBUGGING){
            cat("selected_nodes", selected_nodes, "\n")  
            cat("disease_list", length(disease_list), "\n")  
          }
          incProgress(1, detail = "Evaluating input list...")
          
          output$NodesOfInterest <- renderUI({
            selectInput("InterestingNodes",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
          })
          
          output$NodesOfInterest_totale <- renderUI({
            selectInput("InterestingNodes_totale",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
          })
          
          output$NodesOfInterest_items <- renderUI({
            selectInput("InterestingNodes_items",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
          })
          
          th_p = input$th_slider/100 #threshold per cercare le cliques
          
          W2 = W_ADJ
          diag(W_ADJ) = 0
          
          th_n = 1-th_p
          if(DEBUGGING)
          cat("th ", th_p, "\n")
          
          incProgress(1, detail = "Thresholding...")
          if(DEBUGGING)
          cat("dim(W_ADJ", dim(W_ADJ), "\n")
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
          
          
          incProgress(1, detail = "Removing edges under threshold...")
          
          
          W_ADJ[nano,disease][which(W_ADJ[nano,disease]>0 & W_ADJ[nano,disease]<th_ndis_p)] = 0 #nano disease
          W_ADJ[nano,disease][which(W_ADJ[nano,disease]<0 & W_ADJ[nano,disease]>th_ndis_n)] = 0
          
          W_ADJ[nano,drugs][which(W_ADJ[nano,drugs]>0 & W_ADJ[nano,drugs]<th_nd_p)] = 0 #nano drugs
          W_ADJ[nano,drugs][which(W_ADJ[nano,drugs]<0 & W_ADJ[nano,drugs]>th_nd_n)] = 0
          
          W_ADJ[disease,drugs][which(W_ADJ[disease,drugs]>0 & W_ADJ[disease,drugs]<th_dd_p)] = 0 #disease drugs
          W_ADJ[disease,drugs][which(W_ADJ[disease,drugs]<0 & W_ADJ[disease,drugs]>th_dd_n)] = 0
          
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]>0 & W_ADJ[nano,chemical]<th_nc_p)] = 0 #nano chemical
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]<0 & W_ADJ[nano,chemical]>th_nc_n)] = 0
          
          W_ADJ[drugs,chemical][which(W_ADJ[drugs,chemical]>0 & W_ADJ[drugs,chemical]<th_dd_p)] = 0 #drugs chemical
          W_ADJ[drugs,chemical][which(W_ADJ[drugs,chemical]<0 & W_ADJ[drugs,chemical]>th_dd_n)] = 0
          
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]>0 & W_ADJ[nano,chemical]<th_drc_p)] = 0 #nano chemical
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]<0 & W_ADJ[nano,chemical]>th_drc_n)] = 0
          
          W_ADJ[disease,chemical][which(W_ADJ[disease,chemical]>0 & W_ADJ[disease,chemical]<th_dc_p)] = 0 #disease chemical
          W_ADJ[disease,chemical][which(W_ADJ[disease,chemical]<0 & W_ADJ[disease,chemical]>th_dc_n)] = 0
          
          W_ADJ[nano,nano] = 0
          W_ADJ[drugs,drugs] = 0
          W_ADJ[disease,disease] = 0
          W_ADJ[chemical,chemical] = 0
          
          W_ADJ[selected_nodes,selected_nodes] = W2[selected_nodes,selected_nodes]
          
          incProgress(1, detail = "Creating graph...")
          if(DEBUGGING)
          cat("greating graph \n")
          
          graph_gw = graph.adjacency(W_ADJ,mode="undirected",weighted=TRUE)
          if(DEBUGGING)
          cat(class(graph_gw)," ", ecount(graph_gw), " ", vcount(graph_gw), "\n")
          V(graph_gw)$type = node_type
          graph_gw = igraph::delete.vertices(graph_gw,which(igraph::degree(graph_gw)<1))
          graph_gw = igraph::simplify(graph_gw)
          
          incProgress(1, detail = "Creating subgraph of selected nodes")
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
          
          incProgress(1, detail = "Searching for cliques")
    
          mcl = cliques(graph=graph_s, min=4, max=4)
          
          if(DEBUGGING)
          cat("Nro cliques: ",length(mcl),"\n")
          
          incProgress(1, detail = "Evaluating cliques")
          
          
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
          idx = which(is_good==TRUE)
          
          #plot(cl,vertex.color=V(cl)$color)
          
          good_cliques = mcl[idx]
          
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
          
          incProgress(1, detail = "Building tables")
          
          
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
           if(DEBUGGING){
             cat("HO CREATO IL PROXY su Mi : ",class(proxy),"\n")
           }
          
        })#end progress bar
        
        output$ggplot_totale = renderPlot({
          triple_type = input$plotTripel_total
          
          Mi = MList[[1]]
          for(index in 2:length(MList)){
            Mi = rbind(Mi,MList[[index]])
          }
          
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
          
          Mi_count = Mi_count[order(Mi_count$freq),]
          nElem_c = length(columns_)
          
          d_id = input$InterestingNodes_totale
          d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
          column_index = which(c("disease","nano","drugs","chemical") %in% d_node_type)
          
          indexes_c = 1:(dim(Mi_count)[2]-1)
          indexes_c = indexes_c[-column_index]
          
          
          if(length(columns_)>2){
            nElem = round((dim(Mi_count)[1]) * 0.10,digits = 0)
            Mi_count2 = Mi_count$freq[1:nElem]
            names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep="-")
            
          }
          if(length(columns_) == 2){
            nElem = dim(Mi_count)[1]
            Mi_count2 = Mi_count$freq[1:nElem]
            names(Mi_count2) =  paste(Mi_count[1:nElem,columns_[2:nElem_c]],sep=" ")#,Mi_count[1:nElem,3],sep=" ")
            
          }
          
          index = which(Mi_count[1:nElem,column_index] %in% d_id) 
          
          validate(
            need(length(index)>0, paste("No clique for", d_id ,"with this threshold"))
          )
          
          mar.default = c(5,4,4,2) + 0.1
          par(mar = mar.default+ c(0,15,0,0))
          
          validate(
            need(length(Mi_count2[index])>0, "No correspondence for the selected criteria!")
          )
          
          barplot(as.vector(Mi_count2[index]),horiz = TRUE,col = rainbow(length(index)),main=d_id,
                  names.arg=names(Mi_count2)[index], cex.names=1, las=1)
          
        })
        
        output$ggplot = renderPlot({
          type = input$NetworkPattern #type of clique
          triple_type = input$plotTripel
          
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
          
          Mi_count = Mi_count[order(Mi_count$freq),]
          nElem_c = length(columns_)
          
          d_id = input$InterestingNodes_items
          d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
          column_index = which(c("disease","nano","drugs","chemical") %in% d_node_type)
          
          validate(
            need(input$percentuale != "", "Please select a percentage")
          )
          
if(DEBUGGING){
          cat("INPUT PERCENTUALE CAMBIATA: ",input$percentuale,"\n")
}
          if(length(columns_)>2){
            indexes_c = 1:(dim(Mi_count)[2]-1)
            indexes_c = indexes_c[-column_index]
            nElem = round((dim(Mi_count)[1]) * (input$percentuale)/100,digits = 0)
            Mi_count2 = Mi_count$freq[1:nElem]
            names(Mi_count2) =  paste(Mi_count[1:nElem,indexes_c[1]],Mi_count[1:nElem,indexes_c[2]],sep="-")
            
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
        
        output$Subnetwork_plot = renderForceNetwork({
          if(DEBUGGING){
            cat("Subnetwork plot function\n")
          }
          ADJ_toPlot = ADJ_S
          nano_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% nano)]
          chemical_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% chemical)]
          drug_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% drugs)]
          disease_s = colnames(ADJ_S)[which(colnames(ADJ_S) %in% disease)]
          
          chemMat_index = which(chemMat[,1] %in% chemical_s)
          chemMat = chemMat[chemMat_index,]
          
          d = input$InterestingNodes
          if(DEBUGGING) cat("INTERESTING NODES: ",d,"\n")
          
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
          
          di = which(selected_nodes %in% d)
          if(length(di)>0){
            pd = selected_nodes[-di]
            pd[which(pd %in% V(gto_plot)$name)]->pd
            gto_plot = delete.vertices(gto_plot,pd)
          }
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
        
        
        output$dig_dist = renderPlot({
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
          
          degree_dist = range01(igraph::degree(gto_plot))
          closeness_dist = range01(igraph::closeness(gto_plot,weights = abs(E(gto_plot)$weight)))
          betweenness_dist = range01(igraph::betweenness(gto_plot,weights = abs(E(gto_plot)$weight)))
          eccentricity_dist = range01(igraph::eccentricity(gto_plot))
          
          y_m = max(degree_dist,closeness_dist,betweenness_dist,eccentricity_dist)
          
          names(degree_dist) = V(gto_plot)$name
          
          plot(degree_dist, type="l", col="red", ylim = c(0,y_m), ylab = "",xlab="",lwd= 4)
          lines(closeness_dist, col="green",lwd= 4)
          lines(betweenness_dist, col="yellow",lwd= 4)
          lines(eccentricity_dist, col = "purple",lwd= 4)
          legend(x = "top",legend = c("Degree","Closeness","Betwenness","Eccentricity"),fill=c("red","green","yellow","purple"))
        })
        
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
          pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(items_type))))
    
        })
        
        
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
          for(i in 1:length(type)){
            if(type[i]=="nano"){
              items = c(items,nano_s)
              items_type = c(items_type,"nano")
              items_lengt = c(items_lengt,length(nano_s))
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
          
          idx_g = which(V(g)$node_type %in% "gene")
          gto_plot = igraph::induced.subgraph(graph = g,vids = c(V(g)$name[idx_g],items))
          
          gto_plot = igraph::delete.vertices(gto_plot,items)
          idx_gg = which(V(gto_plot)$name %in% V(g_geni2)$name)
          geni_toPlot = igraph::induced.subgraph(g_geni2,V(gto_plot)$name[idx_gg])
          geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
          
          data_frame = get.data.frame(x = geni_toPlot,what = "both")
          edges = data_frame$edges
          edges$value = round(abs(edges$weight * 10),digits = 0)
          colnames(edges) = c("source","target","weight","value")
          edges$value = edges$value + 0.2
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
          chem_group = input$sub_ChemicalGroup
          
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
          
          idx_g = which(V(g)$node_type %in% "gene")
          gto_plot = igraph::induced.subgraph(graph = g,vids = c(V(g)$name[idx_g],items))
          
          gto_plot = igraph::delete.vertices(gto_plot,items)
          idx_gg = which(V(gto_plot)$name %in% V(g_geni2)$name)
          geni_toPlot = igraph::induced.subgraph(g_geni2,V(gto_plot)$name[idx_gg])
          geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
          
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
          pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(vertices$group))))
          
        })
        
        
      }) #End Listener 1
      
      #ANGELA
      observeEvent(input$Go2, {
        output$info2_1 <- renderUI({
          HTML(info_text)
        }) 
        
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
        
       
        CLIQUE_TYPE = input$clique_type
        if(DEBUGGING)
        cat("CLIQUE_TYPE", CLIQUE_TYPE, "\n")
        
        validate(
          need(input$percentuale_somma != "", "Please select the treshold")
        )
        
        validate(
          need(input$nroCliques != "", "Please select the number of nodes")
        )
        
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
        
        query_th = input$percentuale_somma
        if(DEBUGGING)
        cat("query_th", query_th, "\n")
        nElem_cliques = input$nroCliques
        if(DEBUGGING)
        cat("nElem_cliques", nElem_cliques, "\n")
        
        query_nodes = c(nano_query,drug_query,chemical_query,disease_query)
        
        if(DEBUGGING)
        cat("query_nodes:",query_nodes,"\n")
        
        for(i in query_nodes){
          disease_list[[i]] = i
          selected_nodes = c(selected_nodes,i)
        }
        withProgress(message = 'Progress...', min = 1,max = 9, {
          
          
          incProgress(1, detail = "Evaluating input list...")
          
          output$NodesOfInterest <- renderUI({
            selectInput("InterestingNodes",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
          })
          
          output$NodesOfInterest_totale <- renderUI({
            selectInput("InterestingNodes_totale",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
          })
          
          output$NodesOfInterest_items <- renderUI({
            selectInput("InterestingNodes_items",label = "Node of Interest",multiple = FALSE,choices = disease_list,selected = disease_list[[1]])
          })
          
          th_p = input$th_slider2/100 #threshold per cercare le cliques
          if(DEBUGGING)
          cat("th_p:",th_p,"\n")
          
          W2 = W_ADJ
          diag(W_ADJ) = 0
          
          th_n = 1-th_p
          
          incProgress(1, detail = "Thresholding...")
          
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
          if(DEBUGGING)
          cat("W --> \n")
          
          incProgress(1, detail = "Removing edges under threshold...")
          
          
          W_ADJ[nano,disease][which(W_ADJ[nano,disease]>0 & W_ADJ[nano,disease]<th_ndis_p)] = 0 #nano disease
          W_ADJ[nano,disease][which(W_ADJ[nano,disease]<0 & W_ADJ[nano,disease]>th_ndis_n)] = 0
          
          W_ADJ[nano,drugs][which(W_ADJ[nano,drugs]>0 & W_ADJ[nano,drugs]<th_nd_p)] = 0 #nano drugs
          W_ADJ[nano,drugs][which(W_ADJ[nano,drugs]<0 & W_ADJ[nano,drugs]>th_nd_n)] = 0
          
          W_ADJ[disease,drugs][which(W_ADJ[disease,drugs]>0 & W_ADJ[disease,drugs]<th_dd_p)] = 0 #disease drugs
          W_ADJ[disease,drugs][which(W_ADJ[disease,drugs]<0 & W_ADJ[disease,drugs]>th_dd_n)] = 0
          
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]>0 & W_ADJ[nano,chemical]<th_nc_p)] = 0 #nano chemical
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]<0 & W_ADJ[nano,chemical]>th_nc_n)] = 0
          
          W_ADJ[drugs,chemical][which(W_ADJ[drugs,chemical]>0 & W_ADJ[drugs,chemical]<th_dd_p)] = 0 #drugs chemical
          W_ADJ[drugs,chemical][which(W_ADJ[drugs,chemical]<0 & W_ADJ[drugs,chemical]>th_dd_n)] = 0
          
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]>0 & W_ADJ[nano,chemical]<th_drc_p)] = 0 #nano chemical
          W_ADJ[nano,chemical][which(W_ADJ[nano,chemical]<0 & W_ADJ[nano,chemical]>th_drc_n)] = 0
          
          W_ADJ[disease,chemical][which(W_ADJ[disease,chemical]>0 & W_ADJ[disease,chemical]<th_dc_p)] = 0 #disease chemical
          W_ADJ[disease,chemical][which(W_ADJ[disease,chemical]<0 & W_ADJ[disease,chemical]>th_dc_n)] = 0
          
          W_ADJ[nano,nano] = 0
          W_ADJ[drugs,drugs] = 0
          W_ADJ[disease,disease] = 0
          W_ADJ[chemical,chemical] = 0
          
          W_ADJ[query_nodes,query_nodes] = W2[query_nodes,query_nodes]
          if(DEBUGGING)
          cat("W2 --> \n")
          
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
          graph_s = set.edge.attribute(graph_s,name = "color",index = E(graph_s)[which(sign(E(graph_s)$weight)<0)],value = "darkgreen")
          graph_s = set.edge.attribute(graph_s,name = "color",index = E(graph_s)[which(sign(E(graph_s)$weight)>0)],value = "red")
          
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
            
            Mi_count = Mi_count[order(Mi_count$freq),]
            nElem_c = length(columns_)
            
            d_id = input$InterestingNodes_items
            d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
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
            d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
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
            d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
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
            d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
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
            d_node_type = get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
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
          
          
          plot(g_,vertex.color = V(g_)$color,
               vertex.size = 50,edge.width = abs(E(g_)$weight)+2,
               vertex.label.color = "black")
          legend(x = "bottom",legend = c("Positive Correlation","Negative Correlation"),fill = c("red","darkgreen"))
          
        })
        
        output$Subnetwork_plot = renderForceNetwork({ #in go2
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
        
        
        output$dig_dist = renderPlot({
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
          
          degree_dist = range01(igraph::degree(gto_plot))
          closeness_dist = range01(igraph::closeness(gto_plot,weights = abs(E(gto_plot)$weight)))
          betweenness_dist = range01(igraph::betweenness(gto_plot,weights = abs(E(gto_plot)$weight)))
          eccentricity_dist = range01(igraph::eccentricity(gto_plot))
          
          y_m = max(degree_dist,closeness_dist,betweenness_dist,eccentricity_dist)
          
          names(degree_dist) = V(gto_plot)$name
          
          plot(degree_dist, type="l", col="red", ylim = c(0,y_m), ylab = "",xlab="",lwd= 4)
          lines(closeness_dist, col="green",lwd= 4)
          lines(betweenness_dist, col="yellow",lwd= 4)
          lines(eccentricity_dist, col = "purple",lwd= 4)
          legend(x = "top",legend = c("Degree","Closeness","Betwenness","Eccentricity"),fill=c("red","green","yellow","purple"))
        })
        
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
          pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(items_type))))
          
        })
        
        
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
          
          idx_g = which(V(g)$node_type %in% "gene")
          gto_plot = igraph::induced.subgraph(graph = g,vids = c(V(g)$name[idx_g],items))
          
          gto_plot = igraph::delete.vertices(gto_plot,items)
          idx_gg = which(V(gto_plot)$name %in% V(g_geni2)$name)
          geni_toPlot = igraph::induced.subgraph(g_geni2,V(gto_plot)$name[idx_gg])
          geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
          
          data_frame = get.data.frame(x = geni_toPlot,what = "both")
          edges = data_frame$edges
          edges$value = round(abs(edges$weight * 10),digits = 0)
          colnames(edges) = c("source","target","weight","value")
          edges$value = edges$value + 0.2
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
          chem_group = input$sub_ChemicalGroup
          
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
          
          idx_g = which(V(g)$node_type %in% "gene")
          gto_plot = igraph::induced.subgraph(graph = g,vids = c(V(g)$name[idx_g],items))
          
          gto_plot = igraph::delete.vertices(gto_plot,items)
          idx_gg = which(V(gto_plot)$name %in% V(g_geni2)$name)
          geni_toPlot = igraph::induced.subgraph(g_geni2,V(gto_plot)$name[idx_gg])
          geni_toPlot = igraph::delete.vertices(geni_toPlot,which(igraph::degree(geni_toPlot)<1))
          
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
          pie(slices, labels = lbls, main="Network statistics",col = rainbow(length(table(vertices$group))))
          
        })
        
        
      }) #End Listener 2
      
      #Outside Listener
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

     
      
      output$geneNetwork = renderForceNetwork({
        
        groups_path = input$Patway_g
        validate(
          need(input$Patway_g != "", "Please select a Pathway")
        )
      
        if(DEBUGGING)
        cat("Number of groups ",length(groups_path),"\n")
        
        good_index = which(V(g_geni2)$group %in% groups_path)
        
        geni_toPlot = igraph::induced.subgraph(graph = g_geni2,vids = V(g_geni2)$name[good_index])
        
        
        geni_toPlot = delete.vertices(graph = geni_toPlot,v = which(igraph::degree(geni_toPlot)<1))
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
        
        MyClickScript <- 
          '      d3.select(this).select("circle").transition()
        .duration(750)
        .attr("r", 30)'
        
        MyClickScript2 <- 
          'd3.select(this).select(dLink.target).select("circle").transition().duration(750).attr("r",30),
        d3.select(dLink.target).select("circle").transition().duration(750).attr("r",30)
        '
        
        forceNetwork(Links = edges, Nodes = vertices,
                     Source = "source", Target = "target",
                     Value = "value", NodeID = "name",Nodesize="size",
                     zoom = TRUE,opacity = 0.85,fontSize = 10,Group = "group",
                     legend = TRUE, height = input$gene_width,width =input$gene_height, 
                     clickAction = MyClickScript,charge = -input$gene_repulserad,
                     linkDistance = JS(paste0("function(d){return d.value*",input$gene_length,"}")))
      })
      
      
    }
      
  }) #end event login
}) #end server




