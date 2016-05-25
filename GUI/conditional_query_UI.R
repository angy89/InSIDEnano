conditionl_query_UI = function(input,output){
  output$info2_1 <- renderUI({
    HTML(info_text)
  }) 
  
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
    sliderInput("th_slider2", label = "Strenght of Similarity", min = 1, max = 99, value = 99,step=1)
  })
  
  output$percentuale_somma = renderUI({
    numericInput(inputId = "percentuale_somma",label = "Number of query items connected to the selected nodes based on the threshold",
                 value = 2,min = 1,max = 10,step=1)
  })
  
  output$nroCliques = renderUI({
    numericInput(inputId = "nroCliques",label = "Number of query items that must be selected in the cliques",
                 value = 2,min = 1,max = 10,step=1)
  })
  
  output$clique_type = renderUI({
    selectizeInput("clique_type", label = "Select the type of clique you want to analyze", 
                   choices = list("All" = "ALL","Nano-Drug-Chemical-Disease" = "NDCD",
                                  "Nano-Drug-Disease" = "NDD",
                                  "Nano-Drug-Chemical" = "NDC",
                                  "Drug-Chemical-Disease" = "DCD"),
                   selected="ALL")})
  
  output$Refresh = renderUI({
    actionButton(inputId = "Refresh",label  = "Refresh",icon("refresh", lib = "glyphicon"))
  })
  
  output$Go2 = renderUI({
    actionButton(inputId = "Go2",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
  })
}

conditional_query_UI_refresh = function(input,output){
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
    sliderInput("th_slider2", label = "Strenght of Similarity", min = 1, max = 99, value = 99,step=1)
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
}

conditional_query_UI_set_query_values = function(input,output,nanomaterials,drug,diseases,chemicals,th,n1,n2,type){
  output$nano_input = renderUI({
    selectizeInput('nano_input', label = "Nanomaterials", choices = c("ALL",nano), multiple = TRUE,
                   options = list(create = TRUE),selected = nanomaterials)
  })
  
  output$extimatedTime = renderUI({
    HTML(paste("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size<br/> 
               <strong> Estimated Iteration:"))
  })
  
  output$drug_input = renderUI({
    selectizeInput('drug_input', label = "Drugs", choices = c(unlist(ATC_choice_list),drugs), multiple = TRUE,
                   options = list(create = TRUE),selected = drug)
  })
  
  output$disease_input = renderUI({
    selectizeInput('disease_input', label = "Diseases", choices = c("ALL",disease), multiple = TRUE,
                   options = list(create = TRUE),selected = diseases)
  })
  
  output$chemical_input = renderUI({
    selectizeInput('chemical_input', label = "Chemicals", choices = c(unlist(chemical_choice_list),chemical), multiple = TRUE,
                   options = list(create = TRUE),selected = chemicals)
  })
  
  output$th_slider2 = renderUI({
    sliderInput("th_slider2", label = "Strenght of Similarity", min = 10, max = 99, value = th,step=1)
  })
  
  output$percentuale_somma = renderUI({
    numericInput(inputId = "percentuale_somma",label = "Number of query items connected to the selected nodes based on the threshold",
                 value = n1,min = 1,max = 10,step=1)
  })
  
  output$nroCliques = renderUI({
    numericInput(inputId = "nroCliques",label = "Number of query items that must be selected in the cliques",
                 value = n2,min = 1,max = 10,step=1)
  })
  
  output$Refresh = renderUI({
    actionButton(inputId = "Refresh",label  = "Refresh",icon("refresh", lib = "glyphicon"))
  })
  
  output$Go2 = renderUI({
    actionButton(inputId = "Go2",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
  })
  
#   output$info2_1 <- renderUI({
#     HTML(" ")
#   }) 
  
  output$clique_type = renderUI({
    selectizeInput("clique_type", label = "Cliques type", 
                   choices = list("All" = "ALL","Nano-Drug-Chemical-Disease" = "NDCD",
                                  "Nano-Drug-Disease" = "NDD",
                                  "Nano-Drug-Chemical" = "NDC",
                                  "Drug-Chemical-Disease" = "DCD"),
                   selected=type)})
  }

