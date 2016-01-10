gene_query_UI_part2 = function(input,output,gene_nano,gene_drug,gene_chem,gene_dis){

  output$gene_query_nano_input = renderUI({
    selectizeInput('gene_query_nano_input', label = "Nanomaterials", choices = c("ALL",gene_nano), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  
  output$gene_query_drug_input = renderUI({
    selectizeInput('gene_query_drug_input', label = "Drugs", choices = c(unlist(ATC_choice_list),gene_drug), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  output$gene_query_disease_input = renderUI({
    selectizeInput('gene_query_disease_input', label = "Diseases", choices = c("ALL",gene_dis), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  output$gene_query_chemical_input = renderUI({
    selectizeInput('gene_query_chemical_input', label = "Chemicals", choices = c(unlist(chemical_choice_list),gene_chem), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  output$th_slider3 = renderUI({
    sliderInput("th_slider3", label = "Strenght of Similarity", min = 10, max = 99, value = 99,step=1)
  })
  
  output$gene_query_percentuale_somma = renderUI({
    numericInput(inputId = "gene_query_percentuale_somma",label = "Number of query items connected to the selected nodes based on the threshold",
                 value = 2,min = 1,max = 10,step=1)
  })
  
  output$gene_query_nroCliques = renderUI({
    numericInput(inputId = "gene_query_nroCliques",label = "Number of query items that must be selected in the cliques",
                 value = 2,min = 1,max = 10,step=1)
  })
  
  output$gene_query_clique_type = renderUI({
    selectizeInput("gene_query_clique_type", label = "Select the type of clique you want to analyze", 
                   choices = list("All" = "ALL","Nano-Drug-Chemical-Disease" = "NDCD",
                                  "Nano-Drug-Disease" = "NDD",
                                  "Nano-Drug-Chemical" = "NDC",
                                  "Drug-Chemical-Disease" = "DCD"),
                   selected="ALL")})
  
  output$gene_query_Refresh = renderUI({
    actionButton(inputId = "gene_query_Refresh",label  = "Refresh",icon("refresh", lib = "glyphicon"))
  })
  
  output$Go4 = renderUI({
    actionButton(inputId = "Go4",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
  })
  
  
    
  
}

gene_query_UI = function(input,output,g,genes_input){
  
  genes_input = V(g)[which(V(g)$node_type == "gene")]$name
  genes_input = genes_input[1:100]
  
  if(DEBUGGING) cat("Gene list --> ",genes_input,"\n")
  
  output$gene_input = renderUI({
    selectizeInput('gene_input', label = "Gene", choices = c("ALL",genes_input), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  output$Go3 = renderUI({
    actionButton(inputId = "Go3",label  = "Start",icon("circle-arrow-right", lib = "glyphicon"))
  })
  

  return(genes_input)
  
}

