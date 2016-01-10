
check_login = function(output,loggedIN){
  if(!loggedIN){
    output$checkLOGIN = renderUI({
      HTML("<strong>Please login to proceed with the analisys! <strong/> <br/>")
    })
    
    output$checkLOGIN_free = renderUI({
      HTML("<strong>Please login to proceed with the analisys! <strong/> <br/>")
    })
    
    output$checkLOGIN_gene = renderUI({
      HTML("<strong>Please login to visualize the gene network! <strong/> <br/>")
    })
    
    output$checkLOGIN_geneQuery = renderUI({
      HTML("<strong>Please login to start the analisys! <strong/> <br/>")
    })
    
    output$checkLOGIN_entities = renderUI({
      HTML("<strong>Please login to visualize the network of entities! <strong/> <br/>")
    })
    
    output$extimatedTime = renderUI({
      HTML("How to estimate iteration:  O(n^k * k^2) where n is the number of nodes and k is the clique size.<br/>")
    })
  }else{
    output$checkLOGIN = renderUI({
      HTML(" ")
    })
    
    output$checkLOGIN_free = renderUI({
      HTML(" ")
    })
    
    output$checkLOGIN_gene = renderUI({
      HTML(" ")
    })
    
    output$checkLOGIN_geneQuery = renderUI({
      HTML(" ")
    })
    
    output$checkLOGIN_entities = renderUI({
      HTML(" ")
    })
  }
}

