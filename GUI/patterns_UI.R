output_cliques_pattern = function(input,output,cliques_groups){
  output$NetworkPattern <- renderUI({
    cliques_LL = list()
    
    for(i in 1:length(cliques_groups)){
      cliques_LL[[paste0("Type",i)]] = paste0("M",i)
    }
    selectInput("NetworkPattern",label = "Select Pattern",choices = cliques_LL,selected = "M1")
  })
  
}


building_table_and_proxy = function(input,output,MM_list){
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
  
#   output$
#   selectInput("",label = "Plot of Association Frequencies",
#               choices = list("Disease-Nano-Drug" = 1,
#                              "Disease-Nano-Chemical"= 2,
#                              "Disease-Nano"=3,
#                              "Disease-Drug"=4,
#                              "Disease-Chemical"=5,
#                              "Chemical-Nano" = 6,
#                              "Chemical-Drug" = 7,
#                              "Nano-Drug" = 8),selected = 1)
  
  proxy = dataTableProxy("clique_data_table")
  
  return(proxy)
}