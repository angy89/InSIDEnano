free_query_UI = function(input,output){
  output$disease = renderUI({
    selectizeInput(
      'disease', label = "Select the item to start the analysis", choices = c(disease,nano,drugs,chemical), multiple = TRUE,
      options = list(create = TRUE)
    )
  })
  
  output$th_slider = renderUI({
    sliderInput("th_slider", label = "Strenght of Similarity",
                min = 0, max = 99, value = 99,step=1)
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
  
  output$Check = renderUI({
    actionButton(inputId = "Check",label  = "Check",icon("circle-arrow-right", lib = "glyphicon"))
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
}

free_query_UI_refresh = function(input,output){
  output$disease = renderUI({
    selectizeInput(
      'disease', label = "Select the item to start the analysis", choices = c(disease,nano,drugs,chemical), multiple = TRUE,
      options = list(create = TRUE)
    )
  })
  
  output$th_slider = renderUI({
    sliderInput("th_slider", label = "Strenght of Similarity",
                min = 0, max = 99, value = 99,step=1)
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
  
}

free_query_UI_node_of_interest_output = function(input,output,disease_list){

  output$NodesOfInterest <- renderUI({
    selectInput("NodesOfInterest",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
  })
  
  output$NodesOfInterest_totale <- renderUI({
    selectInput("InterestingNodes_totale",label = "Node of Interest",multiple = TRUE,choices = disease_list,selected = disease_list[[1]])
  })
  
  output$NodesOfInterest_items <- renderUI({
    selectInput("InterestingNodes_items",label = "Node of Interest",multiple = FALSE,choices = disease_list,selected = disease_list)
  })
  
}

free_query_UI_ggplot_totale = function(input,output,MList,graph_gw){
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
    d_node_type = igraph::get.vertex.attribute(graph = graph_gw,name = "type",index = d_id)
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
}



