gene_network_NANO_DISEASE_query = function(input,output,g){
  output$geneNetwork_net3 = renderPlot({
    
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
      cat("Query from ",input$gene_nano_query ,"to: ",input$gene_disease_query,"\n with th -> ",drug_perc/100,"\n v(count) and e(count)",vcount(g)," ",ecount(g),"\n")
    
    geni_toPlot = subgraph_nano_disease(g,from_nano = nano_query,to_disease=disease_query,drug_perc = drug_perc/100)
    
    net3 = geni_toPlot$net
    
    if(DEBUGGING){
      cat("Net3 class --> ",class(net3),"\n")
    }
    
    plot(net3, vertex.cex=(net3 %v% "size")/2, vertex.col=net3 %v% "col" ,label=network.vertex.names(net3),label.pad=0,edge.col = net3 %e% "color")
    legend(x = "topleft",legend = names(table(net3 %v% "node_type")),fill = c("yellow","green","pink","skyblue","red"))  
    
    
  })
  
  #       output$geneNetwork_query = renderForceNetwork({
  #         
  #         nano_query = input$gene_nano_query
  #         disease_query = input$gene_disease_query
  #         drug_perc = input$th_gene_query
  #         
  #         validate(
  #           need(input$gene_nano_query != "", "Please select a nano")
  #         )
  #         
  #         validate(
  #           need(input$gene_disease_query != "", "Please select a disease")
  #         )
  #         
  #         if(DEBUGGING)
  #           cat("Query from ",input$gene_nano_query ,"to: ",input$gene_disease_query,"\n with th -> ",drug_perc/100,"\n v(count) and e(count)",vcount(g)," ",ecount(g),"\n")
  #       
  #         geni_toPlot = subgraph_nano_disease(g,from_nano = nano_query,to_disease=disease_query,drug_perc = drug_perc/100)
  #         geni_toPlot = geni_toPlot$g2
  #         
  #         net3 = geni_toPlot$net
  #         
  #           if(DEBUGGING)
  #             cat("geni_toPlot", class(geni_toPlot), "\n")
  #         
  #           
  #           
  #           if(DEBUGGING)
  #             cat("Final graph ",ecount(geni_toPlot), " ",vcount(geni_toPlot),"\n")
  #           
  #         data_frame = get.data.frame(x = geni_toPlot,what = "both")
  #         edges = data_frame$edges
  #         edges$value = round(abs(edges$weight * 10),digits = 0)
  #         colnames(edges) = c("source","target","weight","value")
  #         
  #         vertices = data_frame$vertices
  #         vertices$size = igraph::degree(geni_toPlot)
  #         colnames(vertices) = c("name","group","type","size")
  #         
  #         
  #         for(i in 1:dim(edges)[1]){
  #           edges[i,"source"] = which(vertices[,"name"] %in% edges[i,"source"]) - 1
  #           edges[i,"target"] = which(vertices[,"name"] %in% edges[i,"target"]) - 1
  #         }
  #         
  #         vertices$name = as.factor(vertices$name)
  #         vertices$group = as.factor(vertices$group)
  #         vertices$size = as.numeric(vertices$size)
  #         vertices$type = as.factor(vertices$type)
  #         
  #         edges$source = as.integer(edges$source)
  #         edges$target  = as.integer(edges$target)
  #         edges$value = as.integer(edges$value)
  #         
  #         MyClickScript <- 
  #           '      d3.select(this).select("circle").transition()
  #         .duration(750)
  #         .attr("r", 30)'
  #         
  #         MyClickScript2 <- 
  #           'd3.select(this).select(dLink.target).select("circle").transition().duration(750).attr("r",30),
  #         d3.select(dLink.target).select("circle").transition().duration(750).attr("r",30)
  #         '
  #         
  #         forceNetwork(Links = edges, Nodes = vertices,
  #                      Source = "source", Target = "target",
  #                      Value = "value", NodeID = "name",Nodesize="size",
  #                      zoom = TRUE,opacity = 0.85,fontSize = 10,Group = "group",
  #                      legend = TRUE, 
  #                      clickAction = MyClickScript,charge = -input$th_gene_query_repulseration,
  #                      linkDistance = JS(paste0("function(d){return d.value*",input$th_gene_query_el,"}")))
  #       })
}