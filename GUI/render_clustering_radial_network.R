render_nano_collapsible_tree = function(input,output,NANO){
  hls.list =  as.collapsible.tree.list(NANO,"NANO")
  output$treeNANO = renderCollapsibleTree(
    collapsibleTree(List = hls.list)
  )
}

render_drugs_collapsible_tree = function(input,output,DRUGS){
  hls.list =  as.collapsible.tree.list(DRUGS,"DRUG")
  output$treeDRUGS = renderCollapsibleTree(
    collapsibleTree(List = hls.list)
  )
}

render_chemical_collapsible_tree = function(input,output,CHEMICAL){
  hls.list =  as.collapsible.tree.list(CHEMICAL,"CHEMICAL")
  output$treeCHEMICAL = renderCollapsibleTree(
    collapsibleTree(List = hls.list)
  )
}

render_disease_collapsible_tree = function(input,output,DISEASE){
  hls.list =  as.collapsible.tree.list(DISEASE,"DISEASE")
  output$treeDISEASE = renderCollapsibleTree(
    collapsibleTree(List = hls.list)
  )
}

render_nano_gene_topTable = function(input,output,nnode,toSave,MYSYMBOL){
  output$node_DT = DT::renderDataTable({
  message("In render_nano_gene_topTable:: nnode ",nnode,"\n")
  
  validate(need(nnode!="","Click on node name"))
  validate(need(nnode!=" ","Click on node name"))
  validate(need(input$GP!=" ","Select a gene percentage"))
  
  input$GP/100 -> percentuale
  
  n_genes = nrow(toSave[[1]]) * percentuale
  n_genes = round(n_genes/2,digits = 0)
  
  idx = which(nano %in% nnode)
  MM = toSave[[idx]]
  MM = cbind(MM,MYSYMBOL)
  logFCpValue = MM[,1] * -log(MM[,5])
  names(logFCpValue) = rownames(MM)
  logFCpValue = sort(logFCpValue,decreasing = TRUE)
  
  length(logFCpValue)->last_id
  second_last_id = (last_id - n_genes) + 1
  logFCpValue=logFCpValue[c(1:n_genes,second_last_id:last_id)]
  
  #dt_toPlot = data.frame(logFCpValue,MM[names(logFCpValue),7])
  dt_toPlot = data.frame(MM[names(logFCpValue),c(1,5,7)])
  
  colnames(dt_toPlot) = c("logFC","Adj.pValue","GeneSymbol")
  dt_toPlot$Adj.pValue = format(dt_toPlot$Adj.pValue,scientific = TRUE,nsmall = 2)
  DT::datatable(data = format(dt_toPlot,digits=3),
                options = list(scrollX=TRUE,scrollY = "400px", 
                               scrollCollapse = TRUE,paging=FALSE,
                               fixedColumns.leftColumns=1),
                escape=FALSE,rownames = FALSE,
                selection = "single") %>% formatStyle(
                  "logFC", target = "row",
                  #red, green
                  backgroundColor = styleInterval(0, c('#FFE1DC', '#E0FFDF'))
                ) 
  })
  
}

render_clustering_radial_network = function(input,output,NANO,CHEMICAL,DISEASE,DRUGS){
  output$nano_radial = renderRadialNetwork({
    as.radialNetwork(NANO) -> rnet
    radialNetwork(rnet)
  })
  
  output$disease_radial = renderRadialNetwork({
    as.radialNetwork(DISEASE) -> rnet
    radialNetwork(rnet,fontSize = 7,height = 4000,width = 4000)
  })
  
  output$drugs_radial = renderRadialNetwork({
    as.radialNetwork(DRUGS) -> rnet
    radialNetwork(rnet,fontSize = 7,height = 4000,width = 4000)
  })
  
  output$chemical_radial = renderRadialNetwork({
    as.radialNetwork(CHEMICAL) -> rnet
    radialNetwork(rnet,fontSize = 4,height = 4000,width = 4000)
  })
}