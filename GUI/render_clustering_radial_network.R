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