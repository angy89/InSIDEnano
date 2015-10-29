remove_edges_between_nodes_by_properties = function(g, name="type",value1 = "nano", value2 = "nano",list_vertices= NULL){
  if(name =="type"){
    if(is.null(list_vertices)){
      list1 = V(g)[which(V(g)$type %in% value1)]$name
      list2 = V(g)[which(V(g)$type %in% value2)]$name
    }else{
      list1 = V(g)[which(V(g)$type %in% value1)]$name
      list2 = list_vertices
    }
  }else{
    warning("This is not a vertex attribute")
    return()
  }
  
  
  for(ni in list1){
    for(nj in  list2){
      if(ni != nj){
        idx = get.edge.ids(graph = g,vp = c(ni,nj))
        if(idx!=0){
          g = delete.edges(graph = g,edges = paste(ni,"|",nj,sep=""))
        }
      }
    }
  }
  return(g)
}

print_informazione = function(g, name="weight",value1 = "nano", value2 = "nano",list_vertices= NULL){
  
  
  if(is.null(list_vertices)){
    list1 = V(g)[which(V(g)$type %in% value1)]$name
    list2 = V(g)[which(V(g)$type %in% value2)]$name
  }else{
    list1 = V(g)[which(V(g)$type %in% value1)]$name
    list2 = list_vertices
  }
  
  
  for(ni in list1){
    for(nj in  list2){
      if(ni != nj){
        idx = get.edge.ids(graph = g,vp = c(ni,nj))
        if(idx!=0){
          w = get.edge.attribute(g,name = name,index = idx)
          cat(ni," ",nj," ",w,"\n")
        }
      }
    }
  }
  #return(g)
}
