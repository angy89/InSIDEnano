subgraph_nano_disease = function(g,from_nano = "MWCNT",to_disease="Asthma",drug_perc = 0.99){

  if(DEBUGGING){
    cat("Inside subgraph_nano_disease function... \n")
    cat("Query from ",from_nano, " to ",to_disease, "\n")
  }
  sp = get.all.shortest.paths(graph = g,from = from_nano,to = to_disease,weights = NA)
  
  if(DEBUGGING){
    cat("shortest paths: ",length(sp),"\n")
  }
  
  genes = c()
  for(i in 1:length(sp$res)){
    genes = c(genes,sp$res[[i]][2])
  }
  
  if(DEBUGGING){
    cat("genes: ",genes,"\n")
  }
  
    
    gsea_pvalue = gsea_res$pvalue[to_disease,-1]
    gsea_stat = gsea_res$statistica[to_disease,-1]


  entrexID = gsub(pattern = "_at",replacement = "",x = rownames(MatBig_cmap))
  symbols = lookUp(entrexID, 'org.Hs.eg', 'SYMBOL')   
  symbols = unlist(symbols)
  rownames(MatBig_cmap) = symbols

  entrexID = gsub(pattern = "_at",replacement = "",x = rownames(MatBig_chemical))
  symbols = lookUp(entrexID, 'org.Hs.eg', 'SYMBOL')   
  symbols = unlist(symbols)
  rownames(MatBig_chemical) = symbols

  entrexID = gsub(pattern = "_at",replacement = "",x = rownames(MatBig_disease))
  symbols = lookUp(entrexID, 'org.Hs.eg', 'SYMBOL')   
  symbols = unlist(symbols)
  rownames(MatBig_disease) = symbols
  
  entrexID = gsub(pattern = "_at",replacement = "",x = rownames(MatBig_nano))
  symbols = lookUp(entrexID, 'org.Hs.eg', 'SYMBOL')   
  symbols = unlist(symbols)
  rownames(MatBig_nano) = symbols

  colSums(abs(MatBig_cmap[names(genes),])) -> somme
  KTd = KTDD[from_nano,colnames(MatBig_cmap)]
  
  rank_ = somme * abs(KTd) * (abs(-log(gsea_pvalue) * sign(gsea_stat)))
  n_genes = quantile(rank_,probs = drug_perc)
  
  names(which(rank_>n_genes)) -> dd
  
  drugs_to_44_genes = MatBig_cmap[names(genes),dd]
  
  genes_nano = list()
  genes_nano_count = c()
  for(i in 1:length(dd)){
    sp = get.all.shortest.paths(graph = g,from = from_nano,to = dd[i],weights = NA)
    vect = c()
    for(j in 1:length(sp$res)){
      vect = c(vect,sp$res[[j]][2])
    }
    genes_nano[[dd[i]]] = setdiff(names(vect),names(genes))
    genes_nano_count = c(genes_nano_count,length(setdiff(names(vect),names(genes))))
  }
  
  genes_disease = list()
  genes_disease_count = c()
  for(i in 1:length(dd)){
    sp = get.all.shortest.paths(graph = g,from = to_disease,to = dd[i],weights = NA)
    vect = c()
    for(j in 1:length(sp$res)){
      vect = c(vect,sp$res[[j]][2])
    }
    genes_disease[[dd[i]]] = setdiff(names(vect),names(genes))
    genes_disease_count = c(genes_disease_count,length(setdiff(names(vect),names(genes))))
  }
  
  symbols = names(genes)
  
  ADJ = matrix(data = 0,nrow = length(genes)+length(dd) + 2,ncol = length(genes)+length(dd) + 2)
  rownames(ADJ) = colnames(ADJ) = c(symbols,dd,from_nano,to_disease)
  
  ADJ[,from_nano] = ADJ[from_nano,] = c(MatBig_nano[names(genes),from_nano],genes_nano_count,0,0)
  ADJ[,to_disease] = ADJ[to_disease,] = c(MatBig_disease[names(genes),to_disease],genes_disease_count,0,0)
  for(d in 1:length(dd)){
    ADJ[,dd[d]] = ADJ[dd[d],] = c(MatBig_cmap[names(genes),dd[d]],rep(0,length(dd)),genes_nano_count[d],genes_disease_count[d])
  }
  
  ADJ2 = ADJ
  ADJ2[ADJ2 == -1] = 2
  ADJ2[which(ADJ2[1:length(genes),to_disease]==1),to_disease] = 3
  
  type_genes = rep("down",length(genes))
  type_genes[which(ADJ2[1:length(genes),from_nano] == 1)] = "up"
  
  ADJ2[which(ADJ2[1:length(genes),from_nano]!=0),from_nano] = 3
  type = c(type_genes,rep("drug",length(dd)),"nano","disease")
  
  ADJ3 = ADJ
  n_links = rowSums(abs(ADJ3[1:length(genes),]))  
  index_to_rem = which(n_links==2)
  tab_2 = table(type[index_to_rem]) 
  
  ADJ3_ = ADJ2[-index_to_rem,-index_to_rem]
  ADJ3__ = cbind(ADJ3_, c(rep(0,(dim(ADJ3_)[1])-2),2,2))
  ADJ3__ = cbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),1,1))
  ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),2,2,0,0))
  ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),1,1,0,0))
  rownames(ADJ3__)[(dim(ADJ3__)[1]-1):dim(ADJ3__)[1]] = c(paste(tab_2["down"],"Genes"),paste(tab_2["up"],"Genes"))
  colnames(ADJ3__)[(dim(ADJ3__)[1]-1):dim(ADJ3__)[1]] = c(paste(tab_2["down"],"Genes"),paste(tab_2["up"],"Genes"))
  
  #rimuovo gli archi tra geni e nano e geni e disease
  ADJ3__[1:length(genes[-index_to_rem]),from_nano] = rep(0,length(genes[-index_to_rem]))
  ADJ3__[from_nano,1:length(genes[-index_to_rem])] = rep(0,length(genes[-index_to_rem]))
  ADJ3__[1:length(genes[-index_to_rem]),to_disease] = rep(0,length(genes[-index_to_rem]))
  ADJ3__[to_disease,1:length(genes[-index_to_rem])] = rep(0,length(genes[-index_to_rem]))
  
  #aggiungo un arco tra nano e drugs rispetto al kendall tau
  ktv = sign(KTDD[from_nano,dd])
  ktv[ktv<0]=2
  ADJ3__[from_nano,dd] = ktv
  ADJ3__[dd, from_nano] = ktv
  
  P = gsea_res$pvalue
  S = gsea_res$statistica
  
  gseav = sign(-log(P[to_disease,dd] * sign(S[to_disease,dd])))
  #gseav[is.na(gseav)]=2
  gseav[gseav<0] = 2
  ADJ3__[to_disease,dd] = gseav
  ADJ3__[dd, to_disease] = gseav
  
  type = c(type_genes[-index_to_rem],rep("drug",length(dd)),"nano","disease","down","up")
  tab = table(type)
  Col = cbind(names(tab),rainbow(length(tab)))
  rownames(Col) = Col[,1]
  Col["disease",2] = "yellow"
  Col["nano",2] = "skyblue"
  Col["up",2] = "red"
  Col["down",2] = "green"
  Col["drug",2] = "pink"
  
  ADJ3__[which(ADJ3__>10)] = ADJ3__[which(ADJ3__>10)]/10
  g2 = graph.adjacency(adjmatrix = ADJ3__,mode = "undirected",weighted = TRUE)
  
  rownames(Col) = names(tab)
  edge_color = c("red","green","black")
  
  
  g2= igraph::set.vertex.attribute(graph = g2,name = "node_type",value = type)
  g2= igraph::set.vertex.attribute(graph = g2,name = "color",value = Col[V(g2)$node_type,2])

  g2= igraph::set.edge.attribute(graph = g2,name = "color",value = edge_color[E(g2)$weight])
  
  labels = edge_color[E(g2)$weight]
  labels[which(labels != "grey")] = ""
  labels[labels %in% "grey"] = paste(E(g2)$weight[labels %in% "grey"],"Genes") 
  g2= igraph::set.edge.attribute(graph = g2,name = "labels",value = labels)

  g2=igraph::delete.vertices(g2, which(igraph::degree(g2)==0))
 
  net3 = network(as.matrix(get.adjacency(g2,names = TRUE,attr = "weight")),directed = FALSE,loops=FALSE,multiple=FALSE)
  net3 %n% "net.name" = "Query network"
  net3 %v% "col" = V(g2)$color
  net3 %v% "name" = V(g2)$name
  net3 %v% "node_type" = V(g2)$node_type
  net3 %v% "size" = igraph::degree(g2)
  
  net3 %e% "weight" = E(g2)$weight
  net3 %e% "color" = E(g2)$color
  net3 %e% "labels" = E(g2)$labels

  return(list(g2=g2,net=net3))
}

group_up_down = function(g2){
  up_index = which(V(g2)$node_type=="up")
  down_index = which(V(g2)$node_type=="down")
  
  node_degree = igraph::degree(g2)
  
}
