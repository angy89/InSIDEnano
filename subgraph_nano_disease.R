subgraph_nano_disease = function(g,from_nano = "TiO2",to_disease="Alzheimer Disease",drug_perc = 0.99){
  library(org.Hs.eg.db)
  library(annotate)
  library(igraph)
  
  sp = get.all.shortest.paths(graph = g,from = from_nano,to = to_disease,weights = NA)
  
  genes = c()
  for(i in 1:length(sp$res)){
    genes = c(genes,sp$res[[i]][2])
  }
  
  #all genes between nano and disease
  
  load("./alzaimer_gsea.RData")
  gsea_stat = toRet$S[,-1]
  gsea_pvalue = toRet$M[,-1]
  
  colSums(abs(MatBig_cmap[names(genes),])) -> somme
  KTd = KTDD[from_nano,colnames(MatBig_cmap)]
  
  rank_ = somme * abs(KTd) * (abs(-log(gsea_pvalue) * sign(gsea_stat)))
  
  #n_genes = quantile(somme,probs = drug_perc)
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
  
  
  entrexID = gsub(pattern = "_at",replacement = "",x = names(genes))
  symbols = lookUp(entrexID, 'org.Hs.eg', 'SYMBOL')   
  symbols = unlist(symbols)
  
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
  #     ADJ3__ = cbind(ADJ3_, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["down"],tab_2["down"]))
  #     ADJ3__ = cbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["up"],tab_2["up"]))
  #     ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["down"],tab_2["down"],0,0))
  #     ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["up"],tab_2["up"],0,0))
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
  
  load("./gsea_drug_disease.RData")
  P = gsea_res$pvalue
  S = gsea_res$statistica
  
  gseav = sign(-log(P[to_disease,dd] * sign(S[to_disease,dd])))
  gseav[gseav<0] = 2
  ADJ3__[to_disease,dd] = gseav
  ADJ3__[dd, to_disease] = gseav
  
  type = c(type_genes[-index_to_rem],rep("drug",length(dd)),"nano","disease","down","up")
  tab = table(type)
  Col = cbind(names(tab),rainbow(length(tab)))
  rownames(Col) = Col[,1]
  Col["disease",2] = "yellow"
  Col["nano",2] = "green"
  Col["up",2] = "red"
  Col["down",2] = "green"
  Col["drug",2] = "purple"
  
  ADJ3__[which(ADJ3__>10)] = ADJ3__[which(ADJ3__>10)]/10
  g2 = graph.adjacency(adjmatrix = ADJ3__,mode = "undirected",weighted = TRUE)
  
  rownames(Col) = names(tab)
  #edge_color = c("red","green","black",rep("grey",200))
  edge_color = c("red","green","black")
  
  plot(g2,vertex.label.cex = 0.6,layout = layout.kamada.kawai,vertex.color = Col[type,2],edge.width = E(g2)$weight,edge.color = edge_color[E(g2)$weight])
  
  g2= set.vertex.attribute(graph = g2,name = "node_type",value = type)
  g2= set.edge.attribute(graph = g2,name = "color",value = edge_color[E(g2)$weight])
  
  labels = edge_color[E(g2)$weight]
  labels[which(labels != "grey")] = ""
  labels[labels %in% "grey"] = paste(E(g2)$weight[labels %in% "grey"],"Genes") 
  g2= set.edge.attribute(graph = g2,name = "labels",value = labels)
  
  write.graph(g2,paste(from_nano,"_",to_disease,"_with_KT_and_Gsea_already_grouped2.graphml",sep=""),"graphml")
  return(g2)
}


subgraph_nano_disease_with_chemical = function(g,from_nano = "ZnO",to_disease="Asthma",drug_perc = 0.99, chem_perc = 0.999){
  library(org.Hs.eg.db)
  library(annotate)
  library(igraph)
  
  sp = get.all.shortest.paths(graph = g,from = from_nano,to = to_disease,weights = NA)
  
  genes = c()
  for(i in 1:length(sp$res)){
    genes = c(genes,sp$res[[i]][2])
  }
  
  #XX = gsub(pattern = "_at",replacement = "",x = genes)
  #library(GOSim)
  #M = getGeneSim(XX)
  #M[is.na(M)] = 0
  #rownames(M) = colnames(M) = paste(rownames(M),"_at",sep="")
  
  #hls = hclust(d = as.dist(M),method = "complete")
  #clusters = cutree(hls,k=10)
  
  #all genes between nano and disease
  
  load("./gsea_drug_disease.RData")
  gsea_stat = gsea_res$statistica[to_disease,-1]
  gsea_pvalue = gsea_res$pvalue[to_disease,-1]
  colSums(abs(MatBig_cmap[names(genes),])) -> somme
  KTd = KTDD[from_nano,colnames(MatBig_cmap)]
  
  rank_ = somme * abs(KTd) * (abs(-log(gsea_pvalue) * sign(gsea_stat)))
  
  n_genes = quantile(rank_,probs = drug_perc)
  names(which(rank_>n_genes)) -> dd
  drugs_to_44_genes = MatBig_cmap[names(genes),dd]
  
  load("./gsea_nano_chem.RData")
  load("./jaccard_chem_dis.RData")
  
  gsea_stat = gsea_res$statistica[,from_nano]
  gsea_pvalue = gsea_res$pvalue[,from_nano]
  colSums(abs(MatBig_chemical[names(genes),names(gsea_stat)])) -> somme_chem
  
  dis_ = rownames(JACC_CHEM_DIS)
  index_dis = which(dis_ %in% to_disease)
  jaccard_val = JACC_CHEM_DIS[index_dis,]
  jaccard_val = jaccard_val[,names(gsea_stat)]
  
  rank_ = somme_chem * colSums(jaccard_val) * (abs(-log(gsea_pvalue) * sign(gsea_stat)))
  
  n_genes = quantile(somme_chem,probs = chem_perc)
  names(which(somme_chem>n_genes)) -> cc
  chem_to_44_genes = MatBig_chemical[names(genes),cc]
  
  #cerco i geni che interagiscono tra i nanomateriali e i farmaci
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
  
  #cerco i geni che interagiscono tra i nanomateriali e i chemical
  genes_nano_chem = list()
  genes_nano_chem_count = c()
  for(i in 1:length(cc)){
    sp = get.all.shortest.paths(graph = g,from = from_nano,to = cc[i],weights = NA)
    vect = c()
    for(j in 1:length(sp$res)){
      vect = c(vect,sp$res[[j]][2])
    }
    genes_nano_chem[[cc[i]]] = setdiff(names(vect),names(genes))
    genes_nano_chem_count = c(genes_nano_chem_count,length(setdiff(names(vect),names(genes))))
  }
  
  
  #cerco i geni che interagiscono tra i farmaci e la malattia
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
  
  #cerco i geni che interagiscono tra i chemical e la malattia
  genes_chem_disease = list()
  genes_disease_chem_count = c()
  for(i in 1:length(cc)){
    sp = get.all.shortest.paths(graph = g,from = to_disease,to = cc[i],weights = NA)
    vect = c()
    for(j in 1:length(sp$res)){
      vect = c(vect,sp$res[[j]][2])
    }
    genes_chem_disease[[cc[i]]] = setdiff(names(vect),names(genes))
    genes_disease_chem_count = c(genes_disease_chem_count,length(setdiff(names(vect),names(genes))))
  }
  
  entrexID = gsub(pattern = "_at",replacement = "",x = names(genes))
  symbols_ = lookUp(entrexID, 'org.Hs.eg', 'SYMBOL')   
  symbols_ = unlist(symbols_)
  
  ADJ = matrix(data = 0,nrow = length(genes)+length(dd) + length(cc) + 2,ncol = length(genes)+length(dd) +length(cc)+ 2)
  rownames(ADJ) = colnames(ADJ) = c(symbols_,dd,cc,from_nano,to_disease)
  
  #riempio la colonna del nanomateriale inserendo le interazioni tra nano e geni
  #e tra ogni nanomateriale e le medicine (geni up e down regolati da entrambi)
  #e tra il nanomateriale e i chemical (geni up o down regolati da entrambi)
  ADJ[,from_nano] = ADJ[from_nano,] = c(MatBig_nano[names(genes),from_nano],genes_nano_count,genes_nano_chem_count,0,0)
  #la stessa cosa la faccio per i disease
  ADJ[,to_disease] = ADJ[to_disease,] = c(MatBig_disease[names(genes),to_disease],genes_disease_count,genes_disease_chem_count,0,0)
  
  #per ogni medicina aggiungo le interazioni di quella medicina con i geni, nessuna connessione con le altre medicine,
  #nessuna connessione con i chemical, il numero di interazioni che ha con il nanomateriale e con la malattia
  for(d in 1:length(dd)){
    ADJ[,dd[d]] = ADJ[dd[d],] = c(MatBig_cmap[names(genes),dd[d]],rep(0,length(dd)),rep(0,length(cc)),
                                  genes_nano_count[d],genes_disease_count[d])
  }
  
  #per ogni medicina aggiungo le interazioni di quella medicina con i geni, nessuna connessione con le altre medicine,
  #nessuna connessione con i chemical, il numero di interazioni che ha con il nanomateriale e con la malattia
  for(ci in 1:length(cc)){
    ADJ[,cc[ci]] = ADJ[cc[ci],] = c(MatBig_chemical[names(genes),cc[ci]],rep(0,length(dd)),rep(0,length(cc)),
                                    genes_nano_chem_count[ci],genes_disease_chem_count[ci])
  }
  
  
  ADJ2 = ADJ
  ADJ2[ADJ2 == -1] = 2 #2 significa downregolato
  ADJ2[which(ADJ2[1:length(genes),to_disease]==1),to_disease] = 3 #voglio che tutti gli archi associati al disease siano grigi
  
  type_genes = rep("down_nano",length(genes))
  type_genes[which(ADJ2[1:length(genes),from_nano] == 1)] = "up_nano"
  
  ADJ2[which(ADJ2[1:length(genes),from_nano]!=0),from_nano] = 3
  type = c(type_genes,rep("drug",length(dd)),rep("chemical",length(cc)),"nano","disease")
  
  ADJ3 = ADJ
  n_links = rowSums(abs(ADJ3[1:length(genes),]))  
  index_to_rem = which(n_links==2)
  tab_2 = table(type[index_to_rem]) 
  
  ADJ3_ = ADJ2[-index_to_rem,-index_to_rem]
  #     ADJ3__ = cbind(ADJ3_, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["down"],tab_2["down"]))
  #     ADJ3__ = cbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["up"],tab_2["up"]))
  #     ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["down"],tab_2["down"],0,0))
  #     ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),tab_2["up"],tab_2["up"],0,0))
  ADJ3__ = cbind(ADJ3_, c(rep(0,(dim(ADJ3_)[1])-2),0,0))
  ADJ3__ = cbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),0,0))
  ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),0,0,0,0))
  ADJ3__ = rbind(ADJ3__, c(rep(0,(dim(ADJ3_)[1])-2),0,0,0,0))
  rownames(ADJ3__)[(dim(ADJ3__)[1]-1):dim(ADJ3__)[1]] = c(paste(tab_2["down_nano"],"Genes"),paste(tab_2["up_nano"],"Genes"))
  colnames(ADJ3__)[(dim(ADJ3__)[1]-1):dim(ADJ3__)[1]] = c(paste(tab_2["down_nano"],"Genes"),paste(tab_2["up_nano"],"Genes"))
  
  #gene_clusters = clusters[-index_to_rem]
  
  type = c(type_genes[-index_to_rem],rep("drug",length(dd)),rep("chemical",length(cc)),"nano","disease","down","up")
  
  node_labels = rep("",length(type))
  node_labels[1:length(type_genes)]=type_genes
  
  toRem = which(colnames(ADJ3__) %in% "NA Genes")
  if(length(toRem)>0){
    ADJ3__ = ADJ3__[-toRem,-toRem]
    type = type[-toRem]
  }
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
  
  #aggiungo arco tra drugs e disease rispetto a gsea
  load("./gsea_drug_disease.RData")
  P = gsea_res$pvalue
  S = gsea_res$statistica
  
  gseav = sign(-log(P[to_disease,dd]) * sign(S[to_disease,dd]))
  gseav[gseav<0] = 2
  ADJ3__[to_disease,dd] = gseav
  ADJ3__[dd, to_disease] = gseav
  
  #aggiungo un arco tra chem e nano rispetto a gseas
  
  load("./gsea_nano_chem.RData")
  P = gsea_res$pvalue
  S = gsea_res$statistica
  
  gseav = sign(-log(P[cc,from_nano]) * sign(S[cc,from_nano]))
  gseav[gseav<0] = 2
  ADJ3__[from_nano,cc] = gseav
  ADJ3__[cc, from_nano] = gseav
  
  #aggiungo un arco tra chem e disease rispetto a jaccard
  
  load("./jaccard_chem_dis.RData")
  dis_ = rownames(JACC_CHEM_DIS)
  index_dis = which(dis_ %in% to_disease)
  jaccard_val = JACC_CHEM_DIS[index_dis,cc]
  #jaccard_val = jaccard_val[,names(gsea_stat)]
  jaccard_val = jaccard_val[1,] + (-1 * jaccard_val[2,])
  jaccard_val = sign(jaccard_val)
  
  jaccard_val[jaccard_val<0] = 2
  ADJ3__[to_disease,cc] = jaccard_val
  ADJ3__[cc, to_disease] = jaccard_val
  
  tab = table(type)
  Col = cbind(names(tab),rainbow(length(tab)))
  rownames(Col) = Col[,1]
  Col["chemical",2]="black"
  Col["disease",2] = "yellow"
  Col["nano",2] = "green"
  Col["up_nano",2] = "red"
  Col["down_nano",2] = "lightblue"
  Col["up",2] = "red"
  Col["down",2] = "lightblue"
  Col["drug",2] = "purple"
  
  
  #ADJ2[which(ADJ2[1:length(genes),from_nano]!=0),from_nano] = 3 #voglio che gli archi associati al nanomateriale siano grigi
  
  ADJ3__[which(ADJ3__>10)] = ADJ3__[which(ADJ3__>10)]/10
  g2 = graph.adjacency(adjmatrix = ADJ3__,mode = "undirected",weighted = TRUE)
  
  rownames(Col) = names(tab)
  #edge_color = c("red","green","black",rep("grey",200))
  edge_color = c("red","green","black")
  
  plot(g2,vertex.label.cex = 0.6,layout = layout.fruchterman.reingold,vertex.color = Col[type,2],edge.width = E(g2)$weight,edge.color = edge_color[E(g2)$weight])
  
  g2= set.vertex.attribute(graph = g2,name = "node_type",value = type)
  g2= set.edge.attribute(graph = g2,name = "color",value = edge_color[E(g2)$weight])
  
  labels = edge_color[E(g2)$weight]
  labels[which(labels != "grey")] = ""
  labels[labels %in% "grey"] = paste(E(g2)$weight[labels %in% "grey"],"Genes") 
  g2= set.edge.attribute(graph = g2,name = "labels",value = labels)
  
  #g2 = set.vertex.attribute(graph = g2,name="gene_cluster",value = node_labels)
  
  write.graph(g2,paste(from_nano,"_",to_disease,"_chem_KT_Gsea_already_grouped.graphml",sep=""),"graphml")
  
  return(g2) 
  #   g2 = graph.adjacency(adjmatrix = ADJ,mode = "undirected",weighted = TRUE)
  #   g3 = graph.adjacency(adjmatrix = ADJ2,mode = "undirected",weighted = TRUE)
  #   
  #   rownames(Col) = names(tab)
  #   edge_color = c("red","green","grey")
  #   
  #   plot(g2,vertex.label.cex = 0.6,layout = layout.circle,vertex.color = Col[type,2],edge.width = E(g3)$weight,edge.color = edge_color[E(g3)$weight])
  #   
  #   g2= set.vertex.attribute(graph = g2,name = "node_type",value = type)
  #   g2= set.edge.attribute(graph = g2,name = "color",value = edge_color[E(g3)$weight])
  #   g2= set.edge.attribute(graph = g2,name = "label",value = edge_color[E(g3)$weight])
  #   
  #   write.graph(g2,paste(from_nano,"_",to_disease,"_with_chem.graphml",sep=""),"graphml")
  #   
  #   return(g2)
}

# shortest_path_subgrap = function(to = "ZnO",from = "Asthma",layout = layout.fruchterman.reingold){
#   cat("Calculating shortest path...\n")
#   sp = get.all.shortest.paths(graph = g,to = to, from = from,weights = NA,output = "both")
#     
#   node_name = c()
#   node_type = c()
#   for(i in 1:length(sp$vpath)){
#     for(j in 1:length(sp$vpath[[i]])){
#       node_name = c(node_name, sp$vpath[[i]][j])
#       node_type = c(node_type,get.vertex.attribute(graph = g,name = "node_type",index = sp$vpath[[i]][j]))
#     }  
#   }
#   
#   connections = list()
#   for(i in 1:length(node_type)){
#     if(node_type[i] =="nano"){
#       connections[[names(node_name[i])]] = nano_degree[names(node_name[i])]
#     }
#     if(node_type[i] =="drug"){
#       connections[[names(node_name[i])]] = drug_degree[names(node_name[i])]
#     }
#     if(node_type[i] =="disease"){
#       connections[[names(node_name[i])]] = disease_degree[names(node_name[i])]
#     }
#     if(node_type[i] =="gene"){
#       connections[[names(node_name[i])]] = gene_degree[[names(node_name[i])]]
#     }
#   }
#   
#   sub_g = induced.subgraph(graph = g, v=names(node_name))
#   ADJ = get.adjacency(sub_g)
#   
#   n = dim(ADJ)[1]
#   plus = 4
#   
#   ADJ = cbind(ADJ,matrix(data=0,nrow = n,ncol=plus))
#   colnames(ADJ)[(n+1):(n+plus)] = c("nDrug","nNano","nGene","nDisease")
#   
#   ADJ = rbind(ADJ,matrix(data=0,nrow = plus,ncol=n+plus))
#   rownames(ADJ)[(n+1):(n+plus)] = c("nDrug","nNano","nGene","nDisease")
#   
#   k = 1
#   for(i in names(node_name)){
#     con = connections[[i]]
#     if(node_type[k] == "nano"){
#       elem = names(con[1])
#       num = con[[1]]
#       ADJ[elem,"nGene"] = num
#       ADJ["nGene",elem] = num
#     }
#     if(node_type[k] == "drug"){
#       elem = names(con[1])
#       num = con[[1]]
#       ADJ[elem,"nGene"] = num
#       ADJ["nGene",elem] = num
#     }
#     if(node_type[k] == "disease"){
#       elem = names(con[1])
#       num = con[[1]]
#       ADJ[elem,"nGene"] = num
#       ADJ["nGene",elem] = num
#     }
#     if(node_type[k] == "gene"){
#       elem = names(con$disease)
#       ADJ[elem,"nDisease"] = ADJ["nDisease",elem] = as.numeric(con$disease)
#       ADJ[elem,"nNano"] = ADJ["nNano",elem] = as.numeric(con$nano)
#       ADJ[elem,"nDrug"] = ADJ["nDrug",elem] = as.numeric(con$drug)
#     }
#     k = k+1
#   }
#   
#   tab = table(c(node_type,rep("info",4)))
#   Col = cbind(names(tab),rainbow(length(tab),start = 0.15,end=0.5))
#   
#   rownames(Col) = names(tab)
#   
#   positioning = c(names(node_name),colnames(ADJ)[(length(node_name)+1):(dim(ADJ)[2])])
#   g2 = graph.adjacency(adjmatrix = ADJ[positioning,positioning],mode = "undirected",weighted = TRUE)
#   labels = E(g2)$weight
#   labels[which(labels == 1)] = ""
#   
#   edge_size = E(g2)$weight
#   edge_size[which(edge_size %in% 10:100)] = 15
#   edge_size[which(edge_size %in% 101:1000)] = 18 
#   edge_size[which(edge_size>1000)] = 23
#   
#   plot(g2, vertex.size = 30, vertex.label.cex = 0.8,
#        vertex.color= Col[c(node_type,rep("info",4)),2],
#        edge.label = labels, edge.width = edge_size, 
#        main=paste("Shortest path between ",from, " and ",to,sep=""),
#        layout = layout)
#   
#   g2 = set.vertex.attribute(graph = g2,name = "color",value = Col[c(node_type,rep("info",4)),2])
#   
#   return(g2)
# }
