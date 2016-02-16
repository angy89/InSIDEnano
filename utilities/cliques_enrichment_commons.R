cliques_enrichment_commons = function(clique_list,g,g_,gene_sets_list,gene_sets_name,th_l,items_list){
  if(DEBUGGING){
    message("In clique_enrichment: clique_list.length: ",length(clique_list),"\n")
    message("In clique_enrichment: vcount(g): ",vcount(g),"\n")
    message("In clique_enrichment: vcount(g_): ",vcount(g_),"\n")
    message("In clique_enrichment: gene_sets_list.length: ",length(gene_sets_list),"\n")
    message("In clique_enrichment: gene_sets_name.length: ",length(gene_sets_name),"\n")
    message("In clique_enrichment: th_l.length: ",length(th_l),"\n")
    message("In clique_enrichment: items_list.length: ",length(items_list),"\n")    
  }    
  
  ncliques = length(clique_list)
  nElem_clique = length(clique_list[[1]])
  n_genes_sets = length(gene_sets_list)
  vds = clique_list[[1]]
  
  sinle_group_length = unlist(lapply(X = gene_sets_list,FUN = function(elem){
    length(elem$genesets)
  }))
  totalSets = sum(sinle_group_length)
  
  gene_sets_row_col_names = c()
  for(i in gene_sets_list){
    gene_sets_row_col_names=c(gene_sets_row_col_names,i$geneset.names)
  }
  
  n_nodes =  nElem_clique  + length(gene_sets_row_col_names) #+ length(genes)
  enriched_ADJ = Matrix(data = 0,nrow = n_nodes,ncol = n_nodes)
  rownames(enriched_ADJ) = colnames(enriched_ADJ) = c(vds,gene_sets_row_col_names)#,genes)
  
  gc()
  
  clique_commons = items_list[[1]]$entrez.genes
  for(j in 2:length(vds)){
    clique_commons = intersect(clique_commons,items_list[[j]]$entrez.genes)
  }
  
  clique_commons = gsub(pattern = "_at",replacement = "",x = clique_commons)
  
  for(gs in 1:n_genes_sets){
    group_j = gene_sets_list[[gs]]$genesets
    for(gj in group_j){
      
    }
  }
  
  for(i in vds){
    il = items_list[[i]]
    ij = unlist(il$jaccard_sets_list)
    idx = which(names(ij) %in% gene_sets_row_col_names)
    enriched_ADJ[i,names(ij)[idx]] = ij[idx]
  }
  
  g_ADJ = igraph::get.adjacency(graph = g_,attr = "weight",sparse = FALSE)
  enriched_ADJ[vds,vds] = g_ADJ[vds,vds]
  
  node_type = rep("cliques",nElem_clique)
  
  for(i in colnames(enriched_ADJ)){
    for(j in 1:n_genes_sets){
      lj = gene_sets_list[[j]]
      if(i %in% lj$geneset.names){
        node_type = c(node_type,unlist(gene_sets_name[j]))
      }
    }
  }
  
  
  index_no_genes = 1:(nElem_clique  + length(gene_sets_row_col_names))
  for(i in 1:nElem_clique){
    val = quantile(enriched_ADJ[vds[i],index_no_genes],th_l[i])
    s_adj = enriched_ADJ[vds[i],index_no_genes]
    s_adj[s_adj<val]=0
    enriched_ADJ[vds[i],index_no_genes] = s_adj
    enriched_ADJ[index_no_genes,vds[i]] = s_adj
    
  }
  
  sums_ = colSums(enriched_ADJ[index_no_genes,index_no_genes])
  toRem = which(sums_==0)
  
  if(length(toRem)>0){
    idx = which(toRem %in% (1:nElem_clique))
    if(length(idx)>0){
      toRem = toRem[-idx]
    }
    enriched_ADJ = enriched_ADJ[-toRem,-toRem]
    node_type = node_type[-toRem]
  }
  
  enriched_ADJ = enriched_ADJ * 10
  enriched_ADJ[vds,vds] = g_ADJ[vds,vds] * 4
  
  enriched_clique_graph = graph.adjacency(enriched_ADJ,mode="undirected",weighted = TRUE)
  V(enriched_clique_graph)$node_type = node_type
  V(enriched_clique_graph)$size = igraph::degree(graph = enriched_clique_graph,v = V(enriched_clique_graph))
  
  enriched_clique_graph = igraph::delete_vertices(graph = enriched_clique_graph,v = which(V(enriched_clique_graph)$size==0))
  
  data_frame = from_igraph_to_data_frame_cluster(enriched_clique_graph,enriched_ADJ)
  
  if(DEBUGGING){
    message("In cliques_enrichment: build dataframe of enriched graph\n")
  }
  return(data_frame)
}

