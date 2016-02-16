clique_enrichment_only_common_genes = function(clique_list,g,g_,gene_sets_list,gene_sets_name,th_l){
  ncliques = length(clique_list)
  nElem_clique = length(clique_list[[1]])
  n_genes_sets = length(gene_sets_list)
  
  sinle_group_length = unlist(lapply(X = gene_sets_list,FUN = function(elem){
    length(elem$genesets)
  }))
  totalSets = sum(sinle_group_length)
  
  Clique_jaccard = list()
  Clique_gene = list()
  
  vds = clique_list[[1]]
  gene_attached = list()
  for(i in vds){ #for each item in the clique
    entrez_at = names(igraph::neighbors(graph = g,v = i))
    gene_attached[[i]] = gsub(pattern = "_at",x = entrez_at,replacement = "") #find all genes connected to the item
    
    jaccard_sets_list = list()
    gene_for_each_sets_list = list()
    
    for(j in 1:n_genes_sets){ #for each grouping (for example c1,c2...)
      group_j = gene_sets_list[[j]]$genesets
      name_j = gene_sets_list[[j]]$geneset.names
      
      sets_list = list()
      jaccard_sets_j = c()
      xx=1
      for(k in group_j){ #for each gene_sets (for example for each KEGG pathway)
        #evaluate the jaccard index between the group and the gene attached to the clique.
        jj = length(intersect(gene_attached[[i]],k))/length(k)
        jaccard_sets_j = c(jaccard_sets_j,jj)
        sets_list[[name_j[[xx]]]]= intersect(gene_attached[[i]],k)
        xx=xx+1
      }
      
      names(jaccard_sets_j) = name_j
      jaccard_sets_list[[j]] = jaccard_sets_j
      gene_for_each_sets_list[[j]] = sets_list
      
    }
    
    Clique_jaccard[[i]] = jaccard_sets_list
    Clique_gene[[i]] = gene_for_each_sets_list
    
  }
  
  gene_sets_row_col_names = c()
  for(i in vds){
    gene_sets_row_col_names=c(gene_sets_row_col_names,names(unlist(Clique_jaccard[[i]][1:n_genes_sets])))
  }
  gene_sets_row_col_names = unique(gene_sets_row_col_names)
  
  gene_row_col_names = c()
  for(i in vds){
    gene_row_col_names=c(gene_row_col_names,unlist(Clique_gene[[i]][1:n_genes_sets]))
  }
  gene_row_col_names = unique(gene_row_col_names)
  
  
  n_nodes =  nElem_clique  + length(gene_sets_row_col_names) + length(gene_row_col_names)
  enriched_ADJ = matrix(0,n_nodes,n_nodes)
  
  rownames(enriched_ADJ) = colnames(enriched_ADJ) = c(vds,gene_sets_row_col_names,gene_row_col_names)
  
  for(i in vds){
    #set edjes between clique items and gene sets
    i_elem = unlist(Clique_jaccard[[i]][1:n_genes_sets])
    enriched_ADJ[i,names(i_elem)] = i_elem 
    enriched_ADJ[names(i_elem),i] = i_elem 
    
    #set edjes between gene sets and genes
    i_gene_sets = Clique_gene[[i]]
    for(j in 1:length(i_gene_sets)){
      gsj = i_gene_sets[[j]]
      for(k in 1:length(gsj)){
        set = names(gsj[k])
        genes = gsj[[k]]
        enriched_ADJ[set,genes] = 1
        enriched_ADJ[genes,set] = 1
        
        #         enriched_ADJ[i,genes] = 1
        #         enriched_ADJ[genes,i] = 1
        
      }
      
    }
    
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
  
  node_type=c(node_type,rep("gene",length(gene_row_col_names)))
  
  #node_type = c(node_type,rep("genes",length(unlist(Clique_jaccard[[1]][1:n_genes_sets]))))
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
  enriched_ADJ[vds,vds] = g_ADJ[vds,vds]
  
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
