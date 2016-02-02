#Function that perform enrichment analisys for a list of cliques
# @param clique_list is the list of cliques
# @gene_sets_list is the list of different gene sets for enrichment analysis

vds = c("Asthma","MWCNT","zomepirac","1-aminopyrene")
gene_sets_list = list(KEGG_file)

cliques_enrichment = function(clique_list,g,gene_sets_list){
  ncliques = length(clique_list)
  n_genes_sets = length(gene_sets_list)
  
  totalSets = unlist(lapply(lapply(X = gene_sets_list,FUN = function(elem){
    length(elem$genesets)
  }),FUN = sum))
  
  Clique_jaccard = list()
  for(i in 1:ncliques){
    vds = clique_list[[1]]
    
    gene_attached = c()
    for(i in vds){
      gene_attached = c(gene_attached,names(igraph::neighbors(graph = g,v = i)))
    }
    gene_attached = unique(gene_attached)
    gene_attached = gsub(pattern = "_at",x = gene_attached,replacement = "")
    
    jaccard_sets_list = list()
    for(j in 1:n_genes_sets){ #for each grouping (for example c1,c2...)
        group_j = gene_sets_list[[j]]$genesets
        name_j = gene_sets_list[[j]]$geneset.names
        
        jaccard_sets_j = c()
        for(k in group_j){ #for each gene_sets (for example for each KEGG pathway)
          #evaluate the jaccard index between the group and the gene attached to the clique.
          jj = length(intersect(gene_attached,k))/length(k)
          jaccard_sets_j = c(jaccard_sets_j,jj)
        }
        
        names(jaccard_sets_j) = name_j
        jaccard_sets_list[[j]] = jaccard_sets_j
    }
    
    Clique_jaccard[[i]] = jaccard_sets_list
  }
  
  n_nodes =  ncliques + n_genes_sets + totalSets
  ADJ = matrix(0,n_nodes,n_nodes)
  
  rownames(ADJ) = colnames(ADJ) = c("Clique","KEGG",gene_sets_list[[1]]$geneset.names)
  for(i in 1:length(jaccard_sets_list)){
    ADJ[i,]
  }
  
  
}