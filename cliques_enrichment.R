#Function that perform enrichment analisys for a list of cliques
# @param clique_list is the list of cliques
# @gene_sets_list is the list of different gene sets for enrichment analysis
library(GSA)

load("/home/INSIdEapp/big_net_with_chemical_up_down80_2_th_30.RData")

c1_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c1.all.v4.0.entrez.gmt")
KEGG_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c2.cp.kegg.v4.0.entrez.gmt")
biocarta_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c2.cp.biocarta.v4.0.entrez.gmt")
reactome_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c2.cp.reactome.v4.0.entrez.gmt")
c3Mir_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c3.mir.v4.0.entrez.gmt")
c3Tft_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c3.tft.v4.0.entrez.gmt")
c4_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c4.all.v4.0.entrez.gmt")
c5_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c5.all.v4.0.entrez.gmt")
c6_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c6.all.v4.0.entrez.gmt")
c7_file = GSA.read.gmt(filename = "/home/aserra/interaction_network/MSigDB/c7.all.v4.0.entrez.gmt")

vds = c("Asthma","MWCNT","zomepirac","1-aminopyrene")
clique_list = list(vds)
gene_sets_list = list(c1_file,KEGG_file,biocarta_file,reactome_file,
                      c3Mir_file,c3Tft_file,c4_file,c5_file,c6_file,c7_file)
gene_sets_name = c("c1","c2_KEGG","c2_biocarta","c2_reactome","c3_miRNA","c3_TFT","c4","c5","c6","c7")
cliques_enrichment(clique_list,g,gene_sets_list,gene_sets_name,th = 0.99)
  

cliques_enrichment = function(clique_list,g,gene_sets_list,gene_sets_name,th = 0.99){
  ncliques = length(clique_list)
  n_genes_sets = length(gene_sets_list)
  
  sinle_group_length = unlist(lapply(X = gene_sets_list,FUN = function(elem){
    length(elem$genesets)
  }))
  totalSets = sum(sinle_group_length)
  
  Clique_jaccard = list()
  for(i in 1:ncliques){
    vds = clique_list[[i]]
    
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
  
  #n_nodes =  ncliques + n_genes_sets + totalSets
  n_nodes =  ncliques  + totalSets
  enriched_ADJ = matrix(0,n_nodes,n_nodes)
  
  col_row_names = c()
  for(i in 1:ncliques){
    col_row_names = c(col_row_names,paste(clique_list[[i]],collapse="_"))
  }
  
#   for(i in 1:n_genes_sets){
#     col_row_names = c(col_row_names,paste("Group ",i,sep=""))
#   }
  
  for(i in 1:ncliques){
    col_row_names=c(col_row_names,names(unlist(Clique_jaccard[[i]][1:n_genes_sets])))
  }
  
  rownames(enriched_ADJ) = colnames(enriched_ADJ) = col_row_names
  
  #enriched_ADJ[1:ncliques,(ncliques+1):(ncliques + n_genes_sets)] = 1
  #enriched_ADJ[(ncliques+1):(ncliques + n_genes_sets),1:ncliques] = 1
  
  for(i in 1:ncliques){
    enriched_ADJ[i,(ncliques +1):n_nodes] = unlist(Clique_jaccard[[i]][1:n_genes_sets])
    enriched_ADJ[(ncliques +1):n_nodes,i] = unlist(Clique_jaccard[[i]][1:n_genes_sets])
    
  }

  
#   group_idx = list()
#   start = ncliques + n_genes_sets + 1
#   for(i in 1:n_genes_sets){
#     end = start + length(Clique_jaccard[[1]][[i]]) - 1  
#     group_idx[[i]] = c(start,end)
#     start = end + 1
#   }
#   i = 1
#   for(j in group_idx){
#      start = j[1]
#      end = j[2]
#       enriched_ADJ[ncliques + i,start:end]=Clique_jaccard[[1]][[i]]  
#       enriched_ADJ[start:end, ncliques + i]=Clique_jaccard[[1]][[i]]  
#       i = i+1
#   }
  
 # myCol = c(rep("red",ncliques),rep("green",n_genes_sets),rep("blue",length(unlist(Clique_jaccard[[1]][1:n_genes_sets]))))

  node_type = rep("cliques",ncliques)
#   for(i in 1:n_genes_sets){
#     node_type = c(node_type,paste("gene_set",i,sep=""))
#   }
  
    for(i in 1:n_genes_sets){
      node_type = c(node_type,rep(gene_sets_name[i],length(Clique_jaccard[[1]][[i]])))
    }
  
  #node_type = c(node_type,rep("genes",length(unlist(Clique_jaccard[[1]][1:n_genes_sets]))))
  val = quantile(enriched_ADJ[1,],th)
  enriched_ADJ[enriched_ADJ<val] = 0
  
  sums_ = colSums(enriched_ADJ)
  toRem = which(sums_==0)
  
  if(length(toRem)>0){
    enriched_ADJ = enriched_ADJ[-toRem,-toRem]
    node_type = node_type[-toRem]
  }
  
  enriched_clique_graph = graph.adjacency(enriched_ADJ,mode="undirected",weighted = TRUE)
  V(enriched_clique_graph)$node_type = node_type
  V(enriched_clique_graph)$size = igraph::degree(graph = enriched_clique_graph,v = V(enriched_clique_graph))
  
 # plot(enriched_clique_graph,vertex.label.cex = 0.01,vertex.size = 2,vertex.color= myCol)
  
  data_frame = from_igraph_to_data_frame_cluster(enriched_clique_graph,enriched_ADJ)
  edges = data_frame$edges
  vertices = data_frame$vertices
  colnames(vertices)[3]="size"
  

  #output$cluster_output<- renderNanoCluster(
    nanocluster(Links = edges, Nodes = vertices,
                Source = "source", Target = "target",cluster_group = "cliques",
                last_level = 3, groups = names(table(vertices$group)),
                Value = "value", NodeID = "name",Nodesize = "size",
                Group = "group",zoom = TRUE,opacity = 0.95,fontSize = 20,
                legend = TRUE,
                charge = -500,
                linkDistance = JS(paste0("function(d){return d.value*",2,
                                         "}")))
  #)
  
}