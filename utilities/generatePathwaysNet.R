load("~/InsideNano/gene_network_KEGG.RData")
load("~/InsideNano/big_net_with_chemical_up_down80_2.RData")
load("~/InsideNano/gene_item_association80.RData")

geneNetwork = g_geni2
genesItemsNetwork = g
library(org.Hs.eg.db)
library(GSA)
library(igraph)
KEGG_file = GSA.read.gmt(filename = "~/InsideNano/c2.cp.kegg.v5.0.entrez.gmt")

KEGG_ADJ = generatePathwaysNet(gene_list,KEGG_file)
save(KEGG_ADJ,file="/home/aserra/InsideNano/KEGG_PATH_ADJ.RData")

generatePathwaysNet = function(gene_list,KEGG_file){
  
  x <- org.Hs.egSYMBOL
  # Get the gene symbol that are mapped to an entrez gene identifiers
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])
  
  genes_id = as.numeric(unlist(lapply(gene_list,FUN = function(item){
    gsub(x = item$gene,pattern = "_at",replacement = "")
  })))
  
  Items_for_KEGG_pathway = lapply(KEGG_file$genesets,FUN = function(KEGG_set){
    KEGG_set_id = as.numeric(KEGG_set)
    KEGG_gene_sublist = gene_list[which(genes_id %in% KEGG_set_id)]
    KEGG_items = unique(unlist(lapply(KEGG_gene_sublist,FUN = function(item){
      n_item = length(item$item)
      if(n_item > 1){
        item$item[2:n_item]
      }
    })))
  })
  
  KEGG_ADJ = matrix(0,length(Items_for_KEGG_pathway),length(Items_for_KEGG_pathway))
  rownames(KEGG_ADJ) = colnames(KEGG_ADJ) = KEGG_file$geneset.names
  
  pb = txtProgressBar(min=1,max = length(Items_for_KEGG_pathway),style = 3)
  for(i in 1:length(Items_for_KEGG_pathway)){
    for(j in 1:length(Items_for_KEGG_pathway)){
      KEGG_ADJ[i,j] = length(intersect(Items_for_KEGG_pathway[[i]],Items_for_KEGG_pathway[[j]]))/length(union(Items_for_KEGG_pathway[[i]],Items_for_KEGG_pathway[[j]]))
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)
  diag(KEGG_ADJ) = 0

  return(KEGG_ADJ)
  
}