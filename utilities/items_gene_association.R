load("~/INSIdEapp/big_net_with_chemical_up_down80_2_th_30.RData")
load("~/INSIdEapp/entities.RData")
#source("~/InsideNano/cliques_enrichment.R")
library(GSA)

c1_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c1.all.v4.0.entrez.gmt")
KEGG_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c2.cp.kegg.v4.0.entrez.gmt")
biocarta_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c2.cp.biocarta.v4.0.entrez.gmt")
reactome_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c2.cp.reactome.v4.0.entrez.gmt")
c3Mir_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c3.mir.v4.0.entrez.gmt")
c3Tft_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c3.tft.v4.0.entrez.gmt")
c4_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c4.all.v4.0.entrez.gmt")
c5_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c5.all.v4.0.entrez.gmt")
c6_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c6.all.v4.0.entrez.gmt")
c7_file = GSA.read.gmt(filename = "~/INSIdEapp/MSigDB/c7.all.v4.0.entrez.gmt")

gene_sets_list = list(c1_file,KEGG_file,biocarta_file,reactome_file,
                      c3Mir_file,c3Tft_file,c4_file,c5_file,c6_file,c7_file)
gene_sets_name = c("c1","c2_KEGG","c2_biocarta","c2_reactome","c3_miRNA","c3_TFT","c4","c5","c6","c7")
n_genes_sets = length(gene_sets_list)

sinle_group_length = unlist(lapply(X = gene_sets_list,FUN = function(elem){
  length(elem$genesets)
}))
totalSets = sum(sinle_group_length)

items = c(nano,drugs,disease,chemical)
items_list = list()
pb_index= 1
pb = txtProgressBar(min = 1,max=length(items),style = 3)

gene_attached = list()

for(i in items){
  gene_names = names(igraph::neighbors(graph = g,v = i))
  
  n_genes = length(gene_names)
  couples = rep(i,n_genes * 2)
  couples[seq(from = 2,to = length(couples),by = 2)] = gene_names
  
  idx = igraph::get.edge.ids(graph = g,vp = couples)
  edge.weigth = igraph::get.edge.attribute(graph = g,name = "weight",index = idx)
  
  gene_attached = gsub(pattern = "_at",x =gene_names,replacement = "") #find all genes connected to the item
    
  jaccard_sets_list = list()
 # gene_for_each_sets_list = list()
    
    for(j in 1:n_genes_sets){ #for each grouping (for example c1,c2...)
      group_j = gene_sets_list[[j]]$genesets
      name_j = gene_sets_list[[j]]$geneset.names
      
  #    sets_list = list()
      jaccard_sets_j = c()
      xx=1
      for(k in group_j){ #for each gene_sets (for example for each KEGG pathway)
        #evaluate the jaccard index between the group and the gene attached to the clique.
        jj = length(intersect(gene_attached,k))/length(k)
        jaccard_sets_j = c(jaccard_sets_j,jj)
   #     sets_list[[name_j[[xx]]]]= intersect(gene_attached,k)
        xx=xx+1
      }
      
      names(jaccard_sets_j) = name_j
      jaccard_sets_list[[j]] = jaccard_sets_j
    #  gene_for_each_sets_list[[j]] = sets_list
      
    }
    
 #entrez.genes.no.at = gsub(pattern = "_at",x =gene_names,replacement = "")
    items_list[[i]] = list(entrez.genes = gene_names,edge.weigth = edge.weigth,
                           jaccard_sets_list=jaccard_sets_list)#,gene_for_each_sets_list=gene_for_each_sets_list)
    
    
  setTxtProgressBar(pb,pb_index)
  pb_index = pb_index+1
}

close(pb)

save(items_list,file = "~/INSIdEapp/items_gene_association_complete_no_genes.RData")

all_genes = V(g)$name[which(V(g)$node_type=="gene")]
all_genes = gsub(pattern = "_at",replacement = "",x = all_genes)

items_gene_list = list()
pb= txtProgressBar(min=1,max=length(c(nano,drugs,chemical,disease)),style=3)
k=1
for(i in c(nano,drugs,chemical,disease)){
  items_gene_list[[i]] = all_genes %in% items_list[[i]]$gene_for_each_sets_list # fare un altro for per ogni set di geni
  setTxtProgressBar(pb,k)
  k = k+1
}
close(pb)
