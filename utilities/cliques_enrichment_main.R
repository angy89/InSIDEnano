#Function that perform enrichment analisys for a list of cliques
# @param clique_list is the list of cliques
# @gene_sets_list is the list of different gene sets for enrichment analysis
library(GSA)
library(igraph)

load("~/INSIdEapp/big_net_with_chemical_up_down80_2_th_30.RData")
load("~/INSIdEapp/items_gene_association_complete_no_genes.RData")
source("~/InsideNano/cliques_enrichment.R")

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

# gc1 = unique(unlist(c1_file$genesets))
# gcKEGG = unique(unlist(KEGG_file$genesets))
# gbiocarta = unique(unlist(biocarta_file$genesets))
# greactome = unique(unlist(reactome_file$genesets))
# gc3Mir = unique(unlist(c3Mir_file$genesets))
# gc3Tft = unique(unlist(c3Tft_file$genesets))
# gc4 = unique(unlist(c4_file$genesets))
# gc5 = unique(unlist(c5_file$genesets))
# gc6 = unique(unlist(c5_file$genesets))
# gc7 = unique(unlist(c7_file$genesets))
# 
# MSIg_genes = union(union(union(union(union(union(union(union(union(gc1,gcKEGG),gbiocarta),greactome),gc3Mir),gc3Tft),gc4),gc5),gc6),gc7)
# 
# 
# c1_list = list()
# for(i in 1:length(c1_file$genesets)){
#   gsi = c1_file$genesets[[i]]
#   c1_list[[c1_file$geneset.names[[i]]]] = MSIg_genes %in% gsi
# }
# 
# KEGG_list = list()
# for(i in 1:length(KEGG_file$genesets)){
#   gsi = KEGG_file$genesets[[i]]
#   KEGG_list[[KEGG_file$geneset.names[[i]]]] = MSIg_genes %in% gsi
# }
# 
# Gene_sets_list = list(c1_list,KEGG_list)
# MSIgDB = list(gene.names = MSIg_genes,gene.tree= Gene_sets_list)

vds = c("Asthma","MWCNT","zomepirac","1-aminopyrene")
#vds = c("Asthma","MWCNT","zomepirac")

th_l = rep(99/100,length(vds))
sets = c(2,3)
M = matrix(1,length(vds),length(vds))
M[1,3] = -1
M[3,1] = -1
M[1,2] = -1
M[2,1] = -1
rownames(M) = colnames(M) = vds
diag(M) = 0
g_ = graph.adjacency(M,mode="undirected",weighted = TRUE)

clique_list = list(vds)
gene_sets_list = list(c1_file,KEGG_file,biocarta_file,reactome_file,
                      c3Mir_file,c3Tft_file,c4_file,c5_file,c6_file,c7_file)
gene_sets_name = c("c1","c2_KEGG","c2_biocarta","c2_reactome","c3_miRNA","c3_TFT","c4","c5","c6","c7")

# gene_sets_list = list(KEGG_file,biocarta_file)
# gene_sets_name = list("KEGG","biocarta")
system.time({
  DF = cliques_enrichment(clique_list,g,g_,gene_sets_list[sets],gene_sets_name[sets],th_l = th_l,items_list)
})

  

edges = DF$edges
vertices = DF$vertices
colnames(vertices)[3]="size"

validate(need(nrow(vertices)>4,"No enrichmet for these sets and threshold"))

nanocluster(Links = edges, Nodes = vertices,
            Source = "source", Target = "target",cluster_group = "cliques",
            last_level = 3, groups = names(table(vertices$group)),
            Value = "value", NodeID = "name",Nodesize = "size",
            Group = "group",zoom = TRUE,opacity = 0.95,fontSize = 20,
            legend = TRUE,
            charge = -500,
            linkDistance = JS(paste0("function(d){return d.value*",2,
                                     "}")))

