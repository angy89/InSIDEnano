source("./my_igraph_function.R",local = FALSE)
source("./subgraph_nano_disease.R",local = FALSE)
source("check_login.R")
source("observer.R")
source("GUI/fenotypic_network_UI.R")
source("GUI/free_query_UI.R")
source("GUI/conditional_query_UI.R")
source("GUI/gene_network_UI.R")
source("GUI/patterns_UI.R")
source("GUI/gene_query_UI.R")
source("GUI/render_clustering_radial_network.R")

source("GUI/query_outputs.R")
source("GUI/conditional_query_outputs.R")
source("GUI/browse_dataset.R")

source("query/save_free_query.R")
source("query/my_cliques_function.R")
source("query/gene_network_NANO_DISEASE_query.R")
source("query/query_utilities.R")
source("query/conditional_query_clique_search.R")

source("query/gene_query_like_conditional.R")
source("query/gene_query.R")
source("query/free_query.R")
source("query/conditional_query.R")

source("query/check_already_existing_conditional_query.R")
source("query/loa_conditional_query.R")
source("login.R")

library(network)

# #for gene_query
# load(paste(APP_PATH,"matrix_gene_disease_01_80.RData",sep=""))
# load(paste(APP_PATH,"matrix_gene_nano_01.RData",sep=""))
# load(paste(APP_PATH,"matrix_gene_cmap_01.RData",sep=""))
# load(paste(APP_PATH,"gsea_drug_disease.RData",sep=""))
# load(paste(APP_PATH,"gsea_drug_disease.RData",sep=""))
# load(paste(APP_PATH,"matrix_gene_chemical_01_inc_dec.RData",sep=""))
# load(paste(APP_PATH,"chemical_disease_gene80.RData",sep=""))
# load(paste(APP_PATH,"degree.RData",sep=""))
# load(paste(APP_PATH,"KTDD_adjacency_red.RData",sep=""))

load(paste(APP_PATH,"entities.RData",sep=""))
load(paste(APP_PATH,"nano_top_table.RData",sep=""))

x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
entrez = gsub(x = rownames(toSave[[1]]),pattern = "_at",replacement = "")
xx[entrez] -> MYSYMBOL
MYSYMBOL = unlist(MYSYMBOL)

load(paste(APP_PATH,"nano_based_clustering.RData",sep=""))

load(paste(APP_PATH,"nano_chemical_disease_drugs_hierarchical_clustering.RData",sep=""))
load(paste(APP_PATH,"join10.RData",sep=""))
join10 = unique(join10)
join10$ATC_lev1 = substr(x = join10$code,start = 1,stop = 1)
load(paste(APP_PATH,"chemicals_classes.RData",sep=""))

selected_nodes = c()
disease_list = list()
#good_cliques = list()
MM_list = list()
th_p = 0.99
info_text = " "
username = "user2015"
password = "helsinki2015"
LOGGED_IN = FALSE
proxy = NULL
DEBUGGING = TRUE
genes_input = c()


ATC_choice_list = list("ALL" = "ALL",
                       "(A) Alimentary tract and metabolism" = "A", 
                       "(C) Cardiovascular system" = "C", 
                       "(D) Dermatologicals" = "D",
                       "(G) Genito-urinary system and sex hormones" = "G",
                       "(H) Systemic hormonal preparations, excluding sex hormones and insulins" = "H",
                       "(J) Antiinfective for systemic use" = "J",
                       "(L) Antineoplastic and immunomodulating agents" = "L",
                       "(M) Musculo-skeletal system" = "M",
                       "(N) Nervous system" = "N",
                       "(P) Antiparasitics products, insecticides and repellent" = "P",
                       "(R) Respiratory system" = "R",
                       "(S) Sensory organs" = "S")

ATC_letter_vector = c("A","C","D","G","H","J","L","M","N","P","R","S")

chemical_choice_list = list("ALL" = "ALL",
                            "Amino acids" = "Amino acids",
                            "Biological factor" = "Biological factor",
                            "Carbohydrates" = "Carbohydrates",
                            "Nucleic Acids" = "Nucleic Acids",
                            "Chemical Actions" = "Chemical Actions",
                            "Complex Mixtures" = "Complex Mixtures",
                            "Enzymes and Coenzymes" = "Enzymes and Coenzymes",
                            "Heterocyclic Compounds" = "Heterocyclic Compounds",
                            "Hormones" = "Hormones",
                            "Inorganic Chemicals" = "Inorganic Chemicals",
                            "Lipids" = "Lipids",
                            "Organic Chemicals" = "Organic Chemicals",
                            "Polycyclic Compounds" = "Polycyclic Compounds")

chemical_group_vector = c("Amino acids","Biological factor","Carbohydrates","Nucleic Acids","Chemical Actions","Complex Mixtures",
                          "Enzymes and Coenzymes","Heterocyclic Compounds","Hormones","Inorganic Chemicals","Lipids",
                          "Organic Chemicals","Polycyclic Compounds")


range01 <- function(x){(x-min(x))/(max(x)-min(x))}

subnetwork_creation <- function(d,gr4,good_cliques){
  
  if(DEBUGGING){
    cat("Disease: ",d, "\n")
    cat("gr4 class: ",class(gr4),"\n")
    cat("vcount(gr4): ",vcount(gr4),"\n")
    cat("ecount(gr4): ",ecount(gr4),"\n")
    cat("Class and length good_cliques", class(good_cliques), length(good_cliques),"\n")
  }
  
  neig = lapply(X = good_cliques,FUN = function(i_cliques){
    v_names = names(i_cliques)
    if(sum(v_names %in% d)>0){
      return(v_names)
    }
  })
  
  vds = unique(unlist(neig))
  if(DEBUGGING)
  cat("length neig: ",length(vds))
  #vds = neighborhood(graph = gr4,order = 1,nodes = d,mindist = 1)
  asthma = induced_subgraph(graph = gr4,vids = vds) 
  if(DEBUGGING)
  cat("vcount(asthma): ",vcount(asthma),"\n")
  if(DEBUGGING)
  cat("ecount(asthma): ",ecount(asthma),"\n")
  
  E(asthma)$weight = abs(E(asthma)$weight)
  max_w = max(E(asthma)$weight)
  E(asthma)$weight = (max_w - E(asthma)$weight) + 1
  
  asthma = minimum.spanning.tree(graph = asthma, weight = E(asthma)$weight + max_w )
  di = which(selected_nodes %in% d)
  pd = selected_nodes[-di]
  pd[which(pd %in% V(asthma)$name)]->pd
  asthma = delete.vertices(asthma,pd)
  #asthma = igraph::delete.vertices(asthma,which(igraph::degree(asthma)<1))
  #asthma = igraph::delete.vertices(asthma,which(igraph::degree(asthma)<=1))
  #asthma = igraph::simplify(asthma)
  
  res = from_igraph_to_data_frame(asthma= asthma)
  
  return(list(edges = res$edges,vertices = res$vertices))
}

from_igraph_to_data_frame <- function(asthma= asthma){
  asthma_data = get.data.frame(x = asthma,what = "both")
  edges = asthma_data$edges
  vertices = asthma_data$vertices
  
  edges$link_color = rep("red",dim(edges)[1])
  edges$link_color[which(edges$weight<0)] = "darkgreen"
  
  if(ecount(asthma)>0){
    edges$weight = round(abs(edges$weight),digits = 0)
    colnames(edges)[1:3] = c("source","target","value")
    edges$value = edges$value + 0.2
    for(i in 1:dim(edges)[1]){
      edges[i,"source"] = which(vertices[,"name"] %in% edges[i,"source"]) - 1
      edges[i,"target"] = which(vertices[,"name"] %in% edges[i,"target"]) - 1
    }
    
    edges$link_color[which(edges$link_color=="red")] = "#FFFFFF"
    edges$link_color[which(edges$link_color=="darkgreen")] = "#FFA07A"
    
    vertices[which(vertices$color=="pink"),"color"] = 1
    vertices[which(vertices$color=="skyblue"),"color"] = 2
    vertices[which(vertices$color=="yellow"),"color"] = 3
    vertices[which(vertices$color=="violet"),"color"] = 4
    vertices$size = igraph::degree(asthma)
    colnames(vertices) = c("name","type","group","size")
    rownames(vertices) = 1:dim(vertices)[1]
    vertices$name = as.factor(vertices$name)
    vertices$group = as.integer(vertices$group)
    vertices$size = as.numeric(vertices$size)
    vertices$type = as.factor(vertices$type)
    
    edges$source = as.integer(edges$source)
    edges$target  = as.integer(edges$target)
    edges$value = as.integer(edges$value)
    edges$link_color = as.factor(edges$link_color)
    
    return(list(edges=edges,vertices = vertices))
  }else{
    cat("No weigthed edges")
  }
  
}



internal_render_plot = function(Mi,gr4,s,proxyList){
  
  if(DEBUGGING){
    cat("In internal render plot. I have proxy: ",class(proxyList),"\n")
    cat("Number of selected rows: ",length(s),"\n")
  }
  
  if(length(s)==0){
    s = 1
    
    dis = gsub(pattern = '</a>',x = Mi[s,"Disease"],replacement = "")
    #dis = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
    dis = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
    dis = gsub(pattern = '">',x = dis,replacement = " ")
    dis = unlist(strsplit(dis, split=" "))
    dis = paste(dis[-which(duplicated(dis))], collapse = ' ')  
    
    nano_ = gsub(pattern = '</a>',x = Mi[s,"Nano"],replacement = "")
    #nano_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
    nano_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
    nano_ = gsub(pattern = '">',x = nano_,replacement = " ")
    nano_ = unlist(strsplit(nano_, split=" "))
    nano_ = paste(nano_[-which(duplicated(nano_))], collapse = ' ')  
    
    drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
    drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
    drug_ = unlist(strsplit(drug_, split=" "))
    drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
    
    chem = gsub(pattern = '</a>',x = Mi[s,"Chemical"],replacement = "")
    #chem = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
    chem = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
    chem = gsub(pattern = '">',x = chem,replacement = " ")
    chem = unlist(strsplit(chem, split=" "))
    chem = paste(chem[-which(duplicated(chem))], collapse = ' ') 

    
    if(DEBUGGING){
      cat("dis: ",dis,"\n")
      cat("nano_: ",nano_,"\n")
      cat("drug_: ",drug_,"\n")
      cat("chem: ",chem,"\n")
    }
    
    g = induced_subgraph(gr4,c(dis,nano_,drug_,chem))
    
    np = which(V(g)$name %in% nano)
    ndrug = which(V(g)$name %in% drugs)
    nchem = which(V(g)$name %in% chemical)
    ndisease = which(V(g)$name %in% disease)
    
    intersect(nchem,ndrug) -> ii
    if(length(ii)>0){
      which(nchem %in% ndrug) -> index_ii
      nchem = nchem[-index_ii]
    }
    
    #V(g)$name=c("Disease","Nano","Drug","Chemical")
    V(g)$name[np] = "Nano" 
    V(g)$name[ndrug] = "Drug" 
    V(g)$name[nchem] = "Chemical" 
    V(g)$name[ndisease] = "Disease" 
    
    V(g)$color[np] = "pink" 
    V(g)$color[ndrug] = "skyblue" 
    V(g)$color[nchem] = "violet" 
    V(g)$color[ndisease] = "yellow" 
    
    E(g)$weight = 4
    return(g)
  }
  if(length(s)>1){
    l = length(s)
    s = s[l]
    selectRows(proxyList, s)
  }
  dis = gsub(pattern = '</a>',x = Mi[s,"Disease"],replacement = "")
  #dis = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
  dis = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
  dis = gsub(pattern = '">',x = dis,replacement = " ")
  dis = unlist(strsplit(dis, split=" "))
  dis = paste(dis[-which(duplicated(dis))], collapse = ' ')  
  
  nano_ = gsub(pattern = '</a>',x = Mi[s,"Nano"],replacement = "")
  #nano_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
  nano_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
  nano_ = gsub(pattern = '">',x = nano_,replacement = " ")
  nano_ = unlist(strsplit(nano_, split=" "))
  nano_ = paste(nano_[-which(duplicated(nano_))], collapse = ' ')  
  
  drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
  drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
  drug_ = unlist(strsplit(drug_, split=" "))
  drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
  
  chem = gsub(pattern = '</a>',x = Mi[s,"Chemical"],replacement = "")
  #chem = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
  chem = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
  chem = gsub(pattern = '">',x = chem,replacement = " ")
  chem = unlist(strsplit(chem, split=" "))
  chem = paste(chem[-which(duplicated(chem))], collapse = ' ') 
  
  
  g = induced_subgraph(gr4,c(dis,nano_,drug_,chem))
  
  np = which(V(g)$name %in% nano)
  ndrug = which(V(g)$name %in% drugs)
  nchem = which(V(g)$name %in% chemical)
  ndisease = which(V(g)$name %in% disease)
  
  
  
  intersect(nchem,ndrug) -> ii
  if(length(ii)>0){
    which(nchem %in% ndrug) -> index_ii
    nchem = nchem[-index_ii]
  }
  
  V(g)$color[np] = "pink" 
  V(g)$color[ndrug] = "skyblue" 
  V(g)$color[nchem] = "violet" 
  V(g)$color[ndisease] = "yellow" 
  
  return(g)    
}

internal_render_plotNDD = function(Mi,gr4,s,proxyList){
  if(length(s)==0){
    s = 1
    dis = gsub(pattern = '</a>',x = Mi[s,"Disease"],replacement = "")
    #dis = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
    dis = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
    dis = gsub(pattern = '">',x = dis,replacement = " ")
    dis = unlist(strsplit(dis, split=" "))
    dis = paste(dis[-which(duplicated(dis))], collapse = ' ')  
    
    nano_ = gsub(pattern = '</a>',x = Mi[s,"Nano"],replacement = "")
    #nano_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
    nano_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
    nano_ = gsub(pattern = '">',x = nano_,replacement = " ")
    nano_ = unlist(strsplit(nano_, split=" "))
    nano_ = paste(nano_[-which(duplicated(nano_))], collapse = ' ')  
    
    drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
    drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
    drug_ = unlist(strsplit(drug_, split=" "))
    drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
    
    g = induced_subgraph(gr4,c(dis,nano_,drug_))
    
    np = which(V(g)$name %in% nano)
    ndrug = which(V(g)$name %in% drugs)
    # nchem = which(V(g)$name %in% chemical)
    ndisease = which(V(g)$name %in% disease)
    
    #     intersect(nchem,ndrug) -> ii
    #     if(length(ii)>0){
    #       which(nchem %in% ndrug) -> index_ii
    #       nchem = nchem[-index_ii]
    #     }
    
    #V(g)$name=c("Disease","Nano","Drug","Chemical")
    V(g)$name[np] = "Nano" 
    V(g)$name[ndrug] = "Drug" 
    #V(g)$name[nchem] = "Chemical" 
    V(g)$name[ndisease] = "Disease" 
    
    V(g)$color[np] = "pink" 
    V(g)$color[ndrug] = "skyblue" 
    # V(g)$color[nchem] = "violet" 
    V(g)$color[ndisease] = "yellow" 
    
    E(g)$weight = 4
    return(g)
  }
  if(length(s)>1){
    l = length(s)
    s = s[l]
    selectRows(proxyList, s)
  }
  dis = gsub(pattern = '</a>',x = Mi[s,"Disease"],replacement = "")
  #dis = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
  dis = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
  dis = gsub(pattern = '">',x = dis,replacement = " ")
  dis = unlist(strsplit(dis, split=" "))
  dis = paste(dis[-which(duplicated(dis))], collapse = ' ')  
  
  nano_ = gsub(pattern = '</a>',x = Mi[s,"Nano"],replacement = "")
  #nano_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
  nano_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
  nano_ = gsub(pattern = '">',x = nano_,replacement = " ")
  nano_ = unlist(strsplit(nano_, split=" "))
  nano_ = paste(nano_[-which(duplicated(nano_))], collapse = ' ')  
  
  drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
  drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
  drug_ = unlist(strsplit(drug_, split=" "))
  drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
  
  g = induced_subgraph(gr4,c(dis,nano_,drug_))
  
  np = which(V(g)$name %in% nano)
  ndrug = which(V(g)$name %in% drugs)
  #  nchem = which(V(g)$name %in% chemical)
  ndisease = which(V(g)$name %in% disease)
  
  #   intersect(nchem,ndrug) -> ii
  #   if(length(ii)>0){
  #     which(nchem %in% ndrug) -> index_ii
  #     nchem = nchem[-index_ii]
  #   }
  
  V(g)$color[np] = "pink" 
  V(g)$color[ndrug] = "skyblue" 
  # V(g)$color[nchem] = "violet" 
  V(g)$color[ndisease] = "yellow" 
  
  return(g)    
}

internal_render_plotNDC = function(Mi,gr4,s,proxyList){
  if(length(s)==0){
    s = 1
    
    nano_ = gsub(pattern = '</a>',x = Mi[s,"Nano"],replacement = "")
    #nano_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
    nano_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
    nano_ = gsub(pattern = '">',x = nano_,replacement = " ")
    nano_ = unlist(strsplit(nano_, split=" "))
    nano_ = paste(nano_[-which(duplicated(nano_))], collapse = ' ')  
    
    drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
    drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
    drug_ = unlist(strsplit(drug_, split=" "))
    drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
    
    chem = gsub(pattern = '</a>',x = Mi[s,"Chemical"],replacement = "")
    #chem = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
    chem = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
    chem = gsub(pattern = '">',x = chem,replacement = " ")
    chem = unlist(strsplit(chem, split=" "))
    chem = paste(chem[-which(duplicated(chem))], collapse = ' ') 
    
    
    np = which(V(g)$name %in% nano)
    ndrug = which(V(g)$name %in% drugs)
    nchem = which(V(g)$name %in% chemical)
    #ndisease = which(V(g)$name %in% disease)
    
    intersect(nchem,ndrug) -> ii
    if(length(ii)>0){
      which(nchem %in% ndrug) -> index_ii
      nchem = nchem[-index_ii]
    }
    
    #V(g)$name=c("Disease","Nano","Drug","Chemical")
    V(g)$name[np] = "Nano" 
    V(g)$name[ndrug] = "Drug" 
    V(g)$name[nchem] = "Chemical" 
    #V(g)$name[ndisease] = "Disease" 
    
    V(g)$color[np] = "pink" 
    V(g)$color[ndrug] = "skyblue" 
    V(g)$color[nchem] = "violet" 
    #V(g)$color[ndisease] = "yellow" 
    
    E(g)$weight = 4
    return(g)
  }
  if(length(s)>1){
    l = length(s)
    s = s[l]
    selectRows(proxyList, s)
  }
  
  nano_ = gsub(pattern = '</a>',x = Mi[s,"Nano"],replacement = "")
  #nano_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
  nano_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = nano_,replacement = "")
  nano_ = gsub(pattern = '">',x = nano_,replacement = " ")
  nano_ = unlist(strsplit(nano_, split=" "))
  nano_ = paste(nano_[-which(duplicated(nano_))], collapse = ' ')  
  
  drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
  drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
  drug_ = unlist(strsplit(drug_, split=" "))
  drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
  
  chem = gsub(pattern = '</a>',x = Mi[s,"Chemical"],replacement = "")
  #chem = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
  chem = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
  chem = gsub(pattern = '">',x = chem,replacement = " ")
  chem = unlist(strsplit(chem, split=" "))
  chem = paste(chem[-which(duplicated(chem))], collapse = ' ') 
  
  g = induced_subgraph(gr4,c(nano_,drug_,chem))
  
  np = which(V(g)$name %in% nano)
  ndrug = which(V(g)$name %in% drugs)
  nchem = which(V(g)$name %in% chemical)
  #ndisease = which(V(g)$name %in% disease)
  
  intersect(nchem,ndrug) -> ii
  if(length(ii)>0){
    which(nchem %in% ndrug) -> index_ii
    nchem = nchem[-index_ii]
  }
  
  #V(g)$name=c("Disease","Nano","Drug","Chemical")
  V(g)$name[np] = "Nano" 
  V(g)$name[ndrug] = "Drug" 
  V(g)$name[nchem] = "Chemical" 
  #V(g)$name[ndisease] = "Disease" 
  
  V(g)$color[np] = "pink" 
  V(g)$color[ndrug] = "skyblue" 
  V(g)$color[nchem] = "violet" 
  #V(g)$color[ndisease] = "yellow" 
  
  return(g)    
}

internal_render_plotDCD = function(Mi,gr4,s,proxyList){
  if(length(s)==0){
    s = 1
    dis = gsub(pattern = '</a>',x = Mi[s,"Disease"],replacement = "")
    #dis = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
    dis = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
    dis = gsub(pattern = '">',x = dis,replacement = " ")
    dis = unlist(strsplit(dis, split=" "))
    dis = paste(dis[-which(duplicated(dis))], collapse = ' ')  
    
    drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
    drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
    drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
    drug_ = unlist(strsplit(drug_, split=" "))
    drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
    
    chem = gsub(pattern = '</a>',x = Mi[s,"Chemical"],replacement = "")
    #chem = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
    chem = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
    chem = gsub(pattern = '">',x = chem,replacement = " ")
    chem = unlist(strsplit(chem, split=" "))
    chem = paste(chem[-which(duplicated(chem))], collapse = ' ') 
    
    
    
    g = induced_subgraph(gr4,c(dis,drug_,chem))
    
    #np = which(V(g)$name %in% nano)
    ndrug = which(V(g)$name %in% drugs)
    nchem = which(V(g)$name %in% chemical)
    ndisease = which(V(g)$name %in% disease)
    
    intersect(nchem,ndrug) -> ii
    if(length(ii)>0){
      which(nchem %in% ndrug) -> index_ii
      nchem = nchem[-index_ii]
    }
    
    #V(g)$name=c("Disease","Nano","Drug","Chemical")
    #V(g)$name[np] = "Nano" 
    V(g)$name[ndrug] = "Drug" 
    V(g)$name[nchem] = "Chemical" 
    V(g)$name[ndisease] = "Disease" 
    
    # V(g)$color[np] = "pink" 
    V(g)$color[ndrug] = "skyblue" 
    V(g)$color[nchem] = "violet" 
    V(g)$color[ndisease] = "yellow" 
    
    E(g)$weight = 4
    return(g)
  }
  if(length(s)>1){
    l = length(s)
    s = s[l]
    selectRows(proxyList, s)
  }
  dis = gsub(pattern = '</a>',x = Mi[s,"Disease"],replacement = "")
  #dis = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
  dis = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = dis,replacement = "")
  dis = gsub(pattern = '">',x = dis,replacement = " ")
  dis = unlist(strsplit(dis, split=" "))
  dis = paste(dis[-which(duplicated(dis))], collapse = ' ')  
  
  drug_ = gsub(pattern = '</a>',x = Mi[s,"Drug"],replacement = "")
  drug_ = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  #drug_ = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = drug_,replacement = "")
  drug_ = gsub(pattern = '">',x = drug_,replacement = " ")
  drug_ = unlist(strsplit(drug_, split=" "))
  drug_ = paste(drug_[-which(duplicated(drug_))], collapse = ' ')  
  
  chem = gsub(pattern = '</a>',x = Mi[s,"Chemical"],replacement = "")
  #chem = gsub(pattern = '<a target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
  chem = gsub(pattern = '<a target=\"_blank\" href=\"https://www.google.com/\\?q=',x = chem,replacement = "")
  chem = gsub(pattern = '">',x = chem,replacement = " ")
  chem = unlist(strsplit(chem, split=" "))
  chem = paste(chem[-which(duplicated(chem))], collapse = ' ') 
  
  
  g = induced_subgraph(gr4,c(dis,drug_,chem))
  
  # np = which(V(g)$name %in% nano)
  ndrug = which(V(g)$name %in% drugs)
  nchem = which(V(g)$name %in% chemical)
  ndisease = which(V(g)$name %in% disease)
  
  intersect(nchem,ndrug) -> ii
  if(length(ii)>0){
    which(nchem %in% ndrug) -> index_ii
    nchem = nchem[-index_ii]
  }
  
  # V(g)$color[np] = "pink" 
  V(g)$color[ndrug] = "skyblue" 
  V(g)$color[nchem] = "violet" 
  V(g)$color[ndisease] = "yellow" 
  
  return(g)    
}

UI_query = function(input,output,nano,drugs,chemical,disease){

  output$nano_cluster_input = renderUI({
    selectizeInput('nano_cluster_input', label = "Nanomaterials", choices = c("ALL",nano), multiple = TRUE,
                   options = list(create = TRUE),selected= "ALL")
  })
  
  output$drug_cluster_input = renderUI({
    selectizeInput('drug_cluster_input', label = "Drugs", choices = c(unlist(ATC_choice_list),drugs), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  output$disease_cluster_input = renderUI({
    selectizeInput('disease_cluster_input', label = "Diseases", choices = c("ALL",disease), multiple = TRUE,
                   options = list(create = TRUE))
  })
  
  output$chemical_cluster_input = renderUI({
    selectizeInput('chemical_cluster_input', label = "Chemicals", choices = c(unlist(chemical_choice_list),chemical), multiple = TRUE,
                   options = list(create = TRUE))
  })
}

conditional_query_nodes_cluster = function(input,output,DEBUGGING,nano,drugs,chemical,disease,chemMat,join10){
  selected_nodes = c()
  disease_list = list()
  
  xx = paste(input$nano_cluster_input,input$drug_cluster_input, input$chemical_cluster_input, input$disease_cluster_input,sep="")
  if(DEBUGGING){
    message("query_utilities::conditional_query_nodes. Concatenazione: ",xx,"\n")
    message("query_utilities::conditional_query_nodes. length(xx): ",length(xx),"\n")
  }
  if(length(xx)==0){
    output$info2_1 <- renderUI({
      HTML("Please insert at least one object for the query!")
    }) 
    validate(need(length(xx)>0, "Please insert at least one object for the query!"))
  }
  
  ALL_TRUE = ("ALL" %in% input$nano_cluster_input) & ("ALL" %in% input$drug_cluster_input) &   
  ("ALL" %in% input$chemical_cluster_input) & ("ALL" %in% input$disease_cluster_input)
  
  nano_query = parse_nano_query_input(input$nano_cluster_input,nano)
  
#   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {nano_query = ",nano_query, "}\n")
#   }
  
  drug_query = parse_drug_query_input(input$drug_cluster_input,drugs,ATC_letter_vector)
  #   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {drug_query before checking= ",drug_query, "}\n")
#   }
  
  
#   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {drug_query = ",drug_query, "}\n")
#   }
#   
  chemical_query = input$chemical_cluster_input
  if(length(chemical_query) != 0){
    if(chemical_query=="ALL"){
      chemical_query = chemical
    }
    else{
      if(chemical_query %in% names(table(chemMat[,2]))){
         index_C = which(chemMat[,2] %in% chemical_query)
         chemical_query = unique(chemMat[index_C,1])
      }
    }
  }
  
#   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {chemical_query = ",chemical_query, "}\n")
#   }
  
  disease_query = input$disease_cluster_input
  if(length(disease_query)!=0){
    if(disease_query=="ALL"){
      disease_query = disease
    }
  }
  
#   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {disease_query = ",disease_query, "}\n")
#   }
  query_nodes = c(nano_query,drug_query,disease_query,chemical_query)
  type_qn = c(rep("nano",length(nano_query)),
              rep("drugs",length(drug_query)),
              rep("disease",length(disease_query)),
              rep("chemical",length(chemical_query)))
  
  for(i in query_nodes){
    disease_list[[i]] = i
    selected_nodes = c(selected_nodes,i)
  }
  
#   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {query_nodes = ",query_nodes, "}\n")
#   }
  
  combination_vect=c(length(nano_query),length(drug_query),length(disease_query),length(chemical_query))
  combination_vect[combination_vect>0]=1
  names(combination_vect)=c("nano","drug","dis","chem")
  
#   if(DEBUGGING){
#     message("query_utilities::conditional_query_nodes {combination_vect = ",combination_vect, "}\n")
#   }
  
  return(list(query_nodes=query_nodes,ALL_TRUE=ALL_TRUE,
              nano_query = nano_query,
              drug_query = drug_query,
              chemical_query = chemical_query,
              disease_query = disease_query,
              combination_vect=combination_vect,
              disease_list=disease_list,selected_nodes=selected_nodes))
}

ADJ_matrix = function(W_ADJ,input,output,nano,drugs,chemical,disease,chemMat,join10,nNano = 29,nNodes = 3866){
  ADJ = matrix(0,length(cluster),length(cluster))
  rownames(ADJ) = colnames(ADJ) = V(graph_gw)$name
  
  for(i in 1:nNano){
    idx = which(cluster==i)
    ADJ[i,idx]=1
    ADJ[idx,i]=1
  }
  
  ADJ = ADJ[colnames(W_ADJ),colnames(W_ADJ)]
  ADJ2 = ADJ
  #ADJ2 = ADJ * as.matrix(W_ADJ)
  ADJ2[1:nNano,1:nNano] = W_ADJ[1:nNano,1:nNano]
  
  CQN = conditional_query_nodes_cluster(input,output,TRUE,nano,drugs,chemical,disease,chemMat,join10)
  selected_nodes = CQN$selected_nodes
  ALL_TRUE = CQN$ALL_TRUE
  #message("ADJ_matrix. selected_nodes: ",selected_nodes,"\n")
  
  if(ALL_TRUE){
    message("Grafo intero \n")
    g_clust = graph.adjacency(adjmatrix = ADJ2,mode = "undirected",weighted = TRUE)
    
    type = find_items_type(items = V(g_clust)$name,nano = nano,drugs = drugs,chemical = chemical,disease = disease)
    type = type$items_type
    V(g_clust)$type = type
    return(list(ADJ2=ADJ2,g_clust=g_clust))
  }else{
    clustering_idx = unique(cluster[selected_nodes])
    sel_n = names(cluster)[cluster %in% clustering_idx]
    sel_n = unique(c(selected_nodes,sel_n)) 
    ADJ2 = ADJ2[sel_n,sel_n]
    nano_selected = sel_n[sel_n %in% nano]
    ADJ2[nano_selected,nano_selected] = W_ADJ[nano_selected,nano_selected]
    
    #message("dim(ADJ2) ",dim(ADJ2),"\n")
    
    g_clust = graph.adjacency(adjmatrix = ADJ2,mode = "undirected",weighted = TRUE)
    idx_n = which(colnames(ADJ2) %in% nano)
    idx_dr = which(colnames(ADJ2) %in% drugs)
    idx_c = which(colnames(ADJ2) %in% chemical)
    idx_di = which(colnames(ADJ2) %in% disease)
    
    V(g_clust)$type = rep("nano",dim(ADJ2)[1])
    V(g_clust)$type[idx_dr] ="drug"
    V(g_clust)$type[idx_c] ="chem"
    V(g_clust)$type[idx_di] ="dise"
    
    #message("V(g_clust)$type  ",V(g_clust)$type ,"\n")
    return(list(ADJ2=ADJ2,g_clust=g_clust))
  }
}

from_igraph_to_data_frame_cluster= function(g_clust,ADJ2){
  
  
#   edges = get.data.frame(x = g_clust,what = "edge")
#   colnames(edges) = c("source","target","value")
#   vertices = data.frame(V(g_clust)$name,V(g_clust)$type)
#   colnames(vertices) = c("name","group")
#   vertices$size = igraph::degree(g_clust)
#   
#   edges$source = match(edges$source,vertices$name) - 1
#   edges$target = match(edges$target,vertices$name) - 1
#   
#   return(list(edges=edges,vertices=vertices))
  
  #edges = get.data.frame(x = g_clust,what = "edge")
  data_ = get.data.frame(x = g_clust,what = "both")
  edges = data_[[2]]
  colnames(edges) = c("source","target","value")
  #vertices = data.frame(V(g_clust)$name,V(g_clust)$type)
  vertices = data_[[1]]
  colnames(vertices) = c("name","group")
  #vertices$size = igraph::degree(g_clust)
  
  edges$source = match(edges$source,vertices$name) - 1
  edges$target = match(edges$target,vertices$name) - 1
  
  return(list(edges=edges,vertices=vertices))
}
