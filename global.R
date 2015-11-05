source("./my_igraph_function.R",local = FALSE)
load(paste(APP_PATH,"entities.RData",sep=""))
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


ATC_choice_list = list("All" = "ALL",
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

chemical_choice_list = list("All" = "ALL",
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