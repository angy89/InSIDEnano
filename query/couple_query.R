range01 <- function(x){(x-min(x))/(max(x)-min(x))}

couple_query3 = function(input,output,disease_list,selected_nodes,ADJ,ADJ01,ADJ01_RANK,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,items_list){
  output$infoFree <- renderUI({HTML(info_free_query_text)}) 
  
  withProgress(message = 'Progress...', min = 1,max = 10, {
    validate(need(input$disease_couple != "", "Please select a Disease"))
    validate(need(length(th_p)!=0,"Please select a th"))
    query_node = input$disease_couple
    
    output$nano_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    
    output$drug_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    output$chemical_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    output$disease_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    
    if(DEBUGGING){
      message("query_node:: ",query_node,"\n")
      message("input$th_slider/100:: ",th_p,"\n")
    }
    
    incProgress(1, detail = "Thresholding...")
    nElem = round(nrow(ADJ01_RANK)*th_p)
    
 
    ADJ_RW = ADJ01_RANK  
    ADJ_RW[ADJ_RW>nElem] = 0
    
#     ADJ_RW01 = ADJ_RW
#     ADJ_RW01[ADJ_RW01>0]=1
#     SADJ_RW01 = ADJ_RW01* t(ADJ_RW01)
    
    # moltiplicando la matrice ADJ_RW01 per la sua trasposta, ottengo l'intersezione del k-vicinato (che è una matrice simmetrica)
    # moltiplicando poi per ADJ risalgo ai pesi originali
   # W_ADJ = ADJ_RW01* t(ADJ_RW01)* sign(ADJ)
    
    incProgress(1, detail = "Removing edges under threshold...")
    
    incProgress(1, detail = "Creating graph...")
    
   # ADJ_RW = ADJ_RW * (ADJ_RW01* t(ADJ_RW01)) 
    neigh = ADJ_RW[which(ADJ_RW[,query_node]!=0),query_node]
    tab = cbind(names(neigh),sign(ADJ[neigh,query_node]),neigh,round(ADJ01[neigh,query_node],2)*sign(ADJ[neigh,query_node]))    
#tab = cbind(names(neigh),W_ADJ[names(neigh),query_node],neigh)
    colnames(tab)=c("name","Weight","Rank","UW")
    rownames(tab) = NULL
    tab = as.data.frame(tab)
    
    incProgress(1, detail = "Creating results...")
    
    nano_t = tab[tab$name %in% nano,]
    drug_t = tab[tab$name %in% drugs,]
    chemical_t = tab[tab$name %in% chemical,]
    disease_t = tab[tab$name %in% disease,]
    
    
    output$boxplot_statistics = renderPlot({
      num_nano = as.numeric(nano_t[,4])
      num_drug = as.numeric(drug_t[,4])
      num_chem = as.numeric(chemical_t[,4])
      num_dise = as.numeric(disease_t[,4])
      
      nn_p = num_nano[num_nano>=0]
      nn_n = num_nano[num_nano<0]
      
      ndr_p = num_drug[num_drug>=0]
      ndr_n = num_drug[num_drug<0]
      
      nc_p = num_chem[num_chem>=0]
      nc_n = num_chem[num_chem<0]
      
      ndi_p = num_dise[num_dise>=0]
      ndi_n = num_dise[num_dise<0]
      
      par(mfrow=c(2,1),oma = c(5,4,0,0) + 0.1,mar=c(1, 1, 1, 12) + 0.1,xpd=TRUE)
      boxplot(nn_p,ndr_p,nc_p,ndi_p,col="red",main="Connections Weight",names=c("Nano","Drugs","Chemical","Disease"),horizontal = FALSE)
      legend(x = "topright",inset=c(-0.2,0),legend = c("Positive","Negative"),fill=c("red","green")) 
      boxplot(nn_n,ndr_n,nc_n,ndi_n,col="green",names=c("Nano","Drugs","Chemical","Disease"),horizontal = FALSE)
      
    })
    
    output$barplot_statistics = renderPlot({
      num_nano = as.numeric(nano_t[,4])
      num_drug = as.numeric(drug_t[,4])
      num_chem = as.numeric(chemical_t[,4])
      num_dise = as.numeric(disease_t[,4])
      
      nn_p = num_nano[num_nano>=0] /length(nano)
      nn_n = num_nano[num_nano<0]/length(nano)
      
      ndr_p = num_drug[num_drug>=0]/length(drugs)
      ndr_n = num_drug[num_drug<0]/length(drugs)
      
      nc_p = num_chem[num_chem>=0]/length(chemical)
      nc_n = num_chem[num_chem<0]/length(chemical)
      
      ndi_p = num_dise[num_dise>=0]/length(disease)
      ndi_n = num_dise[num_dise<0]/length(disease)
      
      vect_p = c(length(nn_p),length(ndr_p),length(nc_p),length(ndi_p))
      vect_n = c(length(nn_n),length(ndr_n),length(nc_n),length(ndi_n))
      par(mfrow=c(2,1),xpd=TRUE,oma = c(5,4,0,0) + 0.1,mar=c(1, 2, 2, 1) + 0.1)
      names(vect_n)= c("Nanos","Drugs","Chemicals","Diseases")
      barplot(vect_p,col="red",main = "Number of positive connections")      
      barplot(vect_n,col="green",main="Number of negative connections")
    })
    
    
    if(nrow(nano_t)>0){
      
      nano_t = as.matrix(nano_t)
      
      for(i in 1:nrow(nano_t)){
        if(as.numeric(nano_t[i,2])>=0){
          nano_t[i,1] = paste('<font color="red"><b>',as.character(nano_t[i,1]),'</b></font>')
          nano_t[i,2] = nano_t[i,2]
        }else{
          nano_t[i,1] = paste('<font color="green"><b>',as.character(nano_t[i,1]),'</b></font>')
          nano_t[i,2] = nano_t[i,2]      
          
        }
      }

      output$nano_couple_table = renderDataTable({
        validate(need(nrow(nano_t)>0,"No connections with nano"))
        
        nano_t = nano_t[,c(1,3)]
        cat("dim(nano_t)--",dim(nano_t),"\n")
        nano_t = nano_t[order(as.numeric(as.vector(nano_t[,2]))),]
        DT::datatable(nano_t,escape=FALSE,selection = "single")
      })
      
      
    }
    
    if(nrow(drug_t)>0){
      
      drug_t = as.matrix(drug_t)
      
      for(i in 1:nrow(drug_t)){
        if(as.numeric(drug_t[i,2])>=0){
          drug_t[i,1] = paste('<font color="red"><b>',as.character(drug_t[i,1]),'</b></font>')
          drug_t[i,2] = drug_t[i,2]
        }else{
          drug_t[i,1] = paste('<font color="green"><b>',as.character(drug_t[i,1]),'</b></font>')
          drug_t[i,2] = drug_t[i,2]                
        }
      }
      

      output$drug_couple_table = renderDataTable({
        validate(need(nrow(drug_t)>0,"No connections with drugs"))
        drug_t = drug_t[,c(1,3)]
        drug_t = drug_t[order(as.numeric(as.vector(drug_t[,2]))),]
        
        DT::datatable(drug_t,escape=FALSE,selection = "single") 
      })
    }
    
    if(nrow(chemical_t)>0){
      chemical_t = as.matrix(chemical_t)

      for(i in 1:nrow(chemical_t)){
        if(as.numeric(chemical_t[i,2])>=0){
          chemical_t[i,1] = paste('<font color="red"><b>',as.character(chemical_t[i,1]),'</b></font>')
          chemical_t[i,2] = chemical_t[i,2]
        }else{
          chemical_t[i,1] = paste('<font color="green"><b>',as.character(chemical_t[i,1]),'</b></font>')
          chemical_t[i,2] = chemical_t[i,2]      
          
        }
      }
      
      output$chemical_couple_table = renderDataTable({
        validate(need(nrow(chemical_t)>0,"No connections with drugs"))
        chemical_t = chemical_t[,c(1,3)]
        chemical_t = chemical_t[order(as.numeric(as.vector(chemical_t[,2]))),]
        
        DT::datatable(chemical_t,escape=FALSE,selection = "single") 
      })
    }
    
    if(nrow(disease_t)>0){
      disease_t = as.matrix(disease_t)
      
      for(i in 1:nrow(disease_t)){
        if(as.numeric(disease_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(disease_t[i,1]),'</b></font>'),disease_t[i,2]))
          disease_t[i,1] = paste('<font color="red"><b>',as.character(disease_t[i,1]),'</b></font>')
          disease_t[i,2] = disease_t[i,2]
        }else{
          disease_t[i,1] = paste('<font color="green"><b>',as.character(disease_t[i,1]),'</b></font>')
          disease_t[i,2] = disease_t[i,2]                
        }
      }
      
      output$disease_couple_table = renderDataTable({
        validate(need(nrow(disease_t)>0,"No connections with drugs"))
        disease_t = disease_t[,c(1,3)]
        disease_t = disease_t[order(as.numeric(as.vector(disease_t[,2]))),]
        
        DT::datatable(disease_t,escape=FALSE,selection = "single") 
      })
    }    
  }) 
}  



couple_query2 = function(input,output,disease_list,selected_nodes,ADJ,ADJ01,ADJ01_RANK,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,items_list){
  output$infoFree <- renderUI({HTML(info_free_query_text)}) 
  
  withProgress(message = 'Progress...', min = 1,max = 10, {
    validate(need(input$disease_couple != "", "Please select a Disease"))
    validate(need(length(th_p)!=0,"Please select a th"))
    query_node = input$disease_couple
    
    output$nano_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    
    output$drug_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    output$chemical_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    output$disease_couple_table = renderDataTable({
      mat = matrix("",1,2)
      colnames(mat) = c("name","rank")
      
      DT::datatable(mat,escape=FALSE,selection = "single")
    })
    
    if(DEBUGGING){
      message("query_node:: ",query_node,"\n")
      message("input$th_slider/100:: ",th_p,"\n")
    }
    
    incProgress(1, detail = "Thresholding...")
    
    #cat("in couple_query:: W_ADJ=ADJ: ",length(which(W_ADJ==0)),"\n")
    
    nElem = round(nrow(ADJ01_RANK)*th_p)
    
    ADJ_RW = ADJ01_RANK  
    ADJ_RW[ADJ_RW>nElem] = 0
    
    ADJ_RW01 = ADJ_RW
    ADJ_RW01[ADJ_RW01>0]=1
    
    SADJ_RW01 = ADJ_RW01* t(ADJ_RW01)
    # moltiplicando la matrice ADJ_RW01 per la sua trasposta, ottengo l'intersezione del k-vicinato (che è una matrice simmetrica)
    # moltiplicando poi per ADJ risalgo ai pesi originali
    W_ADJ = ADJ_RW01* t(ADJ_RW01)* sign(ADJ)
    #W_ADJ = ADJ_RW01* t(ADJ_RW01)* ADJ01
   
    incProgress(1, detail = "Removing edges under threshold...")
   
    incProgress(1, detail = "Creating graph...")
    
    ADJ_RW = ADJ_RW * (ADJ_RW01* t(ADJ_RW01)) 
    neigh = ADJ_RW[which(ADJ_RW[,query_node]!=0),query_node]
    tab = cbind(names(neigh),W_ADJ[names(neigh),query_node],neigh)
    colnames(tab)=c("name","Weight","Rank")
    rownames(tab) = NULL
    tab = as.data.frame(tab)
    
    incProgress(1, detail = "Creating results...")
    
    nano_t = tab[tab$name %in% nano,]
    drug_t = tab[tab$name %in% drugs,]
    chemical_t = tab[tab$name %in% chemical,]
    disease_t = tab[tab$name %in% disease,]
  
    
    output$boxplot_statistics = renderPlot({
      num_nano = as.numeric(nano_t[,2])
      num_drug = as.numeric(drug_t[,2])
      num_chem = as.numeric(chemical_t[,2])
      num_dise = as.numeric(disease_t[,2])
      
      nn_p = num_nano[num_nano>=0]
      nn_n = num_nano[num_nano<0]
      
      ndr_p = num_drug[num_drug>=0]
      ndr_n = num_drug[num_drug<0]
      
      nc_p = num_chem[num_chem>=0]
      nc_n = num_chem[num_chem<0]
      
      ndi_p = num_dise[num_dise>=0]
      ndi_n = num_dise[num_dise<0]
      
      par(mfrow=c(2,1),oma = c(5,4,0,0) + 0.1,mar=c(1, 1, 1, 12) + 0.1,xpd=TRUE)
      boxplot(nn_p,ndr_p,nc_p,ndi_p,col="red",main="Connections Weight",names=c("Nano","Drugs","Chemical","Disease"),horizontal = FALSE)
      legend(x = "topright",inset=c(-0.2,0),legend = c("Positive","Negative"),fill=c("red","green")) 
      boxplot(nn_n,ndr_n,nc_n,ndi_n,col="green",names=c("Nano","Drugs","Chemical","Disease"),horizontal = FALSE)
      
    })
    
    output$barplot_statistics = renderPlot({
      num_nano = as.numeric(nano_t[,2])
      num_drug = as.numeric(drug_t[,2])
      num_chem = as.numeric(chemical_t[,2])
      num_dise = as.numeric(disease_t[,2])
      
      nn_p = num_nano[num_nano>=0] /length(nano)
      nn_n = num_nano[num_nano<0]/length(nano)
      
      ndr_p = num_drug[num_drug>=0]/length(drugs)
      ndr_n = num_drug[num_drug<0]/length(drugs)
      
      nc_p = num_chem[num_chem>=0]/length(chemical)
      nc_n = num_chem[num_chem<0]/length(chemical)
      
      ndi_p = num_dise[num_dise>=0]/length(disease)
      ndi_n = num_dise[num_dise<0]/length(disease)
      
      vect_p = c(length(nn_p),length(ndr_p),length(nc_p),length(ndi_p))
      vect_n = c(length(nn_n),length(ndr_n),length(nc_n),length(ndi_n))
      par(mfrow=c(2,1),xpd=TRUE,oma = c(5,4,0,0) + 0.1,mar=c(1, 2, 2, 1) + 0.1)
      names(vect_n)= c("Nanos","Drugs","Chemicals","Diseases")
      barplot(vect_p,col="red",main = "Number of positive connections")      
      barplot(vect_n,col="green",main="Number of negative connections")
      #       legend(x = "topright",inset=c(-0.2,0),legend = c("Positive","Negative"),fill=c("red","green")) 
      
      #       nn_p = range01(num_nano[num_nano>=0])
      #       nn_n = range01(abs(num_nano[num_nano<0]))*-1
      #       
      #       ndr_p = range01(num_drug[num_drug>=0])
      #       ndr_n = range01(abs(num_drug[num_drug<0])) * -1
      #       
      #       nc_p = range01(num_chem[num_chem>=0])
      #       nc_n = range01(abs(num_chem[num_chem<0])) * -1
      #       
      #       ndi_p = range01(num_dise[num_dise>=0])
      #       ndi_n = range01(abs(num_dise[num_dise<0])) * -1
      
      #       boxplot(nn_p,nn_n,ndr_p,ndr_n,nc_p,nc_n,ndi_p,ndi_n,col=rep(c("red","green"),4),
      #               names=c(rep("Nano",2),rep("Drugs",2),rep("Chemicals",2),rep("Diseases",2)))
      #       legend(x = "topright",legend = c("Positive Connections","Negative Connections"),fill=c("red","green"))      
    })
    
    
    if(nrow(nano_t)>0){
      #       num_nano = as.numeric(nano_t[,2])
      #       
      #       num_nano[num_nano>=0] = num_nano[num_nano>=0]
      #       num_nano[num_nano<0] = num_nano[num_nano<0]
      #       
      #       nano_t[,2] = round(num_nano,2)
      #       nano_t = nano_t[order(nano_t[,2],decreasing = TRUE),]
      #       
      #       datTab = matrix(0,nrow = 1,ncol = 2)
      
      nano_t = as.matrix(nano_t)
      
      for(i in 1:nrow(nano_t)){
        if(as.numeric(nano_t[i,2])>=0){
          #  datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(nano_t[i,1]),'</b></font>'),nano_t[i,2]))
          nano_t[i,1] = paste('<font color="red"><b>',as.character(nano_t[i,1]),'</b></font>')
          nano_t[i,2] = nano_t[i,2]
        }else{
          nano_t[i,1] = paste('<font color="green"><b>',as.character(nano_t[i,1]),'</b></font>')
          nano_t[i,2] = nano_t[i,2]      
          #          datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(nano_t[i,1]),'</b></font>'),nano_t[i,2]))
          
        }
      }
      #datTab = datTab[-1,]
      #colnames(datTab) = c("Name","Weight")
      output$nano_couple_table = renderDataTable({
        validate(need(nrow(nano_t)>0,"No connections with nano"))
        nano_t = nano_t[,c(1,3)]
       # nano_t = as.data.frame(nano_t)
        nano_t = nano_t[order(as.numeric(as.vector(nano_t[,2]))),]
        
        #nano_t[,2] = as.numeric(as.vector(nano_t[,2]))
        DT::datatable(nano_t,escape=FALSE,selection = "single")
      })
      
      
    }
    
    if(nrow(drug_t)>0){
      
      drug_t = as.matrix(drug_t)
      
      #       num_nano = as.numeric(drug_t[,2])
      #       
      #       num_nano[num_nano>=0] = num_nano[num_nano>=0]
      #       num_nano[num_nano<0] = num_nano[num_nano<0]
      #       
      #       drug_t[,2] = round(num_nano,2)
      #       drug_t = drug_t[order(as.numeric(drug_t[,2]),decreasing = TRUE),]
      # 
      #       #datTab = matrix(0,nrow = 1,ncol = 2)
      
      for(i in 1:nrow(drug_t)){
        if(as.numeric(drug_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(drug_t[i,1]),'</b></font>'),drug_t[i,2]))
          drug_t[i,1] = paste('<font color="red"><b>',as.character(drug_t[i,1]),'</b></font>')
          drug_t[i,2] = drug_t[i,2]
        }else{
          drug_t[i,1] = paste('<font color="green"><b>',as.character(drug_t[i,1]),'</b></font>')
          drug_t[i,2] = drug_t[i,2]      
          #datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(drug_t[i,1]),'</b></font>'),drug_t[i,2]))
          
        }
      }
      
      #       datTab = datTab[-1,]
      #       colnames(datTab) = c("Name","Weight")
      output$drug_couple_table = renderDataTable({
        validate(need(nrow(drug_t)>0,"No connections with drugs"))
        drug_t = drug_t[,c(1,3)]
        #drug_t = as.data.frame(drug_t)
        #drug_t[,2] = as.numeric(as.vector(drug_t[,2]))
        drug_t = drug_t[order(as.numeric(as.vector(drug_t[,2]))),]
        
        DT::datatable(drug_t,escape=FALSE,selection = "single") 
      })
    }
    
    if(nrow(chemical_t)>0){
      chemical_t = as.matrix(chemical_t)
      #       num_nano = as.numeric(chemical_t[,2])
      #       
      #       num_nano[num_nano>=0] = num_nano[num_nano>=0]
      #       num_nano[num_nano<0] = num_nano[num_nano<0]
      #       
      #       chemical_t[,2] = round(num_nano,2)
      #       chemical_t = chemical_t[order(chemical_t[,2],decreasing = TRUE),]
      #       
      #       datTab = matrix(0,nrow = 1,ncol = 2)
      
      for(i in 1:nrow(chemical_t)){
        if(as.numeric(chemical_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(chemical_t[i,1]),'</b></font>'),chemical_t[i,2]))
          chemical_t[i,1] = paste('<font color="red"><b>',as.character(chemical_t[i,1]),'</b></font>')
          chemical_t[i,2] = chemical_t[i,2]
        }else{
          chemical_t[i,1] = paste('<font color="green"><b>',as.character(chemical_t[i,1]),'</b></font>')
          chemical_t[i,2] = chemical_t[i,2]      
          #datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(chemical_t[i,1]),'</b></font>'),chemical_t[i,2]))
          
        }
      }
      
      #       datTab = datTab[-1,]
      #       colnames(datTab) = c("Name","Weight")
      output$chemical_couple_table = renderDataTable({
        validate(need(nrow(chemical_t)>0,"No connections with drugs"))
        chemical_t = chemical_t[,c(1,3)]
        #chemical_t = as.data.frame(chemical_t)
        #chemical_t[,2] = as.numeric(as.vector(chemical_t[,2]))
        chemical_t = chemical_t[order(as.numeric(as.vector(chemical_t[,2]))),]
        
        DT::datatable(chemical_t,escape=FALSE,selection = "single") 
      })
    }
    
    if(nrow(disease_t)>0){
      disease_t = as.matrix(disease_t)
      #       num_nano = as.numeric(disease_t[,2])
      #       
      #       num_nano[num_nano>=0] = num_nano[num_nano>=0]
      #       num_nano[num_nano<0] = num_nano[num_nano<0]
      #       
      #       disease_t[,2] = round(num_nano,2)
      #       disease_t = disease_t[order(disease_t[,2],decreasing = TRUE),]
      #       
      #       
      #       datTab = matrix(0,nrow = 1,ncol = 2)
      
      for(i in 1:nrow(disease_t)){
        if(as.numeric(disease_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(disease_t[i,1]),'</b></font>'),disease_t[i,2]))
          disease_t[i,1] = paste('<font color="red"><b>',as.character(disease_t[i,1]),'</b></font>')
          disease_t[i,2] = disease_t[i,2]
        }else{
          disease_t[i,1] = paste('<font color="green"><b>',as.character(disease_t[i,1]),'</b></font>')
          disease_t[i,2] = disease_t[i,2]      
          #datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(disease_t[i,1]),'</b></font>'),disease_t[i,2]))
          
        }
      }
      
      #       datTab = datTab[-1,]
      #       colnames(datTab) = c("Name","Weight")
      output$disease_couple_table = renderDataTable({
        validate(need(nrow(disease_t)>0,"No connections with drugs"))
        disease_t = disease_t[,c(1,3)]
        #disease_t = as.data.frame(disease_t)
        #disease_t[,2] = as.numeric(as.vector(disease_t[,2]))
        disease_t = disease_t[order(as.numeric(as.vector(disease_t[,2]))),]
        
        DT::datatable(disease_t,escape=FALSE,selection = "single") 
      })
    }
    
    #     output$nano_couple_boxplot = renderPlot({
    #       validate(need(nrow(nano_t)>0,"No connections with nano"))
    #       num = as.numeric(nano_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="Nanomaterials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="Nanomaterials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    #     
    #     output$drug_couple_boxplot = renderPlot({
    #       validate(need(nrow(drug_t)>0,"No connections with drug"))
    #       num = as.numeric(drug_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="drugmaterials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="drugmaterials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    #     
    #     output$chemical_couple_boxplot = renderPlot({
    #       validate(need(nrow(chemical_t)>0,"No connections with chemical"))
    #       num = as.numeric(chemical_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="chemicalmaterials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="chemicalmaterials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    #     
    #     output$disease_couple_boxplot = renderPlot({
    #       validate(need(nrow(disease_t)>0,"No connections with disease"))
    #       num = as.numeric(disease_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="diseasematerials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="diseasematerials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    
    
    
  }) 
  
}  
  


couple_query = function(input,output,disease_list,selected_nodes,W_ADJ,th_p = input$th_slider/100,node_type,chemMat,join10,g,g_geni2,items_list){
  output$infoFree <- renderUI({HTML(info_free_query_text)}) 
  
  withProgress(message = 'Progress...', min = 1,max = 10, {
    validate(need(input$disease_couple != "", "Please select a Disease"))
    
    query_node = input$disease_couple
    if(DEBUGGING){
      message("query_node:: ",query_node,"\n")
      message("input$th_slider/100:: ",th_p,"\n")
    }
    
    incProgress(1, detail = "Thresholding...")
   
    
    cat("in couple_query:: W_ADJ=ADJ: ",length(which(W_ADJ==0)),"\n")
#     par(mfrow=c(2,5))
#     hist(W_ADJ[nano,nano])
#     hist(W_ADJ[drugs,drugs])
#     hist(W_ADJ[chemical,chemical])
#     hist(W_ADJ[disease,disease])
#     hist(W_ADJ[nano,drugs])
#     hist(W_ADJ[nano,chemical])
#     hist(W_ADJ[nano,disease])
#     hist(W_ADJ[drugs,chemical])
#     hist(W_ADJ[drugs,disease])
#     hist(W_ADJ[disease,chemical])
    
    THS = find_thresholds(W_ADJ,th_p) #in query_utilities.R
    
    incProgress(1, detail = "Removing edges under threshold...")
    W_ADJ = apply_thresholds(W_ADJ,THS) #in query_utilities.R
    
#     cat("in couple_query:: W_ADJ=ADJ: ",length(which(W_ADJ==0)),"\n")
#     par(mfrow=c(2,5))
#     hist(W_ADJ[nano,nano])
#     hist(W_ADJ[drugs,drugs])
#     hist(W_ADJ[chemical,chemical])
#     hist(W_ADJ[disease,disease])
#     hist(W_ADJ[nano,drugs])
#     hist(W_ADJ[nano,chemical])
#     hist(W_ADJ[nano,disease])
#     hist(W_ADJ[drugs,chemical])
#     hist(W_ADJ[drugs,disease])
#     hist(W_ADJ[disease,chemical])
    
    incProgress(1, detail = "Creating graph...")
    graph_gw = creating_graph(W_ADJ,node_type) #in query_utilities.R
       
    incProgress(1, detail = "Creating results...")
    
    tab <- data.frame(cbind(
      name=V(graph_gw)[neighbors(graph_gw,query_node, mode='out')]$name,
      weight=E(graph_gw)[query_node %->% neighbors(graph_gw,query_node, mode='out')]$weight))

    
    tab[,2] = round(as.numeric(as.vector(tab$weight)),2)
    tab <- tab[order(abs(tab$weight), decreasing=TRUE),]
    
    nano_t = tab[tab$name %in% nano,]
    drug_t = tab[tab$name %in% drugs,]
    chemical_t = tab[tab$name %in% chemical,]
    disease_t = tab[tab$name %in% disease,]
    
    
    output$boxplot_statistics = renderPlot({
      num_nano = as.numeric(nano_t[,2])
      num_drug = as.numeric(drug_t[,2])
      num_chem = as.numeric(chemical_t[,2])
      num_dise = as.numeric(disease_t[,2])
      
      nn_p = num_nano[num_nano>=0]
      nn_n = num_nano[num_nano<0]
      
      ndr_p = num_drug[num_drug>=0]
      ndr_n = num_drug[num_drug<0]
      
      nc_p = num_chem[num_chem>=0]
      nc_n = num_chem[num_chem<0]
      
      ndi_p = num_dise[num_dise>=0]
      ndi_n = num_dise[num_dise<0]
      
      par(mfrow=c(2,1),oma = c(5,4,0,0) + 0.1,mar=c(1, 1, 1, 12) + 0.1,xpd=TRUE)
      boxplot(nn_p,ndr_p,nc_p,ndi_p,col="red",main="Connections Weight",names=c("Nano","Drugs","Chemical","Disease"),horizontal = FALSE)
      legend(x = "topright",inset=c(-0.2,0),legend = c("Positive","Negative"),fill=c("red","green")) 
      boxplot(nn_n,ndr_n,nc_n,ndi_n,col="green",names=c("Nano","Drugs","Chemical","Disease"),horizontal = FALSE)
      
#       legend(x = "topright",inset=c(-0.2,0),legend = c("Positive","Negative"),fill=c("red","green")) 
      
#       nn_p = range01(num_nano[num_nano>=0])
#       nn_n = range01(abs(num_nano[num_nano<0]))*-1
#       
#       ndr_p = range01(num_drug[num_drug>=0])
#       ndr_n = range01(abs(num_drug[num_drug<0])) * -1
#       
#       nc_p = range01(num_chem[num_chem>=0])
#       nc_n = range01(abs(num_chem[num_chem<0])) * -1
#       
#       ndi_p = range01(num_dise[num_dise>=0])
#       ndi_n = range01(abs(num_dise[num_dise<0])) * -1
          
#       boxplot(nn_p,nn_n,ndr_p,ndr_n,nc_p,nc_n,ndi_p,ndi_n,col=rep(c("red","green"),4),
#               names=c(rep("Nano",2),rep("Drugs",2),rep("Chemicals",2),rep("Diseases",2)))
#       legend(x = "topright",legend = c("Positive Connections","Negative Connections"),fill=c("red","green"))      
    })
    
    output$barplot_statistics = renderPlot({
      num_nano = as.numeric(nano_t[,2])
      num_drug = as.numeric(drug_t[,2])
      num_chem = as.numeric(chemical_t[,2])
      num_dise = as.numeric(disease_t[,2])
      
      nn_p = num_nano[num_nano>=0] /length(nano)
      nn_n = num_nano[num_nano<0]/length(nano)
      
      ndr_p = num_drug[num_drug>=0]/length(drugs)
      ndr_n = num_drug[num_drug<0]/length(drugs)
      
      nc_p = num_chem[num_chem>=0]/length(chemical)
      nc_n = num_chem[num_chem<0]/length(chemical)
      
      ndi_p = num_dise[num_dise>=0]/length(disease)
      ndi_n = num_dise[num_dise<0]/length(disease)
      
      vect_p = c(length(nn_p),length(ndr_p),length(nc_p),length(ndi_p))
      vect_n = c(length(nn_n),length(ndr_n),length(nc_n),length(ndi_n))
      par(mfrow=c(2,1),xpd=TRUE,oma = c(5,4,0,0) + 0.1,mar=c(1, 2, 2, 1) + 0.1)
      names(vect_n)= c("Nanos","Drugs","Chemicals","Diseases")
      barplot(vect_p,col="red",main = "Number of positive connections")      
      barplot(vect_n,col="green",main="Number of negative connections")
      #       legend(x = "topright",inset=c(-0.2,0),legend = c("Positive","Negative"),fill=c("red","green")) 
      
      #       nn_p = range01(num_nano[num_nano>=0])
      #       nn_n = range01(abs(num_nano[num_nano<0]))*-1
      #       
      #       ndr_p = range01(num_drug[num_drug>=0])
      #       ndr_n = range01(abs(num_drug[num_drug<0])) * -1
      #       
      #       nc_p = range01(num_chem[num_chem>=0])
      #       nc_n = range01(abs(num_chem[num_chem<0])) * -1
      #       
      #       ndi_p = range01(num_dise[num_dise>=0])
      #       ndi_n = range01(abs(num_dise[num_dise<0])) * -1
      
      #       boxplot(nn_p,nn_n,ndr_p,ndr_n,nc_p,nc_n,ndi_p,ndi_n,col=rep(c("red","green"),4),
      #               names=c(rep("Nano",2),rep("Drugs",2),rep("Chemicals",2),rep("Diseases",2)))
      #       legend(x = "topright",legend = c("Positive Connections","Negative Connections"),fill=c("red","green"))      
    })

    
    if(nrow(nano_t)>0){
#       num_nano = as.numeric(nano_t[,2])
#       
#       num_nano[num_nano>=0] = num_nano[num_nano>=0]
#       num_nano[num_nano<0] = num_nano[num_nano<0]
#       
#       nano_t[,2] = round(num_nano,2)
#       nano_t = nano_t[order(nano_t[,2],decreasing = TRUE),]
#       
#       datTab = matrix(0,nrow = 1,ncol = 2)
      
      nano_t = as.matrix(nano_t)
      
      for(i in 1:nrow(nano_t)){
        if(as.numeric(nano_t[i,2])>=0){
        #  datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(nano_t[i,1]),'</b></font>'),nano_t[i,2]))
           nano_t[i,1] = paste('<font color="red"><b>',as.character(nano_t[i,1]),'</b></font>')
           nano_t[i,2] = nano_t[i,2]
        }else{
           nano_t[i,1] = paste('<font color="green"><b>',as.character(nano_t[i,1]),'</b></font>')
           nano_t[i,2] = nano_t[i,2]      
#          datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(nano_t[i,1]),'</b></font>'),nano_t[i,2]))
          
        }
      }
      #datTab = datTab[-1,]
      #colnames(datTab) = c("Name","Weight")
      output$nano_couple_table = renderDataTable({
        validate(need(nrow(nano_t)>0,"No connections with nano"))
        DT::datatable(nano_t,escape=FALSE,selection = "single")
      })
      
      
    }
    
    if(nrow(drug_t)>0){
      
        drug_t = as.matrix(drug_t)
      
#       num_nano = as.numeric(drug_t[,2])
#       
#       num_nano[num_nano>=0] = num_nano[num_nano>=0]
#       num_nano[num_nano<0] = num_nano[num_nano<0]
#       
#       drug_t[,2] = round(num_nano,2)
#       drug_t = drug_t[order(as.numeric(drug_t[,2]),decreasing = TRUE),]
# 
#       #datTab = matrix(0,nrow = 1,ncol = 2)
      
      for(i in 1:nrow(drug_t)){
        if(as.numeric(drug_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(drug_t[i,1]),'</b></font>'),drug_t[i,2]))
                     drug_t[i,1] = paste('<font color="red"><b>',as.character(drug_t[i,1]),'</b></font>')
                     drug_t[i,2] = drug_t[i,2]
        }else{
                     drug_t[i,1] = paste('<font color="green"><b>',as.character(drug_t[i,1]),'</b></font>')
                     drug_t[i,2] = drug_t[i,2]      
          #datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(drug_t[i,1]),'</b></font>'),drug_t[i,2]))
          
        }
      }
      
#       datTab = datTab[-1,]
#       colnames(datTab) = c("Name","Weight")
      output$drug_couple_table = renderDataTable({
        validate(need(nrow(drug_t)>0,"No connections with drugs"))
        DT::datatable(drug_t,escape=FALSE,selection = "single") 
      })
    }
    
    if(nrow(chemical_t)>0){
      chemical_t = as.matrix(chemical_t)
#       num_nano = as.numeric(chemical_t[,2])
#       
#       num_nano[num_nano>=0] = num_nano[num_nano>=0]
#       num_nano[num_nano<0] = num_nano[num_nano<0]
#       
#       chemical_t[,2] = round(num_nano,2)
#       chemical_t = chemical_t[order(chemical_t[,2],decreasing = TRUE),]
#       
#       datTab = matrix(0,nrow = 1,ncol = 2)
      
      for(i in 1:nrow(chemical_t)){
        if(as.numeric(chemical_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(chemical_t[i,1]),'</b></font>'),chemical_t[i,2]))
                    chemical_t[i,1] = paste('<font color="red"><b>',as.character(chemical_t[i,1]),'</b></font>')
                    chemical_t[i,2] = chemical_t[i,2]
        }else{
                    chemical_t[i,1] = paste('<font color="green"><b>',as.character(chemical_t[i,1]),'</b></font>')
                    chemical_t[i,2] = chemical_t[i,2]      
          #datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(chemical_t[i,1]),'</b></font>'),chemical_t[i,2]))
          
        }
      }
      
#       datTab = datTab[-1,]
#       colnames(datTab) = c("Name","Weight")
      output$chemical_couple_table = renderDataTable({
        validate(need(nrow(chemical_t)>0,"No connections with drugs"))
        DT::datatable(chemical_t,escape=FALSE,selection = "single") 
      })
    }
    
    if(nrow(disease_t)>0){
      disease_t = as.matrix(disease_t)
#       num_nano = as.numeric(disease_t[,2])
#       
#       num_nano[num_nano>=0] = num_nano[num_nano>=0]
#       num_nano[num_nano<0] = num_nano[num_nano<0]
#       
#       disease_t[,2] = round(num_nano,2)
#       disease_t = disease_t[order(disease_t[,2],decreasing = TRUE),]
#       
#       
#       datTab = matrix(0,nrow = 1,ncol = 2)
      
      for(i in 1:nrow(disease_t)){
        if(as.numeric(disease_t[i,2])>=0){
          #datTab = rbind(datTab,c(paste('<font color="red"><b>',as.character(disease_t[i,1]),'</b></font>'),disease_t[i,2]))
                    disease_t[i,1] = paste('<font color="red"><b>',as.character(disease_t[i,1]),'</b></font>')
                    disease_t[i,2] = disease_t[i,2]
        }else{
                    disease_t[i,1] = paste('<font color="green"><b>',as.character(disease_t[i,1]),'</b></font>')
                    disease_t[i,2] = disease_t[i,2]      
          #datTab = rbind(datTab,c(paste('<font color="green"><b>',as.character(disease_t[i,1]),'</b></font>'),disease_t[i,2]))
          
        }
      }
      
#       datTab = datTab[-1,]
#       colnames(datTab) = c("Name","Weight")
      output$disease_couple_table = renderDataTable({
        validate(need(nrow(disease_t)>0,"No connections with drugs"))
        DT::datatable(disease_t,escape=FALSE,selection = "single") 
      })
    }
    
    #     output$nano_couple_boxplot = renderPlot({
    #       validate(need(nrow(nano_t)>0,"No connections with nano"))
    #       num = as.numeric(nano_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="Nanomaterials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="Nanomaterials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    #     
    #     output$drug_couple_boxplot = renderPlot({
    #       validate(need(nrow(drug_t)>0,"No connections with drug"))
    #       num = as.numeric(drug_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="drugmaterials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="drugmaterials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    #     
    #     output$chemical_couple_boxplot = renderPlot({
    #       validate(need(nrow(chemical_t)>0,"No connections with chemical"))
    #       num = as.numeric(chemical_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="chemicalmaterials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="chemicalmaterials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
    #     
    #     output$disease_couple_boxplot = renderPlot({
    #       validate(need(nrow(disease_t)>0,"No connections with disease"))
    #       num = as.numeric(disease_t[,2])
    #       idx_p = which(num>=0)
    #       idx_n = which(num<0)
    #       if(length(idx_n)==0){
    #         boxplot(num[idx_p],main="diseasematerials connections",col = "red",axes=FALSE)
    #         axis(2)
    #         legend(x = "topright",legend = "Positive",fill="red")
    #       }else{
    #         boxplot(num[idx_p],num[idx_n],main="diseasematerials connections",axes=FALSE,col = c("red","green"))
    #         axis(2)
    #         legend(x = "topright",legend = c("Positive","Negative"),fill=c("red","green"))
    #       }
    #     })
 
    

  })
    
}