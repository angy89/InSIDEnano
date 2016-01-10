LOG_CONDITIONAL = matrix("",nrow = 1,ncol = 9)
colnames(LOG_CONDITIONAL) = c("Nano","Drugs","Disease","Chemical","Th","n1","n2","type","File_name")
head(LOG_CONDITIONAL)
LOG_CONDITIONAL=rbind(LOG_CONDITIONAL,c("MWCNT","","Asthma","",100,2,2,"ALL","1.RData"))
LOG_CONDITIONAL=as.data.frame(LOG_CONDITIONAL)

save(LOG_CONDITIONAL,file="LOG.RData")

row_to_find <- data.frame(Nano="MWCNT",Drugs = "",Disease="Asthma",Chemical="",Th=80,n1=2,n2 = 2,type="ALL")
xx = match_df(LOG_CONDITIONAL[,1:8], row_to_find)
row_idx = as.numeric(rownames(xx))
file_name = as.character(LOG_CONDITIONAL[row_idx,9])
