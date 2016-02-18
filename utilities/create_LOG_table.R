create_LOG_table = function(){
  LOG_CONDITIONAL = matrix("",1,9)
  colnames(LOG_CONDITIONAL) = c("Nano","Drugs","Disease","Chemical","Th","n1","n2","Type","File_name")
  LOG_CONDITIONAL = as.data.frame(LOG_CONDITIONAL)
  save(LOG_CONDITIONAL,file = "LOG.RData")
} 

create_LOG_table()
