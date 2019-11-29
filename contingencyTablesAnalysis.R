

getStratContingencyArray <- function(df,sequence_col,disordered_index_list,phosphosite_regex="[ST]P"){
  
  disorderedRegions <- disordered_index_list
  TPphosphoSites <- gregexpr(phosphosite_regex,df[,sequence_col])
  
  stratContingencyArray <- array(dim = c(2,2,nrow(df)))
  for (i in 1:nrow(df)) {
    
    indexProt <- 1:nchar(df[i,sequence_col])
    indexDiso <- disorderedRegions[[i]]
    indexOrd <- setdiff(indexProt,indexDiso)
    indexPhospho <- as.numeric(TPphosphoSites[[i]])
    indexNonPhospho <- setdiff(indexProt,indexPhospho)
    
    ndp <- length(which(indexPhospho %in% indexDiso))
    nop <- length(which(indexPhospho %in% indexOrd))
    ndnp <- length(which(indexNonPhospho %in% indexDiso))
    nonp<- length(which(indexNonPhospho %in% indexOrd))
    stratContingencyArray[1,1,i]<-ndp
    stratContingencyArray[1,2,i]<-nop
    stratContingencyArray[2,1,i]<-ndnp
    stratContingencyArray[2,2,i]<-nonp
    
  }
return(stratContingencyArray)
}



