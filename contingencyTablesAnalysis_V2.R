
# Calculate the contingency table for S and T diso vs struct and Phospho vs Non Phospho

getStratContingencyArray <- function(df,sequence_col,diso_index_col,psites_col){
  
  indexST <- gregexpr("(S|T)",df[,sequence_col])
  
  stratContingencyArray <- array(dim = c(2,2,nrow(df)))
  for (i in 1:nrow(df)) {
    
    
    indexDiso <- as.numeric(strsplit(df[i,diso_index_col],",")[[1]])
    indexPhosphoST <- df[i,psites_col][[1]]
    indexNonPhosphoST <- setdiff(as.numeric(indexST[[i]]),indexPhosphoST)
    
    ndp <- length(which(indexPhosphoST %in% indexDiso))
    nop <- length(which(!indexPhosphoST %in% indexDiso))
    ndnp <- length(which(indexNonPhosphoST %in% indexDiso))
    nonp<- length(which(!indexNonPhosphoST %in% indexDiso))
    # stratContingencyArray[1,1,i]<-ndp
    # stratContingencyArray[1,2,i]<-nop
    # stratContingencyArray[2,1,i]<-ndnp
    # stratContingencyArray[2,2,i]<-nonp
    stratContingencyArray[1,1,i]<-ndp
    stratContingencyArray[1,2,i]<-ndnp
    stratContingencyArray[2,1,i]<-nop
    stratContingencyArray[2,2,i]<-nonp
    # indexProt <- 1:nchar(df[i,sequence_col])
    # indexDiso <- disorderedRegions[[i]]
    # indexOrd <- setdiff(indexProt,indexDiso)
    # indexPhospho <- as.numeric(TPphosphoSites[[i]])
    # indexNonPhospho <- setdiff(indexProt,indexPhospho)
    # 
    # ndp <- length(which(indexPhospho %in% indexDiso))
    # nop <- length(which(indexPhospho %in% indexOrd))
    # ndnp <- length(which(indexNonPhospho %in% indexDiso))
    # nonp<- length(which(indexNonPhospho %in% indexOrd))
    # stratContingencyArray[1,1,i]<-ndp
    # stratContingencyArray[1,2,i]<-nop
    # stratContingencyArray[2,1,i]<-ndnp
    # stratContingencyArray[2,2,i]<-nonp
    
  }
return(stratContingencyArray)
}








# [,1] [,2]
# [1,]    4    0
# [2,]   63   60
# 
# , , 249
# 
# [,1] [,2]
# [1,]    4    0
# [2,]    8   26
# 
# , , 250
# 
# [,1] [,2]
# [1,]    3    0
# [2,]  130   92
