library(dplyr)

getSpotPredDiso <- function(df,output,accession_col){
  predDisoPercAll <- c()
  predStretchDist <- c()
  disoPosList <- list()
  disoSeqList <- list()
  for (i in 1:nrow(df)) {
    mat<-spotDisorderList[[df[i,accession_col]]]
    #select only disordered regions
    mat <- matrix(mat[mat[,3]=="D",],ncol=4)
    # print(df$acc[[i]])
    # print(mat)
    if (nrow(mat)!=0) {
      predStretchDist <- c(predStretchDist,as.numeric(mat[,4]))
      #calculate the number of disoredered AAs
      nDiso<-sum(as.numeric(mat[,4]))
      # print(nDiso)
      proteinLength <- nchar(df$sequence[[i]])
      # print(proteinLength)
      predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
      # print((nDiso/proteinLength)*100)
      disoPos <-numeric()
      disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
      return(disoPos)})
      disoPos <- as.numeric(unlist(disoPos))
      disoPosList[[i]]<-disoPos
      disoSeq <- character()
      for (row in 1:nrow(mat)) {
        startIdx <- as.numeric(mat[row,1])
        endIdx <- as.numeric(mat[row,2])
        disoSeq <- paste(disoSeq,substr(yeastDiso[i,"sequence"],startIdx,endIdx),sep = "")
      }
      disoSeqList[[df[i,accession_col]]]<-disoSeq
      
    } else {
      nDiso<-0
      proteinLength <- nchar(df$sequence[i])
      predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
      
      disoPos <-numeric()
      disoPos <- 0
      disoPosList[[i]]<-disoPos
    }
  }
  if (output=="percentage") {
    return(predDisoPercAll)
  } else if (output=="lenghts") {
    return(predStretchDist)
  } else if (output=="indices") {
    return(disoPosList)
  } else if (output=="sequences") {
    return(disoSeqList)
  } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}
}



holtPsitesTable <- Holt_MS %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
holtCDK1PsitesTable <- Holt_MS_Cdk1_all_phosphosites %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
HoltAllPsites <- merge.data.frame(holtPsitesTable,holtCDK1PsitesTable,by = "Uniprot",all.x = T)

colnames(HoltAllPsites)[2:3] <- c("Psites","Cdk1Psites")
HoltAllPsites$Psites<-strsplit(HoltAllPsites$Psites,",")
HoltAllPsites$Cdk1Psites<-strsplit(HoltAllPsites$Cdk1Psites,",")


spotDisorderRegions <- getSpotPredDiso(Holt_data,accession_col = "Uniprot",output = "indices")

stratContingencyArray <- array(dim = c(2,2,nrow(HoltAllPsites)))
for (i in 1:nrow(HoltAllPsites)) {
    
  # indexProt <- 1:nchar(df[i,sequence_col])
  indexDiso <- spotDisorderRegions[[i]]
  # indexOrd <- setdiff(indexProt,indexDiso)
  indexPhospho <- as.numeric(HoltAllPsites$Cdk1Psites[[i]])
  indexNonPhospho <- setdiff(as.numeric(HoltAllPsites$Psites[[i]]),as.numeric(HoltAllPsites$Cdk1Psites[[i]]))
  
    ndp <- length(which(indexPhospho %in% indexDiso))
    nop <- length(which(!indexPhospho %in% indexDiso))
    ndnp <- length(which(indexNonPhospho %in% indexDiso))
    nonp<- length(which(!indexNonPhospho %in% indexDiso))
    stratContingencyArray[1,1,i]<-ndp
    stratContingencyArray[1,2,i]<-nop
    stratContingencyArray[2,1,i]<-ndnp
    stratContingencyArray[2,2,i]<-nonp
    }
