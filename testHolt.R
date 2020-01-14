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
holtCDK1PsitesConsTable <- Holt_MS_Cdk1_consensus_phosphosites %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
#The contingency table analysis it only make sense in the CDK1 targets if I am comparing T/S phosphorylated in al conditions vs lost T/S phosphorylation upon cdk1 inhibition. Non cdk1 targets don't have lost T/S phosphorylation which makes the first row of the contingency table empty

HoltAllPsites <- merge.data.frame(holtPsitesTable,holtCDK1PsitesTable,by = "Uniprot",all.x = T)
HoltAllPsites <- merge.data.frame(HoltAllPsites,holtCDK1PsitesConsTable,by = "Uniprot",all.x = T)
  
colnames(HoltAllPsites)[2:4] <- c("Psites","Cdk1Psites","Cdk1PsitesCons")
HoltAllPsites$Psites<-strsplit(HoltAllPsites$Psites,",")

HoltAllPsites$Cdk1Psites<-strsplit(HoltAllPsites$Cdk1Psites,",")
# HoltAllPsites$Cdk1Psites[is.na(HoltAllPsites$Cdk1Psites)]<-0

HoltAllPsites$Cdk1PsitesCons<-strsplit(HoltAllPsites$Cdk1PsitesCons,",")
# HoltAllPsites$Cdk1PsitesCons[is.na(HoltAllPsites$Cdk1PsitesCons)]<-0

spotDisorderRegions <- getSpotPredDiso(Holt_data,accession_col = "Uniprot",output = "percentage")


Psites<-numeric()
PsitesCdk1<-numeric()
PsitesCdk1Cons<-numeric()

for (i in 1:nrow(HoltAllPsites)) {
    indexProt <- 1:nchar(Holt_data[i,"sequence"])
    indexDiso <- spotDisorderRegions[[i]]
    
 
    indexPsites <- as.numeric(HoltAllPsites$Psites[[i]])
    print(indexPsites)
    if (!is.na(HoltAllPsites$Cdk1Psites[[i]])) {
      indexPsitesCdk1 <- as.numeric(HoltAllPsites$Cdk1Psites[[i]])
    } else {indexPsitesCdk1 <- as.numeric()}
    if (!is.na(HoltAllPsites$Cdk1PsitesCons[[i]])) {
      indexPsitesCdk1Cons <- as.numeric(HoltAllPsites$Cdk1PsitesCons[[i]])
    } else {indexPsitesCdk1Cons <- as.numeric()}
    print(indexPsitesCdk1)
    print(indexPsitesCdk1Cons)
    cDisoPsites <- length(which(indexPsites %in% indexDiso))/length(indexDiso)
    #when no phsopho is detected we set a 0 value so the length of the index phsopho should be considering the values greater than 0
    cTotalPsites <- length(which(indexPsites>0))/length(indexProt)
    Psites[i] <- (cTotalPsites-cTotalPsites)/cTotalPsites
    
    cDisoPsitesCdk1 <- length(which(indexPsitesCdk1 %in% indexDiso))/length(indexDiso)
    cTotalPsitesCdk1 <- length(which(indexPsitesCdk1>0))/length(indexProt)
    PsitesCdk1[i] <- (cTotalPsitesCdk1-cTotalPsitesCdk1)/cTotalPsitesCdk1
    
    cDisoPsitesCdk1Cons <- length(which(indexPsitesCdk1Cons %in% indexDiso))/length(indexDiso)
    cTotalPsitesCdk1Cons <- length(which(indexPsitesCdk1Cons>0))/length(indexProt)
    PsitesCdk1Cons[i] <- (cTotalPsitesCdk1Cons-cTotalPsitesCdk1Cons)/cTotalPsitesCdk1Cons
    }

# rbind(cbind(Psites,rep("Phosphosite",length(Psites))),cbind(PsitesCdk1,rep("Phosphosite CDK1",length(PsitesCdk1))),cbind(PsitesCdk1Cons,rep("Phosphosite CDK1 Consensus",length(PsitesCdk1Cons))))

# stratContingencyArray <- array(dim = c(2,2,nrow(HoltAllPsites)))
# for (i in 1:nrow(HoltAllPsites)) {
#     
#   # indexProt <- 1:nchar(df[i,sequence_col])
#   indexDiso <- spotDisorderRegions[[i]]
#   # indexOrd <- setdiff(indexProt,indexDiso)
#   indexPhospho <- as.numeric(HoltAllPsites$Cdk1Psites[[i]])
#   indexNonPhospho <- setdiff(as.numeric(HoltAllPsites$Psites[[i]]),as.numeric(HoltAllPsites$Cdk1Psites[[i]]))
#   
#     ndp <- length(which(indexPhospho %in% indexDiso))
#     nop <- length(which(!indexPhospho %in% indexDiso))
#     ndnp <- length(which(indexNonPhospho %in% indexDiso))
#     nonp<- length(which(!indexNonPhospho %in% indexDiso))
#     stratContingencyArray[1,1,i]<-ndp
#     stratContingencyArray[1,2,i]<-nop
#     stratContingencyArray[2,1,i]<-ndnp
#     stratContingencyArray[2,2,i]<-nonp
#     }


