# library(dplyr)
# 
# getSpotPredDiso <- function(df,output,accession_col){
#   predDisoPercAll <- c()
#   predStretchDist <- c()
#   disoPosList <- list()
#   disoSeqList <- list()
#   for (i in 1:nrow(df)) {
#     mat<-spotDisorderList[[df[i,accession_col]]]
#     #select only disordered regions
#     mat <- matrix(mat[mat[,3]=="D",],ncol=4)
#     # print(df$acc[[i]])
#     # print(mat)
#     if (nrow(mat)!=0) {
#       predStretchDist <- c(predStretchDist,as.numeric(mat[,4]))
#       #calculate the number of disoredered AAs
#       nDiso<-sum(as.numeric(mat[,4]))
#       # print(nDiso)
#       proteinLength <- nchar(df$sequence[[i]])
#       # print(proteinLength)
#       predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
#       # print((nDiso/proteinLength)*100)
#       disoPos <-numeric()
#       disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
#       return(disoPos)})
#       disoPos <- as.numeric(unlist(disoPos))
#       disoPosList[[i]]<-disoPos
#       disoSeq <- character()
#       for (row in 1:nrow(mat)) {
#         startIdx <- as.numeric(mat[row,1])
#         endIdx <- as.numeric(mat[row,2])
#         disoSeq <- paste(disoSeq,substr(yeastDiso[i,"sequence"],startIdx,endIdx),sep = "")
#       }
#       disoSeqList[[df[i,accession_col]]]<-disoSeq
#       
#     } else {
#       nDiso<-0
#       proteinLength <- nchar(df$sequence[i])
#       predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
#       
#       disoPos <-numeric()
#       disoPos <- 0
#       disoPosList[[i]]<-disoPos
#     }
#   }
#   if (output=="percentage") {
#     return(predDisoPercAll)
#   } else if (output=="lenghts") {
#     return(predStretchDist)
#   } else if (output=="indices") {
#     return(disoPosList)
#   } else if (output=="sequences") {
#     return(disoSeqList)
#   } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}
# }



# holtPsitesTable <- Holt_MS %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
# holtCDK1PsitesTable <- Holt_MS_Cdk1_all_phosphosites %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
# holtCDK1PsitesConsTable <- Holt_MS_Cdk1_consensus_phosphosites %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
# #The contingency table analysis it only make sense in the CDK1 targets if I am comparing T/S phosphorylated in al conditions vs lost T/S phosphorylation upon cdk1 inhibition. Non cdk1 targets don't have lost T/S phosphorylation which makes the first row of the contingency table empty
# 
# HoltAllPsites <- merge.data.frame(holtPsitesTable,holtCDK1PsitesTable,by = "Uniprot",all.x = T)
# HoltAllPsites <- merge.data.frame(HoltAllPsites,holtCDK1PsitesConsTable,by = "Uniprot",all.x = T)
#   
# colnames(HoltAllPsites)[2:4] <- c("Psites","Cdk1Psites","Cdk1PsitesCons")
# 
# spotDisorderRegions <- getSpotPredDiso(Holt_data,accession_col = "Uniprot",output = "indices")


#Due to the limmited numbers of psites per protein detected by the MS/MS method, stratified Cochram-Mantel-Haezel cannot be used. I'll summarize all the counts in only one contingency table.

STPsitesCT<-array(rep(0,4),dim = c(2,2))
PsitesCdk1CT<-array(rep(0,4),dim = c(2,2))
PsitesCdk1ConsCT<-array(rep(0,4),dim = c(2,2))

for (i in 1:nrow(HoltAllPsites)) {
    indexProt <- 1:nchar(Holt_data[i,"sequence"])
    indexDiso <- spotDisorderRegions[[i]]
    indexPsites <- as.numeric(strsplit(HoltAllPsites$Psites[[i]],",")[[1]])
    
    

    
    indexSTP <-gregexpr("[ST]P",Holt_data[,"sequence"])
    
    
    nSTPsitesDiso <- length(which(as.numeric(indexSTP[[i]]) %in% indexDiso))
    nSTPsitesOrd <- length(which(!(as.numeric(indexSTP[[i]]) %in% indexDiso)))
    nNonSTPsitesDiso <- length(which(setdiff(indexProt,as.numeric(indexSTP[[i]])) %in% indexDiso))
    nNonSTPsitesOrd<- length(which(!(setdiff(indexProt,as.numeric(indexSTP[[i]])) %in% indexDiso)))
    
    # print(nPsitesCdk1Diso)
    # print(nPsitesCdk1Ord)
    # print(nNonPsitesCdk1Diso)
    # print(nNonPsitesCdk1Ord)
    
    STPsitesCT[1,1]<-STPsitesCT[1,1]+nSTPsitesDiso
    STPsitesCT[1,2]<-STPsitesCT[1,2]+nSTPsitesOrd
    STPsitesCT[2,1]<-STPsitesCT[2,1]+nNonSTPsitesDiso
    STPsitesCT[2,2]<-STPsitesCT[2,2]+nNonSTPsitesOrd
    
    
    if (!is.na(HoltAllPsites$Cdk1Psites[[i]][1])) {
      indexPsitesCdk1 <- as.numeric(strsplit(HoltAllPsites$Cdk1Psites[[i]],",")[[1]])
      nPsitesCdk1Diso <- length(which(indexPsitesCdk1 %in% indexDiso))
      nPsitesCdk1Ord <- length(which(!(indexPsitesCdk1 %in% indexDiso)))
      nNonPsitesCdk1Diso <- length(which(setdiff(indexPsites,indexPsitesCdk1) %in% indexDiso))#recalculate
      nNonPsitesCdk1Ord<- length(which(!(setdiff(indexPsites,indexPsitesCdk1) %in% indexDiso)))#recalculate
      # print(nPsitesCdk1Diso+nPsitesCdk1Ord==length(indexPsitesCdk1))
      # print(nPsitesCdk1Diso+nPsitesCdk1Ord+nNonPsitesCdk1Diso+nNonPsitesCdk1Ord==length(indexPsites))
      # print("____")
      # print(indexPsites)
      # print(indexPsitesCdk1)
      # print(indexPsitesCdk1Cons)
      # print(nPsitesCdk1Diso)
      # print(nPsitesCdk1Ord)
      # print(nNonPsitesCdk1Diso)
      # print(nNonPsitesCdk1Ord)
      PsitesCdk1CT[1,1]<-PsitesCdk1CT[1,1]+nPsitesCdk1Diso
      PsitesCdk1CT[1,2]<-PsitesCdk1CT[1,2]+nPsitesCdk1Ord
      PsitesCdk1CT[2,1]<-PsitesCdk1CT[2,1]+nNonPsitesCdk1Diso
      PsitesCdk1CT[2,2]<-PsitesCdk1CT[2,2]+nNonPsitesCdk1Ord
    } 

    
    if (!is.na(HoltAllPsites$Cdk1PsitesCons[[i]][1])) {
      indexPsitesCdk1Cons <- as.numeric(strsplit(HoltAllPsites$Cdk1PsitesCons[[i]],",")[[1]])
      nPsitesCdk1ConsDiso <- length(which(indexPsitesCdk1Cons %in% indexDiso))
      nPsitesCdk1ConsOrd <- length(which(!(indexPsitesCdk1Cons %in% indexDiso)))
      nNonPsitesCdk1ConsDiso <- length(which(setdiff(indexPsites,indexPsitesCdk1Cons) %in% indexDiso))
      nNonPsitesCdk1ConsOrd<- length(which(!(setdiff(indexPsites,indexPsitesCdk1Cons) %in% indexDiso)))
      
      PsitesCdk1ConsCT[1,1]<-PsitesCdk1ConsCT[1,1]+nPsitesCdk1ConsDiso
      PsitesCdk1ConsCT[1,2]<-PsitesCdk1ConsCT[1,2]+nPsitesCdk1ConsOrd
      PsitesCdk1ConsCT[2,1]<-PsitesCdk1ConsCT[2,1]+nNonPsitesCdk1ConsDiso
      PsitesCdk1ConsCT[2,2]<-PsitesCdk1ConsCT[2,2]+nNonPsitesCdk1ConsOrd
    }
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


