library(dplyr)

# THis fuction sholud be added to getDisorder.R . Maybe into a general function on SPOT disorder predictor
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


# ----------------------------------------Data Prep-------------------------------------------------------------
holtPsitesTable <- holt_MS %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
holtCDK1PsitesTable <- holt_MS_Cdk1_all_phosphosites %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
holtCDK1PsitesConsTable <- holt_MS_Cdk1_consensus_phosphosites %>% group_by(Uniprot) %>% summarise(sites=paste(substring(site,2,1000), collapse=","))
#The contingency table analysis it only make sense in the CDK1 targets if I am comparing T/S phosphorylated in al conditions vs lost T/S phosphorylation upon cdk1 inhibition. Non cdk1 targets don't have lost T/S phosphorylation which makes the first row of the contingency table empty

holtAllPsites <- merge.data.frame(holtPsitesTable,holtCDK1PsitesTable,by = "Uniprot",all.x = T)
holtAllPsites <- merge.data.frame(HoltAllPsites,holtCDK1PsitesConsTable,by = "Uniprot",all.x = T)
  
colnames(HoltAllPsites)[2:4] <- c("Psites","Cdk1Psites","Cdk1PsitesCons")
holtAllPsites$Psites<-strsplit(holtAllPsites$Psites,",")

holtAllPsites$Cdk1Psites<-strsplit(holtAllPsites$Cdk1Psites,",")
# HoltAllPsites$Cdk1Psites[is.na(HoltAllPsites$Cdk1Psites)]<-0

holtAllPsites$Cdk1PsitesCons<-strsplit(holtAllPsites$Cdk1PsitesCons,",")
# HoltAllPsites$Cdk1PsitesCons[is.na(HoltAllPsites$Cdk1PsitesCons)]<-0
# ----------------------------------------------------------------------------------------------------------------
# Data prep could be added to the first section of the Index.Rmd



Holt_indexes_data <- Holt_data[,c(1,3,2,13)]
Holt_indexes_data$SPOT<- getSpotPredDiso(Holt_data,accession_col = "Uniprot",output = "indices")
Holt_indexes_data$SPOTperc<- getSpotPredDiso(Holt_data,accession_col = "Uniprot",output = "percentage")
Holt_indexes_data$Lite<- getPredLiteDiso(Holt_data,accession_col = "Uniprot",output = "indices")
Holt_indexes_data$Liteperc<- getPredLiteDiso(Holt_data,accession_col = "Uniprot",output = "percentage")
Holt_indexes_data <- merge.data.frame(Holt_indexes_data,HoltAllPsites,by="Uniprot",all.x = T)
Holt_indexes_data$Psites<-lapply(Holt_indexes_data$Psites,as.numeric)
Holt_indexes_data$Cdk1Psites<-lapply(Holt_indexes_data$Cdk1Psites,as.numeric)
Holt_indexes_data$Cdk1PsitesCons<-lapply(Holt_indexes_data$Cdk1PsitesCons,as.numeric)
Holt_indexes_data$`[S|T]`<-lapply(gregexpr("S|T",Holt_indexes_data$sequence),as.numeric)
Holt_indexes_data$`[S|T]P`<-lapply(gregexpr("[ST]P",Holt_indexes_data$sequence),as.numeric)

Holt_indexes_data_CDK1 <- subset(Holt_indexes_data,target_all=="Cdk1 target")

stratContingencyArray <- array(dim = c(2,2,nrow(Holt_indexes_data_CDK1)))
for (i in 1:nrow(Holt_indexes_data_CDK1)) {

  indexProt <- 1:nchar(Holt_indexes_data_CDK1$sequence[[i]])
  indexDiso <- Holt_indexes_data_CDK1$SPOT[[i]]
  # indexOrd <- setdiff(indexProt,indexDiso)
  indexST_Phospho <- Holt_indexes_data_CDK1$Cdk1Psites[[i]]
  indexST_NON_Phospho <- setdiff(Holt_indexes_data_CDK1$`[S|T]`[[i]],Holt_indexes_data_CDK1$Cdk1Psites[[i]])
  #Check for all the phosphosites fall  into ST positions
  # if(!all.equal(intersect(Holt_indexes_data_CDK1$Cdk1Psites[[i]],Holt_indexes_data_CDK1$`[S|T]`[[i]]), Holt_indexes_data_CDK1$Cdk1Psites[[i]])){stop(paste('row number ',as.character(i),': One or more phosphosite annotated is not a S or T'))}
  #Check passed
  
  Diso_Phos <- length(which(indexST_Phospho %in% indexDiso))
  Ord_Phos <- length(which(!(indexST_Phospho %in% indexDiso)))
  Diso_NONPhos <- length(which(indexST_NON_Phospho %in% indexDiso))
  Ord_NONPhos<- length(which(!(indexST_NON_Phospho %in% indexDiso)))
  stratContingencyArray[1,1,i]<-Diso_Phos
  stratContingencyArray[1,2,i]<-Ord_Phos
  stratContingencyArray[2,1,i]<-Diso_NONPhos
  stratContingencyArray[2,2,i]<-Ord_NONPhos
  }


# Individual fisher tests




Holt_indexes_data_CDK1$FisherOR <- apply(stratContingencyArray, 3,function(x){
  ft<-fisher.test(x)
  return(ft$estimate)
})
Holt_indexes_data_CDK1$FisherPval <- apply(stratContingencyArray, 3,function(x){
  ft<-fisher.test(x)
  return(ft$p.value)
})
Holt_indexes_data_CDK1$FisherAdjPval <- p.adjust(Holt_indexes_data_CDK1$FisherPval,method = "BH")

Holt_indexes_data_CDK1$FisherlogPvalue <- -log10(Holt_indexes_data_CDK1$FisherPval)
Holt_indexes_data_CDK1$FisherlogAdjPvalue <- -log10(Holt_indexes_data_CDK1$FisherAdjPval)
Holt_indexes_data_CDK1$FisherlogOr <- log10(Holt_indexes_data_CDK1$FisherOR)

fisherSpotDF$logOr[is.infinite(fisherSpotDF$logOr) & fisherSpotDF$logOr > 0] <- 2
fisherSpotDF$logOr[is.infinite(fisherSpotDF$logOr) & fisherSpotDF$logOr < 0] <- -2


###Trying to find a way to move the dots out from the logPvalue=0 with random noise(Jitter) in order to asses the amount of points we have there.
set.seed(1234)
fisherSpotDF$logPvalue[fisherSpotDF$logPvalue==0]<- 0.1
fisherSpotDF$logPvalue[fisherSpotDF$logPvalue==0.1]<-abs(jitter(fisherSpotDF$logPvalue[fisherSpotDF$logPvalue==0.1],6))

ggplot(fisherSpotDF, aes(x = logOr, y = logPvalue)) +
  geom_point(aes(color = target)) +
  # scale_color_manual(values = c("", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  scale_y_continuous(trans='log10',limits = c(0.01,10)) +
  xlim(-2, 2) +
  scale_shape_manual()+
  # geom_point(data = subset(fisherDF, log10(ORs) < -2), aes(x = -2, y = -log10(pvalue)), alpha = 0.5,shape=60,size=4) +
  # geom_point(data = subset(fisherDF, log10(ORs) > 2), aes(x = 2, y = -log10(pvalue)), alpha = 0.5,shape=62,size=4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")+
  annotate(geom = "text",label="p-value = 0.05",x=-2,y=-log10(0.05),hjust=0,vjust=-0.8)


# Psites<-numeric()
# PsitesCdk1<-numeric()
# PsitesCdk1Cons<-numeric()
# 
# for (i in 1:nrow(HoltAllPsites)) {
#     indexProt <- 1:nchar(Holt_data[i,"sequence"])
#     indexDiso <- spotDisorderRegions[[i]]
#     
#  
#     indexPsites <- as.numeric(HoltAllPsites$Psites[[i]])
#     print(indexPsites)
#     if (!is.na(HoltAllPsites$Cdk1Psites[[i]])) {
#       indexPsitesCdk1 <- as.numeric(HoltAllPsites$Cdk1Psites[[i]])
#     } else {indexPsitesCdk1 <- as.numeric()}
#     if (!is.na(HoltAllPsites$Cdk1PsitesCons[[i]])) {
#       indexPsitesCdk1Cons <- as.numeric(HoltAllPsites$Cdk1PsitesCons[[i]])
#     } else {indexPsitesCdk1Cons <- as.numeric()}
#     print(indexPsitesCdk1)
#     print(indexPsitesCdk1Cons)
#     cDisoPsites <- length(which(indexPsites %in% indexDiso))/length(indexDiso)
#     #when no phsopho is detected we set a 0 value so the length of the index phsopho should be considering the values greater than 0
#     cTotalPsites <- length(which(indexPsites>0))/length(indexProt)
#     Psites[i] <- (cTotalPsites-cTotalPsites)/cTotalPsites
#     
#     cDisoPsitesCdk1 <- length(which(indexPsitesCdk1 %in% indexDiso))/length(indexDiso)
#     cTotalPsitesCdk1 <- length(which(indexPsitesCdk1>0))/length(indexProt)
#     PsitesCdk1[i] <- (cTotalPsitesCdk1-cTotalPsitesCdk1)/cTotalPsitesCdk1
#     
#     cDisoPsitesCdk1Cons <- length(which(indexPsitesCdk1Cons %in% indexDiso))/length(indexDiso)
#     cTotalPsitesCdk1Cons <- length(which(indexPsitesCdk1Cons>0))/length(indexProt)
#     PsitesCdk1Cons[i] <- (cTotalPsitesCdk1Cons-cTotalPsitesCdk1Cons)/cTotalPsitesCdk1Cons
#     }

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


