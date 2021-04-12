library(progress)

kinase_contingency_sim <- function(data,target_number,rep_number,stratified=F,ser_prop){
  
  pb <- progress_bar$new(
    format = "  processing [:bar] :percent eta: :eta",
    total = rep_number, clear = FALSE, width= 60)
  
  
  output_data <- data.frame(OR=numeric(length = rep_number),
                            pval=numeric(length = rep_number)) 
  OR <- numeric()
  pval <- numeric()
  
  #Select proteins that have at least one S or T
  # data <- subset(data,grepl("[S|T]",x = Sequence))
  # try only in known targetts for CDK
  data <- subset(data,target=="Cdk1 target")
  data_sim_subset <- data[,c("ACC#","Sequence","iupl_disoIndexes","vsl_disoIndexes","spot_disoIndexes","psite_CDK1","target")]
  
  
  # data_sim_subset <- data[,c("ACC#","Sequence","iupl_disoIndexes","vsl_disoIndexes","spot_disoIndexes")]
  
  Ser_index<-list()
  Thr_index<-list()
  
  for (i in 1:nrow(data_sim_subset)) {
    Ser_index[[i]] <- as.numeric(gregexpr("S", as.data.frame(data_sim_subset)[i,"Sequence"])[[1]])
    Thr_index[[i]] <- as.numeric(gregexpr("T", as.data.frame(data_sim_subset)[i,"Sequence"])[[1]])
  }
  
  data_sim_subset$ser_index <- Ser_index
  data_sim_subset$thr_index <- Thr_index
  
  ST_total_number <- sapply(Ser_index, length) + sapply(Thr_index, length)
  
  #Can I select son STPs???????
  SerP_index<-list()
  ThrP_index<-list()
  
  for (i in 1:nrow(data_sim_subset)) {
    SerP_index[[i]] <- as.numeric(gregexpr("SP", as.data.frame(data_sim_subset)[i,"Sequence"])[[1]])
    ThrP_index[[i]] <- as.numeric(gregexpr("TP", as.data.frame(data_sim_subset)[i,"Sequence"])[[1]])
  }
  
  data_sim_subset$serP_index <- SerP_index
  data_sim_subset$thrP_index <- ThrP_index
  
  STP_total_number <- sapply(SerP_index, length) + sapply(ThrP_index, length)
  
  
  for (rep in 1:rep_number) {
    
    pb$tick()
    
    sim_targets_index <- sample(nrow(data_sim_subset),target_number)
    # Maybe I don't have to use a dist function, just resampling the number of psites per target in known kinases 
    # psites_Ser_sampling_number <- rpois(target_number, 3)
    # psites_Thr_sampling_number <- rpois(target_number, 2)
    # psites_Total_sampling_number <- psites_Ser_sampling_number + psites_Thr_sampling_number
    
    # re-sampling the number of CDK targets
    cdk_psite_numbers <- sapply(subset(human_data,target=="Cdk1 target")$psite_CDK1,length)
    sim_psite_numbers <- sample(cdk_psite_numbers,size = target_number)
    psites_Ser_sampling_number <- rbinom(rep(1,target_number),sim_psite_numbers,prob = ser_prop)
    psites_Thr_sampling_number <- sim_psite_numbers-psites_Ser_sampling_number
    
    
    #Check that theres on S or T taken on each target y Check that for every protein the number of random psites generated never is bigger that the number of S and T in the protein
    # while (sum(sapply(data_sim_subset[sim_targets_index,"ser_index"], length) < psites_Ser_sampling_number)>0 | sum(sapply(data_sim_subset[sim_targets_index,"thr_index"], length) < psites_Thr_sampling_number)>0) {
    #   sim_psite_numbers <- sample(cdk_psite_numbers,size = target_number)
    #   psites_Ser_sampling_number <- rbinom(rep(1,target_number),sim_psite_numbers,prob = ser_prop)
    #   psites_Thr_sampling_number <- sim_psite_numbers-psites_Ser_sampling_number
    # }
    
    # Same that above but for when considering STP
    while (sum(sapply(data_sim_subset[sim_targets_index,"serP_index"], length) < psites_Ser_sampling_number)>0 | sum(sapply(data_sim_subset[sim_targets_index,"thrP_index"], length) < psites_Thr_sampling_number)>0) {
      sim_psite_numbers <- sample(cdk_psite_numbers,size = target_number)
      psites_Ser_sampling_number <- rbinom(rep(1,target_number),sim_psite_numbers,prob = ser_prop)
      psites_Thr_sampling_number <- sim_psite_numbers-psites_Ser_sampling_number
    }

    # # Check that for every protein the number of random psites generated never is bigger that the number of S and T in the protein. If that the case, reshuffle psites_Ser_sampling_number and psites_Thr_sampling_number until the condition is TRUE
    # while (sum(sapply(data_sim_subset[sim_targets_index,"ser_index"], length) < psites_Ser_sampling_number)>0) {
    #   psites_Ser_sampling_number <- sample(psites_Ser_sampling_number)
    # }
    # 
    # while (sum(sapply(data_sim_subset[sim_targets_index,"thr_index"], length) < psites_Thr_sampling_number)>0) {
    #   psites_Thr_sampling_number <- sample(psites_Thr_sampling_number)
    # }
    
    stratContingencyArray <- array(dim = c(2,2,target_number))
    for (i in 1:target_number) {
      target_idx <- sim_targets_index[i]
      target_info <- data_sim_subset[target_idx,]
      
      
      all_indexST <- c(target_info[["ser_index"]][[1]],target_info[["thr_index"]][[1]])
      # It's already calculated, but it could have been re calculated as:
      # all_indexST <- gregexpr("(S|T)",target_info[["Sequence"]])[[1]]
      indexDiso <- as.numeric(strsplit(target_info[["iupl_disoIndexes"]],",")[[1]])
      # indexPhosphoST <- c(sample(target_info[["ser_index"]][[1]],psites_Ser_sampling_number[i]),sample(target_info[["thr_index"]][[1]],psites_Thr_sampling_number[i]))
      indexPhosphoST <- c(sample(target_info[["serP_index"]][[1]],psites_Ser_sampling_number[i]),sample(target_info[["thrP_index"]][[1]],psites_Thr_sampling_number[i]))
      # indexPhosphoST <- target_info[["psite_CDK1"]][[1]]
      indexNonPhosphoST <- setdiff(as.numeric(all_indexST),indexPhosphoST)
      
      ndp <- length(which(indexPhosphoST %in% indexDiso))
      nop <- length(which(!indexPhosphoST %in% indexDiso))
      ndnp <- length(which(indexNonPhosphoST %in% indexDiso))
      nonp<- length(which(!indexNonPhosphoST %in% indexDiso))
      
      stratContingencyArray[1,1,i]<-ndp
      stratContingencyArray[1,2,i]<-ndnp
      stratContingencyArray[2,1,i]<-nop
      stratContingencyArray[2,2,i]<-nonp
    }
    
    ContingencyArray <- matrix(c(sum(stratContingencyArray[1,1,]),
                                 sum(stratContingencyArray[1,2,]),
                                 sum(stratContingencyArray[2,1,]),
                                 sum(stratContingencyArray[2,2,])),
                               nrow = 2,byrow = T)
    
    if (stratified) {
      current_test <- mantelhaen.test(stratContingencyArray)
    } else {
      current_test <- fisher.test(ContingencyArray)
    }
    
    OR[rep] <- current_test$estimate
    pval[rep] <- current_test$p.value
      
    

  }

  output_data$OR <- OR
  output_data$pval <- pval
  output_data$adj.pval <- p.adjust(pval,"BH")
  return(output_data)
}





# 
# # Psites not located in S or T
# test<-apply(subset(human_data,target_mapk=="mapk target"), 1, function(x){
#   splittedSeq <- strsplit(x[["Sequence"]],"")[[1]]
#   psites<-splittedSeq[x[["psite_MAPK"]]]
#   output<-x[["psite_MAPK"]][which(psites!="S" & psites!="T")]
#   names(output) <- paste(names(output),x[["ACC#"]][[1]],collapse = "_")
#   return(output)
# })
# names(test) <- paste(names(test),human_data[as.numeric(names(test)),"ACC#"],sep = "_")
# test[sapply(test, function(x){length(x)!=0})]
# 
# 
# #summary of the psites
# test<-apply(subset(human_data,target_mapk=="mapk target"), 1, function(x){
#   splittedSeq <- strsplit(x[["Sequence"]],"")[[1]]
#   psites<-splittedSeq[x[["psite_MAPK"]]]
#   return(psites)
# })
# 
# summary(as.factor(unlist(test)))


# STdiso<-apply(human_data,1,function(x){     
#     ST_index<-as.numeric(gregexpr("[S|T]", x[["Sequence"]])[[1]])
#     indexDiso <- as.numeric(strsplit(x[["iupl_disoIndexes"]],",")[[1]])
#     return(sum(ST_index %in% indexDiso))
#     })

# STdiso<-apply(human_data,1,function(x){
# ST_index<-as.numeric(gregexpr("[S|T]", x[["Sequence"]])[[1]])
# indexDiso <- as.numeric(strsplit(x[["iupl_disoIndexes"]],",")[[1]])
# return(sum(ST_index %in% indexDiso))
# })
# STstruct<-apply(human_data,1,function(x){
# ST_index<-as.numeric(gregexpr("[S|T]", x[["Sequence"]])[[1]])
# indexDiso <- as.numeric(strsplit(x[["iupl_disoIndexes"]],",")[[1]])
# return(sum(!ST_index %in% indexDiso))
# })
# sum(STdiso)
# sum(STstruct)
# disoRes<-apply(human_data,1,function(x){
# indexDiso <- as.numeric(strsplit(x[["iupl_disoIndexes"]],",")[[1]])
# return(length(indexDiso))
# })
# 
# structRes<-apply(human_data,1,function(x){
# totalLength<-nchar(x$Sequence)
# indexDiso <- as.numeric(strsplit(x[["iupl_disoIndexes"]],",")[[1]])
# return(totalLength-length(indexDiso))
# })
# sum(disoRes)/(sum(disoRes)+sum(structRes))
# sum(STdiso)/sum(disoRes)
# 1/0.1727681
# sum(STstruct)/sum(structRes)
# 1/0.1244504
