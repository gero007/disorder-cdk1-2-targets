library(readr)
library(ggplot2)
library(ggpubr)
library(dplyr)

aurk_homologs <- read_delim("utrech/kinases_homologs/aurk_humanHomologs_mapped.tab", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
cdk_homologs <- read_delim("utrech/kinases_homologs/cdk_humanHomologs_mapped.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
dyrk_homologs <- read_delim("utrech/kinases_homologs/dyrk_humanHomologs_mapped.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
mapk_homologs <- read_delim("utrech/kinases_homologs/mapk_humanHomologs_mapped.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
nek_homologs <- read_delim("utrech/kinases_homologs/nek_humanHomologs_mapped.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
plk_homologs <- read_delim("utrech/kinases_homologs/plk_humanHomologs_mapped.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

kinase_homologs <- rbind(aurk_homologs,cdk_homologs,dyrk_homologs,mapk_homologs,nek_homologs,plk_homologs)

kinase_homologs <- kinase_homologs %>% group_by(`ACC#`) %>% summarise_at(c("Xenopus","Kinase"),function(x){paste(x, collapse=",")})

mapped_xenopus <- lapply(kinase_homologs$Xenopus, function(x){ return(unique(strsplit(x,",")[[1]]))})
mapped_xenopus <- lapply(mapped_xenopus, function(x){ return(paste(x,collapse = ","))})
kinase_homologs$Xenopus <- unlist(mapped_xenopus)

# write.table(kinase_homologs,file = "utrech/kinases_homologs/all_Kinases_homologs.tab",quote = F,sep = "\t",row.names = F)
# kinase_homologs$Xenopus[grepl(kinase_homologs$Xenopus,pattern = ",")]
#Fixed the inputs the tables are cleaned of redundancy
# kinase_homologs <- read_delim("utrech/kinases_homologs/all_Kinases_homologs_cleaned.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

#Map plots
##Total
barPlotData <- data.frame(target_number=c(415,659,42,408,43,465),kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),subset=rep("Total",6)) 

##Maping in the proteome

xen_Proteome <- data.frame(target_number=c(398,627,41,381,43,448),kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),subset=rep("xenopus proteome",6))
barPlotData <- rbind(barPlotData,xen_Proteome)

##Maping in the phosphoproteome
phosphoproteome_ids <- xenopus_data$ID
aurk_inPhosphoproteome <- sum(aurk_homologs$Xenopus %in% phosphoproteome_ids)
cdk_inPhosphoproteome <- sum(cdk_homologs$Xenopus %in% phosphoproteome_ids)
dyrk_inPhosphoproteome <- sum(dyrk_homologs$Xenopus %in% phosphoproteome_ids)
mapk_inPhosphoproteome <- sum(mapk_homologs$Xenopus %in% phosphoproteome_ids)
nek_inPhosphoproteome <- sum(nek_homologs$Xenopus %in% phosphoproteome_ids)
plk_inPhosphoproteome <- sum(plk_homologs$Xenopus %in% phosphoproteome_ids)

xen_PhosphoProteome <- data.frame(target_number=c(aurk_inPhosphoproteome,cdk_inPhosphoproteome,dyrk_inPhosphoproteome,mapk_inPhosphoproteome,nek_inPhosphoproteome,plk_inPhosphoproteome ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("xenopus phosphoproteome",6))

barPlotData <- rbind(barPlotData,xen_PhosphoProteome)

##Maping in the dynamic phosphoproteome
dynamicPhosphoproteome_ids <- subset(xenopus_data,ANOVA=="Dynamic")$ID
aurk_inDynamicPhosphoproteome <- sum(aurk_homologs$Xenopus %in% dynamicPhosphoproteome_ids)
cdk_inDynamicPhosphoproteome <- sum(cdk_homologs$Xenopus %in% dynamicPhosphoproteome_ids)
dyrk_inDynamicPhosphoproteome <- sum(dyrk_homologs$Xenopus %in% dynamicPhosphoproteome_ids)
mapk_inDynamicPhosphoproteome <- sum(mapk_homologs$Xenopus %in% dynamicPhosphoproteome_ids)
nek_inDynamicPhosphoproteome <- sum(nek_homologs$Xenopus %in% dynamicPhosphoproteome_ids)
plk_inDynamicPhosphoproteome <- sum(plk_homologs$Xenopus %in% dynamicPhosphoproteome_ids)

xen_DynamicPhosphoProteome <- data.frame(target_number=c(aurk_inDynamicPhosphoproteome,cdk_inDynamicPhosphoproteome,dyrk_inDynamicPhosphoproteome,mapk_inDynamicPhosphoproteome,nek_inDynamicPhosphoproteome,plk_inDynamicPhosphoproteome ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("xenopus dynamic phosphoproteome",6))

barPlotData <- rbind(barPlotData,xen_DynamicPhosphoProteome)

##Maping in the cluster A
clusterA_ids <- subset(xenopus_data,cluster_A)$ID
aurk_inclusterA <- sum(aurk_homologs$Xenopus %in% clusterA_ids)
cdk_inclusterA <- sum(cdk_homologs$Xenopus %in% clusterA_ids)
dyrk_inclusterA <- sum(dyrk_homologs$Xenopus %in% clusterA_ids)
mapk_inclusterA <- sum(mapk_homologs$Xenopus %in% clusterA_ids)
nek_inclusterA <- sum(nek_homologs$Xenopus %in% clusterA_ids)
plk_inclusterA <- sum(plk_homologs$Xenopus %in% clusterA_ids)

xen_clusterA <- data.frame(target_number=c(aurk_inclusterA,cdk_inclusterA,dyrk_inclusterA,mapk_inclusterA,nek_inclusterA,plk_inclusterA ),
                                         kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                                         subset=rep("cluster A",6))

barPlotData <- rbind(barPlotData,xen_clusterA)

##Maping in the cluster B
clusterB_ids <- subset(xenopus_data,cluster_B)$ID
aurk_inclusterB <- sum(aurk_homologs$Xenopus %in% clusterB_ids)
cdk_inclusterB <- sum(cdk_homologs$Xenopus %in% clusterB_ids)
dyrk_inclusterB <- sum(dyrk_homologs$Xenopus %in% clusterB_ids)
mapk_inclusterB <- sum(mapk_homologs$Xenopus %in% clusterB_ids)
nek_inclusterB <- sum(nek_homologs$Xenopus %in% clusterB_ids)
plk_inclusterB <- sum(plk_homologs$Xenopus %in% clusterB_ids)

xen_clusterB <- data.frame(target_number=c(aurk_inclusterB,cdk_inclusterB,dyrk_inclusterB,mapk_inclusterB,nek_inclusterB,plk_inclusterB ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster B",6))

barPlotData <- rbind(barPlotData,xen_clusterB)

##Maping in the cluster C
clusterC_ids <- subset(xenopus_data,cluster_C)$ID
aurk_inclusterC <- sum(aurk_homologs$Xenopus %in% clusterC_ids)
cdk_inclusterC <- sum(cdk_homologs$Xenopus %in% clusterC_ids)
dyrk_inclusterC <- sum(dyrk_homologs$Xenopus %in% clusterC_ids)
mapk_inclusterC <- sum(mapk_homologs$Xenopus %in% clusterC_ids)
nek_inclusterC <- sum(nek_homologs$Xenopus %in% clusterC_ids)
plk_inclusterC <- sum(plk_homologs$Xenopus %in% clusterC_ids)

xen_clusterC <- data.frame(target_number=c(aurk_inclusterC,cdk_inclusterC,dyrk_inclusterC,mapk_inclusterC,nek_inclusterC,plk_inclusterC ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster C",6))

barPlotData <- rbind(barPlotData,xen_clusterC)

##Maping in the cluster D
clusterD_ids <- subset(xenopus_data,cluster_D)$ID
aurk_inclusterD <- sum(aurk_homologs$Xenopus %in% clusterD_ids)
cdk_inclusterD <- sum(cdk_homologs$Xenopus %in% clusterD_ids)
dyrk_inclusterD <- sum(dyrk_homologs$Xenopus %in% clusterD_ids)
mapk_inclusterD <- sum(mapk_homologs$Xenopus %in% clusterD_ids)
nek_inclusterD <- sum(nek_homologs$Xenopus %in% clusterD_ids)
plk_inclusterD <- sum(plk_homologs$Xenopus %in% clusterD_ids)

xen_clusterD <- data.frame(target_number=c(aurk_inclusterD,cdk_inclusterD,dyrk_inclusterD,mapk_inclusterD,nek_inclusterD,plk_inclusterD ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster D",6))

barPlotData <- rbind(barPlotData,xen_clusterD)

# Calculate the proporttions

barPlotData$perc_of_subsets <- c(subset(barPlotData,subset=="Total")$target_number / subset(barPlotData,subset=="Total")$target_number,
                               subset(barPlotData,subset=="xenopus proteome")$target_number / 13576,
                               subset(barPlotData,subset=="xenopus phosphoproteome")$target_number / 1844,
                               subset(barPlotData,subset=="xenopus dynamic phosphoproteome")$target_number / 1032, 
                               subset(barPlotData,subset=="cluster A")$target_number / 254,
                               subset(barPlotData,subset=="cluster B")$target_number / 207,
                               subset(barPlotData,subset=="cluster C")$target_number / 200,
                               subset(barPlotData,subset=="cluster D")$target_number / 150)
                               
barPlotData$perc_of_subsets <- barPlotData$perc_of_subsets*100

barPlotData$perc_of_targets <- c(subset(barPlotData,subset=="Total")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="xenopus proteome")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="xenopus phosphoproteome")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="xenopus dynamic phosphoproteome")$target_number / subset(barPlotData,subset=="Total")$target_number, 
                                 subset(barPlotData,subset=="cluster A")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="cluster B")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="cluster C")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="cluster D")$target_number / subset(barPlotData,subset=="Total")$target_number)

barPlotData$perc_of_targets <- barPlotData$perc_of_targets*100

#Double Norm
barPlotData$doubleNorm <- c(subset(barPlotData,subset=="Total")$target_number / subset(barPlotData,subset=="Total")$target_number,
                                 subset(barPlotData,subset=="xenopus proteome")$target_number / (subset(barPlotData,subset=="Total")$target_number*13576),
                                 subset(barPlotData,subset=="xenopus phosphoproteome")$target_number / (subset(barPlotData,subset=="Total")$target_number*1844),
                                 subset(barPlotData,subset=="xenopus dynamic phosphoproteome")$target_number / (subset(barPlotData,subset=="Total")$target_number*1032), 
                                 subset(barPlotData,subset=="cluster A")$target_number / (subset(barPlotData,subset=="Total")$target_number*254),
                                 subset(barPlotData,subset=="cluster B")$target_number / (subset(barPlotData,subset=="Total")$target_number*207),
                                 subset(barPlotData,subset=="cluster C")$target_number / (subset(barPlotData,subset=="Total")$target_number*200),
                                 subset(barPlotData,subset=="cluster D")$target_number / (subset(barPlotData,subset=="Total")$target_number*150))

barPlotData$doubleNorm <- barPlotData$doubleNorm*10000
#

# PLOTS


p1 <- ggplot(subset(barPlotData,subset=="xenopus proteome")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,100) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Proteome")
p2 <- ggplot(subset(barPlotData,subset=="xenopus phosphoproteome")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + ylim(0,100)+ theme_classic2() + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Phosphoproteome")
p3 <- ggplot(subset(barPlotData,subset=="xenopus dynamic phosphoproteome")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + ylim(0,100)+ theme_classic2() + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Dynamic")
p4 <- ggplot(subset(barPlotData,subset=="cluster A")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster A")
p5 <- ggplot(subset(barPlotData,subset=="cluster B")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster B")
p6 <- ggplot(subset(barPlotData,subset=="cluster C")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster C")
p7 <- ggplot(subset(barPlotData,subset=="cluster D")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster D")

p_clusters <- ggarrange(p4,p5,p6,p7,ncol = 2,nrow = 2,align = "hv")

p_all <- ggarrange(p1,p2,p3,p_clusters,ncol = 4,nrow = 1,common.legend = T)
annotate_figure(p_all,left = text_grob("% of targets", rot = 90,size = 16))





p1 <- ggplot(subset(barPlotData,subset=="xenopus proteome")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,25) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Proteome")
p2 <- ggplot(subset(barPlotData,subset=="xenopus phosphoproteome")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + ylim(0,25)+ theme_classic2() + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Phosphoproteome")
p3 <- ggplot(subset(barPlotData,subset=="xenopus dynamic phosphoproteome")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + ylim(0,25)+ theme_classic2() + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Dynamic")
p4 <- ggplot(subset(barPlotData,subset=="cluster A")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,50) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster A")
p5 <- ggplot(subset(barPlotData,subset=="cluster B")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,50) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster B")
p6 <- ggplot(subset(barPlotData,subset=="cluster C")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,50) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster C")
p7 <- ggplot(subset(barPlotData,subset=="cluster D")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,50) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster D")

p_clusters <- ggarrange(p4,p5,p6,p7,ncol = 2,nrow = 2,align = "hv")

p_all <- ggarrange(p1,p2,p3,p_clusters,ncol = 4,nrow = 1,common.legend = T)
annotate_figure(p_all,left = text_grob("% of subset", rot = 90,size = 16))

p_doubleNorm_1 <- ggplot(subset(barPlotData,subset=="xenopus proteome")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Proteome")
p_doubleNorm_2 <- ggplot(subset(barPlotData,subset=="xenopus phosphoproteome")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + ylim(0,7)+ theme_classic2() + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Phosphoproteome")
p_doubleNorm_3 <- ggplot(subset(barPlotData,subset=="xenopus dynamic phosphoproteome")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + ylim(0,7)+ theme_classic2() + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Dynamic")
p_doubleNorm_4 <- ggplot(subset(barPlotData,subset=="cluster A")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster A")
p_doubleNorm_5 <- ggplot(subset(barPlotData,subset=="cluster B")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster B")
p_doubleNorm_6 <- ggplot(subset(barPlotData,subset=="cluster C")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster C")
p_doubleNorm_7 <- ggplot(subset(barPlotData,subset=="cluster D")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity")+ theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster D")

p_doubleNorm_clusters <- ggarrange(p_doubleNorm_4,p_doubleNorm_5,p_doubleNorm_6,p_doubleNorm_7,ncol = 2,nrow = 2,align = "hv")

p_doubleNorm_all <- ggarrange(p_doubleNorm_1,p_doubleNorm_2,p_doubleNorm_3,p_doubleNorm_clusters,ncol = 4,nrow = 1,common.legend = T)
annotate_figure(p_doubleNorm_all,left = text_grob("targets in cluster \n per 100 targets per 100 elements in cluster ", rot = 90,size = 16))



# ____________________________________________________________________________________________________
#Extracts


cluster1_ids <- readLines("utrech/extract/Cluster1_IDs.txt")
cluster2_ids <- readLines("utrech/extract/Cluster2_IDs.txt")
cluster3_ids <- readLines("utrech/extract/Cluster3_IDs.txt")
cluster4_ids <- readLines("utrech/extract/Cluster4_IDs.txt")
cluster5_ids <- readLines("utrech/extract/Cluster5_IDs.txt")
cluster6_ids <- readLines("utrech/extract/Cluster6_IDs.txt")


#Map plots
##Total
extractBarPlotData <- data.frame(target_number=c(415,659,42,408,43,465),kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),subset=rep("Total",6))

##Maping in the cluster 1
aurk_cluster1 <- sum(aurk_homologs$Xenopus %in% cluster1_ids)
cdk_cluster1 <- sum(cdk_homologs$Xenopus %in% cluster1_ids)
dyrk_cluster1 <- sum(dyrk_homologs$Xenopus %in% cluster1_ids)
mapk_cluster1 <- sum(mapk_homologs$Xenopus %in% cluster1_ids)
nek_cluster1 <- sum(nek_homologs$Xenopus %in% cluster1_ids)
plk_cluster1 <- sum(plk_homologs$Xenopus %in% cluster1_ids)

xen_cluster1 <- data.frame(target_number=c(aurk_cluster1,cdk_cluster1,dyrk_cluster1,mapk_cluster1,nek_cluster1,plk_cluster1 ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster 1",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_cluster1)

##Maping in the cluster 2
aurk_cluster2 <- sum(aurk_homologs$Xenopus %in% cluster2_ids)
cdk_cluster2 <- sum(cdk_homologs$Xenopus %in% cluster2_ids)
dyrk_cluster2 <- sum(dyrk_homologs$Xenopus %in% cluster2_ids)
mapk_cluster2 <- sum(mapk_homologs$Xenopus %in% cluster2_ids)
nek_cluster2 <- sum(nek_homologs$Xenopus %in% cluster2_ids)
plk_cluster2 <- sum(plk_homologs$Xenopus %in% cluster2_ids)

xen_cluster2 <- data.frame(target_number=c(aurk_cluster2,cdk_cluster2,dyrk_cluster2,mapk_cluster2,nek_cluster2,plk_cluster2 ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster 2",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_cluster2)


##Maping in the cluster 3
aurk_cluster3 <- sum(aurk_homologs$Xenopus %in% cluster3_ids)
cdk_cluster3 <- sum(cdk_homologs$Xenopus %in% cluster3_ids)
dyrk_cluster3 <- sum(dyrk_homologs$Xenopus %in% cluster3_ids)
mapk_cluster3 <- sum(mapk_homologs$Xenopus %in% cluster3_ids)
nek_cluster3 <- sum(nek_homologs$Xenopus %in% cluster3_ids)
plk_cluster3 <- sum(plk_homologs$Xenopus %in% cluster3_ids)

xen_cluster3 <- data.frame(target_number=c(aurk_cluster3,cdk_cluster3,dyrk_cluster3,mapk_cluster3,nek_cluster3,plk_cluster3 ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster 3",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_cluster3)

##Maping in the cluster 4
aurk_cluster4 <- sum(aurk_homologs$Xenopus %in% cluster4_ids)
cdk_cluster4 <- sum(cdk_homologs$Xenopus %in% cluster4_ids)
dyrk_cluster4 <- sum(dyrk_homologs$Xenopus %in% cluster4_ids)
mapk_cluster4 <- sum(mapk_homologs$Xenopus %in% cluster4_ids)
nek_cluster4 <- sum(nek_homologs$Xenopus %in% cluster4_ids)
plk_cluster4 <- sum(plk_homologs$Xenopus %in% cluster4_ids)

xen_cluster4 <- data.frame(target_number=c(aurk_cluster4,cdk_cluster4,dyrk_cluster4,mapk_cluster4,nek_cluster4,plk_cluster4 ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster 4",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_cluster4)

##Maping in the cluster 5
aurk_cluster5 <- sum(aurk_homologs$Xenopus %in% cluster5_ids)
cdk_cluster5 <- sum(cdk_homologs$Xenopus %in% cluster5_ids)
dyrk_cluster5 <- sum(dyrk_homologs$Xenopus %in% cluster5_ids)
mapk_cluster5 <- sum(mapk_homologs$Xenopus %in% cluster5_ids)
nek_cluster5 <- sum(nek_homologs$Xenopus %in% cluster5_ids)
plk_cluster5 <- sum(plk_homologs$Xenopus %in% cluster5_ids)

xen_cluster5 <- data.frame(target_number=c(aurk_cluster5,cdk_cluster5,dyrk_cluster5,mapk_cluster5,nek_cluster5,plk_cluster5 ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster 5",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_cluster5)

##Maping in the cluster 6
aurk_cluster6 <- sum(aurk_homologs$Xenopus %in% cluster6_ids)
cdk_cluster6 <- sum(cdk_homologs$Xenopus %in% cluster6_ids)
dyrk_cluster6 <- sum(dyrk_homologs$Xenopus %in% cluster6_ids)
mapk_cluster6 <- sum(mapk_homologs$Xenopus %in% cluster6_ids)
nek_cluster6 <- sum(nek_homologs$Xenopus %in% cluster6_ids)
plk_cluster6 <- sum(plk_homologs$Xenopus %in% cluster6_ids)

xen_cluster6 <- data.frame(target_number=c(aurk_cluster6,cdk_cluster6,dyrk_cluster6,mapk_cluster6,nek_cluster6,plk_cluster6 ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("cluster 6",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_cluster6)

##Maping all mitosis
mitosis_ids <- unique(c(cluster1_ids,cluster2_ids,cluster3_ids))
aurk_mitosis <- sum(aurk_homologs$Xenopus %in% mitosis_ids)
cdk_mitosis <- sum(cdk_homologs$Xenopus %in% mitosis_ids)
dyrk_mitosis <- sum(dyrk_homologs$Xenopus %in% mitosis_ids)
mapk_mitosis <- sum(mapk_homologs$Xenopus %in% mitosis_ids)
nek_mitosis <- sum(nek_homologs$Xenopus %in% mitosis_ids)
plk_mitosis <- sum(plk_homologs$Xenopus %in% mitosis_ids)

xen_mitosis <- data.frame(target_number=c(aurk_mitosis,cdk_mitosis,dyrk_mitosis,mapk_mitosis,nek_mitosis,plk_mitosis ),
                          kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                          subset=rep("mitosis",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_mitosis)


##Maping all interphase
interphase_ids <- unique(c(cluster4_ids,cluster5_ids,cluster6_ids))
aurk_interphase <- sum(aurk_homologs$Xenopus %in% interphase_ids)
cdk_interphase <- sum(cdk_homologs$Xenopus %in% interphase_ids)
dyrk_interphase <- sum(dyrk_homologs$Xenopus %in% interphase_ids)
mapk_interphase <- sum(mapk_homologs$Xenopus %in% interphase_ids)
nek_interphase <- sum(nek_homologs$Xenopus %in% interphase_ids)
plk_interphase <- sum(plk_homologs$Xenopus %in% interphase_ids)

xen_interphase <- data.frame(target_number=c(aurk_interphase,cdk_interphase,dyrk_interphase,mapk_interphase,nek_interphase,plk_interphase ),
                           kinase=c("AURK","CDK","DYRK","MAPK","NEK","PLK"),
                           subset=rep("interphase",6))
extractBarPlotData <- rbind(extractBarPlotData,xen_interphase)



# Calculate the proporttions
extractBarPlotData$perc_of_subsets <- c(subset(extractBarPlotData,subset=="Total")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                 subset(extractBarPlotData,subset=="cluster 1")$target_number / length(cluster1_ids),
                                 subset(extractBarPlotData,subset=="cluster 2")$target_number / length(cluster2_ids),
                                 subset(extractBarPlotData,subset=="cluster 3")$target_number / length(cluster3_ids),
                                 subset(extractBarPlotData,subset=="cluster 4")$target_number / length(cluster4_ids),
                                 subset(extractBarPlotData,subset=="cluster 5")$target_number / length(cluster5_ids),
                                 subset(extractBarPlotData,subset=="cluster 6")$target_number / length(cluster6_ids),
                                 subset(extractBarPlotData,subset=="mitosis")$target_number / length(mitosis_ids),
                                 subset(extractBarPlotData,subset=="interphase")$target_number / length(interphase_ids))

extractBarPlotData$perc_of_subsets <- extractBarPlotData$perc_of_subsets*100

extractBarPlotData$perc_of_targets <- c(subset(extractBarPlotData,subset=="Total")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="cluster 1")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="cluster 2")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="cluster 3")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="cluster 4")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="cluster 5")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="cluster 6")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="mitosis")$target_number / subset(extractBarPlotData,subset=="Total")$target_number,
                                        subset(extractBarPlotData,subset=="interphase")$target_number / subset(extractBarPlotData,subset=="Total")$target_number)

extractBarPlotData$perc_of_targets <- extractBarPlotData$perc_of_targets*100

#Double Norm counts per 100 target per 100 in the cluster
extractBarPlotData$doubleNorm <- c((subset(extractBarPlotData,subset=="Total")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * subset(extractBarPlotData,subset=="Total")$target_number),
                                        (subset(extractBarPlotData,subset=="cluster 1")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * length(cluster1_ids)),
                                        (subset(extractBarPlotData,subset=="cluster 2")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * length(cluster2_ids)),
                                        (subset(extractBarPlotData,subset=="cluster 3")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * length(cluster3_ids)),
                                        (subset(extractBarPlotData,subset=="cluster 4")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * length(cluster4_ids)),
                                        (subset(extractBarPlotData,subset=="cluster 5")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * length(cluster5_ids)),
                                        (subset(extractBarPlotData,subset=="cluster 6")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number * length(cluster6_ids)),
                                        (subset(extractBarPlotData,subset=="mitosis")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number*length(mitosis_ids)),
                                        (subset(extractBarPlotData,subset=="interphase")$target_number*10000) / (subset(extractBarPlotData,subset=="Total")$target_number*length(interphase_ids)))
#
#Plots per clusters
pe_subset_1 <- ggplot(subset(extractBarPlotData,subset=="cluster 1")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 1")
pe_subset_2 <- ggplot(subset(extractBarPlotData,subset=="cluster 2")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 2")
pe_subset_3 <- ggplot(subset(extractBarPlotData,subset=="cluster 3")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 3")
pe_subset_4 <- ggplot(subset(extractBarPlotData,subset=="cluster 4")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 4")
pe_subset_5 <- ggplot(subset(extractBarPlotData,subset=="cluster 5")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 5")
pe_subset_6 <- ggplot(subset(extractBarPlotData,subset=="cluster 6")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 6")

p_subset_all <- ggarrange(pe_subset_1,pe_subset_2,pe_subset_3,pe_subset_4,pe_subset_5,pe_subset_6,ncol = 6,nrow = 1,common.legend = T)
annotate_figure(p_subset_all,left = text_grob("% of subset", rot = 90,size = 16))


pe_target_1 <- ggplot(subset(extractBarPlotData,subset=="cluster 1")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 1")
pe_target_2 <- ggplot(subset(extractBarPlotData,subset=="cluster 2")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 2")
pe_target_3 <- ggplot(subset(extractBarPlotData,subset=="cluster 3")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 3")
pe_target_4 <- ggplot(subset(extractBarPlotData,subset=="cluster 4")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 4")
pe_target_5 <- ggplot(subset(extractBarPlotData,subset=="cluster 5")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 5")
pe_target_6 <- ggplot(subset(extractBarPlotData,subset=="cluster 6")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 6")

p_target_all <- ggarrange(pe_target_1,pe_target_2,pe_target_3,pe_target_4,pe_target_5,pe_target_6,ncol = 6,nrow = 1,common.legend = T)
annotate_figure(p_target_all,left = text_grob("% of targets", rot = 90,size = 16))



pe_double_1 <- ggplot(subset(extractBarPlotData,subset=="cluster 1")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 1")
pe_double_2 <- ggplot(subset(extractBarPlotData,subset=="cluster 2")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 2")
pe_double_3 <- ggplot(subset(extractBarPlotData,subset=="cluster 3")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 3")
pe_double_4 <- ggplot(subset(extractBarPlotData,subset=="cluster 4")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 4")
pe_double_5 <- ggplot(subset(extractBarPlotData,subset=="cluster 5")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 5")
pe_double_6 <- ggplot(subset(extractBarPlotData,subset=="cluster 6")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,20) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Cluster 6")

p_double_all <- ggarrange(pe_double_1,pe_double_2,pe_double_3,pe_double_4,pe_double_5,pe_double_6,ncol = 6,nrow = 1,common.legend = T)
annotate_figure(p_double_all,left = text_grob("targets in cluster \n per 100 targets per 100 elements in cluster", rot = 90,size = 16))


# Plots Mitosis vs Interphase

pe_subset_mito <- ggplot(subset(extractBarPlotData,subset=="mitosis")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,30) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Mitosis")
pe_subset_inter <- ggplot(subset(extractBarPlotData,subset=="interphase")) + geom_bar(aes(x=kinase,y=perc_of_subsets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,30) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Interphase")

p_subset_all <- ggarrange(pe_subset_mito,pe_subset_inter,ncol = 2,nrow = 1,common.legend = T)
annotate_figure(p_subset_all,left = text_grob("% of subset", rot = 90,size = 16))


pe_target_mito <- ggplot(subset(extractBarPlotData,subset=="mitosis")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Mitosis")
pe_target_inter <- ggplot(subset(extractBarPlotData,subset=="interphase")) + geom_bar(aes(x=kinase,y=perc_of_targets,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,40) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Interphase")

p_target_all <- ggarrange(pe_target_mito,pe_target_inter,ncol = 2,nrow = 1,common.legend = T)
annotate_figure(p_target_all,left = text_grob("% of targets", rot = 90,size = 16))



pe_double_mito <- ggplot(subset(extractBarPlotData,subset=="mitosis")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Mitosis")
pe_double_inter <- ggplot(subset(extractBarPlotData,subset=="interphase")) + geom_bar(aes(x=kinase,y=doubleNorm,color=kinase,fill=kinase),stat = "identity") + theme_classic2() + ylim(0,7) + theme(legend.position="none",axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + ggtitle("Interphase")

p_double_all <- ggarrange(pe_double_mito,pe_double_inter,ncol = 2,nrow = 1,common.legend = T)
annotate_figure(p_double_all,left = text_grob("targets in cluster \n per 100 targets per 100 elements in cluster ", rot = 90,size = 16))
