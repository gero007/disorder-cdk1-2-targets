library(readr)
library(ggplot2)

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

write.table(kinase_homologs,file = "utrech/kinases_homologs/all_Kinases_homologs.tab",quote = F,sep = "\t",row.names = F)
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
