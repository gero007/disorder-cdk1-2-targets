
library(forcats)
library(ggpubr)

psap_human <- read_delim("Guido/prediction_RF_human_first_0to100_gero.tab", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

psap_xenopus <- read_delim("Guido/prediction_RF_xenopus_first_0to100_gero.tab", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

#Xenopus
xenopus_psap_data <- cbind(psap_xenopus,subset=rep("Proteome",nrow(psap_xenopus)))

phosphoproteomic_xenopus_psap <- psap_xenopus[psap_xenopus$acc %in% xenopus_data$ID,]
xenopus_psap_data <- rbind(xenopus_psap_data,cbind(phosphoproteomic_xenopus_psap,subset=rep("Phospho-\nproteome",nrow(phosphoproteomic_xenopus_psap))))

dynamic_xenopus_psap <- psap_xenopus[psap_xenopus$acc %in% unique(xenopus_ANOVA_data$Proteins),]
xenopus_psap_data <- rbind(xenopus_psap_data,cbind(dynamic_xenopus_psap,subset=rep(" Dynamic",nrow(dynamic_xenopus_psap))))

clusterA_xenopus_psap <- psap_xenopus[psap_xenopus$acc %in% clusterA,]
xenopus_psap_data <- rbind(xenopus_psap_data,cbind(clusterA_xenopus_psap,subset=rep("Cluster A",nrow(clusterA_xenopus_psap))))

clusterB_xenopus_psap <- psap_xenopus[psap_xenopus$acc %in% clusterB,]
xenopus_psap_data <- rbind(xenopus_psap_data,cbind(clusterB_xenopus_psap,subset=rep("Cluster B",nrow(clusterB_xenopus_psap))))

clusterC_xenopus_psap <- psap_xenopus[psap_xenopus$acc %in% clusterC,]
xenopus_psap_data <- rbind(xenopus_psap_data,cbind(clusterC_xenopus_psap,subset=rep("Cluster C",nrow(clusterC_xenopus_psap))))

clusterD_xenopus_psap <- psap_xenopus[psap_xenopus$acc %in% clusterD,]
xenopus_psap_data <- rbind(xenopus_psap_data,cbind(clusterD_xenopus_psap,subset=rep("Cluster D",nrow(clusterD_xenopus_psap))))

xenopus_psap_data$subset <- forcats::as_factor(xenopus_psap_data$subset)

ggplot(xenopus_psap_data) + geom_violin(aes(x=subset,y=average,fill=subset),trim = F,adjust = .7) + 
  geom_boxplot(aes(x=subset,y=average),width = 0.15) +
  ggpubr::theme_classic2() +
  theme(text = element_text(size=20),legend.position = "none",legend.title = element_blank(),axis.ticks.x = element_blank(),panel.grid.major.y = element_line(colour = "grey",linetype = "dashed")) +
  scale_y_continuous(name = "PSAP score",limits = c(-0.20, 1.2), breaks = c(seq(0, 1, by = 0.2))) +
  scale_x_discrete() + xlab(NULL) +
  scale_fill_manual(values = c("#00adeeff","#ffdd15ff","#f6921eff","#f6921eff","#f6921eff","#f6921eff","#f6921eff"))


#Human
human_psap_data <- cbind(psap_human,subset=rep("Proteome",nrow(psap_human)))

phosphoproteomic_human_psap <- psap_human[psap_human$acc %in% human_data$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(phosphoproteomic_human_psap,subset=rep("Phospho-\nproteome",nrow(phosphoproteomic_human_psap))))

cdk_human_psap <- psap_human[psap_human$acc %in% subset(human_data,target=="Cdk1 target")$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(cdk_human_psap,subset=rep("CDK",nrow(cdk_human_psap))))

mapk_human_psap <- psap_human[psap_human$acc %in% subset(human_data,target_mapk=="mapk target")$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(mapk_human_psap,subset=rep("MAPK",nrow(mapk_human_psap))))

aurk_human_psap <- psap_human[psap_human$acc %in% subset(human_data,target_aurk=="aurk target")$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(aurk_human_psap,subset=rep("AURK",nrow(aurk_human_psap))))

plk_human_psap <- psap_human[psap_human$acc %in% subset(human_data,target_plk=="plk target")$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(plk_human_psap,subset=rep("PLK",nrow(plk_human_psap))))

nek_human_psap <- psap_human[psap_human$acc %in% subset(human_data,target_nek=="nek target")$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(nek_human_psap,subset=rep("NEK",nrow(nek_human_psap))))

dyrk_human_psap <- psap_human[psap_human$acc %in% subset(human_data,target_dyrk=="dyrk target")$`ACC#`,]
human_psap_data <- rbind(human_psap_data,cbind(dyrk_human_psap,subset=rep("DYRK",nrow(dyrk_human_psap))))

human_psap_data$subset <- forcats::as_factor(human_psap_data$subset)

ggplot(human_psap_data) + geom_violin(aes(x=subset,y=average,fill=subset),trim = F,adjust = .7) + 
  geom_boxplot(aes(x=subset,y=average),width = 0.15) +
  ggpubr::theme_classic2() +
  theme(text = element_text(size=20),legend.position = "none",legend.title = element_blank(),axis.ticks.x = element_blank(),panel.grid.major.y = element_line(colour = "grey",linetype = "dashed")) +
  scale_y_continuous(name = "PSAP score",limits = c(-0.20, 1.2), breaks = c(seq(0, 1, by = 0.2))) +
  scale_x_discrete() + xlab(NULL) +
  scale_fill_manual(values = c("#00adeeff","#ffdd15ff",pal_tron()(7)[c(3,4,1,2,5,6)]))
