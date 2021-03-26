library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggsci)
library(scales)
library(reshape2)

source("IUpredScoresPlotGenerator.R")


all_protein_id_phospho <- read_lines("utrech/ids/xenopus_allphosphosites_mpi.txt")
cluster_A <- read_lines("utrech/ids/xenopusclusterA_mpi.txt")
cluster_B <- read_lines("utrech/ids/xenopusclusterB_mpi.txt")
cluster_C <- read_lines("utrech/ids/xenopusclusterC_mpi.txt")
cluster_D <- read_lines("utrech/ids/xenopusclusterD_mpi.txt")
significant_anova <- Reduce(union, list(cluster_A,cluster_B,cluster_C,cluster_D))
human_CDK1targets <- read_lines("utrech/ids/humanCDKtargets_mpi.txt")
xenopus_ANOVA_data <- read_delim("utrech/Xen_phospho_AnovaPhosphosites.txt","\t", escape_double = FALSE, trim_ws = TRUE)
xenopus_extract_ANOVA_data <- read_delim("utrech/Xen_Extracts_phospho_AnovaPhosphosites.txt","\t", escape_double = FALSE, trim_ws = TRUE)

# # Read in pre-calculated file or perform iupred
# if (!file.exists("all-predictions.Rdata")) {
#   all_predictions <-
#     iupred.predict.fastafile("phosphositesbase_MPI.FASTA",
#                              restrict = unique(all_protein_id_phospho
#                                                $Protein),
#                              idPattern = "^[^|]+\\|([^ ]+).*")
#   save(all_predictions, file="all-predictions.Rdata")
# } else {
#   print ("Loading prepared file")
#   load("all-predictions.Rdata")




disoPath <- "predictions/IUpred_run_xenopus/"
file.names <- dir(disoPath, pattern =".iupred")
ids <-unlist(strsplit(x = file.names,split = ".iupred",fixed = T))
all_predictions <- data.frame(ids)
colnames(all_predictions) <- "ID"
all_predictions$ID <- as.character(all_predictions$ID)
positions <- list()
scores <- list()
disordered <- list()

for (n in 1:length(file.names)) {
  aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  aux_table <- as.data.frame(aux_table)
  colnames(aux_table) <- c("POS","RES","IUPRED_SCORE")
  aux_table$IUPRED_DISO <- aux_table$IUPRED_SCORE>=0.5
  all_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
  positions[[n]] <- as.numeric(aux_table$POS)
  scores[[n]] <- as.numeric(aux_table$IUPRED_SCORE)
  # all_predictions[n,"phospho"] <- NULL
  all_predictions[n,"threshold"] <- 0.5
  disordered[[n]] <- which(aux_table$IUPRED_DISO)
}
all_predictions$positions <-positions
all_predictions$IUPredScores <-scores
all_predictions$disordered <-disordered

# remove unannotated proteins
# all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% all_predictions$ID]
all_predictions_phospho <- subset(all_predictions, ID %in% all_protein_id_phospho)
all_predictions_phospho$length <- nchar(all_predictions_phospho$sequence)

# Introduce phosphosites (from all_proteins variable) using the apply function
phosphosites <- read_delim("all_phosphosites_xenopus_nr.tab", "\t", escape_double = FALSE, col_types = cols(`Leading proteins` = col_skip(), Protein = col_skip()), trim_ws = TRUE)

# Add phosphosites from the extracts!!! uncomment if wanted.
phosphosites_extracts <- read_delim("all_phosphosites_xenopus_extracts_nr.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
phosphosites <- rbind(phosphosites,phosphosites_extracts)


phosphosites <- phosphosites %>% dplyr::rename(ID=Proteins,psites=`Positions within proteins`,seqWindow=`Sequence window`,UID=`Unique identifier`) %>% group_by(ID) %>% summarise_at(c("psites","seqWindow","UID"),function(x){paste(x, collapse=",")})
# The data has been opened by excel and some proteins have their names modified to be dates.


phosphosites$psites <- lapply(phosphosites$psites, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
phosphosites$psites_count <- sapply(phosphosites$psites, length)
phosphosites$UID <- lapply(phosphosites$UID, function(x){ return(as.character(strsplit(x,",")[[1]]))})

xenopus_dynamic <- xenopus_ANOVA_data$UID
# Add phosphosites from the extracts!!! uncomment if wanted.
xenopus_dynamic <- c(xenopus_dynamic,xenopus_extract_ANOVA_data$UID)
#

phosphosites$anova_psites <- apply(phosphosites, 1, function(x){x$psites[x$UID %in% xenopus_dynamic]})
all_predictions_phospho <- merge.data.frame(all_predictions_phospho,phosphosites,by = "ID")


# phosphoDiso_obs <- list()
# phosphoDiso_expct_uni <- numeric()
# for (i in 1:nrow(all_predictions_phospho)) {
#   diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
#   phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
#   phosphoDiso_obs[[i]] <- sum(all_predictions_phospho[i,"psites"][[1]] %in% all_predictions_phospho[i,"disordered"][[1]])
# }


all_predictions_phospho <- all_predictions_phospho %>% mutate(anova_sig=case_when(
  ID %in% significant_anova ~ "Dynamic",
  TRUE ~ "Non dynamic"
))


all_predictions_phospho <- all_predictions_phospho %>% mutate(hCDK1target=case_when(
  ID %in% human_CDK1targets ~ "Human CDK1 target",
  TRUE ~ "other"
))



phosphoDiso_ST <- list()
phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
phosphoDiso_expct_prob <- numeric()
for (i in 1:nrow(all_predictions_phospho)) {
  # diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
  # phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
  # Fraction of Ser And Thr that fall in disorder region
  TStotalIndexes <- as.numeric(gregexpr("S|T", all_predictions_phospho[i,"sequence"])[[1]]) 
  TSinDiso_count <- sum(TStotalIndexes %in% all_predictions_phospho[i,"disordered"][[1]])
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_ST[[i]] <- TStotalIndexes
  phosphoDiso_expct_prob[i] <- TSinDiso_fraction
  phosphoDiso_expct[i] <- all_predictions_phospho[i,"psites_count"]*TSinDiso_fraction
  phosphoDiso_obs[i] <- sum(all_predictions_phospho[i,"psites"][[1]] %in% all_predictions_phospho[i,"disordered"][[1]])
}

all_predictions_phospho$ST_residues <- phosphoDiso_ST
all_predictions_phospho$psites_obsv_diso <- phosphoDiso_obs
all_predictions_phospho$psites_expct_diso <- phosphoDiso_expct
all_predictions_phospho$psites_expct_diso_prob <- phosphoDiso_expct_prob




# Binomial test:
# Calculate percentage of disordered region from all proteins


all_predictions_phospho <- as.data.table(all_predictions_phospho)
all_predictions_phospho[,binom := purrr::pmap(.(psites_obsv_diso, psites_count, psites_expct_diso_prob), binom.test, alternative="greater")]
all_predictions_phospho[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
all_predictions_phospho[,binom_q := mapply(p.adjust, binom_p)]
all_predictions_phospho[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]



ggplot(all_predictions_phospho) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig,shape=anova_sig),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.60),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))+
  scale_shape_manual(values = c(16,4,4))




melted_all_predictions_phospho <- all_predictions_phospho %>% melt(id.vars=c("ID","anova_sig","hCDK1target","length"),value.name = "psites_diso",measure.vars=c("psites_expct_diso","psites_obsv_diso"))



# ALL
ggplot(melted_all_predictions_phospho) + 
  geom_boxplot(aes(y=psites_diso,x=variable, fill = variable,color = variable),outlier.shape = NA)+
  ggpubr::theme_classic2()  + 
  theme(text = element_text(size=20),legend.position = "none",axis.ticks.x = element_blank()) +
  geom_segment(aes(x = 1, y = 10.1, xend = 2, yend = 10.1)) + annotate(geom="text", x=1.5, y=10.3, label="***",size=10) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_discrete(labels = c("Expected","Observed")) + xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 11),breaks = c(seq(0, 11, by = 2)),expand = c(0.05,0.05))+ ylab("Phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(7,10)]) +
  scale_fill_manual(values = pal_jco()(10)[c(2,5)])

# ANOVA +
ggplot(subset(melted_all_predictions_phospho,anova_sig == "Dynamic")) + 
  geom_boxplot(aes(y=psites_diso,x=variable, fill = variable,color = variable),outlier.shape = NA)+
  ggpubr::theme_classic2()  + 
  theme(text = element_text(size=20),legend.position = "none",axis.ticks.x = element_blank()) +
  geom_segment(aes(x = 1, y = 10.1, xend = 2, yend = 10.1)) + annotate(geom="text", x=1.5, y=10.3, label="***",size=10) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_discrete(labels = c("Expected","Observed")) + xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 11),breaks = c(seq(0, 11, by = 2)),expand = c(0.05,0.05))+ ylab("Phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(7,10)]) +
  scale_fill_manual(values = pal_jco()(10)[c(2,5)])

# HumanCDK targets in Dynamic
ggplot(subset(melted_all_predictions_phospho,hCDK1target == "Human CDK1 target" & anova_sig == "Dynamic")) + 
  geom_boxplot(aes(y=psites_diso,x=variable, fill = variable,color = variable),outlier.shape = NA)+
  ggpubr::theme_classic2()  + 
  theme(text = element_text(size=20),legend.position = "none",axis.ticks.x = element_blank()) +
  geom_segment(aes(x = 1, y = 10.1, xend = 2, yend = 10.1)) + annotate(geom="text", x=1.5, y=10.3, label="***",size=10) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_discrete(labels = c("Expected","Observed")) + xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 11),breaks = c(seq(0, 11, by = 2)),expand = c(0.05,0.05))+ ylab("Phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(7,10)]) +
  scale_fill_manual(values = pal_jco()(10)[c(2,5)])

plotList <- IUpredScoresPlotGenerator(subset(all_predictions_phospho,anova_sig == "Dynamic"))


# 
# pdf("IUpredScores.pdf",width = 15,height = 3)
#  for (plot in plotList) {
#    print(plot)
#  }
# dev.off()

lilaSetXenopus <- IUpredScoresPlotGenerator(subset(all_predictions_phospho,ID %in% c("coil","npm","ki67","tp53b","nup53","nup98")),sites_col = "psites")

# JM_highly_phospho_mlos <- c("cndd3","dnli1","ube4b","tsc2","at2b1","rptor","caf1b","pcm1","abcf1","gemi5","cq028","tdrkh","tdrd6","ctr9","sf3b1","rbp2","armc9","chsp1","dkc1","eif3a","tacc3")
# 
# plotList_JM_highly_phospho_mlos <- IUpredScoresPlotGenerator(subset(all_predictions_phospho,ID %in% JM_highly_phospho_mlos))

pdf("../exportImages/pdfs/suppFig5/IUpredScores_plotList_JM_highly_phospho_mlos",width = 15,height = 3)
 for (plot in plotList_JM_highly_phospho_mlos) {
   print(plot)
 }
dev.off()


#######################################  HUMAN (Lila subset) ###############################################################


disoPath <- "predictions/IUpred_run_human/"
file.names <- dir(disoPath, pattern =".iupred")
ids <-unlist(strsplit(x = file.names,split = ".iupred",fixed = T))
all_predictions <- data.frame(ids)
colnames(all_predictions) <- "ID"
all_predictions$ID <- as.character(all_predictions$ID)
positions <- list()
scores <- list()
disordered <- list()

for (n in 1:length(file.names)) {
  aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  aux_table <- as.data.frame(aux_table)
  colnames(aux_table) <- c("POS","RES","IUPRED_SCORE")
  aux_table$IUPRED_DISO <- aux_table$IUPRED_SCORE>=0.5
  all_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
  positions[[n]] <- as.numeric(aux_table$POS)
  scores[[n]] <- as.numeric(aux_table$IUPRED_SCORE)
  # all_predictions[n,"phospho"] <- NULL
  # all_predictions[n,"threshold"] <- 0.5
  disordered[[n]] <- which(aux_table$IUPRED_DISO)
}
all_predictions$positions <-positions
all_predictions$IUPredScores <-scores
all_predictions$disordered <-disordered

# remove unannotated proteins
# all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% all_predictions$ID]
all_predictions$length <- nchar(all_predictions$sequence)





human_data <- read_delim("PSP/human_data_curated.tab", 
                         "\t", escape_double = FALSE, col_types = cols(MOD_RSD = col_character()), 
                         trim_ws = TRUE)

human_data$psites_CDK1 <- lapply(human_data$MOD_RSD, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
human_data$target <- rep("Cdk1 target",nrow(human_data))
human_data$target <- as.factor(human_data$target)

human_universe_data <- read_delim("PSP/Phosphorylation_site_dataset", 
                                  "\t", escape_double = FALSE, col_types = cols(HU_CHR_LOC = col_skip(), 
                                                                                SITE_GRP_ID = col_skip(), MW_kD = col_skip(), 
                                                                                DOMAIN = col_skip(), `SITE_+/-7_AA` = col_skip(), 
                                                                                LT_LIT = col_skip(), MS_LIT = col_skip(), 
                                                                                MS_CST = col_skip(), `CST_CAT#` = col_skip()), 
                                  trim_ws = TRUE)

human_universe_data<-rename(human_universe_data,c(`ACC#`=ACC_ID))
# Select human proteins
human_universe_data <- subset(human_universe_data, ORGANISM == "human")
# Remove all information related to isoforms
human_universe_data <- subset(human_universe_data, !grepl("-",`ACC#`))
human_universe_data <- subset(human_universe_data, !grepl(" iso[0-9]",`PROTEIN`))
# Select targets with pS or pT
human_universe_data<-subset(human_universe_data,(substr(MOD_RSD,1,1)=="S"|substr(MOD_RSD,1,1)=="T"))
# format the MOD_RSD column and group by gene/protein/uniprot
human_universe_data <- human_universe_data %>% mutate(MOD_RSD=substr(MOD_RSD,1,nchar(MOD_RSD)-2)) 
human_universe_data <- human_universe_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at("MOD_RSD",function(x){paste(substr(x,2,2000), collapse=",")})
# Generate the psite column, with the vectors containing psites (for the contingency table analysis)
human_universe_data <- as.data.frame(human_universe_data)
human_universe_data$psites <- lapply(human_universe_data$MOD_RSD, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})

# merge the tables and mark CDK1 targets and non CDK1 targets. Only by ACC, protein names in human data are still with  the isoform nomenclature
human_data <- merge.data.frame(x = human_data,y = human_universe_data,by = c("ACC#"),all = T,suffixes = c("_CDK1","_ALL"))
#remove PROTEIN and GENE comumns from x, and rename the y columns
human_data$GENE_CDK1<-NULL
human_data$PROTEIN_CDK1<-NULL
human_data <- rename(human_data,c(PROTEIN=PROTEIN_ALL,GENE=GENE_ALL))
# Adding the target category "Non Cdk1 target"
levels(human_data$target) <- c("Cdk1 target","Non Cdk1 target")
human_data$target[is.na(human_data$target)]<-"Non Cdk1 target"

human_data<-merge.data.frame(human_data,all_predictions,by.x = "ACC#",by.y = "ID")
human_data$psites_count <- sapply(human_data$psites, length)
human_data$psites_CDK1_count <- sapply(human_data$psites_CDK1, length)


phosphoDiso_ST <- list()
phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
phosphoDiso_expct_prob <- numeric()
for (i in 1:nrow(human_data)) {
  # diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
  # phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
  # Fraction of Ser And Thr that fall in disorder region
  TStotalIndexes <- as.numeric(gregexpr("S|T", human_data[i,"Sequence"])[[1]]) 
  TSinDiso_count <- sum(TStotalIndexes %in% human_data[i,"disordered"][[1]])
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_ST[[i]] <- TStotalIndexes
  phosphoDiso_expct_prob[i] <- TSinDiso_fraction
  phosphoDiso_expct[i] <- human_data[i,"psites_count"]*TSinDiso_fraction
  phosphoDiso_obs[i] <- sum(human_data[i,"psites"][[1]] %in% human_data[i,"disordered"][[1]])
}

human_data$ST_residues <- phosphoDiso_ST
human_data$psites_obsv_diso <- phosphoDiso_obs
human_data$psites_expct_diso <- phosphoDiso_expct
human_data$psites_expct_diso_prob <- phosphoDiso_expct_prob


human_data <- as.data.table(human_data)
human_data[,binom := purrr::pmap(.(psites_obsv_diso, psites_count, psites_expct_diso_prob), binom.test, alternative="greater")]
human_data[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
human_data[,binom_q := mapply(p.adjust, binom_p)]
human_data[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]



ggplot(human_data) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig,shape=target),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.60),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))+
  scale_shape_manual(values = c(16,4,4))






# cdkTargetsPlots<- IUpredScoresPlotGenerator(subset(human_data,target=="Cdk1 target"),id_col = "ACC#",sites_col = "psites_CDK1")


lila_ps_mitotic_ProDir <- read_delim("lila_ps_mitotic_ProDir.tab","\t", escape_double = FALSE, col_types = cols(Pro_Directed_Psites = col_character()),trim_ws = TRUE)
lila_ps_mitotic_ProDir <- merge.data.frame(lila_ps_mitotic_ProDir,human_data,all.x = T,by.x = "ID",by.y = "ACC#")
lila_ps_mitotic_ProDir$Pro_Directed_Psites <- lapply(lila_ps_mitotic_ProDir$Pro_Directed_Psites, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})

PhaseSep_mitotic_PorDirPsites<- IUpredScoresPlotGenerator(lila_ps_mitotic_ProDir,id_col = "ID",sites_col = "Pro_Directed_Psites",subset_sites_col = "psites_CDK1")


# ______________________________________________________TEST ALL KINASES________________________________________________________

all_predictions_V2<-merge.data.frame(all_predictions,human_data[,c("ACC#","psite_CDK1","psite_MAPK","psite_aurk","target_all","ST_residues")],by.x = "ID",by.y = "ACC#",all.x = T)



test<-IUpredScoresPlotGenerator_AllKinases(subset(all_predictions_V2,target_all=="CDK & MAPK & AURK"),id_col="ID",sites_CDK="psite_CDK1",sites_MAPK="psite_MAPK",sites_AURK="psite_aurk",sequence_col="Sequence")



