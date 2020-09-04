library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggsci)
library(scales)
library(reshape2)



all_protein_id_phospho <- read_lines("utrech/Xen_phospho_allProteins.txt")
cluster_oscillating <- read_lines("utrech/Xen_phospho_ClusterD.txt")
cluster_interphase <- read_lines("utrech/Xen_phospho_ClusterC.txt")
significant_anova <- read_lines("utrech/Xen_phospho_ANOVApos.txt")


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


# Introduce phosphosites (from all_proteins variable) using the apply function
phosphosites <- phosphosites <- read_delim("all_phosphosites_xenopus_nr.tab", "\t", escape_double = FALSE, col_types = cols(`Leading proteins` = col_skip(), Protein = col_skip()), trim_ws = TRUE)
phosphosites <- phosphosites %>% rename(ID=Proteins,psites=`Positions within proteins`,seqWindow=`Sequence window`,UID=`Unique identifier`) %>% group_by(ID) %>% summarise_at(c("psites","seqWindow","UID"),function(x){paste(x, collapse=",")})
# The data has been opened by excel and some proteins have their names modified to be dates.

phosphosites$psites <- lapply(phosphosites$psites, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
phosphosites$psites_count <- sapply(phosphosites$psites, length)
all_predictions_phospho <- merge.data.frame(all_predictions_phospho,phosphosites,by = "ID")


# phosphoDiso_obs <- list()
# phosphoDiso_expct_uni <- numeric()
# for (i in 1:nrow(all_predictions_phospho)) {
#   diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
#   phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
#   phosphoDiso_obs[[i]] <- sum(all_predictions_phospho[i,"psites"][[1]] %in% all_predictions_phospho[i,"disordered"][[1]])
# }


all_predictions_phospho <- all_predictions_phospho %>% mutate(anova_sig=case_when(
  ID %in% significant_anova ~ "ANOVA Significant",
  TRUE ~ "ANOVA Non significant"
))

all_predictions_phospho <- all_predictions_phospho %>% mutate(cluster=case_when(
  ID %in% cluster_oscillating ~ "Oscillating",
  TRUE ~ "other"
))

# Tested. The results coloring the points with the cluster is too confusing. The statistics will be used
# all_predictions_phospho <- all_predictions_phospho %>% mutate(cluster=case_when(
#   ID %in% cluster_oscillating ~ "Cluster D",
#   ID %in% cluster_interphase ~ "Cluster C",
#   TRUE ~ "other"
# ))

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
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.82),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(3,2,4)])

IUpredScoresPlotGenerator <- function(dataframe){
  plotList <- list()
  for(i in 1:nrow(dataframe)){
    plotTittle <- dataframe[i,"ID"]
    aux_df <-  as.data.frame(cbind(dataframe$positions[[i]],dataframe$IUPredScores[[i]]))
    names(aux_df) <- c("positions","IUPredScores")
    
    plotList[[i]] <- ggplot(aux_df,aes(x=positions,y=IUPredScores)) +
      geom_col(aes(fill=IUPredScores >= 0.5 ,color=IUPredScores >= 0.5)) +
      ggpubr::theme_classic2() + 
      theme(text = element_text(size=19),legend.position = "none") +
      scale_x_continuous(limits = c(0, nrow(aux_df)+1),breaks = c(seq(0, nrow(aux_df)+1, by = 50)),expand = c(0.005,0.005))+ xlab("Positions") +
      scale_y_continuous(limits = c(-0.1, 1),breaks = c(seq(0,1,by=0.25)),expand = c(0.05,0.05)) + ylab("IUpred Scores") +
      scale_fill_manual(values = pal_jco()(10)[c(3,4)]) +
      scale_color_manual(values = pal_jco()(10)[c(3,4)]) +
      annotate("point",x = dataframe$psites[[i]],y=rep(-0.06,length(dataframe$psites[[i]])),shape=21,color=pal_jco()(10)[8],fill=pal_jco()(10)[2],size=6) +
      annotate("point",x = dataframe$ST_residues[[i]],y=rep(0,length(dataframe$ST_residues[[i]])),shape=25,color=pal_jco()(10)[10],fill=pal_jco()(10)[5],size=4) +
      geom_hline(yintercept = 0.5,size=0.5,linetype = "dashed",color = "darkslategrey") + 
      ggtitle(plotTittle)
  }
  names(plotList)<-dataframe$ID
  return(plotList)
}

ggplot(all_predictions_phospho) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.82),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(3,2,4)])

# plotList <- IUpredScoresPlotGenerator(subset(all_predictions_phospho,binom_q < 0.05))
# 
# pdf("IUpredScores.pdf",width = 15,height = 3)
#  for (plot in plotList) {
#    print(plot)
#  }
# dev.off()

#######################################  HUMAN (Lila subset) ###############################################################



# disoPath <- "predictions/IUpred_run_xenopus/"
# file.names <- dir(disoPath, pattern =".iupred")
# ids <-unlist(strsplit(x = file.names,split = ".iupred",fixed = T))
# all_predictions <- data.frame(ids)
# colnames(all_predictions) <- "ID"
# all_predictions$ID <- as.character(all_predictions$ID)
# positions <- list()
# scores <- list()
# disordered <- list()
# 
# for (n in 1:length(file.names)) {
#   aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
#   aux_table <- as.data.frame(aux_table)
#   colnames(aux_table) <- c("POS","RES","IUPRED_SCORE")
#   aux_table$IUPRED_DISO <- aux_table$IUPRED_SCORE>=0.5
#   all_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
#   positions[[n]] <- as.numeric(aux_table$POS)
#   scores[[n]] <- as.numeric(aux_table$IUPRED_SCORE)
#   # all_predictions[n,"phospho"] <- NULL
#   all_predictions[n,"threshold"] <- 0.5
#   disordered[[n]] <- which(aux_table$IUPRED_DISO)
# }
# all_predictions$positions <-positions
# all_predictions$IUPredScores <-scores
# all_predictions$disordered <-disordered
# 
# # remove unannotated proteins
# # all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% all_predictions$ID]
# all_predictions_phospho <- subset(all_predictions, ID %in% all_protein_id_phospho)
# 
# 
# # Introduce phosphosites (from all_proteins variable) using the apply function
# phosphosites <- phosphosites <- read_delim("all_phosphosites_xenopus_nr.tab", "\t", escape_double = FALSE, col_types = cols(`Leading proteins` = col_skip(), Protein = col_skip()), trim_ws = TRUE)
# phosphosites <- phosphosites %>% rename(ID=Proteins,psites=`Positions within proteins`,seqWindow=`Sequence window`,UID=`Unique identifier`) %>% group_by(ID) %>% summarise_at(c("psites","seqWindow","UID"),function(x){paste(x, collapse=",")})
# # The data has been opened by excel and some proteins have their names modified to be dates.
# 
# phosphosites$psites <- lapply(phosphosites$psites, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
# phosphosites$psites_count <- sapply(phosphosites$psites, length)
# all_predictions_phospho <- merge.data.frame(all_predictions_phospho,phosphosites,by = "ID")
