library(readr)
library(dplyr)
library(data.table)
library(ggplot2)


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

all_predictions_phospho <- merge.data.frame(all_predictions_phospho,phosphosites,by = "ID")


phosphoDiso_obs <- list()
phosphoDiso_expct_uni <- numeric()
for (i in 1:nrow(all_predictions_phospho)) {
  diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
  phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
  phosphoDiso_obs[[i]] <- sum(all_predictions_phospho[i,"psites"][[1]] %in% all_predictions_phospho[i,"disordered"][[1]])
}


all_predictions_phospho <- all_predictions_phospho %>% mutate(anova_sig=case_when(
  ID %in% significant_anova ~ "ANOVA Significant",
  TRUE ~ "ANOVA Significant"
))

all_predictions_phospho <- all_predictions_phospho %>% mutate(cluster=case_when(
  ID %in% cluster_oscillating ~ "Cluster D",
  ID %in% cluster_interphase ~ "Cluster C",
  TRUE ~ "other"
))

phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
for (i in 1:nrow(all_predictions_phospho)) {
  # diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
  # phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
  # Fraction of Ser And Thr that fall in disorder region
  TStotalIndexes <- as.numeric(gregexpr("S|T", all_predictions_phospho[i,"sequence"])[[1]]) 
  TSinDiso_count <- sum(TStotalIndexes %in% all_predictions_phospho[i,"disordered"][[1]])
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_expct[[i]] <- length(all_predictions_phospho[i,"psites"][[1]])*TSinDiso_fraction
  phosphoDiso_obs[[i]] <- sum(all_predictions_phospho[i,"psites"][[1]] %in% all_predictions_phospho[i,"disordered"][[1]])
}

all_predictions_phospho$psites_expct_diso <- phosphoDiso_expct
all_predictions_phospho$psites_obsv_diso <- phosphoDiso_obs

ggplot(all_predictions_phospho) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso),size=2,alpha=0.5)+
  geom_abline(color="darkslategrey",slope = 1,size=1,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=17),legend.position = "none") + 
  scale_x_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ ylab("Expected phospho S/T in IDR")


# apply(all_protein_id_phospho, MARGIN=1, FUN=function(x) {
#   all_predictions[[x[["Protein"]]]] <<- setPhospho(all_predictions[[x[["Protein"]]]] , as.numeric(x[["p-site"]]));
#   return(invisible());
# })

# Get the disordered/ordered phospho counts using countPhosphoDisordered
# phospho_relative_to_organization<- lapply(all_predictions, countPhosphoDisordered)
# phospho_relative_to_organization<- lapply(all_predictions, countPhosphoDisordered)

# unlist into data.table
# phospho_relative_to_organization <- rbindlist(phospho_relative_to_organization, id="protein")
phospho_relative_to_organization <- rbindlist(phosphoDiso_obs, id="protein")

# Convenient structure for batch processing
phospho_region <- list(
  all            = phospho_relative_to_organization,
  anova_sig      = phospho_relative_to_organization[protein %in% significant_anova],
  oscillating    = phospho_relative_to_organization[protein %in% cluster_oscillating,],
  exit_metaphase = phospho_relative_to_organization[protein %in% cluster_exit_metaphase,],
  nonsignificant = phospho_relative_to_organization[! protein %in% significant_anova,]
)


# Calculate percentage of disordered region from all proteins
disordered_fraction <- mapply(all_predictions, FUN=function(x) {sum(x$.disordered) / nchar(x$.sequence)})

lapply(phospho_region, function(x) {
  x[,expected_disordered := (disordered + ordered)  * disordered_fraction[protein] ];
  x[,expected_ordered := (disordered + ordered)  * (1-disordered_fraction[protein]) ];
})


long_phospho_region_list <- lapply(phospho_region, FUN=melt, id="protein")
long_phospho_region <- rbindlist(long_phospho_region_list, id="group")

wide_phospho_region <- rbindlist(phospho_region, id="group")

# Binomial test
wide_phospho_region[,disordered_prop := disordered / (disordered + ordered)]
wide_phospho_region[,expect_disordered_prop := expected_disordered / (expected_disordered + expected_ordered)]
wide_phospho_region[,n := disordered + ordered]

wide_phospho_region <- wide_phospho_region[n > 0]

wide_phospho_region[,binom := purrr::pmap(.(disordered, n, expect_disordered_prop), binom.test, alternative="t")]
wide_phospho_region[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
wide_phospho_region[,binom_q := mapply(p.adjust, binom_p)]
wide_phospho_region[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]


ggplot(long_phospho_region, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + facet_grid(. ~ group)

ggplot(wide_phospho_region) + 
  geom_point(aes(x=disordered_prop, y=expect_disordered_prop, color=binom_sig)) + 
  geom_abline(slope=1, color="red") + 
  theme_minimal()

ggplot(wide_phospho_region) +
  geom_point(aes(x=disordered, y=expected_disordered, colour = binom_sig), alpha = 0.5, shape=16) +
  geom_abline(slope=1, color="red") +
  facet_grid(. ~ group ) + 
  scale_colour_manual(values = c("gray","orange","red")) + 
  theme_minimal()




