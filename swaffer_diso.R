library(readr)
library(Biostrings)
library(stringr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
source("getDisorder_V2.R")

swaffer_CDK_psites <- read_delim("swaffer_CDK_psites.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
swaffer_CDK_psites <- as.data.frame(swaffer_CDK_psites)
swaffer_diso <- readAAStringSet("pombe_CDKtargets_Swaffer.fasta",format = "fasta")
swaffer_uniprot <- str_split_fixed(names(swaffer_diso),"\\||_",4)[,2]
swaffer_gene <- str_split_fixed(names(swaffer_diso),"\\||_",4)[,3]
swaffer_diso<-as.data.frame(swaffer_diso)

colnames(swaffer_diso) <- "Sequence"
rownames(swaffer_diso) <- NULL
swaffer_diso$Uniprot <- swaffer_uniprot
swaffer_diso$Gene <-swaffer_gene
swaffer_diso$Length <- nchar(swaffer_diso$Sequence)

swaffer_diso <- getIUpredPredDiso(df = swaffer_diso,disoPath = "predictions/IUpred_run_swaffer/",sequence_col = "Sequence",accession_col = "Uniprot",length_col = "Length")
swaffer_diso$iupl_disoIndexes <- lapply(swaffer_diso$iupl_disoIndexes, function(x){as.numeric(strsplit(x,',')[[1]])})


# Check for every Psite if they are or not in a disordered region

swaffer_CDK_psites$location_iupl <- apply(swaffer_CDK_psites, 1, function(x){as.numeric(x["Positions"]) %in% subset(swaffer_diso,Uniprot==x["UNIPROT ID"])[,"iupl_disoIndexes"][[1]]})
swaffer_CDK_psites <- swaffer_CDK_psites %>% mutate(location_iupl=case_when(
  location_iupl ~ "Disorder",
  TRUE ~ "Structure"
))

ggplot(swaffer_CDK_psites) + 
  geom_point(aes(x=IC50_nMInhib_1,y=IC50_nMInhib_2, colour = `Cell cycle dynamics`,shape=location_iupl),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.82),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  scale_colour_manual(values = pal_jco()(10)[c(3,2,4,1)])+
  scale_shape_manual(values = c(16,4,4)) +
  scale_x_continuous(limits = c(0, 1000),breaks = c(seq(0, nrow(aux_df)+1, by = 100)),expand = c(0.005,0.005))+ xlab("Positions") +
  scale_y_continuous(limits = c(-0.1, 1),breaks = c(seq(0,1,by=0.25)),expand = c(0.05,0.05)) + ylab("Score") +