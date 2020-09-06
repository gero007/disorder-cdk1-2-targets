library(readr)
library(ggseqlogo)
library(Biostrings)
require(ggplot2)
require(ggseqlogo)
library(dplyr)
library(stringr)


phosphoSiteBackground <- read_delim("utrech/Phosphoproteomics_Shotgun_Clustering_Gero/allPhosphosites_background.tab", "\t",
                                    escape_double = FALSE,
                                    col_types = cols(
                                      ANOVA = col_factor(levels = c("+","-"))),
                                    trim_ws = TRUE)
phosphoSiteBackground <- phosphoSiteBackground %>% mutate(type=case_when(
  substr(sequence,start = 17,stop = 17) == "P" ~ "Proline directed",
  TRUE ~ "Non proline directed"
))

phosphoSiteBackground <- phosphoSiteBackground %>% mutate(sequence_red=  substr(sequence,start = 9,stop = 23))

#______________________________________________________________________________________________________________________________________#

phosphoSiteClusters <- read_delim("utrech/Phosphoproteomics_Shotgun_Clustering_Gero/clusterGero.tsb", "\t",
                                  escape_double = FALSE, 
                                  col_types = cols(
                                    Cluster = col_factor(levels = c("cluster1","cluster2", "cluster3", "cluster4")),
                                    sequence = col_character()), 
                                  trim_ws = TRUE)
phosphoSiteClusters <- phosphoSiteClusters %>% mutate(type=case_when(
  substr(sequence,start = 17,stop = 17) == "P" ~ "Proline directed",
  TRUE ~ "Non proline directed"
))

phosphoSiteClusters <- phosphoSiteClusters %>% mutate(sequence_red=  substr(sequence,start = 9,stop = 23))

#Calculate the background estimated probabilities for each aminoacid. All the phosphosites detected have been used for stablishing the background calculations
backgroundSeqs <- Biostrings::AAStringSet(unlist(phosphoSiteBackground$sequence_red))
names(backgroundSeqs)<-phosphoSiteBackground$ID

bkgEstProbs <- colSums(alphabetFrequency(backgroundSeqs))[1:20]/sum(colSums(alphabetFrequency(backgroundSeqs))[1:20])




informationContent_calculator <- function(sequences,backgroundProb){
  if (var(nchar(sequences))!=0){
    stop("Sequences with different lenghts cannot be aligned")
  } else {
    
    windowsLength <- unique(nchar(sequences))
    sequencesCount <- length(sequences)
    aux_matrix <- matrix(as.factor(unlist(strsplit(sequences,""))),byrow = T,nrow=sequencesCount,ncol = windowsLength)
    # not working now. Should return a count sumarrry of each aminoacide
    out_matrix <- matrix(nrow = 20, ncol = windowsLength)
    rownames(out_matrix) <- AA_ALPHABET[1:20]
    
    for (i in 1:ncol(aux_matrix)) {
      
      colCounts <- summary(as.factor(aux_matrix[,i]))
      # colCounts <- str_split_fixed(aux_matrix[,i],":( )*",2)
      # aux_names <- colCounts[,1]
      # colCounts <- as.numeric(trimws(colCounts[,2]))
      # names(colCounts) <- aux_names
      colCounts <- colCounts[(names(colCounts)!="_" & names(colCounts)!="")]
      ObservedProb <- colCounts/sum(colCounts)
      for (aminoacid in row.names(out_matrix)) {
        
        out_matrix[aminoacid,i] <- ObservedProb[aminoacid]*log2(ObservedProb[aminoacid]/backgroundProb[aminoacid])
      
      }
      
    }
  }
  out_matrix[is.na(out_matrix)] <- 0
  norm_out_matrix <- apply(out_matrix, 2, function(x){x/sum(abs(x))})
  return(norm_out_matrix)
}

IC_matrix_list <- list()

# Anova +
IC_matrix_list[["ANOVA positive - All phosphosites"]] <- informationContent_calculator(phosphoSiteClusters$sequence_red,bkgEstProbs)
IC_matrix_list[["ANOVA positive - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,type=="Proline directed")$sequence_red,bkgEstProbs)
IC_matrix_list[["ANOVA positive - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,type=="Non proline directed")$sequence_red,bkgEstProbs)

# Cluster A (cluster3)

IC_matrix_list[["Cluster A - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster3")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster A - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster3" & type=="Proline directed")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster A - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster3" & type=="Non proline directed")$sequence_red,bkgEstProbs)

# Cluster B (cluster4)

IC_matrix_list[["Cluster B - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster4")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster B - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster4" & type=="Proline directed")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster B - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster4" & type=="Non proline directed")$sequence_red,bkgEstProbs)

# Cluster C (cluster1)

IC_matrix_list[["Cluster C - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster1")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster C - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster1" & type=="Proline directed")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster C - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster1" & type=="Non proline directed")$sequence_red,bkgEstProbs)

# Cluster D (cluster2)

IC_matrix_list[["Cluster C - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster2")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster C - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster2" & type=="Proline directed")$sequence_red,bkgEstProbs)
IC_matrix_list[["Cluster C - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster2" & type=="Non proline directed")$sequence_red,bkgEstProbs)

# all phosphosites
IC_matrix_list[["Whole set - All phosphosites"]] <- informationContent_calculator(phosphoSiteBackground$sequence_red,bkgEstProbs)
IC_matrix_list[["Whole set - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteBackground, type=="Proline directed")$sequence_red,bkgEstProbs)
IC_matrix_list[["Whole set - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteBackground,type=="Non proline directed")$sequence_red,bkgEstProbs)

# Plots
  
ggseqlogo(IC_matrix_list,seq_type='aa',method='custom',ncol = 3) + 
  scale_y_continuous(limits = c(-0.5, 1),breaks = c(seq(-0.5, 1, by = 0.5)),expand = c(0,0)) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),labels=c("-7","-6","-5","-4","-3","-2","-1","P","+1","+2","+3","+4","+5","+6","+7")) +
  theme( axis.line.y = element_line(size = 0.5, linetype = "solid"),axis.ticks.y = element_line(size = 0.5, linetype = "solid")) + 
  geom_hline(yintercept=0, linetype="solid",size = 0.5)

kinaseMotifList <- list()

# Plk: [D/N/E/Y]-X-[S/T]-[ ?? /F; no P]-[ ?? /X] 
kinaseMotifList[["Plk"]] <- "^.{5}[D|N|E|Y].[S|T][A|V|I||L|F|W|Y|M]"

# # Plk1: [D/E]-X-[S/T]-[??]-X-[D/E] (https://www.jbc.org/content/278/28/25277.long)
# kinaseMotifList[["Plk1"]] <- "^.{5}[D|E].[S|T][A|V|I||L|F|W|Y|M].[D|E]"

# Aurora: R/K-X-S/T-[ ?? /F; no P] -> [KR].[ST][^P] (https:///pwww.ncbi.nlm.nih.govmc/articles/PMC2999863/#b29)
kinaseMotifList[["Aurora A/B"]] <- "^.{5}[R|K].[S|T][^P]" 

# Nek:  (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6570481/)
# Nek general: [L|M|F|W].[S|T][A|V|I||L|F|W|Y|M][K|R]
kinaseMotifList[["Nek1/2/3/4/5/8/10"]] <- "^.{5}[L|M|F|W].[S|T][A|V|I||L|F|W|Y|M][K|R]"
# Nek6/7/9":[F/L/M/W]-X-R-T
kinaseMotifList[["Nek6/7/9"]] <- "^.{4}[L|M|F][D|E|N][D|E|Y]S[A|V|I||L|F|W|Y|M]"

# Casein kinase 1: D/E-D/E-D/E-X-X-S/T or # pS/pT-X-X-S/T 
kinaseMotifList[["Ck1"]] <- "^((.{2}[E|D][E|D][E|D])|(.{4}[S|T]|)).{2}[S|T]"

# Casein kinase 2: [pS/pT]-{P}-x-[E/D] or [pS/pT]-{P}-x-pS, (https://www.sciencedirect.com/science/article/pii/S2405580815000680)
kinaseMotifList[["Ck2"]] <- "^.{6}[S|T][S|T].[E|D|S]"

# CDC7: [E/D]-[pS/pT]-[E/D]-x-[E/D] 
kinaseMotifList[["Cdc7"]] <- "^((.{7}[S|T][E|D].[E|D])|(.{7}[S|T][S|T][P]))"

# PKA: R-R/K-X-S-??
kinaseMotifList[["Pka"]] <- "^.{4}R[R|K].[S|T][A|V|I||L|F|W|Y|M]"

# Cdk1 full
kinaseMotifList[["Cdk full"]] <- "^.{7}[S|T]P.[K|R]"

# Cddk minimal
kinaseMotifList[["Cdk minimal"]] <- "^.{7}[S|T]P"

# TOTAL
kinaseMotifList[["Total"]] <- "^.{7}[S|T|Y]"



motifCounts_clusterA <- sapply(kinaseMotifList, function(x){
  length(grep(x,subset(phosphoSiteClusters, Cluster=="cluster3")$sequence_red))
})
motifCounts_clusterB <- sapply(kinaseMotifList, function(x){
  length(grep(x,subset(phosphoSiteClusters, Cluster=="cluster4")$sequence_red))
})
motifCounts_clusterC <- sapply(kinaseMotifList, function(x){
  length(grep(x,subset(phosphoSiteClusters, Cluster=="cluster1")$sequence_red))
})
motifCounts_clusterD <- sapply(kinaseMotifList, function(x){
  length(grep(x,subset(phosphoSiteClusters, Cluster=="cluster2")$sequence_red))
})
motifCounts_all <- sapply(kinaseMotifList, function(x){
  length(grep(x,phosphoSiteClusters$sequence_red))
})

motifCounts_table <- rbind(motifCounts_clusterA,motifCounts_clusterB,motifCounts_clusterC,motifCounts_clusterD,motifCounts_all)

row.names(motifCounts_table) <- c("Cluster A","Cluster B","Cluster C","Cluster D","Total")
# colnames(NonProDirected_motifCounts_table) <- c("Plk","Aurora","Nek","Ck1","Unknown","Ck2","Pka")

motifCounts_table
