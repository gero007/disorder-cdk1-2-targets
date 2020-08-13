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



#Calculate the background estimated probabilities for each aminoacid. All the phosphosites detected have been used for stablishing the background calculations
backgroundSeqs <- Biostrings::AAStringSet(unlist(phosphoSiteBackground$sequence))
names(backgroundSeqs)<-phosphoSiteBackground$ID

bkgEstProbs <- colSums(alphabetFrequency(backgroundSeqs))[1:20]/sum(colSums(alphabetFrequency(backgroundSeqs))[1:20])




InformationContent_calculator <- function(sequences,backgroundProb){
  if (var(nchar(sequences))!=0){
    stop("Sequences with different lenghts cannot be aligned")
  } else {
    
    windowsLength <- unique(nchar(sequences))
    sequencesCount <- length(sequences)
    aux_matrix <- matrix(unlist(strsplit(sequences,"")),byrow = T,nrow=sequencesCount,ncol = windowsLength)
    aux_matrix <- summary(aux_matrix,maxsum = 100)
    out_matrix <- matrix(nrow = 20, ncol = windowsLength)
    rownames(out_matrix) <- AA_ALPHABET[1:20]
    
    for (i in 1:ncol(aux_matrix)) {
      
      colCounts <- str_split_fixed(aux_matrix[,i],":( )*",2)
      aux_names <- colCounts[,1]
      colCounts <- as.numeric(trimws(colCounts[,2]))
      names(colCounts) <- aux_names
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

