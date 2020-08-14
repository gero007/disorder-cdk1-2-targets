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




informationContent_calculator <- function(sequences,backgroundProb){
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

IC_matrix_list <- list()

# Total
IC_matrix_list[["All phosphosites"]] <- informationContent_calculator(phosphoSiteClusters$sequence,bkgEstProbs)
IC_matrix_list[["Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,type=="Proline directed")$sequence,bkgEstProbs)
IC_matrix_list[["Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,type=="Non proline directed")$sequence,bkgEstProbs)

# Cluster A (cluster3)

IC_matrix_list[["Cluster A - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster3")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster A - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster3" & type=="Proline directed")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster A - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster3" & type=="Non proline directed")$sequence,bkgEstProbs)

# Cluster B (cluster4)

IC_matrix_list[["Cluster B - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster4")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster B - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster4" & type=="Proline directed")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster B - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster4" & type=="Non proline directed")$sequence,bkgEstProbs)

# Cluster C (cluster1)

IC_matrix_list[["Cluster C - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster1")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster C - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster1" & type=="Proline directed")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster C - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster1" & type=="Non proline directed")$sequence,bkgEstProbs)

# Cluster D (cluster2)

IC_matrix_list[["Cluster D - All phosphosites"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster2")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster D - Proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster2" & type=="Proline directed")$sequence,bkgEstProbs)
IC_matrix_list[["Cluster D - Non proline directed"]] <- informationContent_calculator(subset(phosphoSiteClusters,Cluster=="cluster2" & type=="Non proline directed")$sequence,bkgEstProbs)

# Plots
  
ggseqlogo(IC_matrix_list,seq_type='aa',method='custom',ncol = 3) + 
  scale_y_continuous(limits = c(-0.5, 1),breaks = c(seq(-0.5, 1, by = 0.5)),expand = c(0,0)) + 
  theme( axis.line.y = element_line(size = 0.5, linetype = "solid"),axis.ticks.y = element_line(size = 0.5, linetype = "solid")) + 
  geom_hline(yintercept=0, linetype="solid",size = 0.5)


# Plk: [D/N/E/Y]-X-[S/T]-[ ϕ /F; no P]-[ ϕ /X]
Plk_motif_re <- "^.{13}[D|N|E|Y].[S|T][V|I||L|F|W|Y|M]" 
# Aurora: R/K-X-S/T-[ ϕ /F; no P]
Aurora_motif_re <- "^.{13}[R|K].[S|T][V|I||L|F|W|Y|M]" 
# Nek: [F/L/M]-X-X-S/T-[ ϕ /F; no P]-[R/H/X]
Nek_motif_re <- "^.{12}[F|L|M].{2}[S|T][V|I||L|F|W|Y|M][R|H|X]"  
# Casein kinase 1: D/E-D/E-D/E-X-X-S/T-ϕ
Csk1_motif_re <- "^.{10}[E|D][E|D][E|D].{2}[S|T][A|L|I|V]"
# pS/pT-X-X-S/T-ϕ 
Unknown_motif_re <- "^.{15}[S|T].{2}[S|T][A|L|I|V]"
# Casein kinase 2: S/T-D/E-X-D/E
Csk1_motif_re <- "^.{15}[S|T][E|D].[E|D]"
# PKA: R-R/K-X-S-ϕ
Pka_motif_re <- "^.{15}[S|T][E|D].[E|D]"


phosphoSiteClusters$sequence[grep(Plk_motif_re,phosphoSiteClusters$sequence)]
