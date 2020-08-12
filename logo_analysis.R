library(readr)
library(ggseqlogo)
library(Biostrings)
require(ggplot2)
require(ggseqlogo)



phosphoSiteClusters <- read_delim("utrech/Phosphoproteomics_Shotgun_Clustering_Gero/clusterGero.tsb", 
                                  "\t", escape_double = FALSE, 
                                  col_types = cols(
                                    Cluster = col_factor(levels = c("cluster1","cluster2", "cluster3", "cluster4")),
                                    sequence = col_character()), 
                                  trim_ws = TRUE)

phosphoSiteClusters <- phosphoSiteClusters %>% mutate(type=case_when(
  substr(sequence,start = 17,stop = 17) == "P" ~ "Proline directed",
  TRUE ~ "Non proline directed"
))

seqWindows <- Biostrings::AAStringSet(unlist(phosphoSiteClusters$sequence))
names(seqWindows)<-phosphoSiteClusters$ID

bkgEstProbs <- colSums(alphabetFrequency(seqWindows))[1:20]/sum(colSums(alphabetFrequency(seqWindows))[1:20])


# > which(nchar(phosphoSiteClusters$sequence)!=31)
# [1]  162  445 1010
test<-matrix(unlist(strsplit(phosphoSiteClusters$sequence,"")),byrow = T,nrow=1032,ncol = 31)
PSPlogo_function <- function(){}