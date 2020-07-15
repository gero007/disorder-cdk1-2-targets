library(Biostrings)
library(eulerr)
library(flextable)
library(stringr)
library(DT)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(AnnotationDbi)
library(ggridges)
library(hrbrthemes)
library(ggrepel)
source("getDisorder_V2.R")
source("contingencyTablesAnalysis_V2.R")


xenopus_MPI<-readAAStringSet("utrech/Xenopus_Database_MPI.FASTA",format = "fasta")
xenopusDiso<-as.data.frame(xenopus_MPI)
colnames(xenopusDiso)<-"Sequence"
xenopusDiso$ID <- str_extract(names(xenopus_MPI), regex( "(?<=\\|)(.*)"))
row.names(xenopusDiso) <- NULL
xenopusDiso$Length <-nchar(xenopusDiso$Sequence)

names(iupredJuan)
unique(xenopusDiso$ID)
lenght(unique(xenopusDiso$ID))
length(unique(xenopusDiso$ID))
duplicated(xenopusDiso$ID)
which(duplicated(xenopusDiso$ID))
xenopusDiso$ID[which(duplicated(xenopusDiso$ID))]
all_predictions
all_predictions[[1]]$.positions[all_predictions[[1]]$.disordered]
sapply(all_predictions, function(x){return(x$.positions[x$.disordered])})
iupredJuan <- sapply(all_predictions, function(x){return(x$.positions[x$.disordered])})
# IUpred Yo tiene un elemento mas
iupredYo <- subset(xenopusDiso,ID %in% names(iupredJuan))$iupred_disoIndexes
iupredYo <- list(subset(xenopusDiso,ID %in% names(iupredJuan))$iupred_disoIndexes)
iupredYo <- as.list(subset(xenopusDiso,ID %in% names(iupredJuan))$iupred_disoIndexes)
