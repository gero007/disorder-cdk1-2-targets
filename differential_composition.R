library(Biostrings)
library(dplyr)
library(reshape2)

getDifferentialComposition <- function(df,suffix="_disoSeq",seq_col,outputType="ggplot"){
  
  predictorNames <- c("Disembl 465","Disembl Hot loops","Espritz DisProt","Espritz NMR","Espritz Xray","GlobProt","IUPred long","IUPred short","JRONN","Pfilt","SEG","SPOT","VSL2b")
  diso_seq_cols <- sort(colnames(df)[grep(pattern = suffix,colnames(df))])
  
  differentialComposition <- list()
  for (i in 1:length(diso_seq_cols)) {
    universe_idx <- !is.na(df[,diso_seq_cols[i]])
    # Disorder
    disoSequenceList <- df[universe_idx,diso_seq_cols[i]]
    disoAAstringSet <- Biostrings::AAStringSet(unlist(disoSequenceList))
    disoComposition <- colSums(alphabetFrequency(disoAAstringSet))[1:20]/sum(colSums(alphabetFrequency(disoAAstringSet))[1:20])
    # Universe
    univSequenceList <- df[universe_idx,seq_col]
    univAAstringSet <- Biostrings::AAStringSet(unlist(univSequenceList))
    univComposition <- colSums(alphabetFrequency(univAAstringSet))[1:20]/sum(colSums(alphabetFrequency(univAAstringSet))[1:20])
    
    differentialComposition[[i]] <- (disoComposition-univComposition)/univComposition
    
  }
  names(differentialComposition) <- predictorNames
  outputDF<- tibble::rownames_to_column(as.data.frame(differentialComposition),"AMINOACID")
  outputDF$AMINOACID <- as.factor(outputDF$AMINOACID)
  if(outputType=="ggplot"){
    outputDF <- melt(test,id.vars = "AMINOACID")
    colnames(outputDF)[2:3] <- c("Predictor","Differential.Composition")
    return(outputDF)
  } else if (outputType=="table") {
    return(outputDF)
  } else {stop("The outputType parameter should be 'table' or 'ggplot' ")}
  
  }  