library(Biostrings)
library(stringr)
library(ggplot2)

iupred.predict.fastafile <- function(fastafile, model = "long", restrict = c(), idPattern = "^([^ ]+).*") {
    proteins <- read.proteins(fastafile, idPattern)
    proteins <- filter.proteins(proteins, restrict)
    
    lapply(proteins, FUN = run.iupred.exe, model = model)
  }


iupredplot <- function(x) {
  
  x$disordered <- x$score > 0.5
  
  ggplot(x, aes(x=position, y=score)) +
    geom_col(aes(fill=disordered)) + 
    geom_line() + 
    geom_hline(yintercept = 0.5)
}


read.proteins <- function(fastafile, idPattern) {
  
  proteins <- readAAStringSet(fastafile, use.names = T)
  matches <- str_match(names(proteins), pattern = idPattern)[, 2]
  names(proteins) <- matches
  
  proteins
}

filter.proteins <- function(anAaStringSet, restriction) {
  
  anAaStringSet[names(anAaStringSet) %in% restriction]
}

run.iupred.exe <- function(sequence, model) {
  prediction.line.pattern <- "\\s+(\\d+)\\s(.)\\s+([0-9.]+)"
  
  setwd("iupred")
  writeXStringSet(AAStringSet(sequence), "tempseq.fasta")
  
  result <-
    system2("iupred.exe", c("tempseq.fasta", model), stdout = T )
  
  ## Clean up the program's output
  predictionlines <- grep(pattern = prediction.line.pattern, result)
  result <- result[predictionlines]
  
  extractedResults <-
    apply(
      as.data.frame(result),
      MARGIN = 1,
      FUN = str_match,
      pattern = prediction.line.pattern
    )
  
  extractedResults <- t(extractedResults[2:4,])
  
  extractedResults <- as.data.frame(extractedResults, stringsAsFactors=F)
  names(extractedResults) <- c("position","aa","score")
  extractedResults$score <- as.numeric(extractedResults$score)
  extractedResults$position <- as.numeric(extractedResults$position)
  

  file.remove("tempseq.fasta")
  setwd("..")
  
  return(extractedResults)
}

