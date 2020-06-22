########################################
############## Mobi-DB  ################
########################################

#All the anotations in the prediction are for disordered regions

getPredDiso <- function(df,output){
  predDisoPercAll <- c()
  predStretchDist <- c()
  disoPosList <- list()
  disoSeqList <- list()
  for (i in 1:nrow(df)) {
    
    simpleIndex<-grep("simple",df$disorder$predictors[[i]]$method)
    if (length(simpleIndex)!=0) {
      mat<-df$disorder$predictors[[i]]$regions[[simpleIndex]]
      if (class(mat)[1]=="matrix") {
        predStretchDist <- c(predStretchDist,as.numeric(mat[,2])-as.numeric(mat[,1])+1)
        #calculate the number of disoredered AAs
        nDiso<-sum(as.numeric(mat[,2])-as.numeric(mat[,1])+1)
        proteinLength <- nchar(df$sequence[i])
        predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
        
        disoPos <-numeric()
        disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
        return(disoPos)})
        disoPos <- as.numeric(unlist(disoPos))
        disoPosList[[i]]<-disoPos
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,"sequence"],startIdx,endIdx),sep = "")
        }
        disoSeqList[[df$acc[[i]]]]<-disoSeq
      } else {predDisoPercAll <- c(predDisoPercAll,0)
      disoPosList[[i]]<-0
      disoSeqList[[df$acc[[i]]]]<-""
      }
    } else {predDisoPercAll <- c(predDisoPercAll,NA)
    disoPosList[[i]]<-NA
    }
    # mat<-df$disorder$predictors[[i]]$regions[[simpleIndex]]
    # predStretchDist <- c(predStretchDist,as.numeric(mat[,2])-as.numeric(mat[,1])+1)
    # #calculate the number of disoredered AAs
    # nDiso<-sum(as.numeric(mat[,2])-as.numeric(mat[,1])+1)
    # proteinLength <- nchar(df$sequence[i])
    # predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
    # 
    # disoPos <-numeric()
    # disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
    # return(disoPos)})
    # disoPos <- as.numeric(unlist(disoPos))
    # disoPosList[[i]]<-disoPos
    # disoSeq <- character()
    # for (row in 1:nrow(mat)) {
    #   startIdx <- as.numeric(mat[row,1])
    #   endIdx <- as.numeric(mat[row,2])
    #   disoSeq <- paste(disoSeq,substr(yeastDiso[i,"sequence"],startIdx,endIdx),sep = "")
    # }
    # disoSeqList[[df$acc[[i]]]]<-disoSeq
  }
  if (output=="percentage") {
    return(predDisoPercAll)
  } else if (output=="lenghts") {
    return(predStretchDist)
  } else if (output=="indices") {
    return(disoPosList)
  } else if (output=="sequences") {
    return(disoSeqList)
  } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}
}


getPredLiteDiso <- function(df,output){
  predLiteDisoPercAll <- c()
  predLiteStretchDist <- c()
  disoPosList <- list()
  disoSeqList <- list()
  for (i in 1:nrow(df)) {
    # the regular expression is meant to capture the names of the methods for human and yeast database. one is dash and the other underscore
    liteIndex<-grep("mobidb[_-]lite",df$disorder$predictors[[i]]$method)
    if (length(liteIndex)!=0) {
      mat<-df$disorder$predictors[[i]]$regions[[liteIndex]]
      if (class(mat)[1]=="matrix") {
        predLiteStretchDist <- c(predLiteStretchDist,as.numeric(mat[,2])-as.numeric(mat[,1])+1)
        #calculate the number of disoredered AAs
        nDiso<-sum(as.numeric(mat[,2])-as.numeric(mat[,1])+1)
        proteinLength <- nchar(df$sequence[i])
        predLiteDisoPercAll <- c(predLiteDisoPercAll,(nDiso/proteinLength)*100)
        
        disoPos <-numeric()
        disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
        return(disoPos)})
        disoPos <- as.numeric(unlist(disoPos))
        disoPosList[[i]]<-disoPos
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,"sequence"],startIdx,endIdx),sep = "")
        }
        disoSeqList[[df$acc[[i]]]]<-disoSeq
      } else {predLiteDisoPercAll <- c(predLiteDisoPercAll,0)
        disoPosList[[i]]<-0
        disoSeqList[[df$acc[[i]]]]<-""
      }
    } else {predLiteDisoPercAll <- c(predLiteDisoPercAll,NA)
    disoPosList[[i]]<-NA
    }
  }
  if (output=="percentage") {
    return(predLiteDisoPercAll)
  } else if (output=="lenghts") {
    return(predLiteStretchDist)
  } else if (output=="indices") {
    return(disoPosList)
  } else if (output=="sequences") {
    return(disoSeqList)
  } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}
}



#For the derived methods (full -> the third matrix) the anotations could be Disorder, Structure or Conflict
getDerivedDiso <- function(df,output){
  derivedStretchDist <- c()
  derivedDisoPercAll <- c()
  disoPosList <- list()
  disoSeqList <- list()
  for (i in 1:nrow(df)) {
    if (!is.null(df$disorder$derived[[i]])) {
      fullIndex<-which(df$disorder$derived[[i]]$method=="full")
      mat<-df$disorder$derived[[i]]$regions[[fullIndex]]
      #getting only the disordered rows. Is really important to reconvert to a matrix of ncol=3 just in case I index only one row and I get a vector. the the subindex mat[,2] will throw incorrect number of dimensions
      mat<-matrix(mat[mat[,3]=="D",],ncol=3)
      derivedStretchDist <- c(derivedStretchDist,as.numeric(mat[,2])-as.numeric(mat[,1])+1)
      #calculate the number of disoredered AAs
      nDiso<-sum(as.numeric(mat[,2])-as.numeric(mat[,1])+1)
      proteinLength <- nchar(df$sequence[i])
      derivedDisoPercAll <- c(derivedDisoPercAll,(nDiso/proteinLength)*100)
      # print((nDiso/proteinLength)*100)
      if (nrow(mat)>0) {
        disoPos <-numeric()
        # print(i)
        disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
        return(disoPos)})
        disoPos <- as.numeric(unlist(disoPos))
        disoPosList[[i]]<-disoPos
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,"sequence"],startIdx,endIdx),sep = "")
        }
        disoSeqList[[df$acc[[i]]]]<-disoSeq
      } else {
        disoPosList[[i]]<- 0
        disoSeqList[[df$acc[[i]]]]<-""        
      }

    } else {
      derivedDisoPercAll <- c(derivedDisoPercAll,NA)
      disoPosList[[i]]<-NA
    }
    
  }
  if (output=="percentage") {
    return(derivedDisoPercAll)
  } else if (output=="lenghts") {
    return(derivedStretchDist)
  } else if (output=="indices") {
    return(disoPosList)
  } else if (output=="sequences") {
    return(disoSeqList)
  } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}  
}


#For the db curated assigments (full -> the only matrix) the anotations could be Disorder, Structure or Conflict
getDBDiso <- function(df,output){
  dbDisoPercAll <- c()
  dbStretchDist <- c()
  disoPosList <- list()
  disoSeqList <- list()
  for (i in 1:nrow(df)) {
    if (!is.null(df$disorder$db[[i]])) {
      fullIndex<-which(df$disorder$db[[i]]$method=="full")
      mat<-df$disorder$db[[i]]$regions[[fullIndex]]
      #getting only the disordered rows. Is really important to reconvert to a matrix of ncol=3 just in case I index only one row and I get a vector. the the subindex mat[,2] will throw incorrect number of dimensions
      mat<-matrix(mat[mat[,3]=="D",],ncol=3)
      #calculate the number of disoredered AAs
      dbStretchDist <- c(dbStretchDist,as.numeric(mat[,2])-as.numeric(mat[,1])+1)
      nDiso<-sum(as.numeric(mat[,2])-as.numeric(mat[,1])+1)
      proteinLength <- nchar(df$sequence[i])
      dbDisoPercAll <- c(dbDisoPercAll,(nDiso/proteinLength)*100)
      
      if (nrow(mat)>0) {
        disoPos <-numeric()
        disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
        return(disoPos)})
        disoPos <- as.numeric(unlist(disoPos))
        disoPosList[[i]]<-disoPos
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,"sequence"],startIdx,endIdx),sep = "")
        }
        disoSeqList[[df$acc[[i]]]]<-disoSeq
      } else {
        disoPosList[[i]]<- 0
        disoSeqList[[df$acc[[i]]]]<-""
      }
      
    } else {
      dbDisoPercAll <- c(dbDisoPercAll,NA)
      disoPosList[[i]]<-NA
    }
  }
  if (output=="percentage") {
    return(dbDisoPercAll)
  } else if (output=="lenghts") {
    return(dbStretchDist)
  } else if (output=="indices") {
    return(disoPosList)
  } else if (output=="sequences") {
    return(disoSeqList)
  } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}  
}
# yeastDiso$dbDisoPercAll <- dbDisoPercAll

########################################
################ SPOT ##################
########################################


getSpotDisoMatrix <- function(disoPath){ 
  output_list <- list()
  for (disoFile in dir(path = disoPath,pattern = "*.spotd")) {
    id=strsplit(disoFile,split = "\\.")[[1]][1]
    # print(id)
    # print(paste(disoPath,disoFile,sep = ""))
    auxdf <- read.delim(paste(disoPath,disoFile,sep = ""), header=FALSE, comment.char="#")
    disoRunningLenghts <- rle(as.vector(auxdf$V4))
    start=1
    disoMatrix <- matrix(ncol = 4,nrow = length(disoRunningLenghts$lengths))
    for (i in 1:length(disoRunningLenghts$lengths)) {
      end=start+disoRunningLenghts$lengths[i]-1
      disoMatrix[i,] <- as.character(c(start,end,disoRunningLenghts$values[i],disoRunningLenghts$lengths[i]))
      # print(c(start,end,test$values[i],test$lengths[i]))
      start=end+1
    }
    output_list[[id]] <- disoMatrix
  }
  # print(disoMatrix)
  # print(which(targetDf$acc==id))
  return(output_list)
}


getSpotPredDiso <- function(df,spotDisorderList,output,accession_col){
  predDisoPercAll <- c()
  predStretchDist <- c()
  disoPosList <- list()
  disoSeqList <- list()
  for (i in 1:nrow(df)) {
    mat<-spotDisorderList[[df[i,accession_col]]]
    # print(i)
    # print(accession_col)
    #select only disordered regions
    mat <- matrix(mat[mat[,3]=="D",],ncol=4)
    # print(df$acc[[i]])
    # print(mat)
    if (nrow(mat)!=0) {
      predStretchDist <- c(predStretchDist,as.numeric(mat[,4]))
      #calculate the number of disoredered AAs
      nDiso<-sum(as.numeric(mat[,4]))
      # print(nDiso)
      proteinLength <- nchar(df$sequence[[i]])
      # print(proteinLength)
      predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
      # print((nDiso/proteinLength)*100)
      disoPos <-numeric()
      disoPos <-apply(mat, 1, function(x){disoPos<-c(disoPos,x[1]:x[2])
      return(disoPos)})
      disoPos <- as.numeric(unlist(disoPos))
      disoPosList[[i]]<-disoPos
      disoSeq <- character()
      for (row in 1:nrow(mat)) {
        startIdx <- as.numeric(mat[row,1])
        endIdx <- as.numeric(mat[row,2])
        disoSeq <- paste(disoSeq,substr(df[i,"sequence"],startIdx,endIdx),sep = "")
      }
      disoSeqList[[df[i,accession_col]]]<-disoSeq

    } else {
      nDiso<-0
      proteinLength <- nchar(df$sequence[i])
      predDisoPercAll <- c(predDisoPercAll,(nDiso/proteinLength)*100)
      
      disoPos <-numeric()
      disoPos <- 0
      disoPosList[[i]]<-disoPos
    }
  }
  if (output=="percentage") {
    return(predDisoPercAll)
  } else if (output=="lenghts") {
    return(predStretchDist)
  } else if (output=="indices") {
    return(disoPosList)
  } else if (output=="sequences") {
    return(disoSeqList)
  } else {stop("wrong output specified. Options are 'percentage','lenghts' or 'indices' ")}
}


