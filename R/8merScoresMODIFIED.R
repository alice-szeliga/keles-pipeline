# adds %>% piping operator
# useful for neater code
library(magrittr)

## Assumed file length (without header)
fileLength <- 32896
threshhold <- 0.5


## helper function written by me

getCols3to5 <- function(scoresDT) {
  ncol(scoresDT) %>% min(5, .) %>% seq %>% setdiff(., 1:2) %>% 
                 scoresDT[,.,with=FALSE] %>% return
}

getRangeCols3to5 <-function(Cols3To5) {
  Cols3To5 %>% apply(.,2,range) %>% c(., Cols3To5) %>% return
}

getEscoresInd <- function(RangeCols3To5, threshhold) {
    EscoresInd <- intersect(which(RangeCols3To5[1,] >=  -threshhold), 
                        which(RangeCols3To5[2,] <= threshhold))
    return(EscoresInd)
}

## Helper function to remove first row (header) from tables
## fileLength is the correct file length without a header
removeHeader <- function(filelength, dataTable) {
  if(nrow(dataTable) == (fileLength + 1)) {
    dataTable = dataTable[-1,]
  }
  return(dataTable)
}

## call in construct_summ file is commented out 
CheckEscoresCol <- function(file.name) {

  scoresDT <- fread(file.name, header=FALSE)
  scoresDT <- removeHeader(filelength, scoresDT)

  Cols3To5 <- getCols3to5(scoresDT)
  RangeCols3To5 <- getRangeCols3to5(Cols3To5)
  EscoresInd <- getEscoresInd(RangeCols3To5, threshhold)
  returnValues <- c(EscoresInd, scoresDT, Cols3To5)

  return(returnValues)

## easier version
 # file.name %>% fread(., header=FALSE) 
 #           %>% removeHeader(filelength, .) 
 #            -> scoresDT
 # scoresDT  %>% getRangeCols3to%(.)
 #           %>% getEscoresInd(., threshhold)
 #           %>% c(., scoresDT, Cols3To5)
 #           %>% return(.)
}


##called in construct_summ
CreateEscores<-function(file.name) {

  # might not work saving it like this
  values <- CheckEscoresCol(file.name)
  EscoresInd <- values[1]; scoresDT <- values[2]; Cols3To5 <- values[3]

  Escores <- Cols3To5[,EscoresInd,with=FALSE]
  

  ## doing something to the 1st column of scoresDT
  ## think unlist turns a matrix into a list, but we're working w a column?
  mersSeq1 <- unlist(scoresDT[,1, with=FALSE])
  ## 2nd column of scoresDT
  mersSeq2 <- unlist(scoresDT[,2, with=FALSE])
  mersSeq <- cbind(mersSeq1, mersSeq2)

  matchInd <- mapply(function(subseqs) unlist(apply(as.matrix(subseqs), 1, function(subseq) which(mersSeq==subseq, arr.ind=TRUE)[1,1])), 
                   subseqList, 
                   SIMPLIFY=FALSE)
  EscoresList <- mapply(function(x) Escores[x], matchInd)
  return(EscoresList)
}


# load("/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_human.Rda")
# sum(abs(enrscoresList$FOXN2-Escores[matchInd[[1]]]))

