# adds %>% piping operator
# useful for neater code
library(magrittr)

## Assumed file length (without header)
fileLength <- 32896
threshhold <- 0.5


## helper functions written by me


#' Returns 3rd to 5th columns
#' 
#' \code{getCols3to5} returns the 3rd - 5th columns of a data.table.
#' If the table has less than 5 columns, it returns 3rd - max.
#' If the table has less than 3 columns, it returns a null table.
#'
#' @param DT data.table 
#' @return a data.table of up to 3 columns
getCols3to5 <- function(DT) {
  ncol(DT) %>% min(5, .) %>% seq %>% setdiff(., 1:2) %>%
                 DT[,.,with=FALSE] %>% return
}

#' Get range of a data.table's columns
#' 
#' \code{getRangeCols} finds the minimum and maximum of each column in DT, then
#' returns a data table with the minimum in the first row and the maximum in
#' the second row.
#' 
#' @param DT data.table
#' @return a data.table with 2 rows and the same number of columns
getRangeCols <-function(DT) {
  apply(DT, 2, range) %>% return
}

getEscoresInd <- function(RangeCols3To5, threshhold) {
    EscoresInd <- intersect(which(RangeCols3To5[1,] >=  -threshhold),
                        which(RangeCols3To5[2,] <= threshhold))
    return(EscoresInd)
}

#' Removes header from a data.table
#' 
#' \code{removeHeader} takes the proper number of rows for dt (numRowsNoHeader)
#' and if dt has one extra row, it removes the first row (header) from dt.
#' If dt does not have numRowsNoHeader rows, or numRowsNoHeader + 1 rows,
#' this function returns an error.
removeHeader <- function(numRowsNoHeader, dt) {
  if(nrow(dt) == (numRowsNoHeader + 1)) {
    dt = dt[-1,]
  }
  if(nrow(dt) != numRowsNoHeader) {
    stop("Error, incorrect file length given in removeHeader")
  }
  return(dt)
}

#'
#'
#' \code{CheckEscoresCol}
#'
#' called in construct_summ, but commented out
CheckEscoresCol <- function(file.name) {

  scoresDT <- fread(file.name, header=FALSE)
  scoresDT <- removeHeader(filelength, scoresDT)

  Cols3To5 <- getCols3to5(scoresDT)
  RangeCols3To5 <- getRangeCols(Cols3To5)
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

#'
#'
#'\code{CreateEscores} is used to 
#'
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

