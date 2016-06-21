## removing all existing variables from the work space
### remove this in package
rm(list=ls())
## importing libraries
library(Biostrings)
library(data.table)
library(gtools)


### can remove these in package
## Loading the functions CheckEscoresCol and CreateEscores
#source("/p/keles/CAGI2015/volumeB/Exploratory/8merScores.R")
#source("8merScoresMODIFIED.R")

extract8mers <- function(datadir, sampleFilename, outputLocation) {

  ### fixed location - don't have to change
  ### input from same location every time
  #datadir<-"/p/keles/CAGI2015/volumeB/Data/UniProbe/Human"
  all.result.files<-mixedsort(list.files(path=datadir, full.names=TRUE))
  #all.result.files <- c("Sox4_8mers_v5.txt")


  #load("/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData")
  load("sampleFilename")


  ## looking at 8mers in the sequences
  k <- 8
  ## list of lengths of all the reference sequences
  ## trainingDT needs to already be loaded into the environment
  lenCenter<-mapply(nchar, trainingDT$RefSeq)

  ## returns a list of all 8mer subsequences in each RefSeq
  subseqList <- mapply(function(seq, len)  mapply(function(s) substr(seq, s, s+k-1), 1:(len-k+1)),
                     trainingDT$RefSeq, lenCenter)
  ## labeling the subsequences with the eQTL positions?
  names(subseqList) <- trainingDT$ID


  ## This was already commented out
  #escoresCol<-mapply(CheckEscoresCol, all.result.files)
  #table(escoresCol)

  enrscoresList <- mclapply(as.matrix(all.result.files),
                          function(x) CreateEscores(file.name=x),
                          mc.cores=20)

  ## so it's splitting it into the File name and the Protein name for each file
  fname <- mapply(function(file.name) strsplit(file.name, "/")[[1]][9],
                  all.result.files)
  ## generating 1:length(fname). i dunno why
  names(fname) <- seq(length(fname))

  ## old expression for determining protein names
  #pname<-mapply(function(fname) strsplit(fname, "_")[[1]][1], fname)
  ### NEW EXPRESSION FOR UNIQUE NAMES
  pname = sub("*(_8mers|_contig8mers).*", "", names)
  ## NEW check protein names are unique
  if (anyDuplicated(pname)!=0) {
    print("Error: protein names repeated in extract8mers")
  }

  NameM <- cbind(fname, pname)
  colnames(NameM) <- c("File", "Protein")

  #save(enrscoresList, NameM, file="/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_human_v2.Rda")
  save(enrscoresList, NameM, file="escores8mer_ref_human_v2.Rda")
}
