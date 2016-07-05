## removing all existing variables from the work space
rm(list=ls())
## importing libraries
library(Biostrings)
library(data.table)
library(gtools)

## this script takes a list of DNA sequences
## returns two files - a list of protein names and associated filenames and 
## a list of scores for each subsequence in the proteins
datadir <- "/p/keles/CAGI2015/volumeB/Data/UniProbe/Human"
all.result.files <- mixedsort(list.files(path=datadir, full.names=TRUE))

print("Directory read in")

load("/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData")
#### THIS IS WHAT ACTUALLY CHANGES IN THE PIPELINE
### change to IO location
#load("SampleDT.RData")

## Loading the functions CheckEscoresCol and CreateEscores
#source("/p/keles/CAGI2015/volumeB/Exploratory/8merScores.R")
source("8merScoresMODIFIED.R")

print("8merScores ran")

## looking at 8mers in the sequences
k <- 8
## list of lengths of all the reference sequences
## trainingDT needs to already be loaded into the environment
lenCenter <- mapply(nchar, trainingDT$RefSeq)

## returns a list of all 8mer subsequences in each RefSeq
subseqList <- mapply(function(seq, len)  mapply(function(s) substr(seq, s, s+k-1), 1:(len-k+1)),
                   trainingDT$RefSeq, lenCenter)
## labeling the subsequences with the eQTL positions?
names(subseqList) <- trainingDT$ID

## I don't know why this is commented out. It was like this
#escoresCol<-mapply(CheckEscoresCol, all.result.files)
#table(escoresCol)

## turns all.result.files into a matrix then runs CreateEscores on each filename
## not sure why they're defining a function instead of doing a fn call? maybe you need to?
#enrscoresList <- mclapply(as.matrix(all.result.files), function(x) CreateEscores(file.name=x), mc.cores=20)
## unparalizing it for this old-ass version of R
enrscoresList<-lapply(as.matrix(all.result.files), function(x) CreateEscores(file.name=x))

print("CreateEscores done running")

## so it's splitting it into the File name and the Protein name for each file
fname <- mapply(function(file.name) strsplit(file.name, "/")[[1]][9], all.result.files)
## generating 1:length(fname). i dunno why
names(fname) <- seq(length(fname))                       

#pname<-mapply(function(fname) strsplit(fname, "_")[[1]][1], fname)
### NEW EXPRESSION FOR UNIQUE NAMES
pname = sub("*(_8mers|_contig8mers).*", "", names)
## NEW check protein names are unique
if (anyDuplicated(pname)!=0) {
  print("Protein names repeated")
}

NameM <- cbind(fname, pname)
colnames(NameM) <- c("File", "Protein")                       


#save(enrscoresList, NameM, file="/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_human_v2.Rda")
save(enrscoresList, NameM, file="escores8mer_ref_human_v2.Rda")

