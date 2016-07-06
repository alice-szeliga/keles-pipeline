rm(list=ls())
library(Biostrings)
library(data.table)
library(gtools)

extract_8merscores <- function() {
datadir<-"/p/keles/CAGI2015/volumeB/Data/UniProbe/Human"
all.result.files<-mixedsort(list.files(path=datadir, full.names=TRUE))
load("/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData")
source("/p/keles/CAGI2015/volumeB/Exploratory/8merScores.R")

k<-8
lenCenter<-mapply(nchar, trainingDT$RefSeq)
subseqList<-mapply(function(seq, len)  mapply(function(s) substr(seq, s, s+k-1), 1:(len-k+1)), trainingDT$RefSeq, lenCenter)
names(subseqList)<-trainingDT$ID

#escoresCol<-mapply(CheckEscoresCol, all.result.files)
#table(escoresCol)

enrscoresList<-mclapply(as.matrix(all.result.files), function(x) CreateEscores(file.name=x), mc.cores=20)

fname<-mapply(function(file.name) strsplit(file.name, "/")[[1]][9], all.result.files)
names(fname)<-seq(length(fname))
pname<-mapply(function(fname) strsplit(fname, "_")[[1]][1], fname)
NameM<-cbind(fname, pname)
colnames(NameM)<-c("File", "Protein")

save(enrscoresList, NameM, file="/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_human_v2.Rda")
}
