CheckEscoresCol<-function(file.name) {
scoresDT<-fread(file.name, header=TRUE)
if(nrow(scoresDT)==32895) {
  scoresDT<-fread(file.name, header=FALSE)
} else if(nrow(scoresDT)==32896)
  scoresDT<-fread(file.name, header=TRUE)
  
Cols3To5<-scoresDT[,setdiff(seq(min(5, ncol(scoresDT))), 1:2), with=FALSE]
RangeCols3To5<-apply(Cols3To5, 2, range)
EscoresInd<-intersect(which(RangeCols3To5[1,]>=-0.5), which(RangeCols3To5[2,]<=0.5))
EscoresInd
}

CreateEscores<-function(file.name) {
scoresDT<-fread(file.name, header=TRUE)
if(nrow(scoresDT)==32895) {
  scoresDT<-fread(file.name, header=FALSE)
} else if(nrow(scoresDT)==32896)
  scoresDT<-fread(file.name, header=TRUE)
  
mersSeq1<-unlist(scoresDT[,1, with=FALSE])
mersSeq2<-unlist(scoresDT[,2, with=FALSE])
mersSeq<-cbind(mersSeq1, mersSeq2)

Cols3To5<-scoresDT[,setdiff(seq(min(5, ncol(scoresDT))), 1:2), with=FALSE]
RangeCols3To5<-apply(Cols3To5, 2, range)
EscoresInd<-intersect(which(RangeCols3To5[1,]>=-0.5), which(RangeCols3To5[2,]<=0.5))
Escores<-Cols3To5[,EscoresInd,with=FALSE]

matchInd<-mapply(function(subseqs) unlist(apply(as.matrix(subseqs), 1, function(subseq) which(mersSeq==subseq, arr.ind=TRUE)[1,1])), subseqList, SIMPLIFY=FALSE)
EscoresList<-mapply(function(x) Escores[x], matchInd)
EscoresList
}


#load("/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_human.Rda")
#sum(abs(enrscoresList$FOXN2-Escores[matchInd[[1]]]))
