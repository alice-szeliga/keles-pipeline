#############################
rm(list=ls())
library(data.table)

importHumanMouse <- function()
load("/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_human_v2.Rda")
enrscoresList_human <- enrscoresList
NameM_human <- NameM

load("/p/keles/CAGI2015/volumeB/ProcessedData/escores8mer_ref_mouse_v2.Rda")
enrscoresList_mouse <- enrscoresList
NameM_mouse <- NameM

enrscoresList <- c(enrscoresList_human, enrscoresList_mouse)
NameM <- rbind(NameM_human, NameM_mouse)


#####Add unique protein names#####

ProteinNames <- NameM[,2]
names(enrscoresList) <- ProteinNames
##############

SumEscores <- mapply(function(enrList) mapply(sum, enrList), enrscoresList)
MeanEscores <- mapply(function(enrList) mapply(mean, enrList), enrscoresList)
MaxEscores <- mapply(function(enrList) mapply(max, enrList), enrscoresList)


## change to IO?
## or also being called in other. source both into one script and load this once?
load("/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData")

lenCenter <- mapply(nchar, trainingDT$RefAllele)
table(trainingDT$pos - trainingDT$start + 1)
table(lenCenter)

GetSubseqCoordAbs <- function(subseq.len=21, center.len=lenCenter, center.pos=trainingDT$pos) {
  heading <- ceiling((subseq.len - center.len)/2)
  foot <- floor((subseq.len - center.len)/2)
  start <- center.pos - heading
  end <- center.pos + center.len-1 + foot
  subseq.coord <- cbind(start, end)
  subseq.coord
}

subseq.coord.abs <- GetSubseqCoordAbs(subseq.len=21, center.len=lenCenter, center.pos=trainingDT$pos)
subseq.coord.rel <- subseq.coord.abs - matrix(trainingDT$start, nrow(trainingDT), 2) + 1
subseq.coord.rel.list <- lapply(seq_len(nrow(subseq.coord.rel)), function(i) subseq.coord.rel[i,])

enrscoresSubseq21List <- mapply(function(enrscores) t(mapply(function(enrscore, rel.pos) enrscore[rel.pos[1]:rel.pos[2]], enrscores, subseq.coord.rel.list)), enrscoresList, SIMPLIFY=FALSE)

SumSubseq21Escores <- mapply(rowSums, enrscoresSubseq21List)
MeanSubseq21Escores <- mapply(rowMeans, enrscoresSubseq21List)
MaxSubseq21Escores <- mapply(function(enrM) apply(enrM, 1, max), enrscoresSubseq21List)

escoresList <- vector("list")
escoresList[[1]] <- SumEscores
escoresList[[2]] <- MeanEscores
escoresList[[3]] <- MaxEscores
escoresList[[4]] <- SumSubseq21Escores
escoresList[[5]] <- MeanSubseq21Escores
escoresList[[6]] <- MaxSubseq21Escores
names(escoresList) <- c("sum", "mean", "max", "sum_bp21", "mean_bp21", "max_bp21")

## the final output
save(escoresList, file="/p/keles/CAGI2015/volumeB/ProcessedData/Features8mersEscores_ref_v2.Rda")

