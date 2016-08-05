dataDir <- "/p/keles/CAGI2015/volumeB/Data/After_Nov5/"
outDir <- "/p/keles/CAGI2015/volumeB/Fimo/"
factorName <- "Training_alt_thres0.001"
fastaFile <- "alt_SS.v1.fasta"

## Changed fimo threshold from 0.0001 to 0.001

#####################################
#                                   #
# Run FIMO:REF                      #
#                                   #
#####################################


dir.fimo <- paste(outDir, factorName, "/",   sep = "")
pwmFiles <- c("/p/keles/EnhancerPred/volumeD/PWM-Scans/PWM/encode_motifs_for_fimo.txt", "/p/keles/EnhancerPred/volumeD/PWM-Scans/PWM/factorbook_PWMs.txt",
"/p/keles/CAGI2015/volumeA/Jaspar_Motifs/jaspar_motifs_for_fimo.txt")

pwmNames <- c("ENCODE", "FactorBook", "Jaspar")

setwd(outDir)
for (k in 1:length(pwmNames)){
    PWM <- pwmFiles[k]
    dir.out.PWMSpecific <- paste(dir.fimo, pwmNames[k], "/",  sep = "")
    fimo.command <- paste("/p/keles/MEME/meme_4.9.0/src/fimo --verbosity 1 --thresh 0.001 -oc ", dir.out.PWMSpecific, "  ", PWM, "  ", fastaFile, "&",  sep = "")
    system(fimo.command)
}

#####################################
#                                   #
# Summarize FIMO w FactorBook: REF  #
#                                   #
#####################################

## FactorBook
PWM <- scan("/p/keles/EnhancerPred/volumeD/PWM-Scans/PWM/factorbook_PWMs.txt", what = "", sep = "\n")
TF.names <- apply(as.matrix(PWM[grep("MOTIF", PWM)]), 1, function(x){strsplit(x, "MOTIF ")[[1]][2]})
fimo.dir <- paste(outDir, factorName, "/FactorBook", "/",   sep = "")

setwd(fimo.dir)
d <- read.table("fimo.txt")
d1 <- split(d,  d$V2)
    
out.1 <- matrix(1, length(d1), length(TF.names))
colnames(out.1) <- TF.names
rownames(out.1) <- names(d1)

out.2 <- matrix(0, length(d1), length(TF.names))
colnames(out.2) <- TF.names
rownames(out.2) <- names(d1)

out.3 <- matrix(0, length(d1), length(TF.names))
colnames(out.3) <- TF.names
rownames(out.3) <- names(d1)

out.4 <- matrix(0, length(d1), length(TF.names))
colnames(out.4) <- TF.names
rownames(out.4) <- names(d1)

for (k in 1:length(d1)){
    d2 <- d1[[k]]
    d3 <- split(d2, d2$V1)
    d4 <- lapply(d3, function(x){if (nrow(x) >0){return(min(x$V7))} else return(1)}) #return p-val
    d5 <- lapply(d3, function(x){if (nrow(x) >0){return(max(x$V6))} else return(0)}) #return score
    d6 <- lapply(d3, function(x){if (nrow(x) >0){return(length(x$V7))} else return(0)}) #return number of occurrences that make the p-val cut off
    d7 <- lapply(d3, function(x){if (nrow(x) >0){return(sum(x$V6))} else return(0)}) #return sum of the scores
    
    out.1[names(d1)[k], match(names(d4), TF.names)] <- unlist(d4) #matching names of the TFs so that we have correct TFs in each columns
    out.2[names(d1)[k], match(names(d5), TF.names)] <- unlist(d5) #matching names of the TFs so that we have correct TFs in each columns
    out.3[names(d1)[k], match(names(d6), TF.names)] <- unlist(d6) #matching names of the TFs so that we have correct TFs in each columns
    out.4[names(d1)[k], match(names(d7), TF.names)] <- unlist(d7) #matching names of the TFs so that we have correct TFs in each columns

}
cat(factorName, " ", dim(out.1), dim(out.2), dim(out.3), dim(out.4),"\n")

write.table(out.1, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "FactorBook", "_minPval.TXT", sep = ""))
save(out.1, file = paste(outDir, factorName,  "_", "FactorBook", "_minPval.RData", sep = ""))

write.table(out.2, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "FactorBook", "_maxScore.TXT", sep = ""))
save(out.2, file = paste(outDir, factorName,  "_", "FactorBook", "_maxScore.RData", sep = ""))

write.table(out.3, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "FactorBook", "_numOcc.TXT", sep = ""))
save(out.3, file = paste(outDir, factorName,  "_", "FactorBook", "_numOcc.RData", sep = ""))

write.table(out.4, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "FactorBook", "_sumScore.TXT", sep = ""))
save(out.4, file = paste(outDir, factorName,  "_", "FactorBook", "_sumScore.RData", sep = ""))

#####################################
#                                   #
# Summarize FIMO w ENCODE: REF      #
#                                   #
#####################################
PWM <- scan("/p/keles/EnhancerPred/volumeD/PWM-Scans/PWM/encode_motifs_for_fimo.txt", what = "", sep = "\n")
TF.names <- apply(as.matrix(PWM[grep("MOTIF", PWM)]), 1, function(x){strsplit(x, "MOTIF ")[[1]][2]})
fimo.dir <- paste(outDir, factorName, "/ENCODE", "/",   sep = "")

setwd(fimo.dir)
d <- read.csv("fimo.txt", sep = "\t")
d1 <- split(d,  d$sequence.name)
    
out.1 <- matrix(1, length(d1), length(TF.names))
colnames(out.1) <- TF.names
rownames(out.1) <- names(d1)

out.2 <- matrix(0, length(d1), length(TF.names))
colnames(out.2) <- TF.names
rownames(out.2) <- names(d1)

out.3 <- matrix(0, length(d1), length(TF.names))
colnames(out.3) <- TF.names
rownames(out.3) <- names(d1)

out.4 <- matrix(0, length(d1), length(TF.names))
colnames(out.4) <- TF.names
rownames(out.4) <- names(d1)


for (k in 1:length(d1)){
    d2 <- d1[[k]]
    d3 <- split(d2, d2$X.pattern.name)
    d4 <- lapply(d3, function(x){if (nrow(x) >0){return(min(x$p.value))} else return(1)}) #return p-val
    d5 <- lapply(d3, function(x){if (nrow(x) >0){return(max(x$score))} else return(0)}) #return score
    d6 <- lapply(d3, function(x){if (nrow(x) >0){return(length(x$p.value))} else return(0)}) #return score
    d7 <- lapply(d3, function(x){if (nrow(x) >0){return(sum(x$score))} else return(0)}) #return score
    
    out.1[names(d1)[k], match(names(d4), TF.names)] <- unlist(d4) #matching names of the TFs so that we have correct TFs in each columns
    out.2[names(d1)[k], match(names(d5), TF.names)] <- unlist(d5) #matching names of the TFs so that we have correct TFs in each columns
    out.3[names(d1)[k], match(names(d6), TF.names)] <- unlist(d6) #matching names of the TFs so that we have correct TFs in each columns
    out.4[names(d1)[k], match(names(d7), TF.names)] <- unlist(d7) #matching names of the TFs so that we have correct TFs in each columns
}



cat(factorName, " ", dim(out.1), dim(out.2), dim(out.3), dim(out.4),"\n")

write.table(out.1, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "ENCODE", "_minPval.TXT", sep = ""))
save(out.1,  file = paste(outDir, factorName,  "_", "ENCODE", "_minPval.RData", sep = ""))

write.table(out.2, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "ENCODE", "_maxScore.TXT", sep = ""))
save(out.2,  file = paste(outDir, factorName,  "_", "ENCODE", "_maxScore.RData", sep = ""))

write.table(out.3, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "ENCODE", "_numOcc.TXT", sep = ""))
save(out.3, file = paste(outDir, factorName,  "_", "ENCODE", "_numOcc.RData", sep = ""))

write.table(out.4, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "ENCODE", "_sumScore.TXT", sep = ""))
save(out.4, file = paste(outDir, factorName,  "_", "ENCODE", "_sumScore.RData", sep = ""))


#####################################
#                                   #
# Summarize FIMO w Jaspar: REF      #
#                                   #
#####################################

## Jaspar
PWM <- scan("/p/keles/CAGI2015/volumeA/Jaspar_Motifs/jaspar_motifs_for_fimo.txt", what = "", sep = "\n")
TF.names <- apply(as.matrix(PWM[grep("MOTIF", PWM)]), 1, function(x){strsplit(x, "MOTIF ")[[1]][2]})
fimo.dir <- paste(outDir, factorName, "/Jaspar", "/",   sep = "")

setwd(fimo.dir)
d <- read.table("fimo.txt")
d1 <- split(d,  d$V2)

out.1 <- matrix(1, length(d1), length(TF.names))
colnames(out.1) <- TF.names
rownames(out.1) <- names(d1)

out.2 <- matrix(0, length(d1), length(TF.names))
colnames(out.2) <- TF.names
rownames(out.2) <- names(d1)

out.3 <- matrix(0, length(d1), length(TF.names))
colnames(out.3) <- TF.names
rownames(out.3) <- names(d1)

out.4 <- matrix(0, length(d1), length(TF.names))
colnames(out.4) <- TF.names
rownames(out.4) <- names(d1)


for (k in 1:length(d1)){
    d2 <- d1[[k]]
    d3 <- split(d2, d2$V1)
    d4 <- lapply(d3, function(x){if (nrow(x) >0){return(min(x$V7))} else return(1)}) #return p-val
    d5 <- lapply(d3, function(x){if (nrow(x) >0){return(max(x$V6))} else return(0)}) #return score
    d6 <- lapply(d3, function(x){if (nrow(x) >0){return(length(x$V7))} else return(0)}) #return number of occurrences that make the p-val cut off
    d7 <- lapply(d3, function(x){if (nrow(x) >0){return(sum(x$V6))} else return(0)}) #return sum of the scores
    
    out.1[names(d1)[k], match(names(d4), TF.names)] <- unlist(d4) #matching names of the TFs so that we have correct TFs in each columns
    out.2[names(d1)[k], match(names(d5), TF.names)] <- unlist(d5) #matching names of the TFs so that we have correct TFs in each columns
    out.3[names(d1)[k], match(names(d6), TF.names)] <- unlist(d6) #matching names of the TFs so that we have correct TFs in each columns
    out.4[names(d1)[k], match(names(d7), TF.names)] <- unlist(d7) #matching names of the TFs so that we have correct TFs in each columns

 }


cat(factorName, " ", dim(out.1), dim(out.2), dim(out.3), dim(out.4),"\n")

write.table(out.1, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "Jaspar", "_minPval.TXT", sep = ""))
save(out.1, file = paste(outDir, factorName,  "_", "Jaspar", "_minPval.RData", sep = ""))

write.table(out.2, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "Jaspar", "_maxScore.TXT", sep = ""))
save(out.2, file = paste(outDir, factorName,  "_", "Jaspar", "_maxScore.RData", sep = ""))

write.table(out.3, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "Jaspar", "_numOcc.TXT", sep = ""))
save(out.3, file = paste(outDir, factorName,  "_", "Jaspar", "_numOcc.RData", sep = ""))

write.table(out.4, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, file = paste(outDir, factorName,  "_", "Jaspar", "_sumScore.TXT", sep = ""))
save(out.4, file = paste(outDir, factorName,  "_", "Jaspar", "_sumScore.RData", sep = ""))




