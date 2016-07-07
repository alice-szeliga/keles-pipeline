## fimo threshold
fimoThreshold <- 0.001

##run this on REF and on ALT

#' Creates directory name for factor
#'
#' \code{getFactorName} creates a directory name for a directory within
#'   fimo to save to. 
#'  
#' Relies on previously defined fimoThreshold
#' Ex: "Training_ref_thres0.001"
#' @param refOrAlt: string of "ref" or "alt"
getFactorName <- function(refOrAlt) {
  return(paste("Training_", refOrAlt, "_thres", fimoThreshold,   sep = ""))
}

#' Creates filename for FASTA file
#'
#' \code{getFactorName} creates filename
#'  
#' Ex: "ref_SS.v1.fasta"
#' @param refOrAlt: string of "ref" or "alt"
getFastaName <- function (refOrAlt) {
  return(paste(refOrAlt, "_SS.v1.fasta", sep = ""))
}

#' Run Fimo
#' 
#' Sends a system command to run FIMO
#' 
#' @param memeFimoDir: string, path of fimo directory for MEME
#' @param outDir: string, location to save to
#' @param pwmFiles: list of strings, locations of PWM for ENCODE, Factorbook, and Jaspar
#' @param factorName: string, created by getFactorName
#' @param fastaFile: string, created by getFastaName
runFimo <- function(memeFimoDir = "/p/keles/MEME/meme_4.9.0/src/fimo",
                    outDir,
                    pwmFiles = c("/p/keles/EnhancerPred/volumeD/PWM-Scans/PWM/encode_motifs_for_fimo.txt", 
                                 "/p/keles/EnhancerPred/volumeD/PWM-Scans/PWM/factorbook_PWMs.txt",
                                 "/p/keles/CAGI2015/volumeA/Jaspar_Motifs/jaspar_motifs_for_fimo.txt"),
                    factorName,
                    fastaFile) {
  dir.fimo <- paste(outDir, factorName, "/",   sep = "")

  pwmNames <- c("ENCODE", "FactorBook", "Jaspar")
  
  #previousDir <- getwd()  ##potentially reset directory afterwards
  setwd(outDir)
  for (k in 1:length(pwmNames)){
    PWM <- pwmFiles[k]
    dir.out.PWMSpecific <- paste(dir.fimo, pwmNames[k], "/",  sep = "")
    fimo.command <- paste(memeFimoDIr, " --verbosity 1 --thresh ", fimoThreshold, " -oc ", 
                          dir.out.PWMSpecific, "  ", PWM, "  ", fastaFile, "&",  sep = "")
    system(fimo.command)
  }
  #setwd(previousDir)
}

#' Create TF output
#' 
#' 
#' 
#' @param TF.names
#' @param d1
#' @return
createOutput <- function(TF.names, d1) {
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
  
  return(list(out.1, out.2, out.3, out.4))
}

#' Saves output to directory
#'
#'
#'
#'
saveOutput <- function(functionName, factorName, out.1, out.2, out.3, out.4) {
  cat(factorName, " ", dim(out.1), dim(out.2), dim(out.3), dim(out.4),"\n")
  
  write.table(out.1, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, 
              file = paste(outDir, factorName,  "_", functionName, "_minPval.TXT", sep = ""))
  save(out.1, 
       file = paste(outDir, factorName,  "_", functionName, "_minPval.RData", sep = ""))
  
  write.table(out.2, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, 
              file = paste(outDir, factorName,  "_", functionName, "_maxScore.TXT", sep = ""))
  save(out.2, 
       file = paste(outDir, factorName,  "_", functionName, "_maxScore.RData", sep = ""))
  
  write.table(out.3, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, 
              file = paste(outDir, factorName,  "_", functionName, "_numOcc.TXT", sep = ""))
  save(out.3, 
       file = paste(outDir, factorName,  "_", functionName, "_numOcc.RData", sep = ""))
  
  write.table(out.4, quote = F, sep = "\t", col.names = TRUE, row.names = TRUE, 
              file = paste(outDir, factorName,  "_", functionName, "_sumScore.TXT", sep = ""))
  save(out.4, 
       file = paste(outDir, factorName,  "_", functionName, "_sumScore.RData", sep = ""))
  
  
}

#' Summarize FIMO with FactorBook
#' 
#' 
#' 
#' dslkj
summarizeFactorBook <- function(factorBookPWM) {
  PWM <- scan(factorBookPWM, what = "", sep = "\n")
  TF.names <- apply(as.matrix(PWM[grep("MOTIF", PWM)]), 1, function(x){strsplit(x, "MOTIF ")[[1]][2]})
  fimo.dir <- paste(outDir, factorName, "/FactorBook", "/",   sep = "")
  
  setwd(fimo.dir)
  d <- read.table("fimo.txt")
  d1 <- split(d,  d$V2)
  
  output <- createOutput(TF.names, d1)
  out.1 <- output[[1]]; out.2 <- output[[2]]
  out.3 <- output[[3]]; out.4 <- output[[4]]
  
  ### this bit might be redundant
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
  
  saveOutput("FactorBook", factorName, out.1, out.2, out.3, out.4)
  
}

#' Summarize FIMO with ENCODE
#' 
#' 
#' 
#' 
#' lkj
summarizeEncode <- function(encodePWM) {
  PWM <- scan(encodePWM, what = "", sep = "\n")
  TF.names <- apply(as.matrix(PWM[grep("MOTIF", PWM)]), 1, function(x){strsplit(x, "MOTIF ")[[1]][2]})
  fimo.dir <- paste(outDir, factorName, "/ENCODE", "/",   sep = "")
  
  setwd(fimo.dir)
  d <- read.csv("fimo.txt", sep = "\t")
  d1 <- split(d,  d$sequence.name)
  
  
  output <- createOutput(TF.names, d1)
  out.1 <- output[[1]]; out.2 <- output[[2]]
  out.3 <- output[[3]]; out.4 <- output[[4]]
  
  
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
  
  saveOutput("ENCODE", factorName, out.1, out.2, out.3, out.4)
  
}

#' Summarize FIMO with Jaspar
#' 
#' d
#'
#' ddf
summarizeJaspar <- function(jasparPWM) {
  PWM <- scan(jasparPWM, what = "", sep = "\n")
  TF.names <- apply(as.matrix(PWM[grep("MOTIF", PWM)]), 1, function(x){strsplit(x, "MOTIF ")[[1]][2]})
  fimo.dir <- paste(outDir, factorName, "/Jaspar", "/",   sep = "")
  
  setwd(fimo.dir)
  d <- read.table("fimo.txt")
  d1 <- split(d,  d$V2)
  
  output <- createOutput(TF.names, d1)
  out.1 <- output[[1]]; out.2 <- output[[2]]
  out.3 <- output[[3]]; out.4 <- output[[4]]
  
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
  
  saveOutput("Jaspar", factorName, out.1, out.2, out.3, out.4)
  
}

#' Runs FactorBook, ENCODE and Jaspar on FIMO
#' 
#' 
#'
#' @param refOrAlt: string, "ref" or "alt"
#' @param memeFimoDir: string, 
#' @param outDir:
#' @param pwmFiles: list of strings, locations of PWMs for ENCODE, Factorbook,
#'   and Jaspar
run_fimo_and_summarize <- function () {
  factorName <- getFactorName(refOrAlt)
  fastaFile <- getFastaName(refOrAlt)
  runFimo()
  factorBookPWM <- pwmFiles[1]; encodePWM <- pwmFiles[2]; jasparPWM <- pwmFiles[3]
  summarizeFactorBook(factorBookPWM)
  summarizeEncode(encodePWM)
  summarizeJaspar(jasparPWM)
}