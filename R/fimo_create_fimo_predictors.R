#' Creates a score list
#'
#' \code{makeScoreList} creates and gives column names to a list which will
#'   contain the minPval, maxScore, numOcc, and sumScore for Factorbook, 
#'   Jaspar, and Encode scores for the given eQTLs, as well as the Deepbind
#'   score.
#'
#' @return scoreList: a list which will contain the FactorBook, Jaspar, Encode,
#'   and Deepbind scores.
makeScoreList <- function() {
  scoreList <- vector("list", 13)
  names(scoreList) <- c("fb_minPval", "fb_maxScore", "fb_numOcc", "fb_sumScore",
                        "en_minPval", "en_maxScore", "en_numOcc", "en_sumScore",
                        "ja_minPval", "ja_maxScore", "ja_numOcc", "ja_sumScore", 
                        "db_score")
  return(scoreList)
}

#' Load in FIMO summary data
#'
#' \code{loadFimoData} loads in the .Rdata files created by 
#'   fimo_run_fimo_and_summarize.R
#'
#' @param refOrAlt
#' @param functionName
#' @param fimoThreshold
#' @return Null. Loads 4 .Rdata files into the environment
loadFimoData <- function(refOrAlt, functionName, fimoThreshold) {
  for (name in c("minPval", "maxScore", "sumScore", "numOcc")) {
    load(paste("Training_", refOrAlt, "_thres", fimoThreshold, 
               "_", functionName, "_", name, ".RData", sep = ""))
  }
}

#' Adds a score to the score list
#'
#' \code{addToScoreList} adds a FactorBook or Jaspar score to the score list.
#'   Specifically, minPval, maxScore, numOcc, or sumScore.
#'
#' @param trainingD: a table of eQTL causal SNPs
#' @param fimoS: a score for FIMO
#' @param scoreName: string, name indicating analysis type and type of score
#'   Ex. "fb_minPval"
#' @param scoreList: a vector containing all scores for FactorBook, ENCODE,
#'   Jaspar, and Deepbhind 
addToScoreList <- function(trainingD, fimoS, scoreName, scoreList) {
  # split the fimo score rownames into ??
  oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
  fimoS.ordered <- matrix(1, nrow(trainingD), ncol(fimoS))
  colnames(fimoS.ordered) <- colnames(fimoS)
  fimoS.ordered[oIndex, ] <- fimoS
  # add the ordered fimo score to the scoreName's column of the list of scores
  scoreList[[scoreName]] <- fimoS.ordered

  return(scoreList)
}


#' Adds a score to the score list
#'
#' \code{addAllToScoreList} adds the minPval, maxScore, numOcc, and sumScore to
#'   the scoreslist for either FactBook or Jaspar
#'
#' @param trainingD: a table of eQTL causal SNPs
#' @param fimoS: TYPE?, a score for FIMO
#' @param scoreName: list of strings, which are indicating analysis type and type of score
#'   Ex. "fb_minPval"
#' @param scoreList: a vector containing all scores for FactorBook, ENCODE,
#'   Jaspar, and Deepbhind 

addAllToScoreList <- function (trainingD, fimoS, scoreNames, scoreList) {
  # renaming the data which has been loaded in
  fimo.minPval <- out.1
  fimo.maxScore <- out.2
  fimo.numOcc <- out.3
  fimo.sumScore <- out.4
  
  # creating names to label each score with - specific to this analysis
  for (score in scoreNames) {
    scoreList <- addToScoreList(trainingD, fimoS, scoreName, scoreList)
  }
  return(scoreList)
}

#'
#'
#' \code{getCorrelation} performs a simple correlation between a fimo score and
#'    the ????? from the eQTL data given.
#'
#' @param fimoS.ordered, ???
#' @param trainingD: a table of eQTL data
getCorrelation <- function(fimoS.ordered, trainingD) {
  cor(fimoS.ordered, trainingD$C.A.log2FC)
}

#' Perform simple correlations
#' 
#' \code{createAllCorrelations} will take the scores for either FactBook,
#'   Jaspar, or ENCODE and will perform correlations for the minimum P-value,
#'   the maximum score, the sum of the scores, and the number of occurences.
#' 
#' @param scoreNames: list of strings, contains names of the analysis and the
#'    different scores. 
#'    Ex. c("fb_minPval", "fb_maxScore", "fb_numOcc", "fb_sumScore")
#' @param scoreList:
#' @return Null. Prints out.
createAllCorrelations <- function(scoreList, scoreNames) {
  # Getting correlation for the minimum p-value
  cat("minPVal\n")
  fimoS.ordered <- scoreList[[scoreNames[1]]]
  fimoS.ordered.binary <- apply(fimoS.ordered, 2, function(x){res <- x; res[x>=0.05]<-0; res[x <0.05] <-1; return(res)}) #PWM match crude
  corRes.binary <- apply(fimoS.ordered.binary, 2, function(x){ifelse(length(unique(x)) > 1, wilcox.test(trainingD$C.A.log2FC~x)$p.value,1)})
  #order(corRes.binary, decreasing = FALSE)
  summary(-log10(corRes.binary))
  
  # Getting correlation for the maximum score
  cat("maxScore\n")
  fimoS.ordered <- scoreList[[scoreNames[2]]]
  corRes <- apply(fimoS.ordered, 2, getCorrelation(fimoS.ordered, trainingD))
  #order(abs(corRes), decreasing = TRUE)
  summary(abs(corRes))
  
  # Getting correlation for the sum of the scores
  cat("sumScore\n")
  fimoS.ordered <- scoreList[[scoreNames[3]]]
  corRes <- apply(fimoS.ordered, 2, getCorrelation(fimoS.ordered, trainingD))
  #order(abs(corRes), decreasing = TRUE)
  #sort(abs(corRes), decreasing = TRUE)[1:10]
  summary(abs(corRes))
  
  # Getting correlation for the number of occurences
  cat("numOcc\n")
  fimoS.ordered <- scoreList[[scoreNames[4]]]
  corRes <- apply(fimoS.ordered, 2, getCorrelation(fimoS.ordered, trainingD))
  #order(abs(corRes), decreasing = TRUE)
  #sort(abs(corRes), decreasing = TRUE)[1:10]
  summary(abs(corRes))
}

#' Exploratory analysis with FactorBook PWM scores
#'
#' \code{factorBookAnalysis} 
#'
#' @param dataDir
#' @param fimoDir
#' @param refOrAlt
#' @param trainingD
#' @param fimoS
#' @param scoreList
#' @return 
factorBookAnalysis <- function(dataDir, fimoDir) {
  # First load in the eQTL data
  setwd("/p/keles/CAGI2015/volumeB/Data/")
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
 
  # Load in the FIMO data from previous script
  #setwd("/Users/keles/Research_Projects/CAGI2015/Fimo/")
  setwd(fimoDir)
  loadFimoData(refOrAlt, "FactorBook", fimoThreshold)
  
  # Create a specific list of names for the 4 scores 
  # Ex. fb_minPval, etc.
  scoreNames <- paste("fb", scoreTypes, sep = "_")
  
  # Add all 4 scores to the output list
  scoreList <- addAllToScoreList(trainingD, fimoS, scoreNames, scoreList)
  
  # Create correlations for all 4 scores
  createAllCorrelations(scoreList, scoreNames)
  
}

#' Exploratory analysis with Deepbind scores
#'
#' \code{deepbindAnalysis} 
#'
#' @param 
#' @return 
deepbindAnalysis <- function() {
  setwd("/p/keles/CAGI2015/volumeB/DeepBind/deepbind")
  db <- read.table("alt_SS.v1_deepbind.txt", header = TRUE)
  
  scoreList[["db_score"]] <- db
  
  ## Plot the correlations for the Deepbind data?
  setwd("/Users/keles/Research_Projects/CAGI2015/Data/After_Nov5")
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
  setwd("/Users/keles/Research_Projects/CAGI2015/DeepBind/")
  db <- read.table("Training_alt_3044_deepbind_927.txt", header = TRUE)
  corRes <- apply(db, 2, getCorrelation(fimoS.ordered, trainingD))
  index <- order(abs(corRes), decreasing = TRUE)[1:50]
 
  db.sub <- db[, index]
  colnames(db.sub) <- paste("db_", 1:50, sep = "")
  
  dd <- cbind(trainingD, db.sub)
  ddd <- dd
  ddd$grp <- paste(ddd$Regulatory_Hit, ddd$emVar_Hit, sep = "")
  ggplot(ddd, aes(x = db_2, y =C.A.log2FC)) + geom_smooth(method='lm',formula=y~x) + geom_point(aes(colour=grp)) + facet_wrap( ~ grp, nrow = 1)
  

}

#' Exploratory analysis with ENCODE PWM scores
#'
#' \code{deepbindAnalysis} 
#'
#' @param 
#' @return 
encodeAnalysis <- function(fimoDir, trainingD) {
  setwd("/p/keles/CAGI2015/volumeA/Data/After_Nov5")
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
  
  scoreNames <- paste("en", scoreTypes, sep = "_")
  
  setwd(fimoDir)
  loadFimoData(refOrAlt, "ENCODE", fimoThreshold)
  
  scoreList <- addAllToScoreList(trainingD, fimoS, scoreNames, scoreList)
  
  createAllCorrelations(scoreNames)
  
  
}

#' Exploratory analysis with Jaspar PWM scores
#'
#' \code{jasparAnalysis} 
#'
#' @param 
#' @return 
jasparAnalysis <- function(dataDir, fimoDir, trainingD) {
  setwd("/p/keles/CAGI2015/volumeB/Data/")
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
  
  scoreNames <- paste("ja", scoreTypes, sep = "_")
  
  setwd(fimoDir)
  loadFimoData(refOrAlt, "Jaspar", fimoThreshold)
  
  createAllCorrelations(scoreList, scoreNames)
}


#' Output all scores
#' 
#' \code{outputAllScores} saves the FactorBook, Deepbind, ENCODE, and Jaspar
#'   scores.
#' 
#' 
outputAllScores <- function(scoreList) {
  setwd("/p/keles/CAGI2015/volumeB/Predictors/")
  save(scoreList, file = "scoreList_Training_alt_thres0.001_pwm.RData")
}

#' Run all functions in here
#' 
#' 
#' 
#' @param trainingDLocation: string, location of the training D
create_fimo_predictors <- function(trainingDLocation = "/p/keles/CAGI2015/volumeA/Data/After_Nov5/4-eQTL-causal_SNPs_sample_v2.txt", 
                                   dataDir = "/p/keles/CAGI2015/volumeA/Data/After_Nov5",
                                   fimoDir = "/p/keles/CAGI2015/volumeB/Fimo") {
  trainingD <- read.table(trainingDLocation, header = T)
  scoreList <- makeScoreList()
  
  scoreTypes <- c("minPval", "maxScore", "numOcc", "sumScore")
  
  factorBookAnalysis(scoreList, trainingD)
  deepbindAnalysis()
  encodeAnalysis()
  jasparAnalysis()
  outputAllScores(scoreList)
}

