#'
#'
#' \code{loadData} loads in the .Rdata files created by 
#'   fimo_run_fimo_and_summarize.R
#'
#' @param refOrAlt
#' @param functionName
#' @param fimoThreshold
#' @return Null. Loads in 4 .Rdata files
loadData <- function(refOrAlt, functionName, fimoThreshold) {
  for (name in c("minPval", "maxScore", "sumScore", "numOcc")) {
    load(paste("Training_", refOrAlt, "_thres", fimoThreshold, 
               "_", functionName, "_", name, ".RData", sep = ""))
  }
}

#' Exploratory analysis with FactorBook PWM scores
#'
#' \code{factorBookAnalysis} 
#'
#' @param 
#' @return 
factorBookAnalysis <- function(fimoDir) {
  setwd(dataDir)
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
  
  #setwd("/Users/keles/Research_Projects/CAGI2015/Fimo/")
  setwd(fimoDir)
  loadData(refOrAlt, "FactorBook", fimoThreshold)
  
  
}

#' Exploratory analysis with Deepbind scores
#'
#' \code{deepbindAnalysis} 
#'
#' @param 
#' @return 
deepbindAnalysis <- function() {
  setwd(dataDir)
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
  
}

#' Exploratory analysis with ENCODE PWM scores
#'
#' \code{deepbindAnalysis} 
#'
#' @param 
#' @return 
encodeAnalysis <- function(fimoDir) {
  setwd(fimoDir)
  loadData(refOrAlt, "ENCODE", fimoThreshold)
}

#' Exploratory analysis with Jaspar PWM scores
#'
#' \code{jasparAnalysis} 
#'
#' @param 
#' @return 
jasparAnalysis <- function(dataDir, fimoDir) {
  setwd(dataDir)
  trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
  
  setwd(fimoDir)
  loadData(refOrAlt, "Jaspar", fimoThreshold)
}


#' Output all scores
#' 
#' \code{outputAllScores} saves the FactorBook, Deepbind, ENCODE, and Jaspar
#'   scores.
#' 
#' 
outputAllScores <- function() {
  
}

#' Run all functions in here
#' 
#' 
#' 
#' 
create_fimo_predictors <- function(dataDir = "/p/keles/CAGI2015/volumeA/Data/After_Nov5",
                                   fimoDir = "/p/keles/CAGI2015/volumeB/Fimo") {
  factorBookAnalysis()
  deepbindAnalysis()
  encodeAnalysis()
  jasparAnalysis()
  outputAllScores()
}

