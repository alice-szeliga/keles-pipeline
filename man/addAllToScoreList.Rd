% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fimo_create_fimo_predictors.R
\name{addAllToScoreList}
\alias{addAllToScoreList}
\title{Adds a score to the score list}
\arguments{
\item{trainingD:}{a table of eQTL causal SNPs}

\item{fimoS:}{TYPE?, a score for FIMO}

\item{scoreName:}{list of strings, which are indicating analysis type and type of score
Ex. "fb_minPval"}

\item{scoreList:}{a vector containing all scores for FactorBook, ENCODE,
Jaspar, and Deepbhind}
}
\description{
\code{addAllToScoreList} adds the minPval, maxScore, numOcc, and sumScore to
  the scoreslist for either FactBook or Jaspar
}

