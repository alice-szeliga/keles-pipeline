% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fimo_create_fimo_predictors.R
\name{createAllCorrelations}
\alias{createAllCorrelations}
\title{Perform simple correlations}
\arguments{
\item{scoreNames:}{list of strings, contains names of the analysis and the
different scores. 
Ex. c("fb_minPval", "fb_maxScore", "fb_numOcc", "fb_sumScore")}
}
\value{
Null. Prints out.
}
\description{
\code{createAllCorrelations} will take the scores for either FactBook,
  Jaspar, or ENCODE and will perform correlations for the minimum P-value,
  the maximum score, the sum of the scores, and the number of occurences.
}

