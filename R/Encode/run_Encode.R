### delete later
library(GenomicRanges)
library(rtracklayer)
library(data.table)


#' ENCODE data 
#' 
#' Main function for ChIP, DNase, and DNA methylation from ENCODE
#' 
#' Depends on run_Uniprobe_scores.R, extract_encodePeaksSS.v3.R, 
#'   DFToGranges.R, and construct_overlapping.encodepeaks.narrow1_SS.v3.R
#' Default locations set for Keles' lab.
#' The ENCODE data are peak/bed files obtained from encodeproject.org
#' @param chipLocation: string, location of ENCODE ChIP-seq metadata
#' @param dnaseLocation: string, location of ENCODE DNase-seq metadata
#' @param dnamethLocation: string, location of ENCODE DNA methylation metadata
#' @param peakDir: string, directory of ENCODE peaks
#' @param inputDTFile: string, location of eQTL data to analyze
#'   Data set from CAGI 2015 challenge is /ProcessedData/SampleDT.Rdata
#' @param outputDir: string, directory to save output in
#' @return Null. Data saved to disk.
run_encode <- function(chipLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/ChIP-seq/metadata.tsv",
                       dnaseLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/DNase-seq/metadata.tsv",
                       dnamethLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/DNAmethylation/metadata.tsv"
                       peakDir = "/p/keles/CAGI2015/volumeB/Data/ENCODE/peaks",
                       inputDTLocation = "/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData",
                       outputDir = "/p/keles/CAGI2015/volumeB/ProcessedData") {
  
  # takes ChIP, DNase, and DNA methylation data from ENCODE data
  extract_encodepeaks(chipLocation, dnaseLocation, dnamethLocation, peakDir, outputDir)
  
  # then finds the overlapping peaks
  # returns no output - saves to disk at outputDir/FeaturesEncodePeaksNarrow1_v3.Rda
  construct_overlapping_encodepeaks()  
}