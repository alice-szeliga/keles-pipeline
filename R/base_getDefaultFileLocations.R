#' Get file locations of Data
#' 
#' \code{getDefaultFileLocations} sets and loads file paths of important data.
#' 
#' These locations are only stored here in the package. This is the only place
#'   that default locations need to be changed. Default locations set for 
#'   Keles' lab's server.
#' The ENCODE data are peak/bed files obtained from encodeproject.org
#' @param uniprobeDirectory: string, location of /Uniprobe scores
#'   Should have /Human and /Mouse subdirectories
#' @param chipLocation: string, location of ENCODE ChIP-seq metadata
#' @param dnaseLocation: string, location of ENCODE DNase-seq metadata
#' @param dnamethLocation: string, location of ENCODE DNA methylation metadata
#' @param peakDir: string, directory of ENCODE peaks
#' @param inputDTFile: string, location of eQTL data to analyze
#'   Data set from CAGI 2015 challenge is /ProcessedData/SampleDT.Rdata
#' @param outputDir: string, directory to save output in
#' @return A list 
#' 
#' @export
#' @name getDefaultFileLocations
getDefaultFileLocations <- function(
  uniprobeDirectory = "/p/keles/CAGI2015/volumeB/Data/UniProbe", 
  chipLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/ChIP-seq/metadata.tsv",
  dnaseLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/DNase-seq/metadata.tsv",
  dnamethLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/DNAmethylation/metadata.tsv",
  peakDir = "/p/keles/CAGI2015/volumeB/Data/ENCODE/peaks",
  inputDTLocation = "/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData",
  outputDir = "/p/keles/CAGI2015/volumeB/ProcessedData") {
  return(list(uniprobeDirectory, 
              chipLocation, 
              dnaseLocation, 
              dnamethLocation, 
              peakDir, 
              inputDTLocation, 
              outputDir))
}