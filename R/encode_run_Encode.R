
#' ENCODE data 
#' 
#' \code{run_encode} is the main function for interpreting ChIP, DNase, and DNA methylation data
#'
#' Relies on output from  uniprobe_run_Uniprobe_scores.R.
#' Depends on encode_extract_encodePeaks.R, encode_DFToGranges.R, and 
#' encode_construct_overlapping.encodepeaks.narrow1_SS.R.
#' Default locations set for Keles' lab and are read in from
#'   base_getDefaultFileLocations.R. 
#' @return Null. Data saved to disk.
#' 
#' @export
#' @name run_encode
run_encode <- function() {
  fileLocations <- getDefaultFileLocations()
  chipLocation <- fileLocations[[2]]; dnaseLocation <- fileLocations[[3]]
  dnamethLocation <- fileLocations[[4]]; peakDir <- fileLocations[[5]]
  inputDTLocation <- fileLocations[[6]]; outputDir <- fileLocations[[7]]
  
  # takes ChIP, DNase, and DNA methylation data from ENCODE data
  # returns no output
  # saves to disk at outputDir/encodePeaks_v3.Rda
  extract_encodepeaks(chipLocation, dnaseLocation, dnamethLocation, peakDir, outputDir)
  
  # then finds the overlapping peaks
  # returns no output
  # saves to disk at outputDir/FeaturesEncodePeaksNarrow1_v3.Rda
  construct_overlapping_encodepeaks(inputDTLocation, outputDir)  
}