#' Computing Uniprobe 8mer scores
#' 
#' \code{run_Uniprobe_scores} is the main function for computing Uniprobe 8mer
#'   scores
#' 
#' Depends on uniprobe_8merScores.R, uniprobe_extract_8merscores_SS.R, and 
#'   uniprobe_construct_summary.8merscores_ref_SS.R
#' Default locations set for Keles' lab and are read in from
#'   base_getDefaultFileLocations.R.
#' @return Null. Data saved to disk.
run_Uniprobe_scores <- function() {
  fileLocations <- getDefaultFileLocations()
  uniprobeDirectory <- fileLocations[[1]]
  inputDTLocation <- fileLocations[[6]]
  outputDir <- fileLocations[[7]]
  
  # creating human 8mers - saved to outputDir/escores8mer_ref_human_v2.Rda"
  extract8mers(uniprobeDirectory, inputDTLocation, 
               outputDir, organismName = "human")
  
  # creating mouse 8mers - saved to outputDir/escores8mer_ref_mouse_v2.Rda"
  extract8mers(uniprobeDirectory, inputDTLocation,
               outputDir, organismName = "mouse")
  
  # final output saved to disk 
  #   - in outputDir/Features8mersEscores_ref_v2.Rda
  constructSummary8merScores(outputDir, inputDTLocation)
}

