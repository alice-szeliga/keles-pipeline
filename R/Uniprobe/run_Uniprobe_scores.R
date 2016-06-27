# just here for testing purposes - delete for package
source('~/keles/R/eqtl.pipeline/R/8merScoresMODIFIED.R')
source('~/keles/R/eqtl.pipeline/R/extract_8merscores_ref_SS_MODIFIED2.v2.R')
source('~/keles/R/eqtl.pipeline/R/construct_summary.8merscores_ref_SS_MODIFIED.v2.R')

#' Computing Uniprobe 8mer scores
#' 
#' Main function for computing Uniprobe 8mer scores
#' 
#' Depends on 8merScores.R, extract_8merscores.R, and 
#'   construct_summary.8merscores_ref.R
#' Default locations set for Keles' lab.
#' @param uniprobeDirectory: string, location of /Uniprobe scores
#'   Should have /Human and /Mouse subdirectories
#' @param inputDTFile: string, location of eQTL data to analyze
#'   Data set from CAGI 2015 challenge is /ProcessedData/SampleDT.Rdata
#' @param outputDir: string, directory to save results to
#' @return Null. Data saved to disk.
run_Uniprobe_scores <- function(uniprobeDirectory = "/p/keles/CAGI2015/volumeB/Data/UniProbe", 
                                inputDTLocation = "/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData",
                                outputDir = "/p/keles/CAGI2015/volumeB/ProcessedData") {
  
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

