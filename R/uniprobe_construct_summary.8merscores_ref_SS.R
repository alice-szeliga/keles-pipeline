#' Imports human and mouse data
#'
#' \code{createEnrscoresList} imports the human and mouse .Rda files created
#'   in extract_8merscores and creates and labels a list of scores.
#'
#' Helper function for constructSummary8merScores
#' @export
#' @name createEnrscoresList
createEnrscoresList <- function(outputDir) {
  humanEscores <- paste(outputDir, "/escores8mer_ref_human_v2.Rda", sep = "")
  load(humanEscores)
  enrscoresList_human <- enrscoresList
  NameM_human <- NameM
  
  mouseEscores <- paste(outputDir, "/escores8mer_ref_mouse_v2.Rda", sep = "")
  load(mouseEscores)
  enrscoresList_mouse <- enrscoresList
  NameM_mouse <- NameM
  
  enrscoresList <- c(enrscoresList_human, enrscoresList_mouse)
  NameM <- rbind(NameM_human, NameM_mouse)
  
  #####Add unique protein names#####
  ProteinNames <- NameM[,2]
  names(enrscoresList) <- ProteinNames

  enrscoresList
}

#' Get absolute coordinates for the subsequences
#' 
#' \code{GetSubseqCoordAbs} takes something and makes something
#' 
#' #' Helper function for createSubseq21List()
#' @param subseq.len: int, 
#' @param center.len: 
#' @param center.pos: variable of positions from inputDTFile
#' @name GetSubseqCoordAbs
GetSubseqCoordAbs <- function(subseq.len, center.len, center.pos) {
  heading      <- ceiling( (subseq.len - center.len) / 2)
  foot         <- floor( (subseq.len - center.len) / 2)
  start        <- center.pos - heading
  end          <- center.pos + center.len - 1 + foot
  subseq.coord <- cbind(start, end)
  subseq.coord
}

#' Create a list from the DT
#' 
#' \code{createSubseq21List} reads in a DT of eQTL data and produces a list
#' 
#' Helper function for constructSummary8merScores
#' @param inputDTFile: string, location of DT
#' @name createSubseq21List
createSubseq21List <- function(inputDTFile) {
  # contains trainingDT
  load(inputDTFile)
  
  lenCenter <- mapply(nchar, trainingDT$RefAllele)
  table(trainingDT$pos - trainingDT$start + 1)
  table(lenCenter)
  
  subseq.coord.abs <- GetSubseqCoordAbs(subseq.len=21, 
                                        center.len=lenCenter, 
                                        center.pos=trainingDT$pos)
  subseq.coord.rel <- subseq.coord.abs - matrix(trainingDT$start, nrow(trainingDT), 2) + 1
  subseq.coord.rel.list <- lapply(seq_len(nrow(subseq.coord.rel)), 
                                  function(i) subseq.coord.rel[i,])
  
  enrscoresSubseq21List <- mapply(function(enrscores) 
    t(mapply(function(enrscore, rel.pos) 
      enrscore[rel.pos[1]:rel.pos[2]], 
      enrscores, 
      subseq.coord.rel.list)),
    enrscoresList, SIMPLIFY=FALSE)
  
  enrscoresSubseq21List
}

#' Summarizes results of 8mer scores
#'
#' \code{constructSummary8merScores} takes the 8mer files from extract_8merscores.R
#'  and creates a summary 
#'
#' @param inputDTFile: string, location of eQTL data table
#' @param outputDir: string, file location to save results
#' @return Null. Saves results to disk in outputDir
#' 
#' @export
#' @name constructSummary8merScores
constructSummary8merScores <- function(outputDir, inputDTFile) {
  
  enrscoresList <- createEnrscoresList(outputDir)
  
  SumEscores <- mapply(function(enrList) mapply(sum, enrList), enrscoresList)
  MeanEscores <- mapply(function(enrList) mapply(mean, enrList), enrscoresList)
  MaxEscores <- mapply(function(enrList) mapply(max, enrList), enrscoresList)
  
  enrscoresSubseq21List <- createSubseq21List(inputDTFile)
  
  SumSubseq21Escores <- mapply(rowSums, enrscoresSubseq21List)
  MeanSubseq21Escores <- mapply(rowMeans, enrscoresSubseq21List)
  MaxSubseq21Escores <- mapply(function(enrM) apply(enrM, 1, max), 
                               enrscoresSubseq21List)
  
  escoresList <- vector("list")
  escoresList[[1]] <- SumEscores
  escoresList[[2]] <- MeanEscores
  escoresList[[3]] <- MaxEscores
  escoresList[[4]] <- SumSubseq21Escores
  escoresList[[5]] <- MeanSubseq21Escores
  escoresList[[6]] <- MaxSubseq21Escores
  names(escoresList) <- c("sum", "mean", "max", "sum_bp21", 
                          "mean_bp21", "max_bp21")
  
  ## the final output
  outputFile <- paste(outputDir, "/Features8mersEscores_ref_v2.Rda", sep = "")
  save(escoresList, file=outputFile)
}
