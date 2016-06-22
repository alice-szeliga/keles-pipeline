#### remove this??
rm(list=ls())

#' Extract all 8mers for given eQTL sequences
#' 
#' \code{extract8mers} Computes enrichment scores for all 8mers on the 150bp sequences 
#'   given, output to drive
#' 
#' Calls helper functions from 8merScores.R
#' @param datadir: string, directory where Uniprobe scores have been saved
#' @param sampleDTLocation: string, file path of input data (in .Rdata format)
#' @param outputDir: string, directory to save output to 
#' @param organismName: string, "human" or "mouse
#' 
#' @return Returns null, saves output in outputDir/escores8mer_ref_human_v2.Rda

#' 
extract8mers <- function(uniprobeDir, inputDTFile, outputDir, organismName) {
  # getting correct subdirectory name from the input /Uniprobe location
  if (organismName == "human") {
    uniprobeDir = paste(uniprobeDir, "/Human", sep = "")
  }
  else if (organismName == "mouse") {
    uniprobeDir = paste(uniprobeDir, "/Mouse", sep = "")
  }
  else {stop("Invalid organism name. Correct options: 'human' or 'mouse' ")}
  
  # importing Uniprobe data
  all.result.files<-mixedsort(list.files(path=uniprobeDir, full.names=TRUE))
  # importing input data
  load(inputDTFile)

  # examining 8mers
  k <- 8
  ## list of lengths of all the reference sequences
  ## trainingDT needs to already be loaded into the environment
  lenCenter<-mapply(nchar, trainingDT$RefSeq)

  ## returns a list of all 8mer subsequences in each RefSeq
  subseqList <- mapply(function(seq, len)  
                       mapply(function(s) substr(seq, s, s+k-1), 
                              1:(len-k+1)),
                       trainingDT$RefSeq, lenCenter)
  ## labeling the subsequences with the eQTL positions?
  names(subseqList) <- trainingDT$ID


  ## This was already commented out AS
  #escoresCol<-mapply(CheckEscoresCol, all.result.files)
  #table(escoresCol)

  enrscoresList <- mclapply(as.matrix(all.result.files),
                          function(x) CreateEscores(file.name=x),
                          mc.cores=20)

  # splits each file.name into the File name and the Protein name
  fname <- mapply(function(file.name) strsplit(file.name, "/")[[1]][9],
                  all.result.files)
  ## generating 1:length(fname). i dunno why
  names(fname) <- seq(length(fname))

  # Creates unique names for each protein entered
  pname = sub("*(_8mers|_contig8mers).*", "", names)
  # Checks if any names were repeated
  if (anyDuplicated(pname)!=0) {
    print("Error: protein names repeated in extract8mers")
  }

  NameM <- cbind(fname, pname)
  colnames(NameM) <- c("File", "Protein")

  # saving in directory outputDir with a specific filename
  outputName <- paste("escores8mer_ref_", organismName, "_v2.Rda", sep = "")
  outputFile <- paste(outputDir, outputName, sep = "/")
  save(enrscoresList, NameM, file=outputFile)
}
