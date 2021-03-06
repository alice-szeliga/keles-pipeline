#' Capitalizes words in a string
#'
#' \code{simpleCap} capitalizes the first letter of every word in a string
#'
#' Helper function for readChip, readDnase, and readDnameth
#' @param x, a string with words separated by spaces
#' 
#' @export
#' @name simpleCap
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#' read in ChIP data from ENCODE
#'
#' \code{readChip} reads in broad and narrow peak files from ENCODE ChiP metadata
#'
#' @param chipLocation: string, file location of ENCODE ChIP metadata.tsv
#' @return data.table containing the ChIP info
#' 
#' @export
#' @name readChip
readChip <- function(chipLocation) {
  metachip.all<-fread(chipLocation, header=TRUE)
  
  metachip<-metachip.all[which(metachip.all[,7, with=FALSE]=="GM12878"),]
  setnames(metachip, names(metachip), sapply(sapply(names(metachip), simpleCap), function(x) gsub(" ", "", x, fixed = TRUE)))
  
  #table(metachip[OutputType=="peaks",FileFormat]) #check File format for peaks
  #    bed broadPeak    bed narrowPeak  bigBed broadPeak bigBed narrowPeak 
  #              108                62               108                62 
  
  #peak files with ENCODE Consortium as a lab are seen from Richard Myers from ENCODE website.
  metachip_allpeaks<-metachip[OutputType=="peaks"]
  metachip_allpeaks[Lab %in% "ENCODE Consortium", Lab:="Richard Myers, HAIB"]
  
  metachip_narrowPeak<-metachip_allpeaks[FileFormat=="bed narrowPeak"]
  metachip_broadPeak<-metachip_allpeaks[FileFormat=="bed broadPeak"]
  setkey(metachip_narrowPeak, ExperimentTarget, Lab)
  setkey(metachip_broadPeak, ExperimentTarget, Lab)
  metachip_peaks<-rbind(metachip_narrowPeak, metachip_broadPeak)
  return(metachip_peaks)
  
  #length(unlist(unique(metachip_peaks$ExperimentTarget)))
  #99
  
  #table(table(unlist(metachip_peaks$ExperimentTarget), unlist(metachip_peaks$Lab))) #check the number of peakfiles for an experiment target from  a lab
  #  0   1   2   3   4 
  #679  65  42   3   3 
  #dim(metachip_peaks) #170*41
}


#' read in Dnase data from ENCODE
#'
#' \code{readChip} reads in narrow peak files from ENCODE Dnase metadata
#'
#' @param chipLocation: string, file location of ENCODE Dnase metadata.tsv
#' @return data.table containing the Dnase info
#' 
#' @export
#' @name readDnase
readDnase <- function(dnaseLocation) {
  ##DNase##
  metadnase.all<-fread(dnaseLocation, header=TRUE)
  metadnase<-metadnase.all[which(metadnase.all[,7, with=FALSE]=="GM12878"),]
  setnames(metadnase, names(metadnase), sapply(sapply(names(metadnase), simpleCap), function(x) gsub(" ", "", x, fixed = TRUE)))
  metadnase[,ExperimentTarget:=rep("DNase", nrow(metadnase))]
  
  metadnase_allpeaks<-metadnase[OutputType=="peaks"]
  metadnase_peaks<-metadnase_allpeaks[FileFormat=="bed narrowPeak"]
  setkey(metadnase_peaks, ExperimentTarget, Lab)
  return(metadnase_peaks)
}

#Crawford has 5 replicates; 5 fastq files for 3 replicates (error?); 5 bam files; 1 peak (1 bigBed, 1 bed.gz)
##########

#' read in DNA methylation data from ENCODE
#'
#' \code{readChip} reads in files from ENCODE DNA 
#'   methylation metadata
#'
#' @param chipLocation: string, file location of ENCODE DNA methylation 
#'   metadata.tsv
#' @return data.table containing the DNA methylation info
#' 
#' @export
#' @name readDnameth
readDnameth <- function(dnamethLocation) {
  ##DNAmethylation##
  metamethyl.all<-fread(dnamethLocation, header=TRUE)
  metamethyl<-metamethyl.all[which(metamethyl.all[,7, with=FALSE]=="GM12878"),]
  setnames(metamethyl, 
           names(metamethyl), 
           sapply(sapply(names(metamethyl), simpleCap), 
                  function(x) gsub(" ", "", x, fixed = TRUE)))
  metamethyl[,ExperimentTarget:=rep("DNAmethylation", nrow(metamethyl))]
  
  metamethyl_allpeaks<-metamethyl[OutputType=="methylation state at CpG"]
  metamethyl_peaks<-metamethyl_allpeaks[FileFormat=="bed bedMethyl"]
  setkey(metamethyl_peaks, ExperimentTarget, Lab)
  return(metamethyl_peaks)
}

#' Combines all peak files from ENCODE
#'
#' \code{readChip} reads in ChIP seq, DNase seq, DNAmethylation and combines
#'   the peaks into one file
#'
#' @param chipLocation: string, file location of ENCODE ChIP metadata.tsv
#' @param dnaseLocation: string, file location of ENCODE DNase metadata.tsv
#' @param dnamethLocation: string, file location of ENCODE DNA methylation 
#'   metadata.tsv
#' @return data.table containing all info from ChIP, DNase, and DNA methylation
#' 
#' @export
#' @name getMetaPeaks
getMetaPeaks <- function(chipLocation, dnaseLocation, dnamethLocation) {
  metachip_peaks <- readChip(chipLocation)
  metadnase_peaks <- readDnase(dnaseLocation)
  metamethyl_peaks <- readDnameth(dnamethLocation)
  meta_peaks<-rbind(metachip_peaks, metadnase_peaks, metamethyl_peaks)
  return(meta_peaks)
}


#create peak file names to put them under a directory automatically for chipseq; 
# manually for dnaseseq and DNA methylation


##not commented out by AS
#unique(meta_peaks$FileFormat)
#[1] "bed narrowPeak" "bed broadPeak"  "bed bedMethyl"

#' Get name of Bed files for all peaks
#'
#' \code{getBedNames} gets a data.ta______
#'
#' @param file_type: string of file name
#' @param meta_peaks: data.table of the ENCODE metadata
#' @return a data.table of strings of all the bed file names
#' 
#' @export
#' @name getBedNames
getBedNames <- function(file_type, meta_peaks) {
  mapply(function(name) paste(name, ".bed", sep=""), 
         meta_peaks$FileAccession[grep(file_type, meta_peaks$FileFormat)])
}

#' Read in a bed file
#'
#' \code{getBedDF} reads in a bed file and saves it as a data.table
#'
#' @param peakname: string, filename
#' @return a data.table of ___
#' helper fn called by getBedList
#' 
#' @export
#' @name getBedDF
getBedDF <- function(peakname) {
  bed.df<-read.table(peakname, header=FALSE)
  bed.dt<-as.data.table(bed.df)
  setkey(bed.dt, V1, V2)
  bed.dt
}

#' Turns bed files into data.tables
#'
#' \code{getBedList} gets a list of data.tables containing bed file info
#'
#' @param bedNames: a list(?) of .bed file names
#' @return a list of data.tables
#' helper fn called by getAllLists
#' 
#' @export
#' @name getBedList
getBedList <- function(bedNames) {
  apply(as.matrix(bedNames), 
        1, 
        getBedDF)
}

#' Gets normal, unique, and narrowpeak bed files
#' 
#' \code{getAllLists} separates the metadata of peaks into narrow peaks, broad
#'   peaks, and methylation data, then returns a variety of lists of bed files,
#'   including the peaks unique for genomic ranges.
#'
#' @param peakname: 
#' @return a data.table of ___
#' helper fn called by getBedList
#' 
#' @export
#' @name getAllLists
getAllLists <- function(meta_peaks) {
  bednames.narrowpeak <- getBedNames("bed narrowPeak", meta_peaks)
  bednames.broadpeak <- getBedNames("bed broadPeak", meta_peaks)
  bednames.methyl <- getBedNames("bed bedMethyl", meta_peaks)
  
  #BROADPEAK LIST: 108
  bed.broadpeak.list<-getBedList(bednames.broadpeak)
  #METHYLATION LIST: 65
  bed.methyl.list<- getBedList(bednames.methyl)
  ###Handling narrowpeak
  #NARROWPEAK LIST: 4
  bed.narrowpeak.list0<-getBedList(bednames.narrowpeak)
  
  #files have no summit called in the 10th column; refer to ENCODE narrowPeak in https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  #Just for investigation, not for data manipulation
  bed.narrowpeakIndBlank<-which(mapply(function(dt) length(which(dt$V10==-1)), bed.narrowpeak.list0)!=0)
  #names(bed.narrowpeakIndBlank)
  # [1] "ENCFF001HHQ" "ENCFF001HHR" "ENCFF001EXA" "ENCFF001EXB" "ENCFF001EXJ"
  # [6] "ENCFF001EXK" "ENCFF001EXU" "ENCFF001EXV" "ENCFF001CUF" "ENCFF001CUG"
  meta_narrowpeakBlank<-meta_peaks[which(is.element(meta_peaks$FileAccession, names(bed.narrowpeakIndBlank))),]
  
  #files have summit called in the 10th column
  #modified the peaks to be a length of 151bps around the summit.
  bed.narrowpeakIndSummit<-which(mapply(function(dt) length(which(dt$V10==-1)), bed.narrowpeak.list0)==0)
  bed.narrowpeak.list1<-bed.narrowpeak.list0
  bed.narrowpeak.list1[bed.narrowpeakIndSummit]<-mapply(function(bed) {se=apply(cbind(bed$V2-5,bed$V2+bed$V10-75), 1, max); ee=apply(cbind(bed$V3+5, bed$V2+bed$V10+76), 1, min); bed$V2=se; bed$V3=ee; bed; }, bed.narrowpeak.list0[bed.narrowpeakIndSummit], SIMPLIFY=FALSE)
  
  
  
  #### existing comments ####
  #Ignore this at this point because the comparison across the three concludes to have narrowpeak.list0 or narrowpeak.list1.
  #bed.narrowpeak.list2<-bed.narrowpeak.list0
  #Test to see my categories are correct. Correct! No 5 found; 
  #t0<-mapply(function(bed) {t(apply(cbind(bed$V2, bed$V3, bed$V10), 1, function(triple) {if((triple[2]-triple[1])<=151) 1 else {if(triple[3]>=70) {if(triple[2]-(triple[1]+triple[3])>=71) 2 else 3} else {if(triple[2]-(triple[1]+triple[3])>=71) 4 else 5}}}))}, bed.narrowpeak.list0[bed.narrowpeakIndSummit], SIMPLIFY=FALSE)
  
  #bed.narrowpeak.list2[bed.narrowpeakIndSummit]<-mapply(function(bed) {newse=t(apply(cbind(bed$V2, bed$V3, bed$V10), 1, function(triple) {se1=c(triple[1]+triple[3]-75, triple[1]+triple[3]+76); se2=c(triple[1]-5, triple[1]-5+151); se3=c(triple[2]+5-151,triple[2]+5); if((triple[2]-triple[1])<=151) triple[1:2] else {if(triple[3]>=70) {if(triple[2]-(triple[1]+triple[3])>=71) se1 else se3} else {if(triple[2]-(triple[1]+triple[3])>=71) se2 else c(0,0)}}})); bed$V2=newse[,1]; bed$V3=newse[,2]; bed}, bed.narrowpeak.list0[bed.narrowpeakIndSummit], SIMPLIFY=FALSE)
  
  #Comparison across narrowpeak.list0, narrowpeak.list1, narrowpeak.list2
  #rbind(t(mapply(function(bed) mean(bed$V3-bed$V2), bed.narrowpeak.list0[bed.narrowpeakIndSummit])), t(mapply(function(bed) mean(bed$V3-bed$V2), bed.narrowpeak.list1[bed.narrowpeakIndSummit])), t(mapply(function(bed) mean(bed$V3-bed$V2), bed.narrowpeak.list2[bed.narrowpeakIndSummit])))
  #Not much difference between narrowpeak.list1, narrowpeak.list2
  #Only go with narrowpeak.list0 and narrowpeak.list1
  
  #See their distribution; their maximum is 151.
  #mapply(function(bed) summary(bed$V3-bed$V2), bed.narrowpeak.list1[bed.narrowpeakIndBlank])
  #mapply(function(bed) summary(bed$V3-bed$V2), bed.narrowpeak.list1[bed.narrowpeakIndSummit])
  ###########################
  
  #peak list without peak adjustment for narrowPeak
  bed.list.raw<-c(bed.narrowpeak.list0, bed.broadpeak.list, bed.methyl.list) 
  #peak list with peak adjustment for narrowPeak
  bed.list.narrow1<-c(bed.narrowpeak.list1, bed.broadpeak.list, bed.methyl.list) 
  
  metaInd<-apply(as.matrix(names(bed.list.raw)), 1, function(bedname) grep(bedname, meta_peaks$FileAccession))
  TargetNames<-gsub("-human", "", meta_peaks$ExperimentTarget[metaInd])
  LabNames<-mapply(function(x) sub(",", "", x[2]), strsplit(meta_peaks$Lab[metaInd], split=' ', fixed=TRUE))
  FeatureNames<-mapply(function(target, lab) paste(target, "_", lab, sep=""), TargetNames, LabNames)
  UniqueFeatureNames<-unique(FeatureNames)
  
  names(bed.list.raw)<-FeatureNames
  names(bed.list.narrow1)<-FeatureNames
  
  #source("/p/keles/CAGI2015/volumeB/Exploratory/DFToGranges.R")
  #### calling helper functions from DFToGranges.R - not needed in package
  #GenomicRanges object for peak list from bed.list.raw
  uniquepeaks.raw<-dt2uniquegranges(bed.list.raw) 
  #GenomicRanges object for peak list from bed.list.narrow1
  uniquepeaks.narrow1<-dt2uniquegranges(bed.list.narrow1) 
  return(list(meta_peaks, bed.list.raw, bed.list.narrow1, uniquepeaks.raw, uniquepeaks.narrow1))
}

#metachip_narrowPeak[grep("H3", unlist(metachip_narrowPeak[, 16, with=FALSE])), FileAccession]


#' Extracts peaks from ENCODE data
#'
#' \code{extract_encode peaks} takes ENCODE data and outputs a list of peaks
#'   with metadata, as well as lists of bed files and unique peaks.
#'
#' @param chipLocation: string, file location of ENCODE ChIP metadata.tsv
#' @param dnaseLocation: string, file location of ENCODE DNase metadata.tsv
#' @param dnamethLocation: string, file location of ENCODE DNA methylation 
#'   metadata.tsv
#' @param peakdir: string, path of directory containing peak files in .bed format
#' @param outdir: string, path of directory to save output
#' @return NULL, but saves to disk an .Rda file of
#'   meta_peaks: a data.table of all the peak metadata
#'   bed.list.raw: a list of the raw bed file names 
#'   bed.list.narrow1: a list of all the narrowpeak files
#'   uniquepeaks.raw: a list of the raw unique peaks
#'   uniquepeaks.narrow1: a list of the unique narrowpek files
#' This is what is actually called by run_Encode.R
#' 
#' @export
#' @name extract_encodepeaks
extract_encodepeaks <- function(chipLocation, dnaseLocation, 
                                dnamethLocation, peakdir, outdir) {
  meta_peaks <- getMetaPeaks(chipLocation, dnaseLocation, dnamethLocation)
  setwd(peakdir)
  allOutput <- getAllLists(meta_peaks)
  meta_peaks <- allOutput[[1]]; bed.list.raw <- allOutput[[2]]
  bed.list.narrow1 <- allOutput[[3]]; uniquepeaks.raw <- allOutput[[4]]
  uniquepeaks.narrow1 <- allOutput[[5]]
  
  previouswd <- getwd()
  setwd(outdir)
  save(meta_peaks, bed.list.raw, bed.list.narrow1, uniquepeaks.raw, uniquepeaks.narrow1, file="encodePeaks_v3.Rda")
  ### RESET WORKING DIR
  setwd(previouswd)
}
