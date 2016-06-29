#' Converts data.table to Granges
#'
#' \code{dt2granges} converts a data.table object to Granges
#'
#' modified from a function on http://davetang.org/muse/2015/02/04/bed-granges/
#' @param dt: a data.table object
#' @returns Granges object containing same information
dt2granges <- function(dt) { #
  df <- as.data.frame(dt)
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]

  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }

  library("GenomicRanges")

  if(length(df)==3){
    # start position increased by 1 bc/ bed files start from 0, 
    #   and end right after ranges. 
    # i.e. the end position is not a peak region.
    gr <- with(df, GRanges(chr, IRanges(start+1, end))) 
  } 
  else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start+1, end), id=id))
  } 
  else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start+1, end), id=id, score=score))
  } 
  else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start+1, end), id=id, score=score, 
                           strand=strand))
  }
  return(gr)
}

#' Converts data.table to UniqueGranges
#'
#' \code{dt2granges} converts a data.table object to unique Granges
#'
#' @param dt: a data.table object
#' @returns Granges object containing unique peaks
dt2uniquegranges <- function(bed.list) {
  peaks.granges.list0 <- mapply(dt2granges, bed.list)
  peaks.granges.list <- lapply(peaks.granges.list0, reduce)

  #FeatureNames[which(abs(as.numeric(mapply(summary, peaks.granges.list0)[1,])-as.numeric(mapply(summary, peaks.granges.list)[1,]))!=0)] #Find annotations with reduce function working, i.e. they have duplicated ranges.
  #                CTCF                 CTCF                 EZH2 
  #    "CTCF_Bernstein"          "CTCF_Iyer"     "EZH2_Bernstein" 
  #               H2AFZ              H3K27ac             H3K27me3 
  #   "H2AFZ_Bernstein"  "H3K27ac_Bernstein" "H3K27me3_Bernstein" 
  #            H3K36me3              H3K4me1              H3K4me2 
  #"H3K36me3_Bernstein"  "H3K4me1_Bernstein"  "H3K4me2_Bernstein" 
  #             H3K4me3             H3K79me2               H3K9ac 
  # "H3K4me3_Bernstein" "H3K79me2_Bernstein"   "H3K9ac_Bernstein" 
  #             H3K9me3             H4K20me1 
  # "H3K9me3_Bernstein" "H4K20me1_Bernstein" 

  TargetLabInd <- as.numeric(factor(names(bed.list), 
                                    levels=unique(names(bed.list))))
  unique.peaks.granges.list <- vector("list", 
                                      length=length(unique(TargetLabInd)))
  UniqueTargetLabInd <- which(table(TargetLabInd) == 1)
  unique.peaks.granges.list[UniqueTargetLabInd] 
        <- peaks.granges.list[which(is.element(TargetLabInd, 
                                                UniqueTargetLabInd))]
  DuplicateTargetLabInd <- which(table(TargetLabInd) > 1)
  unique.peaks.granges.list[DuplicateTargetLabInd] 
        <- apply(as.matrix(DuplicateTargetLabInd), 
                 1, 
                 function(ind) {dup.ind<-which(is.element(TargetLabInd, ind)); 
                                dup.granges.list<-peaks.granges.list[dup.ind]; 
                                Reduce(union, dup.granges.list)})
  names(unique.peaks.granges.list) <- unique(names(bed.list))
  unique.peaks.granges.list
}
