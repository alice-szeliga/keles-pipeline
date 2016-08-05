#' Calls whole model upon eQTL data
#'
#' \code{runAll} calls all predictors on the eQTL data given to determine the 
#'   regulatory/expression-modulating status of each sequence variant.
#' 
#' The predictors used are FIMO PWM scores from ENCODE, JASAPAR, and 
#'   FACTORBOOK, as well as binding affinity scores from Uniprobe and Deepbind,
#'   and genome annotations from ENCODE for DNase, TF occupancy, histone 
#'   modifications, and methylation.
#' Sequence variants can be SNPs or small insertions/deletions.
#' The model requires MPRA data as a ground truth.
#' @export
#' @name runAll
runModel <- function() {
  run_fimo()
  runDeepbind()
  run_Uniprobe_scores()
  run_encode()
}