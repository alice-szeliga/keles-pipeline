#' Runs FactorBook, ENCODE and Jaspar on FIMO
#'
#' \code{run_Fimo} calculates the PWM for FactorBook, ENCODE, and Jaspar, then
#'   summarizes the scores and creates predictors for the set of eQTL data.
#'
#' @export
#' @name run_Fimo
run_Fimo <- function() {
  run_fimo_and_summarize(refOrAlt, memeFimoDir, outDir, pwmFiles) 
  create_fimo_predictors(trainingDLocation, dataDir, fimoDir) 
}