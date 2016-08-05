#' Calls Deepbind
#'
#' \code{run_Deepbind} submits 2 external calls to the command line to Deepbind
#' 
#' @export
#' @name run_Deepbind
run_Deepbind <- function() {
  deepBindDir = "/p/keles/CAGI2015/volumeB/SK/Fimo"
  commandAlt = paste("deepbind all_ids_from_deepbind.ids \< ", deepBindDir "ref_SS.v1.fasta \> ref_SS.v1_deepbind.txt", sep = "")
  print(commandAlt)
  commandRef = "deepbind all_ids_from_deepbind.ids < ", deepBindDir, "alt_SS.v1.fasta > alt_SS.v1_deepbind.txt", sep = "")

  system(commandAlt)
  system(commandRef)
}
