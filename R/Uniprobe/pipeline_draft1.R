## makefile in shell script

#!/path ##change to where it should be

##getting input??
## this is Rscript, not bash
## name it inputFilenames.R
{{
cat("Enter locations of: 1) 2) 3): \n")
numArgs = 3
args <- commandArgs(TRUE)
if (length(args) != numArgs) {
  cat("Please supply the appropriate number of filenames.\n")
  quit()
}
file1 <- args[1]
}}

Rscript --vanilla inputFilenames.R 

## Filenames should be Input??
## dependencies


all: 8merScores.R extract_8merscores_ref_SS.v2.R construct_summary.8merscores_ref_SS.v2.R

8merScores.R:
  Rscript #inputs - scoresDT

extract_8merscores_ref_SS.v2.R: 8merScores.R
  Rscript #filename inputs here?

construct_summary.8merscores_ref_SS.v2.R: 8merScores.R extract_8merscores_ref_SS.v2.R
  Rscript #filename inputs?

#clean rule
clean: 
  rm -f #where will output be saved?

# run extract on human

# run extract on mouse
source()


# run construct on previous 2 outputs
