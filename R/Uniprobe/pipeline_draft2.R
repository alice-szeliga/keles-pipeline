# # input filenames on command line
# cat("Enter file locations of: \n
#     1) Data table (.Rdata)
#     (following are optional)
#     2) UniProbe Directory (human)
#     3) UniProbe Directory (mouse)
#     4) 8mer scores \n")
# numArgs = 4
# args <- commandArgs(TRUE)
# 
# ## input in script from command line
# # if (length(args) == 0) {
# #   cat("Please enter", numArgs, "filenames.\n")
# #   input <- readLines(con="stdin",1)
# #   filenames <- strsplit(input, " ")
# # }
# 
# if (length(args) != numArgs) {
#   cat("Please supply the appropriate number of filenames:", numArgs, "\n")
#   quit()
# }
# inputDT <- args[1]
# uniprobeHumanDir <- args[2]
# uniprobeMouseDir <- args[3]
# 8merScores <- args[4]

# do this here???
rm(list=ls())

# change it so a fn does the work so you can call it w file name?
source("extract_8merscores_ref_SS.v2.R")

# call fn on mouse
extract8mers(uniprobeMouseDir, samples)

# call fn on human
extract8mers(uniprobeHumanDir, samples)

# maybe also make this a function?
source("construct_summary.8merscores_ref_SS.v2.R")

# output was saved in the previous function


