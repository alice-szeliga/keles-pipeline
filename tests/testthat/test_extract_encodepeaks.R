#source("C:/Users/aszeliga/Documents/keles/eqtl.pipeline/R/encode_extract_encodepeaks.R")

## DNase file of only 500 lines
tinyDnaseLocation <- "C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyDNase.tsv"
## correct output from previous version of code - metadnase_peaks
tinyDnaseactual <- load("C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyDNase.Rda")

## DNA methylation file of only 500 lines
tinyDnamethLocation <- "C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyDNAmeth.tsv"
## correct output from previous version of code - metamethyl_peaks
tinyDnamethactual <- load("C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyDNAmeth.Rda")

## ChIP file of only 500 lines
tinyChIPLocation <- "C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyChIP.tsv"
## correct output from previous version of code - metachip_peaks
tinyChIPactual <- load("C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyChIP.Rda")

## data.table file of only 100 eQTLs
tinyDTLocation <- "C:/Users/aszeliga/Documents/keles/eqtl.pipeline/data/test_data/tinyDT.RData"

test_that("readChip fn gives correct output for a small dataset", {
  tinyChIPout <- readChip(tinyChIPLocation)
  expect_true(all.equal(tinyChIPout, metachip_peaks))
})

test_that("readDnase fn gives correct output for a small dataset", {

  tinyDnaseout <- readDnase(tinyDnaseLocation)
  expect_true(all.equal(tinyDnaseout, metadnase_peaks))
})

test_that("readDnameth fn gives correct output for a small dataset", {
  tinyDnamethout <- readDnameth(tinyDnamethLocation)
  expect_true(all.equal(tinyDnamethout, metamethyl_peaks))
})


test_that("getAllLists() works", {
  allOutput <- getAllLists(getMetaPeaks(tinyChIPLocation, 
                                      tinyDnaseLocation, 
                                      tinyDnamethLocation))
  meta_peaks <- allOutput[[1]]; bed.list.raw <- allOutput[[2]]
  bed.list.narrow1 <- allOutput[[3]]; uniquepeaks.raw <- allOutput[[4]]
  uniquepeaks.narrow1 <- allOutput[[5]]
  for (i in 1:5) {
    expect_true(all.equal(allOutput[[i]], allResults[[i]]))
  }
})

test_that("extract_encodepeaks() works. This fn calls all others in the file", {
  
})