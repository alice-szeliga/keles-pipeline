test_that("Does getFactorName create correct names?", {
  expect_equal("Training_ref_thres0.001", getFactorName("ref"))
})


test_that("", {
  expect_equal("ref_SS.v1.fasta", getFastaName("ref"))
})