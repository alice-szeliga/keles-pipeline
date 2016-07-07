
test_that("loadData works correctly", {
  ##oh crap it loads, not pastes whoops. fix this later
  x <- loadData("ref", "FactorBook", 0.001)
  x1 <- load("Training_ref_thres0.001_FactorBook_minPval.RData")
  x2 <- load("Training_ref_thres0.001_FactorBook_maxScore.RData")
  x3 <- load("Training_ref_thres0.001_FactorBook_sumScore.RData")
  x4 <- load("Training_ref_thres0.001_FactorBook_numOcc.RData")
  expect_equal(x[[1]], x1)
  expect_equal(x[[2]], x2)
  expect_equal(x[[3]], x3)
  expect_equal(x[[4]], x4)
})