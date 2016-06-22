context("8mer scores")

# data tables used for testing
dt0 = data.table()
dt2 = data.table(1,2)
dt3 = data.table(1,2,3)
dt4 = data.table(1,2,3,4)
dt5 = data.table(1,2,3,4,5)
dt6 = data.table(1,2,3,4,5,6)
dtMixed = data.table(c(1,2,3,4), c(5,4,3,2), c(10,20,0,5))

dt6Header = data.table(c(1, 1), c(2, 2), c(3, 3), c(4, 4), c(5, 5), c(6, 6))
scoresDT = data.table(1:10, 2:11, 3:12, 4:13, 5:14, 6:15)



test_that("Helper fn getCols3to5 works for data.table with different # cols", {

  expect_equal(getCols3to5(dt0), data.table())
  expect_equal(getCols3to5(dt2), data.table())
  expect_equal(getCols3to5(dt3), data.table(V3 = 3))
  expect_equal(getCols3to5(dt4), data.table(V3 = 3, V4 = 4))
  expect_equal(getCols3to5(dt5), data.table(V3 = 3, V4 = 4, V5 = 5))
  expect_equal(getCols3to5(dt6), data.table(V3 = 3, V4 = 4, V5 = 5))
})

test_that("Piping operator %>% works correctly in getCols3to5", {
  expect_equal(
    ncol(scoresDT) %>% min(5, .) %>% seq %>% setdiff(., 1:2) %>% scoresDT[,.,with=FALSE],
    scoresDT[,setdiff(seq(min(5, ncol(scoresDT))), 1:2), with=FALSE])
})

test_that("Helper fn getRangeCols works for data.table with different values", {
  expect_equal(getRangeCols(dtMixed), apply(dtMixed, 2, range))
})

test_that("removeHeader does remove first row if file too long", {
  expect_equal(dt6, removeHeader(1, dt6Header))
})

test_that("removeHeader does not remove first row if file right length", {
  expect_equal(dt6, removeHeader(1, dt6))
})
