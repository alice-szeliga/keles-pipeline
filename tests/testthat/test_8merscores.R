context("8mer scores")

test_that("Helper fn getCols3to5 works for data.table with different # cols", {
  dt0 = data.table()
  dt2 = data.table(1,2)
  dt3 = data.table(1,2,3)
  dt4 = data.table(1,2,3,4)
  dt5 = data.table(1,2,3,4,5)
  dt6 = data.table(1,2,3,4,5,6)
  expect_equal(getCols3to5(dt0), data.table())
  expect_equal(getCols3to5(dt2), data.table())
  expect_equal(getCols3to5(dt3), data.table(V3 = 3))
  expect_equal(getCols3to5(dt4), data.table(V3 = 3, V4 = 4))
  expect_equal(getCols3to5(dt5), data.table(V3 = 3, V4 = 4, V5 = 5))
  expect_equal(getCols3to5(dt6), data.table(V3 = 3, V4 = 4, V5 = 5))
})

test_that("Helper fn getRangeCols works for data.table with different values", {
  dtA = data.table(seq(4), seq(5,2), c(10,20,0,5))
  expect_equal()
})
