context("extract_8merscores")


test_that("extract8mers catches invalid organismName", {
  
  expect_that(extract8mers(organismName = "rat"), 
              "Invalid organism name. Correct options: 'human' or 'mouse'")
})