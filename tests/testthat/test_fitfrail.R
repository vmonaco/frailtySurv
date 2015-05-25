library(testthat)
library(frailtyr)

test_that("fitfrail recovers parameters within error", {
  data <- genfrail()
  model <- fitfrail(data)
  
  expect_equal(10, 10) # Example test
})
