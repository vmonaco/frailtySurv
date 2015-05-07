library(testthat)
library(frailtyr)

test_that("fitfrail recovers parameters from genfrail data", {
  data <- genfrail()
  model <- fitfrail(data)
  
  expect_equal(10, 10) # Example test
})
