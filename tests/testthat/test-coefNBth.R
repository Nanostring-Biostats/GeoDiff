test_that("coefNBth produces desired results from output of fitNBthDE", {

  #### Specs for coefNBth


  ## 1. when fullpara=TRUE, the output parameters should be regression coefficients, threshold and r in a list.
  ##    Both threshold and r are positive.
  ## 2. when fullpara=FALSE, the output parameters should be regression coefficients only in a list

  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions


  data(NBthDEmod2)
  coeffull <- coefNBth(NBthDEmod2, fullpara=TRUE)
  coefreg <- coefNBth(NBthDEmod2, fullpara=FALSE)

  ## 1. when fullpara=TRUE, the output parameters should be regression coefficients, threshold and r in a list.
  ##    Both threshold and r are positive.


  expect_true(all(c(colnames(NBthDEmod2$X), c("r", "threshold")) == rownames(coeffull$estimate)))
  expect_true(all(coeffull$estimate["r", ] > 0))
  expect_true(all(coeffull$estimate["threshold", ] > 0))

  ## 2. when fullpara=FALSE, the output parameters should be regression coefficients only in a list

  expect_true(all(colnames(NBthDEmod2$X) == rownames(coefreg$estimate)))
})




test_that("coefNBth produces desired results from output of fitNBthmDE", {

  #### Specs for coefNBth


  ## 1. when fullpara=TRUE, the output parameters should be regression coefficients, threshold and r in a list.
  ##    Both threshold and r are positive.
  ## 2. when fullpara=FALSE, the output parameters should be regression coefficients only in a list

  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)

  data(NBthmDEmod2)
  coefmfull <- coefNBth(NBthmDEmod2, fullpara=TRUE)
  coefmreg <- coefNBth(NBthmDEmod2, fullpara=FALSE)

  ## 1. when fullpara=TRUE, the output parameters should be regression coefficients, threshold and r in a list.
  ##    Both threshold and r are positive.


  expect_true(all(c(colnames(NBthmDEmod2$X), c("r", "threshold")) == rownames(coefmfull$estimate)))
  expect_true(all(coefmfull$estimate["r", ] > 0))
  expect_true(all(coefmfull$estimate["threshold", ] > 0))

  ## 2. when fullpara=FALSE, the output parameters should be regression coefficients only in a list

  expect_true(all(colnames(NBthmDEmod2$X) == rownames(coefmreg$estimate)))
})
