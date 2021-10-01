test_that("contrastNBth produces desired results from output of fitNBthDE", {

  #### Specs for contrastNBth


  ## 1. The function takes in a DE model as an input from fitNBthDE or fitNBthmDE
  ## 2. The user input test:statistical test, choose from c("two-sided", ">", "<")
  ## 3. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1


  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions
  data(NBthDEmod2) #

  coeffull <- contrastNBth(NBthDEmod2)
  coeftest <- contrastNBth(NBthDEmod2, method=matrix(c(0,1), 2, 1), baseline=0)

  coeftestupper <- contrastNBth(NBthDEmod2, method=matrix(c(0,1), 2, 1), baseline=0, test = ">")

  coeftestlower <- contrastNBth(NBthDEmod2, method=matrix(c(0,1), 2, 1), baseline=0, test = "<")



  expect_true(all(coeffull$estimate[2,] == coeftest$estimate[1,]))

  ## 3. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1

  expect_equal(unname(coeftestlower$p_value+coeftestupper$p_value),
               rep(1, length(coeftestlower$p_value)))


})




test_that("contrastNBth produces desired results from output of fitNBthmDE", {

  #### Specs for contrastNBth


  ## 1. The function takes in a DE model as an input from fitNBthDE or fitNBthmDE
  ## 2. The user input test:statistical test, choose from c("two-sided", ">", "<")
  ## 3. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1


  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions
  data(NBthmDEmod2) #

  coefmfull <- contrastNBth(NBthmDEmod2)
  coefmtest <- contrastNBth(NBthmDEmod2, method=matrix(c(0,1), 2, 1), baseline=0)

  coefmtestupper <- contrastNBth(NBthmDEmod2, method=matrix(c(0,1), 2, 1), baseline=0, test = ">")

  coefmtestlower <- contrastNBth(NBthmDEmod2, method=matrix(c(0,1), 2, 1), baseline=0, test = "<")



  expect_true(all(coefmfull$estimate[2,] == coefmtest$estimate[1,]))

  ## 3. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1

  expect_equal(unname(coefmtestlower$p_value+coefmtestupper$p_value),
               rep(1, length(coefmtestlower$p_value)))


})
