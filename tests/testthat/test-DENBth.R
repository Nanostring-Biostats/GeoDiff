test_that("DENBth produces desired results from output of fitNBthDE and coefNBth", {

  #### Specs for DENBth

  ## 1. In the output data.frame, the DE table should not have NA if NAto1=TRUE
  ## 2. In the output data.frame, the DE table has a column padj when adjp=TRUE

  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions
  data(NBthDEmod2)
  coeffull <- coefNBth(NBthDEmod2, fullpara=TRUE)



  DEtab1 <- DENBth(coeffull, variable = "group")
  DEtab2 <- DENBth(coeffull, variable = "group", NAto1=FALSE)
  DEtab3 <- DENBth(coeffull, variable = "group", padj=FALSE)
  ## 1. In the output data.frame, the DE table should not have NA if NAto1=TRUE
  if(sum(is.na(DEtab2$pvalue))>0)
    expect_equal(DEtab1$pvalue[is.na(DEtab2$pvalue)],rep(1, sum(is.na(DEtab2$pvalue))))
  ## 2. In the output data.frame, the DE table has a column padj when adjp=TRUE
  expect_equal(setdiff(colnames(DEtab1), colnames(DEtab3)), "adjp")

})




test_that("DENBth produces desired results from output of fitNBthmDE and coefNBth", {

  #### Specs for DENBth

  ## 1. In the output data.frame, the DE table should not have NA if NAto1=TRUE
  ## 2. In the output data.frame, the DE table has a column padj when adjp=TRUE

  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions
  data(NBthmDEmod2)
  coefmfull <- coefNBth(NBthmDEmod2, fullpara=TRUE)



  DEmtab1 <- DENBth(coefmfull, variable = "regiontubule")
  DEmtab2 <- DENBth(coefmfull, variable = "regiontubule", NAto1=FALSE)
  DEmtab3 <- DENBth(coefmfull, variable = "regiontubule", padj=FALSE)
  ## 1. In the output data.frame, the DE table should not have NA if NAto1=TRUE
  if(sum(is.na(DEmtab2$pvalue))>0)
    expect_equal(DEmtab1$pvalue[is.na(DEmtab2$pvalue)],rep(1, sum(is.na(DEmtab2$pvalue))))
  ## 2. In the output data.frame, the DE table has a column padj when adjp=TRUE
  expect_equal(setdiff(colnames(DEmtab1), colnames(DEmtab3)), "adjp")

})
