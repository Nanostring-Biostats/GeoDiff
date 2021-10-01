
test_that("QuanRange produces desired results, CTA", {

    #### Specs for QuanRange
    # 1 The function outputs a GeoMx S4 class
    #   with length same as length of sample IDs (rownames) in phenoData for each probs input.
    #   The colname is the input prob.


    library(dplyr)
    ### Initializing CTA objects before running tests
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    withr::local_dir(tmp_dir)
    # Run data through (required) upstream functions
    data("demoData") # for tests of structure of demoData itself, see test-scoretest.R
    NSGMS <- demoData # NanoStringGeoMxSet => NSGMS
    # Calling QuanRange this early should error since fitPoisBG(_sp) hasn't been run yet.
    expect_error(
        QuanRange(NSGMS, split = FALSE, probs = c(0.75, 0.8, 0.9, 0.95)),
        "Please run \`fitPoisBG\` first."
    )
    # Two datasets to analyze: without (NSGMS) and with (NSGMS_sp) slide groupings
    # Estimate Poisson background sample-feature factor mode
    NSGMS <- fitPoisBG(NSGMS)
    NSGMS <- aggreprobe(NSGMS, use = "cor")

    # Estimate Poisson background sample-feature factor model for multiple slides
    NSGMS_sp <- fitPoisBG(NSGMS, groupvar = "slide name", size_scale = "sum")
    # Calling QuanRange on NSGMS with split=TRUE should error
    expect_error(
        QuanRange(NSGMS, split = TRUE, probs = c(0.75, 0.8, 0.9, 0.95)),
        "Please run `fitPoisBG` first with `groupvar`."
    )

    # Case 1: "single" slide consideration
    test_probs <- c(0.75, 0.8, 0.9, 0.95)
    case1 <- QuanRange(NSGMS,
        split = FALSE,
        probs = test_probs
    )

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case1, "NanoStringGeoMxSet"))
    #   with length same as length of sample IDs (rownames) in phenoData
    #   for each probs input.
    for (p in test_probs) {
        p_len <- length(na.omit(pData(case1)[, which(colnames(pData(case1)) == p)])) # makes sure no NAs are present
        expect_true(length(row.names(pData(case1))) == p_len)
    }
    #   The colname is the input prob.
    expect_true(all(test_probs %in% colnames(pData(case1))))

    # Case 2: "multiple" slides consideration
    case2 <- QuanRange(NSGMS_sp,
        split = TRUE,
        probs = test_probs
    )

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case2, "NanoStringGeoMxSet"))
    #   with length same as length of sample IDs (rownames) in phenoData
    #   for each probs input.
    unique_ids <- unique(pData(NSGMS_sp)$`slide name`) # pull out the unique ids

    for (p in test_probs) {
        p_len <- length(na.omit(pData(case2)[, which(colnames(pData(case2)) == p)])) # makes sure no NAs are present
        expect_true(length(row.names(pData(case2))) == p_len)
    }
    #   The colname is the input prob.
    expect_true(all(test_probs %in% colnames(pData(case2))))

    # Expect that quantRange values are differnt between single and multiple cases
    expect_false(all(pData(case1)$`0.8` == pData(case2)$`0.8`))
})


## 2 It returns an error without running fitPoisBG.
test_that("It returns an error without running fitPoisBG.", {
    expect_error(
        QuanRange(demoData, probs = c(0.75, 0.8, 0.9, 0.95)),
        "Please run `fitPoisBG` first"
    )
    expect_error(
        QuanRange(demoData, probs = c(0.75, 0.8, 0.9, 0.95), split = TRUE),
        "Please run `fitPoisBG` first"
    )
})


## 3 It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.
test_that("It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.", {
    res <- fitPoisBG(demoData, size_scale = "first")
    expect_error(
        QuanRange(res, split = TRUE, probs = c(0.75, 0.8, 0.9, 0.95)),
        "Please run `fitPoisBG` first with `groupvar`"
    )
})
