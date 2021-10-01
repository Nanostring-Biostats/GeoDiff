# reqs:
# - matrix of lower tail probability in assayData slot called lowtail_prob : lowtail_prob
# - matrix of upper tail probability in assayData slot called uptail_prob: uptail_prob
# - the dispersion parameter in experimentData named disper: disper
# - matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called low_outlier : low_outlier
# - matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called up_outlier : up_outlier

# Create temporary directory that will get destroyed after this block is executed.
tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
# Change to the temporary directory (will set back to getwd() once block is executed.)
# load example data:
data("demoData")

#### test single-slide version: ----------------------------------------------
# run method
res <- fitPoisBG(demoData, size_scale = "first")
res <- aggreprobe(res, use = "cor")

# 1 When generate_ppplot is TRUE, a figure is generated, when FALSE, no figure is generated
res <- diagPoisBG(res, generate_ppplot = FALSE)

# 2 when padj=FALSE, each element of sum of lowtail_prob and uptail_prob in the assay slot named lowtail_prob and uptail_prob equals to 1
# test basic structure:
test_that("single-slide diagPoisBG  returns a well-formatted geomxset object", {
    withr::local_dir(tmp_dir)
    expect_true(class(res) == "NanoStringGeoMxSet")

    # matrix of lower tail probability in assayData slot called lowtail_prob : lowtail_prob
    expect_true(is.matrix(assayDataElement(res, "lowtail_prob")))

    # matrix of upper tail probability in assayData slot called uptail_prob: uptail_prob
    expect_true(is.matrix(assayDataElement(res, "uptail_prob")))

    # the dispersion parameter in experimentData named disper: disper
    expect_true(notes(res)$disper > 0)
    expect_equal(length(notes(res)$disper), 1)

    # matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called low_outlier : low_outlier
    expect_true(is.matrix(assayDataElement(res, "low_outlier")))
    expect_true(length(setdiff(assayDataElement(res, "low_outlier"), c(0, 1))) == 0)

    # matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called up_outlier : up_outlier
    expect_true(is.matrix(assayDataElement(res, "up_outlier")))
    expect_true(length(setdiff(assayDataElement(res, "up_outlier"), c(0, 1))) == 0)
})



#### same tests but for multi-slide version: --------------------------------
# run method
res <- fitPoisBG(demoData, size_scale = "first", groupvar = "slide name")
res <- diagPoisBG(res, split = TRUE)

# test basic structure:
test_that("multi-slide diagPoisBG returns a well-formatted geomxset object", {
    withr::local_dir(tmp_dir)
    expect_true(class(res) == "NanoStringGeoMxSet")

    # matrix of lower tail probability in assayData slot called lowtail_prob : lowtail_prob
    expect_true(is.matrix(assayDataElement(res, "lowtail_prob")))

    # matrix of upper tail probability in assayData slot called uptail_prob: uptail_prob
    expect_true(is.matrix(assayDataElement(res, "uptail_prob")))

    # the dispersion parameter in experimentData named disper: disper
    expect_true(notes(res)$disper > 0)
    expect_equal(length(notes(res)$disper), 1)

    # matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called low_outlier : low_outlier
    expect_true(is.matrix(assayDataElement(res, "low_outlier")))
    expect_true(length(setdiff(assayDataElement(res, "low_outlier"), c(0, 1))) == 0)

    # matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called up_outlier : up_outlier
    expect_true(is.matrix(assayDataElement(res, "up_outlier")))
    expect_true(length(setdiff(assayDataElement(res, "up_outlier"), c(0, 1))) == 0)
})


## test that values haven't changed from June 2021 initial release:
test_that("diagPoisBG is stable", {
    withr::local_dir(tmp_dir)
    res <- fitPoisBG(demoData, size_scale = "first")
    res <- diagPoisBG(res, generate_ppplot = FALSE)
    expect_equal(notes(res)$disper, 1.323681, tol = 1e-5)
})

## 3 It returns an error without running fitPoisBG.
test_that("It returns an error without running fitPoisBG.", {
    withr::local_dir(tmp_dir)
    expect_error(
        diagPoisBG(demoData, generate_ppplot = FALSE),
        "Please run `fitPoisBG` first"
    )
    expect_error(
        diagPoisBG(demoData, generate_ppplot = FALSE, split = TRUE),
        "Please run `fitPoisBG` first"
    )
})


## 4 It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.
test_that("It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.", {
    withr::local_dir(tmp_dir)
    res <- fitPoisBG(demoData, size_scale = "first")
    expect_error(
        diagPoisBG(res, split = TRUE, generate_ppplot = FALSE),
        "Please run `fitPoisBG` first with `groupvar`"
    )
})
