# specs:
# When groupvar is not provided,
# - the function outputs a GeoMx S4 class with length same as length of negative probes, featfact, in featureData. The value is NA for non-negative probes.
# - If size_scale="first", sizefact[1]=1
# - If size_scale="sum", sum(sizefact)=1

# load example data:
data("demoData")

# run function:
res <- fitPoisBG(demoData, size_scale = "first")

## test basic structure:
# spec: "The function outputs a GeoMx S4 class with length same as length of ROIs, sizefact, in phenoData.
#       The value is NA for non-negative probes."
test_that("fitPoisBG returns a well-formatted geomxset object", {
    expect_true(class(res) == "NanoStringGeoMxSet")
    expect_true(nrow(phenoData(res)) == nrow(phenoData(demoData)))
    expect_true(nrow(featureData(res)) == nrow(featureData(demoData)))
    # 1 The function outputs a GeoMx S4 class with length same as length of ROIs, sizefact, in phenoData.
    expect_true("sizefact" %in% varLabels(res))
    expect_false(any(is.na(res[["sizefact"]])))
    # 2 The function outputs a GeoMx S4 class with length same as length of negative probes, featfact, in featureData. The value is NA for non-negative probes
    expect_true(all(is.na(fData(res)$featfact[fData(res)$featfact == "Endogenous"])))
    expect_identical(!is.na(fData(res)$featfact), fData(res)$CodeClass == "Negative")
})

## test size factors are correct:
test_that("sizefact is correct", {
    # Spec 3: If size_scale="first", sizefact[1]=1
    res <- fitPoisBG(demoData, size_scale = "first")
    expect_equal(pData(res)$sizefact[1], 1, tol = 1e-5)

    # Spec 4: If size_scale="sum", sum(sizefact)=1
    res <- fitPoisBG(demoData, size_scale = "sum")
    expect_equal(sum(pData(res)$sizefact), 1, tol = 1e-5)
})


## test that values haven't changed from June 2021 initial release:
test_that("fitPoisBG is stable", {
    res <- fitPoisBG(demoData, size_scale = "sum")
    expect_equal(pData(res)$sizefact[c(1, 10, 50)], c(0.011925349, 0.013114282, 0.007818129), tol = 1e-5)
    expect_equal(fData(res)$featfact[c(7932, 7933, 7934)], c(319, 380, 299), tol = 1e-2)
})

# specs
# When groupvar is provided but not found in the phenodata,
# 1. The function returns an error saying this groupvar is not found in the S4 object.

# load example data:
res <- fitPoisBG(demoData, size_scale = "first")

test_that("The function returns an error", {
    expect_error(
        fitPoisBG(demoData, size_scale = "first", groupvar = "slidename"),
        "is not found in the S4 object"
    )
})

# specs
# When groupvar is provided and found in the phenodata, but the it only has one unique value,
# 1. The function returns a warning message saying that the groupvar has only one value

# load example data:
res <- fitPoisBG(demoData, size_scale = "first")

test_that("The function returns a warning messag", {
    expect_warning(
        fitPoisBG(demoData, size_scale = "first", groupvar = "segment"),
        "has only one value"
    )
})

# specs:
# When groupvar is provided and found in the phenodata with more than one unique value,
# 1 When groupvar is provided, the function outputs a GeoMx S4 class with length same as length of ROIs, sizefact, in phenoData.
# 2 The function outputs a GeoMx S4 class with length same as length of negative probes, featfact, in featureData.
#   The value is NA for non-negative probes.
# 3 If size_scale="first", sizefact[1]=1
# 4 If size_scale="sum", sum(sizefact)=1

# load example data:
res <- fitPoisBG(demoData, size_scale = "first")

# run function:
res <- fitPoisBG(demoData, size_scale = "first", groupvar = "slide name")

## test basic structure:
# - The function outputs a GeoMx S4 class with length same as length of ROIs, sizefact, in phenoData.
# - The function outputs a GeoMx S4 class with length same as length of negative probes, featfact, in featureData for each unique slide value.
# - the group variable name for slide id of ROIs in experimentData: fitPoisBG_sp_var
test_that("fitPoisBG returns a well-formatted geomxset object", {
    expect_true(class(res) == "NanoStringGeoMxSet")
    # 1 The function outputs a GeoMx S4 class with length same as length of ROIs, sizefact, in phenoData
    expect_true(nrow(phenoData(res)) == nrow(phenoData(demoData)))
    expect_true(nrow(featureData(res)) == nrow(featureData(demoData)))
    expect_true(all(is.na(fData(res)$featfact[fData(res)$featfact != "Negative"])))
    # 2 The function outputs a GeoMx S4 class with length same as length of negative probes, featfact, in featureData.
    #   The value is NA for non-negative probes.
    for (i in grep("featfact_", fvarLabels(res))) {
        expect_identical(!is.na(fData(res)[, i]), fData(res)$CodeClass == "Negative")
    }
    expect_identical(
        fvarLabels(res)[grep("featfact_", fvarLabels(res))],
        paste0("featfact_", unique(demoData[["slide name"]]))
    )
    expect_equal(notes(res)$fitPoisBG_sp_var, "slide name")
})

## test size factors are correct:
test_that("sizefact is correct", {
    # 3 Spec: If size_scale="first", sizefact[1]=1
    res <- fitPoisBG(demoData, size_scale = "first", groupvar = "slide name")
    expect_equal(pData(res)$sizefact[1], 1, tol = 1e-5)

    # 4 Spec: If size_scale="sum", sum(sizefact)=1
    res <- fitPoisBG(demoData, size_scale = "sum", groupvar = "slide name")
    expect_equal(sum(pData(res)$sizefact), 1, tol = 1e-5)
})


## test that values haven't changed from June 2021 initial release:
test_that("fitPoisBG is stable", {
    res <- fitPoisBG(demoData, size_scale = "sum", groupvar = "slide name")
    expect_equal(pData(res)$sizefact[c(1, 10, 50)], c(0.011925349, 0.013114282, 0.007818129), tol = 1e-5)
    expect_equal(fData(res)$"featfact_6panel-old-slide1 (PTL-10891)"[c(7932:7934)], c(313.8457, 438.5242, 356.8383), tol = 1e-2)
})
