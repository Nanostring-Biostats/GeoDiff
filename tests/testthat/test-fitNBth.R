test_that("fitNBth produces desired results, CTA", {

    #### Specs for fitNBth
    # 1 Without providing values for features_high, sizefact_BG, threshold_start, the function returns the same value
    # 2 The function outputs a GeoMx S4 class with
    #   para0 in the experimentData as NA.
    # 3 The function outputs a GeoMx S4 class
    #   with para, a matrix of estimated parameters, in the featureData.
    #   This matrix has feature_high_fitNBth in columns(same as features_high)
    #   and parameters(signal, r) in columns.
    # 4 The function outputs sizefact_fitNBth in the phenoData, which is
    #   positive, same length as sizefact_BG
    # 5 The function outputs threshold in the experimentData. When
    #   threshold_fix=TRUE, threshold in the output is the same as
    #   threshold_start.

    library(dplyr)
    ### Initializing CTA objects before running tests
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    withr::local_dir(tmp_dir)
    # Run data through (required) upstream functions
    data("demoData") # for tests of structure of demoData itself, see test-scoretest.R

    # susbet samples
    demoData <- demoData[, c(1:5, 33:37)]
    set.seed(98070)
    NSGMS <- fitPoisBG(demoData, groupvar = "slide name", size_scale = "sum")
    NSGMS <- aggreprobe(NSGMS, use = "cor")
    NSGMS <- BGScoreTest(NSGMS, split = TRUE)
    # Negative and Non-Negative facets:
    NSGMS_neg <- NSGMS[which(fData(NSGMS)$CodeClass == "Negative"), ]
    NSGMS_pos <- NSGMS[-which(fData(NSGMS)$CodeClass == "Negative"), ]
    # feature factors per groupvar (i.e., slide name)
    scores_sp <- fData(NSGMS_pos)[, grep("scores_", fvarLabels(NSGMS_pos))]
    # scaling factors
    features_high <- apply(scores_sp, 2, function(x){
        ((x > quantile(x, probs = 0.4)) & (x < quantile(x, probs = 0.95)))
    }) |> (function(x) apply(x, 1, function(y) all(y)))() |>
        which() |>
        names()
    featfact_sp <- fData(NSGMS_neg)[, grep("featfact_", fvarLabels(NSGMS_neg))]
    thmean <- (featfact_sp |> colMeans())[1] # picks the first slide name's value

    ### Case 1: run examplar function from vignette and check each specification
    case1 <- fitNBth(NSGMS,
        features_high = features_high,
        sizefact_BG = NSGMS_neg$sizefact_sp,
        threshold_start = thmean,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = FALSE, # default but calling it explicitly here
        tol = 1e-8
    )

    # expect same results without specifying the values.
    set.seed(123)
    case1_df <- fitNBth(NSGMS,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = FALSE, # default but calling it explicitly here
        tol = 1e-8
    )
    # 1 Without providing values for features_high, sizefact_BG, threshold_start, the function returns the same value
    # expect same results without specifying the values.
    test_that("expect same results without specifying the values.", {
        expect_true(all.equal(case1, case1_df))
    })

    # 2 The function outputs a GeoMx S4 class...
    expect_true(inherits(case1, "NanoStringGeoMxSet"))
    # ...with para0 in the experimentData as 'NA'.
    expect_false(is.na(notes(case1)$para0)) # not NA
    expect_true(notes(case1)$para0 == "NA") # "NA"

    # 3 The function outputs a GeoMx S4 class...
    # (tested above)
    #   with para, a matrix of estimated parameters, in the featureData.
    expect_true("para" %in% colnames(fData(case1)))
    expect_true(inherits(fData(case1)$para, "matrix"))
    #   This matrix has feature_high_fitNBth in columns(same as features_high)
    to_test <- fData(case1)$para[fData(case1)$feature_high_fitNBth == 1, ]
    expect_true(nrow(to_test) == length(features_high))
    expect_true(all(row.names(to_test) == features_high))
    #   and parameters(signal, r) in columns.
    expect_true(all(c("signal", "r") == colnames(fData(case1)$para)))

    # 4 The function outputs sizefact_fitNBth in the phenoData, which is
    #   positive, same length as sizefact_BG
    expect_true("sizefact_fitNBth" %in% colnames(pData(case1)))
    expect_true(all(pData(case1)$sizefact_fitNBth >= 0)) # 0 is positive
    expect_true(length(pData(case1)$sizefact_fitNBth) == length(NSGMS_neg$sizefact))

    # 5 The function outputs threshold in the experimentData.
    expect_false("threshold" %in% names(notes(NSGMS))) # no threshold in experimentalData originally
    expect_true("threshold" %in% names(notes(case1)))
    # expect_true(is.numeric(notes(case1)$threshold)) # not strictly a spec
    #   When threshold_fix=TRUE, threshold in the output is the same as
    #   threshold_start.
    case1_1 <- fitNBth(NSGMS,
        features_high = features_high,
        sizefact_BG = NSGMS_neg$sizefact,
        threshold_start = thmean,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = TRUE,
        tol = 1e-8
    )
    expect_true(notes(case1_1)$threshold == thmean)
})

test_that("fitNBth produces desired results, WTA", {

    ### Same overall workflow as above but with the WTA kidney dataset

    library(dplyr)
    ### Initializing WTA objects before running tests
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    withr::local_dir(tmp_dir)
    # Run data through (required) upstream functions
    data("kidney")
    set.seed(98070)
    NSGMS <- kidney[, kidney$`slide name` %in% c("disease1B", "disease2B")]
    NSGMS <- NSGMS[, c(1:5, 11:15)]
    NSGMS <- fitPoisBG(NSGMS, groupvar = "slide name", size_scale = "sum")
    all0probeidx <- which(rowSums(exprs(NSGMS))==0)
    NSGMS <- NSGMS[-all0probeidx, ]
    NSGMS <- aggreprobe(NSGMS, use = "cor")
    # Negative and Non-Negative facets:
    NSGMS_neg <- NSGMS[which(fData(NSGMS)$CodeClass == "Negative"), ]
    NSGMS_pos <- NSGMS[-which(fData(NSGMS)$CodeClass == "Negative"), ]
    # feature factors per groupvar (i.e., slide name)
    featfact_sp <- fData(NSGMS_neg)[, grep("featfact_", fvarLabels(NSGMS_neg))]
    # scaling factors
    posdat <- Biobase::exprs(NSGMS_pos)
    gene_sum <- rowSums(posdat)
    features_high <- ((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95))) |>
        which() |>
        names()
    set.seed(123)
    genes_high <- sample(features_high, 1500) # subset
    thmean <- 1 * (featfact_sp |> colMeans())[1] # picks the first slide name's value

    ### Case 1: run examplar function from vignette and check each specification
    set.seed(123)
    case1 <- fitNBth(NSGMS,
        features_high = genes_high,
        sizefact_BG = NSGMS_neg$sizefact_sp,
        threshold_start = thmean,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = FALSE, # default but calling it explicitly here
        tol = 1e-8
    )
    # expect same results without specifying the values.
    set.seed(123)
    case1_df <- fitNBth(NSGMS,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = FALSE, # default but calling it explicitly here
        tol = 1e-8
    )
    # 1 Without providing values for features_high, sizefact_BG, threshold_start, the function returns the same value
    # expect same results without specifying the values.
    test_that("expect same results without specifying the values.", {
        expect_true(all.equal(case1, case1_df))
    })


    # 2 The function outputs a GeoMx S4 class...
    expect_true(inherits(case1, "NanoStringGeoMxSet"))
    # ...with para0 in the experimentData as 'NA'.
    expect_false(is.na(notes(case1)$para0)) # not NA
    expect_true(notes(case1)$para0 == "NA") # "NA"

    # 3 The function outputs a GeoMx S4 class...
    # (tested above)
    #   with para, a matrix of estimated parameters, in the featureData.
    expect_true("para" %in% colnames(fData(case1)))
    expect_true(inherits(fData(case1)$para, "matrix"))
    #   This matrix has feature_high_fitNBth in columns(same as features_high)
    to_test <- fData(case1)$para[fData(case1)$feature_high_fitNBth == 1, ]
    expect_true(nrow(to_test) == length(genes_high))
    expect_true(all(sort(row.names(to_test)) == sort(genes_high)))
    #   and parameters(signal, r) in columns.
    expect_true(all(c("signal", "r") == colnames(fData(case1)$para)))

    # 4 The function outputs sizefact_fitNBth in the phenoData, which is
    #   positive, same length as sizefact_BG
    expect_true("sizefact_fitNBth" %in% colnames(pData(case1)))
    expect_true(all(pData(case1)$sizefact_fitNBth >= 0)) # 0 is positive
    expect_true(length(pData(case1)$sizefact_fitNBth) == length(NSGMS_neg$sizefact))

    # 5 The function outputs threshold in the experimentData.
    expect_false("threshold" %in% names(notes(NSGMS))) # no threshold in experimentalData originally
    expect_true("threshold" %in% names(notes(case1)))
    # expect_true(is.numeric(notes(case1)$threshold)) # not strictly a spec
    #   When threshold_fix=TRUE, threshold in the output is the same as
    #   threshold_start.
    case1_1 <- fitNBth(NSGMS,
        features_high = genes_high,
        sizefact_BG = NSGMS_neg$sizefact_sp,
        threshold_start = thmean,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = TRUE,
        tol = 1e-8
    )
    expect_true(notes(case1_1)$threshold == thmean)
})
