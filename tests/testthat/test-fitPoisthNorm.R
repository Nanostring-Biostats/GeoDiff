### There are two main functions being testing here
### fitPoisthNorm for single slide and multiple slides


test_that("fitPoisthNorm when split = FALSE produces desired results", {

    #### Specs for fitPoisthNorm
    # 1 Without providing values for ROIs_high, features_high, features_all, sizefact_start, sizefact_BG, the function returns the same value
    # 2 user need to set iterations =2 for now
    # 3 The function outputs a GeoMx S4 class with para0_norm, matrix of estimated parameters, in the featureData.
    # This matrix para0_norm has the following structure:
    # 3.1) 1 row for each feature (row name). If a feature is not in features_high, all columns will be NA.
    # 3.2) n+1 columns labeled var1, var2, ..., var<n>, var<n+1> where n is the length of ROIs_high elements.
    # 3.3) the n columns will have log2 expression (if feature is in features_high) or NA (otherwise).
    # 3.4) the n+1th column contains the threshold for each feature in features_high and NA otherwise.
    # 4 The function outputs a GeoMx S4 class with para_norm, matrix of estimated parameters, in the featureData. This matrix para_norm has the following structure:
    # 4.1) 1 row for each feature (row name) which is equal to the length of features_all
    # 4.2) n+1 columns labeled var1, var2, ..., var<n>, var<n+1> where n is the length of ROIs_high elements.
    # 4.3) the n columns will have log2 expression.
    # 4.4) the n+1th column contains the threshold for each feature in features_all.
    # 5 The function outputs a column called conv0 in featureData, with values in [NA, 0, 1] and length of 0s and 1s are the same as the length of features_high.
    # 6 The function outputs a column called conv in featureData, length same as features_all, and has values [NA, 0, 1]. The length of NAs equals the number of negative probes.



    # 7 when confac=0 and prior_type="contrast", preci1_norm value in
    #   experimentData will be single value repeated over an n-by-n matrix (where n is the length of ROI_high).
    #   This single value is equivalent, within 10 digits, to preci1con/n^2.

    library(dplyr)
    ### Initializing CTA objects before running tests
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    withr::local_dir(tmp_dir)
    # Run data through (required) upstream functions
    data("demoData") # for tests of structure of demoData itself, see test-scoretest.R
    set.seed(98070)
    demoData <- demoData[, c(1:5, 33:37)]
    NSGMS <- fitPoisBG(demoData) # "single" slide
    NSGMS <- aggreprobe(NSGMS, use = "cor")
    NSGMS <- BGScoreTest(NSGMS)

    # Negative and Non-Negative facets:
    NSGMS_neg <- NSGMS[which(fData(NSGMS)$CodeClass == "Negative"), ]
    NSGMS_pos <- NSGMS[-which(fData(NSGMS)$CodeClass == "Negative"), ]
    # feature factors
    featfact <- fData(NSGMS_neg)[, grep("featfact", fvarLabels(NSGMS_neg))]
    # scaling factors
    posdat <- Biobase::exprs(NSGMS_pos)
    sc1_scores <- fData(NSGMS_pos)[, "scores"]
    names(sc1_scores) <- fData(NSGMS_pos)[, "TargetName"]

    features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) & (sc1_scores < quantile(sc1_scores, probs = 0.95))) |>
        which() |>
        names()

    backmean <- mean(featfact)
    # Negative Binomial threshold model
    NSGMS <- fitNBth(NSGMS,
        features_high = features_high,
        sizefact_BG = NSGMS_neg$sizefact,
        threshold_start = backmean,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        threshold_fix = FALSE, # default but calling it explicitly here
        tol = 1e-8
    )
    high_ROIs <- sampleNames(NSGMS)[which((quantile(fData(NSGMS)[["para"]][, 1],
        probs = 0.90, na.rm = TRUE
    ) - notes(NSGMS)[["threshold"]]) * NSGMS$sizefact_fitNBth > 2)]

    features_all <- featureNames(NSGMS_pos)

    # Case 1: Similar to vignette but only treated as a single slide & iterations = 3
    # 2 User need to set iterations =2 for now
    expect_error(
        fitPoisthNorm(
            object = NSGMS,
            split = FALSE,
            iterations = 3, # note: not the default of 2
            ROIs_high = high_ROIs,
            features_high = features_high,
            features_all = features_all,
            sizefact_start = NSGMS[, high_ROIs][["sizefact_fitNBth"]],
            sizefact_BG = NSGMS[, high_ROIs][["sizefact"]],
            threshold_mean = backmean,
            preci2 = 10000,
            prior_type = "contrast",
            confac = 1, # default but called explicitly here for the test
            preci1con = 1 / 25
        ),
        "Only iterations=2 is allowed"
    )

    set.seed(98070)
    case1 <- fitPoisthNorm(
        object = NSGMS,
        split = FALSE,
        iterations = 2,
        ROIs_high = high_ROIs,
        features_high = features_high,
        features_all = features_all,
        sizefact_start = NSGMS[, high_ROIs][["sizefact_fitNBth"]],
        sizefact_BG = NSGMS[, high_ROIs][["sizefact"]],
        threshold_mean = backmean,
        sizescalebythreshold = TRUE,
        preci2 = 10000,
        prior_type = "contrast",
        confac = 1, # default but called explicitly here for the test
        preci1con = 1 / 25
    )

    set.seed(98070)
    NSGMS <- BGScoreTest(NSGMS)
    case1_df <- fitPoisthNorm(
        object = NSGMS,
        split = FALSE,
        iterations = 2,
        sizescalebythreshold = TRUE,
        preci2 = 10000,
        prior_type = "contrast",
        confac = 1,
        preci1con = 1 / 25
    )

    # 1 Without providing values for ROIs_high, features_high, features_all, sizefact_start, sizefact_BG, the function returns the same value  # expect same results without specifying the values.
    test_that("expect same results without specifying the values.", {
        expect_true(all.equal(case1, case1_df))
    })


    # 3 The function outputs a GeoMx S4 class...
    expect_true(inherits(case1, "NanoStringGeoMxSet"))
    # ...with para0_norm, matrix of estimated parameters, in the featureData.
    expect_true("para0_norm" %in% colnames(fData(case1)))
    expect_true(inherits(fData(case1)$para0_norm, "matrix"))
    # This matrix para0_norm has the following structure:
    # 3.1) 1 row for each feature (row name). If a feature is not in features_high, all columns will be NA.
    para0_norm <- fData(case1)$para0_norm
    expect_true(nrow(para0_norm) == dim(case1)[1])
    expect_true(all(is.na(para0_norm[setdiff(row.names(para0_norm), features_high), ]))) # i.e., get non-features_high rows and check if all are NAs
    # 3.2) n+1 columns labeled var1, var2, ..., var<n>, var<n+1> where n is the length of ROIs_high elements.
    expect_true(all(colnames(para0_norm) == paste0("var", 1:(length(high_ROIs) + 1))))
    # 3.3) the n columns will have log2 expression (if feature is in features_high) or NA (otherwise).
    expect_false(any(is.na(para0_norm[features_high, 1:length(high_ROIs)])))
    expect_true(all(is.na(para0_norm[setdiff(row.names(para0_norm), features_high), 1:length(high_ROIs)])))
    # 3.4) the n+1th column contains the threshold for each feature in features_high and NA otherwise.
    expect_false(any(is.na(para0_norm[features_high, (length(high_ROIs) + 1)])))
    expect_true(all(is.na(para0_norm[setdiff(row.names(para0_norm), features_high), (length(high_ROIs) + 1)])))

    # 4 The function outputs a GeoMx S4 class...
    # (tested above)
    # with para_norm, matrix of estimated parameters, in the featureData.
    expect_true("para_norm" %in% colnames(fData(case1)))
    expect_true(inherits(fData(case1)$para_norm, "matrix"))
    # This matrix para_norm has the following structure:
    # 4.1) 1 row for each feature (row name) which is equal to the length of features_all
    para_norm <- fData(case1)$para_norm
    expect_true(nrow(para_norm) == dim(case1)[1])
    # 4.2) n+1 columns labeled var1, var2, ..., var<n>, var<n+1> where n is the length of ROIs_high elements.
    expect_true(all(colnames(para_norm) == paste0("var", 1:(length(high_ROIs) + 1))))
    # 4.3) the n columns will have log2 expression.
    expect_false(any(is.na(para_norm[features_all, 1:length(high_ROIs)])))
    # 4.4) the n+1th column contains the threshold for each feature in features_all.
    expect_false(any(is.na(para_norm[features_all, (length(high_ROIs) + 1)])))

    # 5 The function outputs a column called conv0 in featureData, with values in [NA, 0] and length of 0s are the same as the length of features_high.
    expect_true("conv0" %in% colnames(fData(case1)))
    unique_values <- as.character(unique(fData(case1)$conv0))
    expect_true(all(unique_values %in% c("0", "1", NA)))
    to_test <- fData(case1)[fData(case1)$feature_high == 1, ]
    expect_true(nrow(to_test) == length(features_high))
    length_na <- length(which(is.na(fData(case1)$conv0)))
    expect_true(length_na == length(features_all) - length(features_high) + nrow(NSGMS_neg))

    # 6 The function outputs a column called conv in featureData, length same as features_all, and has values [NA, 0, 1]. The length of NAs equals the number of negative probes.
    expect_true("conv" %in% colnames(fData(case1)))
    expect_true(length(fData(case1)$conv) == length(fData(case1)$features_all))
    unique_values <- as.character(unique(fData(case1)$conv))
    expect_true(all(unique_values %in% c("0", "1", NA)))
    expect_true(any(is.na(unique_values)))
    length_na <- length(which(is.na(fData(case1)$conv)))
    expect_true(length_na == nrow(NSGMS_neg))



    # 7 when confac=0 and prior_type="contrast", preci1_norm value in
    #   experimentData will be single value repeated over an n-by-n matrix (where n is the length of ROI_high).
    #   This single value is equivalent, within 10 digits, to preci1con/n^2.
    case1_1 <- fitPoisthNorm(
        object = NSGMS,
        split = FALSE,
        iterations = 2,
        ROIs_high = high_ROIs,
        features_high = features_high,
        features_all = features_all,
        sizefact_start = NSGMS[, high_ROIs][["sizefact_fitNBth"]],
        sizefact_BG = NSGMS[, high_ROIs][["sizefact"]],
        sizescalebythreshold = TRUE,
        threshold_mean = backmean,
        preci2 = 10000,
        prior_type = "contrast",
        confac = 0,
        preci1con = 1 / 25
    )


    expect_false("preci1_norm" %in% names(notes(NSGMS))) # original object
    expect_true("preci1_norm" %in% names(notes(case1))) # case1
    expect_true("preci1_norm" %in% names(notes(case1_1))) # case1

    preci1_norm <- notes(case1_1)$preci1_norm
    expect_true(inherits(preci1_norm, "matrix"))
    expect_true(nrow(preci1_norm) == length(high_ROIs))
    expect_true(ncol(preci1_norm) == length(high_ROIs))
    the_testing_value <- (1 / 25) / (length(high_ROIs)^2)
    expect_true(length(unique(as.numeric(preci1_norm))) == 1)
    expect_true(round(unique(as.numeric(preci1_norm)), 10) == round(the_testing_value, 10))

    # 8 It returns an error without running fitPoisBG.
    test_that("It returns an error without running fitPoisBG.", {
        expect_error(
            fitPoisthNorm(
                object = demoData,
                split = FALSE,
                iterations = 2,
                threshold_mean = 1,
                preci2 = 10000,
                prior_type = "contrast",
                confac = 0,
                preci1con = 1 / 25
            ),
            "Please run `fitPoisBG` first."
        )
        expect_error(
            fitPoisthNorm(
                object = demoData,
                split = TRUE,
                iterations = 2,
                threshold_mean = 1,
                preci2 = 10000,
                prior_type = "contrast",
                confac = 0,
                preci1con = 1 / 25
            ),
            "Please run `fitPoisBG` first with `groupvar`."
        )
        demoData_fit <- fitPoisBG(demoData, split = FALSE)
        expect_error(
            fitPoisthNorm(
                object = demoData_fit,
                split = FALSE,
                iterations = 2,
                threshold_mean = 1,
                preci2 = 10000,
                prior_type = "contrast",
                confac = 0,
                preci1con = 1 / 25
            ),
            "Please run `fitNBth` first."
        )
    })


    # 9 It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.
    test_that("It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.", {
        expect_error(
            fitPoisthNorm(
                object = NSGMS,
                split = TRUE,
                iterations = 2,
                ROIs_high = high_ROIs,
                features_high = features_high,
                features_all = features_all,
                sizefact_start = NSGMS[, high_ROIs][["sizefact_fitNBth"]],
                sizefact_BG = NSGMS[, high_ROIs][["sizefact"]],
                threshold_mean = backmean,
                preci2 = 10000,
                prior_type = "contrast",
                confac = 0,
                preci1con = 1 / 25
            ),
            "Please run `fitPoisBG` first with `groupvar`"
        )
    })
})

test_that("fitPoisthNorm when split = TRUE produces desired results", {

    #### Specs for fitPoisthNorm_sp
    # 1. Given a GeoMx S4 object, fitPoisthNorm_sp runs the Poisson model-based
    #    normalization and log2 transformation on each element in "groupvar" individually. As such,
    #    the results for a given grouping/facet of the data should match the fitPoisthNorm
    #    results when an object is subset down to a single slide. Specifically, the following
    #    should be true:
    #   1.1 For a given element in groupvar, the corresponding column in the "threshold0" matrix, which is within featureData, should be identical to the single-patient case's fetureData's para0_norm[,n+1]th column.
    #   1.2 For a given element in groupvar, the corresponding column in the "threshold" matrix, which is within featureData, should be identical to the single-patient case's fetureData's para_norm[,n+1]th column.
    #   1.3 For a given element in groupvar, the normalized matrix "normmat0_sp", in the assayData slot, should be identical to that element's "normmat0" matrix, also in the assayData slot, for all samples within that element. In other words, the matrix within the "single slide" results (normmat0) should be a subset of the "multiple slide" results (normmat_sp).
    #   1.4 For a given element in groupvar, the normaized matrix "normmat_sp", in the assayData slot, should be identical to that element's "normmat" matrix, also in the assayData slot, for all samples within that element. In other words, the matrix within the "single slide" results (normmat) should be a subset of the "multiple slide" results (normmat_sp).
    #   1.5 For a given element in groupvar, the vector of sizefact, located in phenoData's sizefact_norm column, is identical to that element's sizefact_norm vector from running fitPoisthNorm (i.e., single grouping case).
    #   1.6 For a given element in groupvar, the vector of sizefact0, located in phenoData's sizefact_norm column, is identical to that element's sizefact_norm0 vector from running fitPoisthNorm (i.e., single grouping case).


    # First, create an NanoStringGeoMxSet object as in the fitPoisthNorm test but split into individual
    # objects based on "Subject ID" (i.e., groupvar).
    library(dplyr)
    ### Initializing CTA objects before running tests
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    withr::local_dir(tmp_dir)
    # Run data through (required) upstream functions
    data("demoData") # for tests of structure of demoData itself, see test-scoretest.R
    set.seed(98070)
    demoData <- demoData[, c(1:5, 33:37)]
    NSGMS_sp <- demoData
    # Estimate Poisson background sample-feature factor model for multiple slides (sp)
    NSGMS_sp <- fitPoisBG(NSGMS_sp, groupvar = "slide name", size_scale = "sum")
    NSGMS_sp <- aggreprobe(NSGMS_sp, use = "cor")
    NSGMS_sp <- BGScoreTest(NSGMS_sp, split = TRUE)

    # Negative and Non-Negative facets of each:
    NSGMS_sp_neg <- NSGMS_sp[which(fData(NSGMS_sp)$CodeClass == "Negative"), ]
    NSGMS_sp_pos <- NSGMS_sp[-which(fData(NSGMS_sp)$CodeClass == "Negative"), ]

    # feature factors, features, and backmeans:
    featfact_sp <- fData(NSGMS_sp_neg)[, grep("featfact_", fvarLabels(NSGMS_sp_neg))]
    posdat_sp <- Biobase::exprs(NSGMS_sp_pos)
    sc1_scores_sp <- fData(NSGMS_sp_pos)[, grep("scores", colnames(fData(NSGMS_sp_pos)))]
    rownames(sc1_scores_sp) <- fData(NSGMS_sp_pos)[, "TargetName"]
    features_high_sp <- apply(sc1_scores_sp, 2, function(x) {
        ((x > quantile(x, probs = 0.4)) & (x < quantile(x, probs = 0.95)))
    })
    features_high_sp <- names(which(apply(features_high_sp, 1, all)))
    backmean_sp <- (featfact_sp |> colMeans())[1]

    # Fit NB
    NSGMS_sp <- fitNBth(NSGMS_sp,
        features_high = features_high_sp,
        sizefact_BG = NSGMS_sp_neg$sizefact,
        threshold_start = backmean_sp,
        iterations = 5,
        start_para = c(200, 1),
        lower_sizefact = 0,
        lower_threshold = 100,
        tol = 1e-8
    )

    # Get high ROIs and all features for each:
    message("The Workflow_CTA_demoData_S4.Rmd file doesn't consider multiple slide case for determining high ROIs. Will subset data by slide name manually!")
    unique_pts <- unique(NSGMS_sp$`slide name`)
    high_ROIs_sp <- sampleNames(NSGMS_sp)[which(NSGMS_sp$sizefact_fitNBth * mean(colMeans(fData(NSGMS_sp)[, grep("featfact_", colnames(fData(NSGMS_sp)))], na.rm = TRUE)) > 2)]
    features_all_sp <- featureNames(NSGMS_sp_pos)

    # Now run fitPoisthNorm_sp on NSGMS_sp
    set.seed(98070)
    case_sp <- fitPoisthNorm(
        object = NSGMS_sp,
        split = TRUE,
        ROIs_high = high_ROIs_sp,
        features_high = features_high_sp,
        features_all = features_all_sp,
        sizefact_start = NSGMS_sp[, high_ROIs_sp][["sizefact_fitNBth"]],
        sizefact_BG = NSGMS_sp[, high_ROIs_sp][["sizefact_sp"]],
        sizescalebythreshold = TRUE,
        threshold_mean = backmean_sp,
        preci2 = 10000,
        prior_type = "contrast",
        preci1con = 1 / 25
    )

    set.seed(98070)
    case_sp_df <- fitPoisthNorm(
        object = NSGMS_sp,
        split = TRUE,
        iterations = 2,
        sizescalebythreshold = TRUE,
        preci2 = 10000,
        prior_type = "contrast",
        confac = 1, # default but called explicitly here for the test
        preci1con = 1 / 25
    )

    # 1 Without providing values for features_high, sizefact_BG, threshold_start, the function returns the same value
    # expect same results without specifying the values.
    test_that("expect same results without specifying the values.", {
        expect_true(all.equal(case_sp, case_sp_df))
    })
    # Subset NSGMS_sp into just the first patient and run fitPoisthNorm_sp separately.
    p <- unique_pts[1]
    NSGMS_pt1 <- NSGMS_sp[, which(phenoData(NSGMS_sp)[["slide name"]] == p)]
    high_ROIs_pt1 <- sampleNames(NSGMS_pt1)[which(NSGMS_pt1$sizefact_fitNBth * mean(fData(NSGMS_pt1)[, which(colnames(fData(NSGMS_pt1)) == paste0("featfact_", p))], na.rm = TRUE) > 2)]

    set.seed(98070)
    case_pt1 <- fitPoisthNorm(
        object = NSGMS_pt1,
        split = FALSE,
        iterations = 2,
        ROIs_high = high_ROIs_pt1,
        features_high = features_high_sp,
        features_all = features_all_sp,
        sizefact_start = NSGMS_pt1[, high_ROIs_pt1][["sizefact_fitNBth"]],
        sizefact_BG = NSGMS_pt1[, high_ROIs_pt1][["sizefact_sp"]],
        threshold_mean = backmean_sp,
        sizescalebythreshold = TRUE,
        preci2 = 10000,
        prior_type = "contrast",
        confac = 1, # default but called explicitly here for the test
        preci1con = 1 / 25
    )

    ## Run tests to compare
    # 2.1 For a given element in groupvar, the corresponding column in the "threshold0" matrix, which is within featureData,
    #     should be identical to the single-patient case's fetureData's para0_norm[,n+1]th column.
    threshold0_sp_mat <- fData(case_sp)$threshold0 # pull out the matrix from fData(case_sp)
    sp_compare <- threshold0_sp_mat[, which(colnames(threshold0_sp_mat) == unique_pts[1])] # i.e., pt1 pulled out of sp
    pt1_compare <- fData(case_pt1)$para0_norm[, length(high_ROIs_pt1) + 1]
    expect_true(identical(sp_compare, pt1_compare))
    rm(threshold0_sp_mat)
    rm(sp_compare)
    rm(pt1_compare)

    # 2.2 For a given element in groupvar, the corresponding column in the "threshold" matrix, which is within featureData,
    #     should be identical to the single-patient case's fetureData's para_norm[,n+1]th column.
    threshold_sp_mat <- fData(case_sp)$threshold # pull out the matrix from fData(case_sp)
    sp_compare <- threshold_sp_mat[, which(colnames(threshold_sp_mat) == unique_pts[1])] # i.e., pt1 pulled out of sp
    pt1_compare <- fData(case_pt1)$para_norm[, length(high_ROIs_pt1) + 1]
    expect_true(identical(sp_compare, pt1_compare))
    rm(threshold_sp_mat)
    rm(sp_compare)
    rm(pt1_compare)

    # 2.3 For a given element in groupvar, the normalized matrix "normmat0_sp", in the assayData slot, should be identical to
    #     that element's "normmat0" matrix, also in the assayData slot, for all samples within that element. In other words,
    #     the matrix within the "single slide" results (normmat0) should be a subset of the "multiple slide" results (normmat_sp).
    sp_compare <- assayData(case_sp)$normmat0_sp[, high_ROIs_pt1]
    pt1_compare <- assayData(case_pt1)$normmat0[, high_ROIs_pt1]
    expect_true(identical(sp_compare, pt1_compare))
    rm(sp_compare)
    rm(pt1_compare)

    # 2.4 For a given element in groupvar, the normaized matrix "normmat_sp", in the assayData slot, should be identical to
    #     that element's "normmat" matrix, also in the assayData slot, for all samples within that element. In other words,
    #     the matrix within the "single slide" results (normmat) should be a subset of the "multiple slide" results (normmat_sp).
    sp_compare <- assayData(case_sp)$normmat_sp[, high_ROIs_pt1]
    pt1_compare <- assayData(case_pt1)$normmat[, high_ROIs_pt1]
    expect_true(identical(sp_compare, pt1_compare))
    rm(sp_compare)
    rm(pt1_compare)

    # 2.5 For a given element in groupvar, the vector of sizefact, located in phenoData's sizefact_norm column, is identical to
    #     that element's sizefact_norm vector from running fitPoisthNorm (i.e., single grouping case).
    sp_compare <- pData(
        case_sp[, which(pData(case_sp)[["slide name"]] == p)]
    )$sizefact_norm
    pt1_compare <- pData(case_pt1)$sizefact_norm
    expect_true(identical(sp_compare, pt1_compare))
    rm(sp_compare)
    rm(pt1_compare)

    # 2.6 For a given element in groupvar, the vector of sizefact0, located in phenoData's sizefact_norm column, is identical to
    #     that element's sizefact_norm0 vector from running fitPoisthNorm (i.e., single grouping case).
    sp_compare <- pData(
        case_sp[, which(pData(case_sp)[["slide name"]] == p)]
    )$sizefact0_norm
    pt1_compare <- pData(case_pt1)$sizefact0_norm
    expect_true(identical(sp_compare, pt1_compare))
    rm(sp_compare)
    rm(pt1_compare)
})
