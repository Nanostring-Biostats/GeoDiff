### There are two main functions being testing here
### BGScoreTest for single slide and  multiple slides

test_that("BGScoreTest for single slide produces desired results", {

    # Desired results occurs when:
    # 1 The function outputs a GeoMx S4 class with p values in featureData with length same as length of targets. The p value is NA for negative probes.
    # 2 The function outputs a GeoMx S4 class with score values in featureData with length same as length of targets. The score value is NA for negative probes.
    # 3 All p values are between 0 and 1 (inclusive) for non-negative features.
    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    # 6 The order of pvalues is the same as scores.

    # Preamble/load example data
    data("demoData")
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    # Change to the temporary directory (will set back to getwd() once block is executed.)
    withr::local_dir(tmp_dir)
    expect_true(inherits(demoData, "NanoStringGeoMxSet"))

    # First run tests to ensure that the input data is as expected
    expect_true(nrow(pData(demoData)) == 88)
    expect_true(ncol(pData(demoData)) == 6)
    expect_true(nrow(pData(protocolData(demoData))) == 88)
    expect_true(ncol(pData(protocolData(demoData))) == 21)
    expect_true(inherits(assayData(demoData)[["exprs"]], "matrix"))
    expect_true(nrow(assayData(demoData)[["exprs"]]) == 8707)
    expect_true(ncol(assayData(demoData)[["exprs"]]) == 88)

    # Next estimate Poisson background sample-feature factor model
    set.seed(98070)
    demoData <- fitPoisBG(demoData)
    demoData <- aggreprobe(demoData, use = "cor")
    # Case 1: adjustment factor 5, no outlier removal, no prior
    case1 <- BGScoreTest(demoData,
        adj = 5,
        removeoutlier = FALSE, useprior = FALSE
    )

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case1, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    expect_false("pvalues" %in% colnames(fData(demoData))) # original does not have 'pvalues'
    expect_true("pvalues" %in% colnames(fData(case1))) # new object does have 'pvalues'
    # ...with length same as length of targets.
    expect_true(length(featureNames(case1)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case1_negatives <- case1[which(fData(case1)$CodeClass == "Negative"), ]
    expect_true(all(is.na(fData(case1_negatives)$pvalues)))

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false("scores" %in% colnames(fData(demoData))) # original does not have 'scores'
    expect_true("scores" %in% colnames(fData(case1))) # new object does have 'scores'
    # ...with length same as length of targets.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(all(is.na(fData(case1_negatives)$scores)))

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case1_positives <- case1[fData(case1)[["Negative"]] == FALSE, ] # non-negative features (endog + control)
    expect_true(all(fData(case1_positives)$pvalues >= 0 | fData(case1_positives)$pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case1)$pvalues))) == nrow(case1_positives))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case1)$scores))) == nrow(case1_positives))

    # 6 The order of pvalues is the same as scores.
    expect_identical(
        pnorm(fData(case1_positives)[["scores"]], lower.tail = FALSE),
        fData(case1_positives)[["pvalues"]]
    )

    # Case 2: adjustment factor 5, outlier removal, with prior
    # This runs the same tests as above but with different parameters.
    case2 <- BGScoreTest(demoData,
        adj = 5,
        removeoutlier = TRUE, useprior = TRUE
    )

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case2, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    expect_false("pvalues" %in% colnames(fData(demoData))) # original does not have 'pvalues'
    expect_true("pvalues" %in% colnames(fData(case2))) # new object does have 'pvalues'
    # ...with length same as length of targets.
    expect_true(length(featureNames(case2)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case2_negatives <- case2[which(fData(case2)$CodeClass == "Negative"), ]
    expect_true(all(is.na(fData(case2_negatives)$pvalues)))

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false("scores" %in% colnames(fData(demoData))) # original does not have 'scores'
    expect_true("scores" %in% colnames(fData(case2))) # new object does have 'scores'
    # ...with length same as length of targets.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(all(is.na(fData(case2_negatives)$scores)))

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case2_positives <- case1[fData(case2)[["Negative"]] == FALSE, ] # positive features
    expect_true(all(fData(case2_positives)$pvalues >= 0 | fData(case2_positives)$pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case2)$pvalues))) == nrow(case2_positives))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case2)$scores))) == nrow(case2_positives))

    # 6 The order of pvalues is the same as scores.
    expect_identical(
        pnorm(fData(case2_positives)[["scores"]], lower.tail = FALSE),
        fData(case2_positives)[["pvalues"]]
    )

    # Case 3: adjustment factor 5, outlier removal, without prior
    # This runs the same tests as above but with different parameters.
    case3 <- BGScoreTest(demoData,
        adj = 5,
        removeoutlier = TRUE, useprior = FALSE
    )

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case3, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    expect_false("pvalues" %in% colnames(fData(demoData))) # original does not have 'pvalues'
    expect_true("pvalues" %in% colnames(fData(case3))) # new object does have 'pvalues'
    # ...with length same as length of targets.
    expect_true(length(featureNames(case3)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case3_negatives <- case3[which(fData(case3)$CodeClass == "Negative"), ]
    expect_true(all(is.na(fData(case3_negatives)$pvalues)))

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false("scores" %in% colnames(fData(demoData))) # original does not have 'scores'
    expect_true("scores" %in% colnames(fData(case3))) # new object does have 'scores'
    # ...with length same as length of targets.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(all(is.na(fData(case3_negatives)$scores)))

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case3_positives <- case1[fData(case3)[["Negative"]] == FALSE, ] # positive features
    expect_true(all(fData(case3_positives)$pvalues >= 0 | fData(case3_positives)$pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case3)$pvalues))) == nrow(case3_positives))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case3)$scores))) == nrow(case3_positives))

    # 6 The order of pvalues is the same as scores.
    expect_identical(
        pnorm(fData(case3_positives)[["scores"]], lower.tail = FALSE),
        fData(case3_positives)[["pvalues"]]
    )

    # Case 4: adjustment factor 5, no outlier removal, with prior
    # This runs the same tests as above but with different parameters.
    case4 <- BGScoreTest(demoData,
        adj = 5,
        removeoutlier = FALSE, useprior = TRUE
    )

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case4, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    expect_false("pvalues" %in% colnames(fData(demoData))) # original does not have 'pvalues'
    expect_true("pvalues" %in% colnames(fData(case4))) # new object does have 'pvalues'
    # ...with length same as length of targets.
    expect_true(length(featureNames(case4)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case4_negatives <- case4[which(fData(case4)$CodeClass == "Negative"), ]
    expect_true(all(is.na(fData(case4_negatives)$pvalues)))

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false("scores" %in% colnames(fData(demoData))) # original does not have 'scores'
    expect_true("scores" %in% colnames(fData(case4))) # new object does have 'scores'
    # ...with length same as length of targets.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(all(is.na(fData(case4_negatives)$scores)))

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case4_positives <- case1[fData(case4)[["Negative"]] == FALSE, ] # positive features
    expect_true(all(fData(case4_positives)$pvalues >= 0 | fData(case4_positives)$pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case4)$pvalues))) == nrow(case4_positives))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    expect_true(length(which(!is.na(fData(case4)$scores))) == nrow(case4_positives))

    # 6 The order of pvalues is the same as scores.
    expect_identical(
        pnorm(fData(case4_positives)[["scores"]], lower.tail = FALSE),
        fData(case4_positives)[["pvalues"]]
    )

    # four different settings should yield different results if outliers are present.
    expect_false(identical(fData(case1)[["pvalues"]], fData(case2)[["pvalues"]]))
    expect_false(identical(fData(case1)[["pvalues"]], fData(case3)[["pvalues"]]))
    expect_false(identical(fData(case1)[["pvalues"]], fData(case4)[["pvalues"]]))
    expect_false(identical(fData(case2)[["pvalues"]], fData(case3)[["pvalues"]]))
    expect_false(identical(fData(case2)[["pvalues"]], fData(case4)[["pvalues"]]))
    expect_false(identical(fData(case3)[["pvalues"]], fData(case4)[["pvalues"]]))
})

test_that("BGScoreTest for multiple slides produces desired results", {

    # Desired results occurs when:
    # 1 The function outputs a GeoMx S4 class with p values in featureData with length same as length of targets for each unique id value. The p value is NA for negative probes.
    # 2 The function outputs a GeoMx S4 class with score values in featureData with length same as length of targets for each unique id value. The score value is NA for negative probes.
    # 3 All p values are between 0 and 1 (inclusive) for non-negative features.
    # 4 The order of each column of pvalues is the same as each column of scores for each unique id value.

    # Preamble/load example data
    # Create temporary directory that will get destroyed after this block is executed.
    tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
    # Change to the temporary directory (will set back to getwd() once block is executed.)
    withr::local_dir(tmp_dir)
    data("demoData") # input structure checked above and not repeated here
    # Estimate Poisson background sample-feature factor model for multiple slides
    set.seed(98070)
    demoData <- fitPoisBG(demoData, groupvar = "slide name", size_scale = "sum")
    demoData <- diagPoisBG(demoData, split = TRUE)
    demoData <- aggreprobe(demoData, use = "cor")

    # Case 1: adjustment factor 5, no prior, no outlier removal
    case1 <- BGScoreTest(demoData, split = TRUE, adj = 5, useprior = FALSE, removeoutlier = FALSE)

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case1, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    unique_ids <- unique(pData(demoData)$`slide name`) # pull out the unique ids
    expect_false(all(paste0("pvalues_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("pvalues_", unique_ids) %in% colnames(fData(case1))))
    # ...with length same as length of targets for each unique id value.
    expect_true(length(featureNames(case1)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case1_negatives <- case1[which(fData(case1)$CodeClass == "Negative"), ]
    expect_true(
        all(is.na(fData(case1_negatives)[, grepl("pvalues_", colnames(fData(case1_negatives)))]))
    )

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false(all(paste0("scores_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("scores_", unique_ids) %in% colnames(fData(case1))))
    # ...with length same as length of targets for each unique id value.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(
        all(is.na(fData(case1_negatives)[, grepl("scores_", colnames(fData(case1_negatives)))]))
    )

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case1_positives <- case1[fData(case1)[["Negative"]] == FALSE, ]
    pos_pvalues <- as.numeric(as.matrix(fData(case1_positives)[, grepl("pvalues_", colnames(fData(case1_positives)))]))
    expect_true(all(pos_pvalues >= 0 | pos_pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    compute_length_non_nas_per_column <- function(df) {
        return(
            as.numeric(apply(df, 2, function(x) {
                length(which(!is.na(x)))
            }))
        )
    }
    to_test <- compute_length_non_nas_per_column(
        df = fData(case1)[, grepl("pvalues_", colnames(fData(case1)))]
    )
    expect_true(all(to_test %in% nrow(case1_positives)))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    to_test <- compute_length_non_nas_per_column(
        df = fData(case1)[, grepl("scores_", colnames(fData(case1)))]
    )
    expect_true(all(to_test %in% nrow(case1_positives)))

    # 6 The order of pvalues is the same as scores.
    # This will loop through the different slides (i.e., IDs).
    for (id in unique_ids) {
        expect_identical(
            pnorm(fData(case1_positives)[[paste0("scores_", id)]], lower.tail = FALSE),
            fData(case1_positives)[[paste0("pvalues_", id)]]
        )
    }

    # Case 2: adjustment factor 5, no prior, outlier removal
    case2 <- BGScoreTest(demoData, split = TRUE, adj = 5, useprior = FALSE, removeoutlier = TRUE)

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case2, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    unique_ids <- unique(pData(demoData)$`slide name`) # pull out the unique ids
    expect_false(all(paste0("pvalues_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("pvalues_", unique_ids) %in% colnames(fData(case2))))
    # ...with length same as length of targets for each unique id value.
    expect_true(length(featureNames(case2)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case2_negatives <- case2[which(fData(case2)$CodeClass == "Negative"), ]
    expect_true(
        all(is.na(fData(case2_negatives)[, grepl("pvalues_", colnames(fData(case2_negatives)))]))
    )

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false(all(paste0("scores_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("scores_", unique_ids) %in% colnames(fData(case2))))
    # ...with length same as length of targets for each unique id value.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(
        all(is.na(fData(case2_negatives)[, grepl("scores_", colnames(fData(case2_negatives)))]))
    )

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case2_positives <- case2[fData(case2)[["Negative"]] == FALSE, ]
    pos_pvalues <- as.numeric(as.matrix(fData(case2_positives)[, grepl("pvalues_", colnames(fData(case2_positives)))]))
    expect_true(all(pos_pvalues >= 0 | pos_pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    compute_length_non_nas_per_column <- function(df) {
        return(
            as.numeric(apply(df, 2, function(x) {
                length(which(!is.na(x)))
            }))
        )
    }
    to_test <- compute_length_non_nas_per_column(
        df = fData(case2)[, grepl("pvalues_", colnames(fData(case2)))]
    )
    expect_true(all(to_test %in% nrow(case2_positives)))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    to_test <- compute_length_non_nas_per_column(
        df = fData(case2)[, grepl("scores_", colnames(fData(case2)))]
    )
    expect_true(all(to_test %in% nrow(case2_positives)))

    # 6 The order of pvalues is the same as scores.
    # This will loop through the different slides (i.e., IDs).
    for (id in unique_ids) {
        expect_identical(
            pnorm(fData(case2_positives)[[paste0("scores_", id)]], lower.tail = FALSE),
            fData(case2_positives)[[paste0("pvalues_", id)]]
        )
    }


    # Case 3: adjustment factor 5, with prior, no outlier removal
    case3 <- BGScoreTest(demoData, split = TRUE, adj = 5, useprior = TRUE, removeoutlier = FALSE)

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case3, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    unique_ids <- unique(pData(demoData)$`slide name`) # pull out the unique ids
    expect_false(all(paste0("pvalues_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("pvalues_", unique_ids) %in% colnames(fData(case3))))
    # ...with length same as length of targets for each unique id value.
    expect_true(length(featureNames(case3)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case3_negatives <- case3[which(fData(case3)$CodeClass == "Negative"), ]
    expect_true(
        all(is.na(fData(case3_negatives)[, grepl("pvalues_", colnames(fData(case3_negatives)))]))
    )

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false(all(paste0("scores_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("scores_", unique_ids) %in% colnames(fData(case3))))
    # ...with length same as length of targets for each unique id value.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(
        all(is.na(fData(case3_negatives)[, grepl("scores_", colnames(fData(case3_negatives)))]))
    )

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case3_positives <- case3[fData(case3)[["Negative"]] == FALSE, ]
    pos_pvalues <- as.numeric(as.matrix(fData(case3_positives)[, grepl("pvalues_", colnames(fData(case3_positives)))]))
    expect_true(all(pos_pvalues >= 0 | pos_pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    compute_length_non_nas_per_column <- function(df) {
        return(
            as.numeric(apply(df, 2, function(x) {
                length(which(!is.na(x)))
            }))
        )
    }
    to_test <- compute_length_non_nas_per_column(
        df = fData(case3)[, grepl("pvalues_", colnames(fData(case3)))]
    )
    expect_true(all(to_test %in% nrow(case3_positives)))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    to_test <- compute_length_non_nas_per_column(
        df = fData(case3)[, grepl("scores_", colnames(fData(case3)))]
    )
    expect_true(all(to_test %in% nrow(case3_positives)))

    # 6 The order of pvalues is the same as scores.
    # This will loop through the different slides (i.e., IDs).
    for (id in unique_ids) {
        expect_identical(
            pnorm(fData(case3_positives)[[paste0("scores_", id)]], lower.tail = FALSE),
            fData(case3_positives)[[paste0("pvalues_", id)]]
        )
    }

    # Case 4: adjustment factor 5, no prior, outlier removal
    case4 <- BGScoreTest(demoData, split = TRUE, adj = 5, useprior = TRUE, removeoutlier = TRUE)

    # 1 The function outputs a GeoMx S4 class...
    expect_true(inherits(case4, "NanoStringGeoMxSet"))
    # ...with p values in featureData...
    unique_ids <- unique(pData(demoData)$`slide name`) # pull out the unique ids
    expect_false(all(paste0("pvalues_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("pvalues_", unique_ids) %in% colnames(fData(case4))))
    # ...with length same as length of targets for each unique id value.
    expect_true(length(featureNames(case4)) == length(featureNames(demoData)))
    # The p value is NA for negative probes.
    case4_negatives <- case4[which(fData(case4)$CodeClass == "Negative"), ]
    expect_true(
        all(is.na(fData(case4_negatives)[, grepl("pvalues_", colnames(fData(case4_negatives)))]))
    )

    # 2 The function outputs a GeoMx S4 class...
    # (testing above)...
    # ...with score values in featureData...
    expect_false(all(paste0("scores_", unique_ids) %in% colnames(fData(demoData)))) # original does not have 'pvalues_<id>'
    expect_true(all(paste0("scores_", unique_ids) %in% colnames(fData(case4))))
    # ...with length same as length of targets for each unique id value.
    # (tested above)
    # The score value is NA for negative probes.
    expect_true(
        all(is.na(fData(case4_negatives)[, grepl("scores_", colnames(fData(case4_negatives)))]))
    )

    # 3 All p values are between 0 and 1 (inclusive)
    # for non-negative features.
    case4_positives <- case4[fData(case4)[["Negative"]] == FALSE, ]
    pos_pvalues <- as.numeric(as.matrix(fData(case4_positives)[, grepl("pvalues_", colnames(fData(case4_positives)))]))
    expect_true(all(pos_pvalues >= 0 | pos_pvalues <= 1))

    # 4 The length of non-NA p values is equal to the number of non-negative probes.
    compute_length_non_nas_per_column <- function(df) {
        return(
            as.numeric(apply(df, 2, function(x) {
                length(which(!is.na(x)))
            }))
        )
    }
    to_test <- compute_length_non_nas_per_column(
        df = fData(case4)[, grepl("pvalues_", colnames(fData(case4)))]
    )
    expect_true(all(to_test %in% nrow(case4_positives)))

    # 5 The length of non-NA scores values is equal to the number of non-negative probes.
    to_test <- compute_length_non_nas_per_column(
        df = fData(case4)[, grepl("scores_", colnames(fData(case4)))]
    )
    expect_true(all(to_test %in% nrow(case4_positives)))

    # 6 The order of pvalues is the same as scores.
    # This will loop through the different slides (i.e., IDs).
    for (id in unique_ids) {
        expect_identical(
            pnorm(fData(case4_positives)[[paste0("scores_", id)]], lower.tail = FALSE),
            fData(case4_positives)[[paste0("pvalues_", id)]]
        )
    }

    # four different settings should yield different results if outliers are present.
    pvar_names <- fvarLabels(case1)[grepl("pvalues_", fvarLabels(case1))]

    expect_false(identical(fData(case1)[pvar_names], fData(case2)[pvar_names]))
    expect_false(identical(fData(case1)[pvar_names], fData(case3)[pvar_names]))
    expect_false(identical(fData(case1)[pvar_names], fData(case4)[pvar_names]))
    expect_false(identical(fData(case2)[pvar_names], fData(case3)[pvar_names]))
    expect_false(identical(fData(case2)[pvar_names], fData(case4)[pvar_names]))
    expect_false(identical(fData(case3)[pvar_names], fData(case4)[pvar_names]))
})


## 7 It returns an error without running fitPoisBG.
test_that("It returns an error without running fitPoisBG.", {
    # Preamble/load example data
    data("kidney")
    all0probeidx <- which(rowSums(exprs(kidney))==0)
    kidney <- kidney[-all0probeidx, ]
    kidney <- aggreprobe(kidney, use = "cor")
    expect_error(
        BGScoreTest(kidney),
        "Please run `fitPoisBG` first"
    )
    expect_error(
        BGScoreTest(kidney, split = TRUE),
        "Please run `fitPoisBG` first"
    )
    data("demoData")
    expect_error(
        expect_warning(BGScoreTest(demoData),
                       "No `probenum` is found."),
        "Please run `fitPoisBG` first."
    )
    expect_error(
        expect_warning(BGScoreTest(demoData, split = TRUE),
                       "No `probenum` is found."),
        "Please run `fitPoisBG` first with `groupvar`."
    )
})


## 8 It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.
test_that("It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.", {
    # Preamble/load example data
    data("kidney")
    all0probeidx <- which(rowSums(exprs(kidney))==0)
    kidney <- kidney[-all0probeidx, ]
    kidney <- aggreprobe(kidney, use = "cor")
    res <- fitPoisBG(kidney, size_scale = "first")
    expect_error(
        BGScoreTest(res, split = TRUE),
        "Please run `fitPoisBG` first with `groupvar`"
    )

    # Preamble/load example data
    data("demoData")
    res <- fitPoisBG(demoData, size_scale = "first")
    expect_error(
        expect_warning(BGScoreTest(res, split = TRUE),
                       "No `probenum` is found."),
        "Please run `fitPoisBG` first with `groupvar`."
    )
})
