#### Specs for aggreprobe:
#1. The function shall aggregate the probes depending on the argument provided.
#2. The function returns a GeoMxSet object when given a GeoMxSet object as input.
#3. The expression matrix and the target feature data available prior to collapsing shall match after collapsing except TargetName.
#4. TargetName for negative probes shall be updated to probe IDs after collapsing.
#5. For the non-negative probes, the subset of probes selected by the use method will be aggregated by sum into one target count.
#6. The resulting object shall have the same size of feature names as negative probe names plus non-negative target names.

# load example data:
data("demoData")
demoNeg <- demoData[which(fData(demoData)$CodeClass == "Negative"), ]

#Create background model
demoData <- fitPoisBG(demoData, size_scale = "sum")

# Perform aggregations
scoreAgg <- aggreprobe(demoData, split=FALSE, use="score")
corAgg <- aggreprobe(demoData, split=FALSE, use="cor")
bothAgg <- aggreprobe(demoData, split=FALSE, use="both")

# Spec 1: The function shall aggregate the probes depending on the argument provided.
test_that("Aggreprobe function aggregates the probes by use method", {
    expect_false(all(exprs(scoreAgg) == exprs(corAgg)))
    expect_false(all(exprs(corAgg) == exprs(bothAgg)))
    expect_false(all(exprs(scoreAgg) == exprs(bothAgg)))
})

aggdObjs <- list(score=scoreAgg, cor=corAgg, both=bothAgg)
for (used in names(aggdObjs)) {
    aggd <- aggdObjs[[used]]
    negAggd <- aggd[which(fData(aggd)$CodeClass == "Negative"), ]
    negAggd <- negAggd[featureNames(demoNeg), ]
    # Spec 2: The function returns a GeoMxSet object when given a GeoMxSet
    #         object as input.
    test_that(paste0(used, ": aggreprobe function aggregates by use method"), {
        expect_true(inherits(aggd, "NanoStringGeoMxSet"))
    })
    # Spec 3: For negative probes, the expression matrix and the target
    #         feature data available prior to collapsing shall match after
    #         collapsing except TargetName.
    test_that(paste0(used, ": neg probes are not collapsed"), {
        targetLabels <- intersect(fvarLabels(demoNeg), fvarLabels(negAggd))
        targetLabels <- targetLabels[targetLabels != "TargetName"]
        expect_equal(exprs(negAggd), exprs(demoNeg))
        expect_equal(
            fData(negAggd)[, targetLabels], fData(demoNeg)[, targetLabels])
    })
    # Spec 4: TargetName for negative probes shall be updated to probe IDs
    #         after collapsing.
    test_that(paste0(used, ": neg probe target names replaced by probe IDs"), {
        expect_equal(fData(negAggd)[["TargetName"]], fData(demoNeg)[["RTS_ID"]])
    })
    # Spec 5: For the non-negative probes, selected by the use method will be
    #         aggregated by sum into one target count.
    test_that(paste0(used, "aggregation is by sum"), {
        endoAggd <- aggd[which(fData(aggd)$CodeClass != "Negative"), ]
        testNames <- sample(featureNames(endoAggd), 50)
        for (testName in testNames) {
            probeNames <- unlist(fData(endoAggd)[testName, "proberemained"])
            expectedCount <-
                apply(exprs(demoData)[probeNames, , drop=FALSE], 2, sum)
            expect_equal(exprs(endoAggd)[testName, ], expectedCount)
        }
    })
    # Spec 6: The resulting object shall have the same size of feature names
    #         as negative probe names plus non-negative target names.
    test_that(paste0(used, "feature length is correct"), {
        endoObj <- demoData[which(fData(demoData)$CodeClass != "Negative"), ]
        featLen <-
            length(unique(fData(endoObj)[["TargetName"]])) + dim(demoNeg)[1L]
        expect_equal(dim(aggd)[1L], featLen)
    })
    # Spec 7: Single probe targets shall be returned without aggregation
    test_that(paste0(used, "WTA data (single probe) is not aggregated"), {
        data("kidney")
        all0probeidx <- which(rowSums(exprs(kidney))==0)
        kidney <- kidney[-all0probeidx, ]
        kidney <- fitPoisBG(kidney, size_scale = "sum")
        aggdKid <- aggreprobe(kidney, split=FALSE, use = used)
        endoKid <- kidney[which(fData(kidney)$CodeClass != "Negative"), ]
        negKid <- kidney[which(fData(kidney)$CodeClass == "Negative"), ]
        kidLen <-
            length(unique(fData(endoKid)[["TargetName"]])) + dim(negKid)[1L]
        expect_equal(dim(aggdKid)[1L], kidLen)
        endoAggdKid <- aggdKid[which(fData(aggdKid)$CodeClass != "Negative"), ]
        expect_true(all(exprs(endoAggdKid)[fData(endoKid)$TargetName, ] == exprs(endoKid)))
    })
}
