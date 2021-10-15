#' Negative Binomial threshold model
#'
#' Estimate the signal size factor for features above the background
#'
#' @param object a valid GeoMx S4 object
#' @param features_high subset of features which are well above the background
#' @param sizefact_BG size factors for the background
#' @param sizefact_start initial value for size factors
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param threshold_start initial value for threshold
#' @param threshold_fix whether to fix the threshold, default=FALSE
#' @param tol tolerance to determine convergence, default=1e-3
#' @param iterations maximum iterations to be run, default=5
#' @param start_para starting values for parameter estimation, default=c(threshold_start, 1)
#' @param lower_sizefact lower limit for sizefact, default=0
#' @param lower_threshold lower limit for threshold
#' @param ... additional argument list that might be used
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase experimentData
#' @importFrom Biobase featureNames
#' @importFrom Biobase annotation
#'
#'
#' @return a valid GeoMx S4 object
#' \itemize{
#'   \item para0 = "NA", in experimentData
#'   \item para, estimated parameters, "signal" "r" in rows and features in columns, in featureData
#'   \item sizefact, estimated size factor, in phenoData
#'   \item preci1 = "NA", in experimentData
#'   \item conv0 = "NA", in experimentData
#'   \item conv = "NA", in experimentData
#'   \item Im = "NA", in experimentData
#'   \item features_high, a vector of indicators, in featureData (0: No; 1: Yes; NA: not included in features_high)
#'   \item features_all = "NA", in experimentData
#'   \item threshold, estimated threshold, when threshold_fix, equals to threshold_start, in experimentData
#' }
#'
#' @export
#' @docType methods
#' @rdname fitNBth-methods
#' @examples
#'
#' library(Biobase)
#' library(dplyr)
#' data(demoData)
#' demoData <- fitPoisBG(demoData, size_scale = "sum")
#' demoData <- aggreprobe(demoData, use = "cor")
#' demoData <- BGScoreTest(demoData)
#' thmean <- 1 * mean(fData(demoData)$featfact, na.rm = TRUE)
#' demo_pos <- demoData[which(!fData(demoData)$CodeClass == "Negative"), ]
#' demo_neg <- demoData[which(fData(demoData)$CodeClass == "Negative"), ]
#' sc1_scores <- fData(demo_pos)[, "scores"]
#' names(sc1_scores) <- fData(demo_pos)[, "TargetName"]
#' features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) &
#'    (sc1_scores < quantile(sc1_scores, probs = 0.95))) |>
#'     which() |>
#'     names()
#' set.seed(123)
#' features_high <- sample(features_high, 100)
#' demoData <- fitNBth(demoData,
#'                     features_high = features_high,
#'                     sizefact_BG = demo_neg$sizefact,
#'                     threshold_start = thmean,
#'                     iterations = 5,
#'                     start_para = c(200, 1),
#'                     lower_sizefact = 0,
#'                     lower_threshold = 100,
#'                     tol = 1e-8)
#'

setGeneric("fitNBth",
    signature = c("object"),
    function(object, ...) standardGeneric("fitNBth")
)

#' @rdname fitNBth-methods
#' @aliases fitNBth,NanoStringGeoMxSet-method
setMethod(
    "fitNBth", "NanoStringGeoMxSet",
    function(object,
    features_high = NULL,
    sizefact_BG = NULL, sizefact_start = sizefact_BG,
    size_scale = c("sum", "first"),
    threshold_start = NULL, threshold_fix = FALSE,
    tol = 1e-3, iterations = 5,
    start_para = c(threshold_start, 0.5),
    lower_sizefact = 1e-2, lower_threshold = 1e-2) {
        # check on tol
        tol <- as.double(tol)
        stopifnot(length(tol) == 1)
        stopifnot(tol >= 0)


        fDat <- Biobase::fData(object)
        pDat <- Biobase::pData(object)

        posdat <- object[-which(fDat$CodeClass == "Negative"), ]
        countmat <- Biobase::exprs(posdat)

        fDatNeg <- fDat[which(fDat$CodeClass == "Negative"), ]

        # calculate probenum for the dataset
        if ("probenum" %in% fvarLabels(posdat)) {
            probenum <- fData(posdat)[["probenum"]]
        } else {
            stop("No `probenum` is found. Run `aggreprobe` first.")
        }
        names(probenum) <- rownames(fData(posdat))

        # only calculate thmean if any of the two params are missing
        if (any(c(is.null(sizefact_BG), is.null(threshold_start)))) {
            if (isFALSE(split)) {
                # single slide
                if (!("sizefact" %in% varLabels(object))) {
                    stop("Please run `fitPoisBG` first.")
                } else {
                    # calculate the thmean for WTA or CTA data
                    thmean <- mean(fDatNeg[["featfact"]])
                }
            } else {
                # multiple slides
                if (!("sizefact_sp" %in% varLabels(object))) {
                    stop("Please run `fitPoisBG` first with `groupvar`.")
                } else {
                    # calculate the thmean for WTA or CTA data
                    thmean <- colMeans(fDatNeg[, grep("featfact_", fvarLabels(object))])[1]
                }
            }
        }

        # setting default value for sizefact_BG
        if (is.null(sizefact_BG)) {
            if (isFALSE(split)) {
                # single slide
                sizefact_BG <- pDat[["sizefact"]]
            } else {
                # multiple slides
                sizefact_BG <- pDat[["sizefact_sp"]]
            }
            sizefact_start <- sizefact_BG
        }

        # setting default value for sizefact_start
        if (is.null(threshold_start)) {
            threshold_start <- thmean
        }

        # setting default value for features_high
        if (is.null(features_high)) {
            gene_sum <- rowSums(countmat)

            if (any(grepl("WTA", toupper(Biobase::annotation(object))))) {
                features_high <- names(which(((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95)))))
                features_high <- sample(features_high, 1500)
            } else if (any(grepl("CTA", toupper(Biobase::annotation(object))))) {
                if ( !any(grepl("scores", colnames(fDat))) ) {
                    stop("Please run `BGScoreTest` first. If you run `BGScoreTest` before, please specify `split = TRUE` for multiple slides.")
                } else {
                    if ( any(grepl("scores_", colnames(fDat))) ){
                        # fit the model with multiple slides
                        sc1_scores <- fData(posdat)[, grep("scores_", fvarLabels(posdat))]
                        rownames(sc1_scores) <- fData(posdat)[, "TargetName"]
                        features_high <- apply(sc1_scores, 2, function(x){
                            ((x > quantile(x, probs = 0.4)) & (x < quantile(x, probs = 0.95)))
                        })
                        features_high <- names(which(apply(features_high, 1, all)))

                    } else {
                        sc1_scores <- fData(posdat)[, "scores"]
                        names(sc1_scores) <- fData(posdat)[, "TargetName"]
                        features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) & (sc1_scores < quantile(sc1_scores, probs = 0.95)))
                        features_high <- names(which(features_high))

                    }
                }

            } else {
                stop("No information is found to determine the data type (CTA or WTA).")
            }
        }

        result <- fitNBth(
            object = countmat[features_high, ],
            probenum = probenum,
            features_high = features_high,
            sizefact_BG = sizefact_BG,
            sizefact_start = sizefact_start,
            size_scale = size_scale,
            threshold_start = threshold_start,
            threshold_fix = threshold_fix,
            tol = tol,
            iterations = iterations,
            start_para = start_para,
            lower_sizefact = lower_sizefact,
            lower_threshold = lower_threshold
        )

        # update para0
        Biobase::notes(object)$para0 <- ifelse(is.na(result$para0), "NA", result$para0)

        # update para
        Biobase::fData(object)[["para"]] <- matrix(NA,
            nrow = nrow(object), ncol = 2,
            dimnames = list(featureNames(object), c("signal", "r"))
        )
        Biobase::fData(object)[["para"]][colnames(result$para), ] <- t(result$para)

        # update sizefact
        object[["sizefact_fitNBth"]] <- result$sizefact

        # update preci1
        Biobase::notes(object)$preci1 <- ifelse(is.na(result$preci1), "NA", result$preci1)

        # update conv0
        Biobase::notes(object)$conv0 <- ifelse(is.na(result$conv0), "NA", result$conv0)

        # update conv
        Biobase::notes(object)$conv <- ifelse(is.na(result$conv), "NA", result$conv)

        # update Im
        Biobase::notes(object)$Im <- ifelse(is.na(result$Im), "NA", result$Im)

        # update features_high
        Biobase::fData(object)[["feature_high_fitNBth"]] <- 0
        Biobase::fData(object)[["feature_high_fitNBth"]][match(result$features_high, Biobase::featureNames(object), nomatch = 0)] <- 1

        # update features_all
        Biobase::notes(object)$features_all <- ifelse(is.na(result$features_all), "NA", result$features_all)

        # update threshold
        Biobase::notes(object)$threshold <- ifelse(is.na(result$threshold), "NA", result$threshold)

        return(object)
    }
)

#' Negative Binomial threshold model
#'
#' Estimate the signal size factor for features above the background
#'
#' @param object count matrix with features in rows and samples in columns
#' @param probenum a vector of numbers of probes in each gene
#' @param features_high subset of features which are well above the background
#' @param sizefact_BG size factors for the background
#' @param sizefact_start initial value for size factors
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param threshold_start initial value for threshold
#' @param threshold_fix whether to fix the threshold, default=FALSE
#' @param tol tolerance to determine convergence, default=1e-3
#' @param iterations maximum iterations to be run, default=5
#' @param start_para starting values for parameter estimation, default=c(threshold_start, 1)
#' @param lower_sizefact lower limit for sizefact, default=0
#' @param lower_threshold lower limit for threshold
#'
#' @return a list of following items, some items are place holders = NA
#' \itemize{
#'   \item para0 = NA,
#'   \item para, estimated parameters, "signal" "r" in rows and features in columns
#'   \item sizefact, estimated size factor
#'   \item preci1 = NA
#'   \item conv0 = NA
#'   \item conv = NA
#'   \item Im = NA
#'   \item features_high = features_high
#'   \item features_all = NA
#'   \item threshold, estimated threshold, when threshold_fix, equals to threshold_start
#' }
#'
#' @rdname fitNBth-methods
#' @aliases fitNBth,matrix-method
setMethod(
    "fitNBth", "matrix",
    function(object, features_high, probenum, sizefact_BG, sizefact_start = sizefact_BG, size_scale = c("sum", "first"), threshold_start, threshold_fix = FALSE, tol = 1e-3, iterations = 5,
    start_para = c(threshold_start, 1), lower_sizefact = 0, lower_threshold = threshold_start / 5) {
        size_scale <- match.arg(size_scale)
        sizefact0 <- sizefact <- sizefact_start
        threshold <- threshold_start
        # mat <- matrix(1, nrow(object), 1)
        if (is.null(names(probenum))) names(probenum) <- rownames(object)
        for (iter in seq_len(iterations)) {
            para <- NBth_paraopt(object[features_high, ], probenum[features_high], sizefact, sizefact_BG, threshold, start = start_para)
            # result <- mleprobeNBall(object[,features_high], mat, sizefact0, sizefact,
            #                         matrix(0,1,1), threshold, 0,
            #                         c(rep(0,ncol(mat)), 1, threshold), 0)
            features_NA <- features_high[unique(which(is.na(para), arr.ind = TRUE)[, 2])]
            features_remain <- setdiff(features_high, features_NA)


            for (i in seq_len(length(sizefact))) {
                fun <- NBth_scalenll(object[features_remain, i], probenum[features_remain], t(para[1, features_remain]), t(para[2, features_remain]), sizefact_BG[i], threshold)
                sizefact[i] <- optim(c(sizefact_start[i]), fun, lower = c(lower_sizefact), method = "L-BFGS-B")$par
            }

            if (size_scale == "first") {
                scale_fac <- sizefact[1]
            } else if (size_scale == "sum") {
                scale_fac <- sum(sizefact)
            }

            sizefact <- sizefact / scale_fac


            if (!threshold_fix) {
                fun1 <- NBth_thnll(object[features_remain, ], probenum[features_remain], sizefact, sizefact_BG, scale_fac * t(para[1, features_remain]), t(para[2, features_remain]))

                threshold <- optim(c(threshold_start), fun1, lower = c(lower_threshold), method = "L-BFGS-B")$par
            }




            message(sprintf("Iteration = %s, squared error = %s", iter, sum((sizefact - sizefact0)^2)))

            if (sum((sizefact - sizefact0)^2) < tol) break

            sizefact0 <- sizefact
        }
        message("Model converged.")

        rownames(para) <- c("signal", "r")

        return(list(
            para0 = NA,
            para = para,
            sizefact = sizefact,
            preci1 = NA,
            conv0 = NA,
            conv = NA,
            Im = NA,
            features_high = features_high,
            features_all = NA,
            threshold = threshold
        ))
    }
)
