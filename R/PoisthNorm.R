 #' Poisson threshold model based normalization-log2 transformation for single slide or for multiple slides
#'
#' Poisson threshold model based normalization-log2 transformation for single slide or for multiple slides
#'
#' @param object a valid GeoMx S4 object
#' @param split indicator variable on whether it is for multiple slides (Yes, TRUE; No, FALSE)
#' @param ROIs_high ROIs with high expressions defined based on featfact and featfact
#' @param features_high subset of features which are well above the background
#' @param features_all full feature vector to apply the normalization on
#' @param sizefact_start initial value for size factors
#' @param sizefact_BG size factor for background
#' @param threshold_mean average threshold level
#' @param preci2 precision for threshold, default=10000
#' @param iterations iteration number, default=2,
#'                   the first iteration using the features_high to construct the prior for parameters then refit the model on all features.
#'                   precision matrix for threshold: preci2
#' @param prior_type prior type for preci1, "equal" or "contrast", default="contrast"
#' @param sizefactrec XXXX, default = TRUE
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param sizescalebythreshold XXXX, default = FALSE
#' @param covrob whether to use robust covariance in calculating the prior precision matrix 1, default = FALSE
#' @param preci1con The user input constant term in specifying precision matrix 1, default=1/25
#' @param cutoff term in calculating precision matrix 1, default=15
#' @param confac The user input factor for contrast in precision matrix 1, default=1
#' @param calhes The user input whether to calculate hessian: calhes, default=FALSE
#' @param ... additional argument list that might be used
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase assayDataElement
#'
#'
#' @return if split is FALSE, a valid GeoMx S4 object including the following items
#' \itemize{
#'   \item para0_norm, matrix of estimated parameters for iter=1, features in columns and parameters(log2 expression, threshold) in rows, in featureData.
#'   \item para_norm, matrix of estimated parameters for iter=2, features in columns and parameters(log2 expression, threshold) in rows, in featureData.
#'   \item normmat0, matrix of log2 expression for iter=1, features in columns and log2 expression in rows, in assay slot.
#'   \item normmat, matrix of log2 expression for iter=2, features in columns and log2 expression in rows, in assay lot.
#'   \item sizefact_norm, estimated sizefact, in phenoData.
#'   \item sizefact0_norm, estimated sizefact in iter=1, in phenoData.
#'   \item preci1, precision matrix 1, in experimentData.
#'   \item conv0, vector of convergence for iter=1, 0 converged, 1 not converged, in featureData
#'   \item conv, vector of convergence for iter=2, 0 converged, 1 not converged, in featureData
#'   \item features_high, same as the input features_high, in featureData
#'   \item features_all, same as the input features_all, in featureData
#' }
#'
#' if split is TRUE, a valid GeoMx S4 object with the following items appended.
#' \itemize{
#'   \item threshold0, matrix of estimated threshold for iter=1, features in columns and threshold for different slides in rows, in featureData.
#'   \item threshold, matrix of estimated threshold for iter=2, features in columns and threshold for different slides in rows, in featureData.
#'   \item normmat0_sp, matrix of log2 expression for iter=1, features in columns and log2 expression in rows, in assay slot.
#'   \item normmat_sp, matrix of log2 expression for iter=2, features in columns and log2 expression in rows, in assay slot.
#'   \item sizefact_norm_sp, estimated sizefact, in phenoData
#'   \item sizefact0_norm_sp, estimated sizefact in iter=1, in phenoData
#'   \item preci1, precision matrix 1, in experimentData
#'   \item conv0_sp_XX, vector of convergence for each unique slide value for iter=1, 0 converged, 1 not converged, in featureData for each unique slide.
#'   \item conv_sp_XX, vector of convergence for each unique slide value for iter=2, 0 converged, 1 not converged, in featureData for each unique slide.
#'   \item features_high_sp, same as the input features_high, in featureData.
#'   \item features_all_sp, same as the input features_all, in featureData.
#' }
#'
#' @importFrom Biobase sampleNames
#' @importFrom Biobase featureNames
#' @importFrom Biobase annotation
#'
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
#' ROIs_high <- sampleNames(demoData)[which((quantile(fData(demoData)[["para"]][, 1],
#'                                                    probs = 0.90, na.rm = TRUE) -
#'          notes(demoData)[["threshold"]]) * demoData$sizefact_fitNBth > 2)]
#' features_all <- rownames(demo_pos)
#' thmean <- mean(fData(demo_neg)[["featfact"]])
#' demoData <- fitPoisthNorm(
#'     object = demoData,
#'     split = FALSE,
#'     ROIs_high = ROIs_high,
#'     features_high = features_high,
#'     features_all = features_all,
#'     sizefact_start = demoData[, ROIs_high][["sizefact_fitNBth"]],
#'     sizefact_BG = demoData[, ROIs_high][["sizefact"]],
#'     threshold_mean = thmean,
#'     preci2 = 10000,
#'     prior_type = "contrast",
#'     covrob = FALSE,
#'     preci1con = 1 / 25
#' )
#'
#' @export
#' @docType methods
#' @rdname fitPoisthNorm-methods

setGeneric("fitPoisthNorm",
    signature = c("object"),
    function(object, ...) standardGeneric("fitPoisthNorm")
)

#' @rdname fitPoisthNorm-methods
#' @aliases fitPoisthNorm,NanoStringGeoMxSet-method
setMethod(
    "fitPoisthNorm", "NanoStringGeoMxSet",
    function(object, split = FALSE, ROIs_high = NULL, features_high = NULL,
    features_all = NULL, sizefact_start = NULL, sizefact_BG = NULL,
    threshold_mean = NULL, preci2=10000, iterations = 2, prior_type = c("contrast", "equal"),
    sizefactrec = TRUE, size_scale = c("sum", "first"), sizescalebythreshold = FALSE,
    covrob = FALSE, preci1con = 1 / 25, cutoff = 15, confac = 1, calhes = FALSE) {
        # calculate backmean
        # setting default values for sizefact_BG
        fDat <- Biobase::fData(object)
        fDatNeg <- fDat[which(fDat$CodeClass == "Negative"), ]

        # only calculate backmean if any of the three params are missing
        if (any(c(is.null(sizefact_BG), is.null(sizefact_start), is.null(ROIs_high), is.null(threshold_mean)))) {
            if (isFALSE(split)) {
                # single slide
                if (!("sizefact" %in% varLabels(object))) {
                    stop("Please run `fitPoisBG` first.")
                } else {
                    # calculate the backmean for WTA or CTA data
                    thmean <- mean(fDatNeg[["featfact"]])
                }
            } else {
                # multiple slides
                if (!("sizefact_sp" %in% varLabels(object))) {
                    stop("Please run `fitPoisBG` first with `groupvar`.")
                } else {
                    # calculate the backmean for WTA or CTA data
                    thmean <- colMeans(fDatNeg[, grep("featfact_", fvarLabels(object))])[1]
                }
            }
        }



        # setting default value for ROIs_high
        if (is.null(ROIs_high)) {
            if (!("sizefact_fitNBth" %in% varLabels(object))) {
                stop("Please run `fitNBth` first.")
            } else {
                ROIs_high <- Biobase::sampleNames(object)[which((quantile(fData(object)[["para"]][, 1],
                    probs = 0.90, na.rm = TRUE
                ) - notes(object)[["threshold"]]) * object$sizefact_fitNBth > 2)]
            }
        }

        object_high <- object[, ROIs_high]

        fDat <- Biobase::fData(object_high)
        pDat <- Biobase::pData(object_high)

        posdat <- object_high[-which(fDat$CodeClass == "Negative"), ]
        countmat <- Biobase::exprs(posdat)

        # calculate probenum for the dataset
        if ("probenum" %in% fvarLabels(posdat)) {
            probenum <- fData(posdat)[["probenum"]]
        } else {
            stop("No `probenum` is found. Run `aggreprobe` first.")
        }
        names(probenum) <- rownames(fData(posdat))

        # setting default value for sizefact_BG
        if (is.null(sizefact_BG)) {
            if (isFALSE(split)) {
                # single slide
                sizefact_BG <- pDat[["sizefact"]]
            } else {
                # multiple slides
                sizefact_BG <- pDat[["sizefact_sp"]]
            }
        }


        # setting default value for sizefact_start
        if (is.null(sizefact_start)) {
            if (!("sizefact_fitNBth" %in% colnames(pDat))) {
                stop("Please run `fitNBth` first.")
            } else {
                sizefact_start <- pDat[["sizefact_fitNBth"]]
            }
        }


        # setting default value for features_all
        if (is.null(features_all)) {
            gene_sum <- rowSums(countmat)

            if (any(grepl("WTA", toupper(Biobase::annotation(object))))) {
                features_all <- rownames(countmat)
            } else if (any(grepl("CTA", toupper(Biobase::annotation(object))))) {
                features_all <- rownames(countmat)
            } else {
                stop("No information is found to determine the data type (CTA or WTA).")
            }
        }

        # setting default value for sizefact_start
        if (is.null(threshold_mean)) {
            threshold_mean <- thmean
        }

        # setting default value for features_high
        if (is.null(features_high)) {
            gene_sum <- rowSums(countmat)

            if (any(grepl("WTA", toupper(Biobase::annotation(object))))) {
                features_high <- names(which(((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95)))))
                features_high <- sort(sample(features_high, 1500))
            } else if (any(grepl("CTA", toupper(Biobase::annotation(object))))) {
                if (any(grepl("scores", fvarLabels(object)))) {
                    if (split == TRUE) {
                        sc1_scores <- fData(object)[-which(fData(object)$Negative), grepl("scores_", fvarLabels(object))]
                        rownames(sc1_scores) <- fData(object)[-which(fData(object)$Negative), "TargetName"]
                        features_high <- apply(sc1_scores, 2, function(x) {
                            ((x > quantile(x, probs = 0.4)) & (x < quantile(x, probs = 0.95)))
                        })
                        features_high <- names(which(apply(features_high, 1, all)))
                    } else {
                        sc1_scores <- fData(object)[-which(fData(object)$Negative), "scores"]
                        names(sc1_scores) <- fData(object)[-which(fData(object)$Negative), "TargetName"]
                        features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) & (sc1_scores < quantile(sc1_scores, probs = 0.95)))
                        features_high <- names(which(features_high))
                    }
                } else {
                    stop("Please run score test first using `BGScoreTest`.")
                }
            } else {
                stop("No information is found to determine the data type (CTA or WTA).")
            }
        }


        if (isFALSE(split)) {
            result <- fitPoisthNorm(
                object = countmat,
                probenum = probenum,
                features_high = features_high,
                features_all = features_all,
                sizefact_start = sizefact_start,
                sizefact_BG = sizefact_BG,
                threshold_mean = threshold_mean,
                preci2 = preci2,
                iterations = iterations,
                prior_type = prior_type,
                sizefactrec = sizefactrec,
                size_scale = size_scale,
                sizescalebythreshold = sizescalebythreshold,
                covrob = covrob,
                preci1con = preci1con,
                cutoff = cutoff,
                confac = confac,
                calhes = calhes
            )

            # para0
            Biobase::fData(object)[["para0_norm"]] <- matrix(NA,
                nrow = nrow(object), ncol = nrow(result$para0),
                dimnames = list(Biobase::featureNames(object), paste0("var", seq_len((nrow(result$para0)))))
            )
            Biobase::fData(object)[["para0_norm"]][colnames(result$para0), ] <- t(result$para0)

            # para
            Biobase::fData(object)[["para_norm"]] <- matrix(NA,
                nrow = nrow(object), ncol = nrow(result$para),
                dimnames = list(Biobase::featureNames(object), paste0("var", seq_len((nrow(result$para)))))
            )
            Biobase::fData(object)[["para_norm"]][colnames(result$para), ] <- t(result$para)

            # normmat0
            normmat0_mat <- matrix(NA,
                nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object)
            )
            normmat0_mat[colnames(result$normmat0), ROIs_high] <- t(result$normmat0)
            Biobase::assayDataElement(object, "normmat0") <- normmat0_mat

            # normmat
            normmat_mat <- matrix(NA,
                nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object)
            )
            normmat_mat[colnames(result$normmat), ROIs_high] <- t(result$normmat)
            Biobase::assayDataElement(object, "normmat") <- normmat_mat

            # sizefact
            Biobase::pData(object)[ROIs_high, "sizefact_norm"] <- result$sizefact
            # sizefac0
            Biobase::pData(object)[ROIs_high, "sizefact0_norm"] <- result$sizefact0

            # preci1
            Biobase::notes(object)$preci1_norm <- result$preci1

            # skipping appending Im0 and Im from result for now.

            # conv0
            Biobase::fData(object)[features_high, "conv0"] <- result$conv0
            # conv
            Biobase::fData(object)[features_all, "conv"] <- result$conv

            # features_high
            Biobase::fData(object)[features_high, "features_high"] <- 1
            # features_all
            Biobase::fData(object)[features_all, "features_all"] <- 1
        } else {
            if (is.null(Biobase::notes(object)$fitPoisBG_sp_var)) {
                stop("Please run `fitPoisBG` first with `groupvar`.")
            } else {
                idvar <- Biobase::notes(object)$fitPoisBG_sp_var
                id <- Biobase::pData(object)[[idvar]][match(ROIs_high, colnames(object))]
                message(sprintf("The results are based on stored `groupvar`, %s", idvar))
            }

            result <- fitPoisthNorm_sp(
                object = countmat,
                probenum = probenum,
                features_high = features_high,
                features_all = features_all,
                sizefact_start = sizefact_start,
                sizefact_BG = sizefact_BG,
                threshold_mean = threshold_mean,
                preci2 = preci2,
                id = id,
                iterations = iterations,
                prior_type = prior_type,
                sizefactrec = sizefactrec,
                size_scale = size_scale,
                sizescalebythreshold = sizescalebythreshold,
                covrob = covrob,
                preci1con = preci1con,
                cutoff = cutoff,
                confac = confac
            )

            # threshold0
            Biobase::fData(object)[["threshold0"]] <- matrix(NA,
                nrow = nrow(object), ncol = nrow(result$threshold0),
                dimnames = list(Biobase::featureNames(object), rownames(result$threshold0))
            )
            Biobase::fData(object)[["threshold0"]][colnames(result$threshold0), ] <- t(result$threshold0)

            # threshold
            Biobase::fData(object)[["threshold"]] <- matrix(NA,
                nrow = nrow(object), ncol = nrow(result$threshold),
                dimnames = list(Biobase::featureNames(object), rownames(result$threshold))
            )
            Biobase::fData(object)[["threshold"]][colnames(result$threshold), ] <- t(result$threshold)

            # normmat0
            normmat0_mat <- matrix(NA,
                nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object)
            )
            normmat0_mat[colnames(result$normmat0), ROIs_high] <- t(result$normmat0)
            Biobase::assayDataElement(object, "normmat0_sp") <- normmat0_mat

            # normmat
            normmat_mat <- matrix(NA,
                nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object)
            )
            normmat_mat[colnames(result$normmat), ROIs_high] <- t(result$normmat)
            Biobase::assayDataElement(object, "normmat_sp") <- normmat_mat

            # sizefact
            Biobase::pData(object)[ROIs_high, "sizefact_norm_sp"] <- result$sizefact
            # sizefact0
            Biobase::pData(object)[ROIs_high, "sizefact0_norm_sp"] <- result$sizefact0

            # preci1
            Biobase::notes(object)$preci1_norm_sp <- result$preci1

            # skipping appending Im0 and Im from result for now.

            # conv0
            for (index in names(result$conv0)) {
                Biobase::fData(object)[features_high, paste0("conv0_sp_", index)] <- result$conv0[[index]]
            }
            # conv
            for (index in names(result$conv0)) {
                Biobase::fData(object)[features_all, paste0("conv_sp_", index)] <- result$conv[[index]]
            }
            # features_high
            Biobase::fData(object)[features_high, "features_high_sp"] <- 1
            # features_all
            Biobase::fData(object)[features_all, "features_all_sp"] <- 1
        }

        return(object)
    }
)

#'
#' Poisson model based normalization and log2 transformation with threshold
#'
#' @param object count matrix with features in rows and samples in columns
#' @param probenum a vector of numbers of probes in each gene
#' @param features_high subset of features which are well above the background
#' @param features_all full feature vector to apply the normalization on
#' @param sizefact_start initial value for size factors
#' @param sizefact_BG size factor for background
#' @param threshold_mean average threshold level
#' @param preci2 precision for threshold, default=10000
#' @param iterations iteration number, default=2,
#'                   the first iteration using the features_high to construct the prior for parameters then refit the model on all features.
#'                   precision matrix for threshold: preci2
#' @param prior_type prior type for preci1, "equal" or "contrast", default="contrast"
#' @param sizefactrec XXXX, default = TRUE
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param sizescalebythreshold XXXX, default = FALSE
#' @param covrob whether to use robust covariance in calculating the prior precision matrix 1, default = FALSE
#' @param preci1con The user input constant term in specifying precision matrix 1, default=1/25
#' @param cutoff term in calculating precision matrix 1, default=15
#' @param confac The user input factor for contrast in precision matrix 1, default=1
#' @param calhes The user input whether to calculate hessian: calhes, default=FALSE
#'
#' @importFrom robust covRob
#'
#' @return a list of following items
#' \itemize{
#'   \item para0, matrix of estimated parameters for iter=1, features in columns and parameters(log2 expression, threshold) in rows.
#'   \item para, matrix of estimated parameters for iter=2, features in columns and parameters(log2 expression, threshold) in rows.
#'   \item normmat0, matrix of log2 expression for iter=1, features in columns and log2 expression in rows.
#'   \item normmat, matrix of log2 expression for iter=2, features in columns and log2 expression in rows.
#'   \item sizefact, estimated sizefact
#'   \item sizefact0, estimated sizefact in iter=1
#'   \item preci1, precision matrix 1
#'   \item Im0, Information matrix of parameters in iter=1
#'   \item Im, Information matrix of parameters in iter=2
#'   \item conv0, vector of convergence for iter=1, 0 converged, 1 not converged
#'   \item conv, vector of convergence for iter=2, 0 converged, 1 not converged
#'   \item features_high, same as the input features_high
#'   \item features_all, same as the input features_all
#' }
#'
#' @rdname fitPoisthNorm-methods
#' @aliases fitPoisthNorm,matrix-method

setMethod(
    "fitPoisthNorm", "matrix",
    function(object, probenum = rep(1, NROW(object)), features_high, features_all,
    sizefact_start, sizefact_BG, threshold_mean, preci2=10000, iterations = 2,
    prior_type = c("contrast", "equal"), sizefactrec = TRUE, size_scale = c("sum", "first"),
    sizescalebythreshold = FALSE, covrob = FALSE, preci1con = 1 / 25, cutoff = 15, confac = 1, calhes = FALSE) {
        if (iterations != 2) {
            stop("Only iterations=2 is allowed")
        }

        if (is.null(names(probenum))) names(probenum) <- rownames(object)

        sizefact0 <- sizefact <- sizefact_start


        n_para <- ncol(object)

        X <- diag(1, n_para)

        prior_type <- match.arg(prior_type)


        if (prior_type == "equal") {
            preci1 <- preci1con * t(X) %*%  diag(1, n_para) %*% X
        } else if (prior_type == "contrast") {
            contrvec <- t(rep(1 / n_para, n_para)) %*% X

            preci1 <- (preci1con) * t(contrvec) %*% contrvec
        }



        for (iter in seq_len(iterations)) {
            if (iter == 1) {
                modfit <- PoisthNorm_paraOptall(t(object[features_high, ]), X, sizefact_BG, sizefact, preci1, threshold_mean * probenum, preci2, sizescalebythreshold, calhes & (iter == iterations))
                Im <- modfit$hes
                Im0 <- Im
                para <- modfit$par
                colnames(para) <- features_high
                para0 <- para
                conv0 <- modfit$conv
            } else {
                modfit <- PoisthNorm_paraOptall(t(object[features_all, ]), X, sizefact_BG, sizefact, preci1, threshold_mean * probenum, preci2, sizescalebythreshold, calhes & (iter == iterations))
                Im <- modfit$hes

                para <- modfit$par
                colnames(para) <- features_all
                conv <- modfit$conv
            }


            message("probe finished")

            if (iter == 1) {
                features_remain <- names(which(colMeans(abs(para[seq_len(n_para), , drop = FALSE])) < cutoff))
                if (prior_type == "equal") {
                    if (covrob) {
                        cov_mat <- robust::covRob(t(para[seq_len(n_para), features_remain]), na.action = na.omit)$cov
                    } else {
                        cov_mat <- cov(t(para[seq_len(n_para), features_remain]), use = "pairwise.complete.obs")
                    }

                    preci1 <- solve(cov_mat)
                } else if (prior_type == "contrast") {
                    contrmat <- cbind(rep(1, n_para - 1), -diag(1, (n_para - 1)))
                    para_EB <- para[seq_len(n_para), features_remain]

                    contrpara <- contrmat %*% para_EB


                    if (covrob) {
                        cov_mat <- robust::covRob(t(contrpara), na.action = na.omit)$cov
                    } else {
                        cov_mat <- cov(t(contrpara), use = "pairwise.complete.obs")
                    }




                    preci1 <- t(contrmat) %*% solve(cov_mat) %*% contrmat

                    #   avevec <- rep(1, n_para)/n_para

                    preci1 <- confac * preci1 + (preci1con) * t(contrvec) %*% contrvec
                }
            }

            if (sizefactrec) {
                size_scale <- match.arg(size_scale)
                features_remain <- names(which(colMeans(abs(para[seq_len(n_para), , drop = FALSE])) < cutoff))




                for (i in seq_len(length(sizefact))) {
                    fun <- PoisthNorm_scalenll(X[i, ], object[features_remain, i], probenum[features_remain], para[seq_len(n_para), features_remain], sizefact_BG[i], para[n_para + 1, features_remain], sizescalebythreshold, threshold_mean)
                    sizefact[i] <- optim(c(sizefact[i]), fun, lower = c(0), method = "L-BFGS-B")$par
                }
                if (size_scale == "first") {
                    scale_fac <- sizefact[1]
                } else if (size_scale == "sum") {
                    scale_fac <- sum(sizefact)
                }

                sizefact <- sizefact / scale_fac

                message(sprintf("Iteration = %s, squared error = %s", iter, sum((sizefact - sizefact0)^2)))

                if (iter == 1) {
                    sizefact0 <- sizefact
                }
            }
        }
        normmat0 <- X %*% para0[seq_len(n_para), ]
        normmat <- X %*% para[seq_len(n_para), ]
        message("Model converged.")
        return(list(
            para0 = para0,
            para = para,
            normmat0 = normmat0,
            normmat = normmat,
            sizefact = sizefact,
            sizefact0 = sizefact0,
            preci1 = preci1,
            Im0 = Im0,
            Im = Im,
            conv0 = conv0,
            conv = conv,
            features_high = features_high,
            features_all = features_all
        ))
    }
)


#' Poisson threshold model based normalization-log2 transformation for multiple slides
#'
#' Poisson threshold model based normalization-log2 transformation for multiple slides
#'
#' @param object count matrix with features in rows and samples in columns
#' @param probenum a vector of numbers of probes in each gene
#' @param features_high subset of features which are well above the background
#' @param features_all full feature vector to apply the normalization on
#' @param sizefact_start initial value for size factors
#' @param sizefact_BG size factor for background
#' @param threshold_mean average threshold level
#' @param preci2 precision for threshold, default=10000
#' @param id character vector of slide name of each sample
#' @param iterations iteration number, default=2,
#'                   the first iteration using the features_high to construct the prior for parameters then refit the model on all features.
#'                   precision matrix for threshold: preci2
#' @param prior_type prior type for preci1, "equal" or "contrast", default="contrast"
#' @param sizefactrec XXXX, default = TRUE
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param sizescalebythreshold XXXX, default = FALSE
#' @param covrob whether to use robust covariance in calculating the prior precision matrix 1, default = FALSE
#' @param preci1con The user input constant term in specifying precision matrix 1, default=1/25
#' @param cutoff term in calculating precision matrix 1, default=15
#' @param confac The user input factor for contrast in precision matrix 1, default=1
#' @param ... additional argument list that might be used
#'
#' @return a list of following items
#' \itemize{
#'   \item threshold0, matrix of estimated threshold for iter=1, features in columns and threshold for different slides in rows.
#'   \item threshold, matrix of estimated threshold for iter=2, features in columns and threshold for different slides in rows.
#'   \item normmat0, matrix of log2 expression for iter=1, features in columns and log2 expression in rows.
#'   \item normmat, matrix of log2 expression for iter=2, features in columns and log2 expression in rows.
#'   \item sizefact, estimated sizefact
#'   \item sizefact0, estimated sizefact in iter=1
#'   \item preci1, precision matrix 1
#'   \item Im0, Information matrix in iter=1
#'   \item Im, Information matrix in iter=2
#'   \item conv0, vector of convergence for iter=1, 0 converged, 1 not converged
#'   \item conv, vector of convergence for iter=2, 0 converged, 1 not converged
#'   \item features_high, same as the input features_high
#'   \item features_all, same as the input features_all
#' }
#' @docType methods
#' @rdname fitPoisthNorm_sp-methods

setGeneric("fitPoisthNorm_sp",
    signature = c("object"),
    function(object, ...) standardGeneric("fitPoisthNorm_sp")
)

#' @rdname fitPoisthNorm_sp-methods
#' @aliases fitPoisthNorm_sp,matrix-method
setMethod(
    "fitPoisthNorm_sp", "matrix",
    function(object, probenum, features_high,
    features_all = colnames(object), sizefact_start, sizefact_BG,
    threshold_mean, preci2=10000, id, iterations = 2, prior_type = c("contrast", "equal"),
    sizefactrec = TRUE, size_scale = c("sum", "first"), sizescalebythreshold = FALSE,
    covrob = FALSE, preci1con = 1 / 25, cutoff = 15, confac = 1) {
        uniid <- unique(as.character(id))

        loc <- lapply(uniid, function(x) which(id == x))
        names(loc) <- uniid
        normmod_ls <- lapply(uniid, function(x) {
            fitPoisthNorm(object[, id == x],
                probenum = probenum, features_high, features_all,
                sizefact_start[id == x], sizefact_BG[id == x], threshold_mean, preci2, iterations, prior_type = prior_type, sizefactrec = sizefactrec, size_scale = size_scale, sizescalebythreshold = sizescalebythreshold,
                covrob = covrob, preci1con = preci1con, cutoff = cutoff, confac = confac
            )
        })
        names(normmod_ls) <- uniid

        preci1 <- list()
        threshold0 <- matrix(NA, length(uniid), length(features_high))
        threshold <- matrix(NA, length(uniid), length(features_all))
        rownames(threshold) <- rownames(threshold0) <- uniid
        normmat0 <- matrix(NA, length(id), length(features_high))
        normmat <- matrix(NA, length(id), length(features_all))

        colnames(normmat) <- colnames(threshold) <- features_all
        colnames(normmat0) <- colnames(threshold0) <- features_high

        sizefact <- sizefact0 <- sizefact_BG
        conv0 <- list()
        conv <- list()
        Im0 <- list()
        Im <- list()

        for (idname in uniid) {
            threshold0[idname, ] <- normmod_ls[[idname]]$para0[length(loc[[idname]]) + 1, ]
            threshold[idname, ] <- normmod_ls[[idname]]$para[length(loc[[idname]]) + 1, ]

            normmat0[loc[[idname]], ] <- normmod_ls[[idname]]$normmat0
            normmat[loc[[idname]], ] <- normmod_ls[[idname]]$normmat
            sizefact[loc[[idname]]] <- normmod_ls[[idname]]$sizefact
            sizefact0[loc[[idname]]] <- normmod_ls[[idname]]$sizefact0

            Im0[[idname]] <- normmod_ls[[idname]]$Im0
            Im[[idname]] <- normmod_ls[[idname]]$Im


            conv0[[idname]] <- normmod_ls[[idname]]$conv0
            conv[[idname]] <- normmod_ls[[idname]]$conv

            preci1[[idname]] <- normmod_ls[[idname]]$preci1
        }

        return(list(
            threshold0 = threshold0,
            threshold = threshold,
            normmat0 = normmat0,
            normmat = normmat,
            sizefact = sizefact,
            sizefact0 = sizefact0,
            preci1 = preci1,
            Im0 = Im0,
            Im = Im,
            conv0 = conv0,
            conv = conv,
            features_high = features_high,
            features_all = features_all
        ))
    }
)
