#' Negative Binomial threshold mixed model for differential expression analysis
#'
#' Negative Binomial threshold mixed model for differential expression analysis
#'
#' @param form model formula
#' @param object count matrix with features in rows and samples in columns
#' @param split indicator variable on whether it is for multiple slides (Yes, TRUE; No, FALSE)
#' @param ROIs_high ROIs with high expressions defined based on featfact and featfact
#' @param features_all full list of features
#' @param sizefact size factor
#' @param sizefact_BG size factor for background
#' @param preci1 precision matrix for regression coefficients
#' @param threshold_mean average background level
#' @param preci2 precision for the background, default=10000
#' @param sizescalebythreshold XXX, default=FALSE
#' @param controlRandom list of random effect control parameters, default=list()
#' @param ... additional argument list that might be used
#'
#' @return a list with parameter estimation
#' #' \itemize{
#'   \item X, design matrix for fixed effect
#'   \item Z, design matrix for random effect
#'   \item rt, random effect terms
#'   \item para0, =NA
#'   \item para, estimated parameters, including regression coefficients, r and threshold in rows and features in columns
#'   \item sizefact, same as input sizefact
#'   \item sizefact0, NA
#'   \item preci1, input precision matrix for regression coefficients
#'   \item Im0, NA
#'   \item Im, Information matrix of parameters
#'   \item conv0, NA
#'   \item conv, vector of convergence, 0 converged, 1 not converged
#'   \item features_high, NA
#'   \item features_all, same as the input features_all
#'   \item theta, list of estimated random effect parameters
#'   \item MAP random effect
#' }
#'
#' @importFrom Biobase sampleNames
#' @importFrom Biobase annotation
#'
#' @examples
#'
#' library(Biobase)
#' library(dplyr)
#' data(demoData)
#' demoData <- demoData[, c(1:5, 33:37)]
#' demoData <- fitPoisBG(demoData, size_scale = "sum")
#' demoData <- aggreprobe(demoData, use = "cor")
#' demoData <- BGScoreTest(demoData)
#' demoData$slidename <- substr(demoData[["slide name"]], 12, 17)
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
#' demoData <- fitNBth(demoData,
#'                     features_high = features_high,
#'                     sizefact_BG = demo_neg$sizefact,
#'                     threshold_start = thmean,
#'                     iterations = 5,
#'                     start_para = c(200, 1),
#'                     lower_sizefact = 0,
#'                     lower_threshold = 100,
#'                     tol = 1e-8)
#' ROIs_high <- sampleNames(demoData)[which(demoData$sizefact_fitNBth * thmean > 2)]
#' features_all <- rownames(demo_pos)
#'
#' pData(demoData)$group <- c(rep(1, 5), rep(2, 5))
#'
#'
#' NBthDEmod2 <- fitNBthDE(form = ~group,
#'                      split = FALSE,
#'                      object = demoData,
#'                      ROIs_high = ROIs_high,
#'                      features_high = features_high,
#'                      features_all = features_all,
#'                      sizefact_start = demoData[, ROIs_high][['sizefact_fitNBth']],
#'                      sizefact_BG = demoData[, ROIs_high][['sizefact']],
#'                      threshold_mean = notes(demoData)[["threshold"]],
#'                      preci2=10000,
#'                      prior_type="contrast",
#'                      covrob=FALSE,
#'                      preci1con=1/25,
#'                      sizescalebythreshold=TRUE)
#'
#' set.seed(123)
#' NBthmDEmod1 <- fitNBthmDE(
#'     form = ~ group + (1 | `slide name`),
#'     split = FALSE,
#'     object = demoData,
#'     ROIs_high = ROIs_high,
#'     features_all = features_all[1:5],
#'     sizefact = demoData[, ROIs_high][["sizefact_fitNBth"]],
#'     sizefact_BG = demoData[, ROIs_high][["sizefact"]],
#'     preci1=NBthDEmod2$preci1,
#'     threshold_mean = thmean,
#'     preci2=10000,
#'     sizescale = TRUE,
#'     controlRandom=list(nu=12, nmh_e=400, thin_e=60))
#'
#' @export
#' @docType methods
#' @rdname fitNBthmDE-methods
setGeneric("fitNBthmDE",
    signature = c("object"),
    function(object, ...) standardGeneric("fitNBthmDE")
)

#' @rdname fitNBthmDE-methods
#' @aliases fitNBthmDE,NanoStringGeoMxSet-method
setMethod(
    "fitNBthmDE", "NanoStringGeoMxSet",
    function(object, form, split, ROIs_high = NULL,
    features_all = NULL, sizefact = NULL, sizefact_BG = NULL,
    preci1, threshold_mean = NULL,
    preci2=10000, sizescalebythreshold = FALSE, controlRandom = list()) {
        fDat <- Biobase::fData(object)
        pDat <- Biobase::pData(object)

        posdat <- object[-which(fDat$CodeClass == "Negative"), ]
        countmat <- Biobase::exprs(posdat)

        fDatNeg <- fDat[which(fDat$CodeClass == "Negative"), ]

        # only calculate backmean if any of the three params are missing
        if (any(c(is.null(sizefact_BG), is.null(sizefact), is.null(ROIs_high), is.null(threshold_mean)))) {
            if (isFALSE(split)) {
                # single slide
                if (!("sizefact" %in% varLabels(object))) {
                    stop("Please run `fitPoisBG` first.")
                } else {
                    # calculate the backmean for WTA or CTA data
                    featfact_mean <- mean(fDatNeg[["featfact"]])
                }
            } else {
                # multiple slides
                if (!("sizefact_sp" %in% varLabels(object))) {
                    stop("Please run `fitPoisBG` first with `groupvar`.")
                } else {
                    # calculate the backmean for WTA or CTA data
                    featfact_mean <- colMeans(fDatNeg[, grep("featfact_", fvarLabels(object))])[1]
                }
            }
        }



        # calculate probenum for the dataset
        if ("probenum" %in% fvarLabels(posdat)) {
            probenum <- fData(posdat)[["probenum"]]
        } else {
            stop("No `probenum` is found. Run `aggreprobe` first.")
        }
        names(probenum) <- rownames(fData(posdat))

        # extract annot from object
        annot <- Biobase::pData(object)

        # setting default value for ROIs_high
        if (is.null(ROIs_high)) {
            if (!("sizefact_fitNBth" %in% varLabels(object))) {
                stop("Please run `fitNBth` first.")
            } else {
                # estimate values for ROIs_high
                ROIs_high <- Biobase::sampleNames(object)[which((quantile(Biobase::fData(object)[["para"]][, 1], probs = 0.90, na.rm = TRUE) -
                    Biobase::notes(object)[["threshold"]]) * object$sizefact_fitNBth > 2)]
            }
        }

        # setting default value for sizefact_BG
        if (is.null(sizefact_BG)) {
            if (isFALSE(split)) {
                # single slide
                sizefact_BG <- pDat[ROIs_high, "sizefact"]
            } else {
                # multiple slides
                sizefact_BG <- pDat[ROIs_high, ][["sizefact_sp"]]
            }
            names(sizefact_BG) <- rownames(pDat[ROIs_high, ])
        }


        # setting default value for sizefact
        if (is.null(sizefact)) {
            if (!("sizefact_fitNBth" %in% colnames(pDat))) {
                stop("Please run `fitNBth` first.")
            } else {
                sizefact <- pDat[ROIs_high, ][["sizefact_fitNBth"]]
            }
            names(sizefact) <- rownames(pDat[ROIs_high, ])
        }

        # setting default value for features_high
        if (is.null(features_all)) {
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
            features_all <- features_high[seq_len(5)]
        }

        result <- fitNBthmDE(
            form = form,
            annot = annot[ROIs_high, ],
            object = countmat[, ROIs_high],
            probenum = probenum,
            features_all = features_all,
            sizefact = sizefact,
            sizefact_BG = sizefact_BG,
            preci1 = preci1,
            threshold_mean = threshold_mean,
            preci2 = preci2,
            sizescalebythreshold = sizescalebythreshold,
            controlRandom = controlRandom
        )
    }
)

#' Negative Binomial threshold mixed model for differential expression analysis
#'
#' Negative Binomial threshold mixed model for differential expression analysis
#'
#' @param form model formula
#' @param annot annotations files with variables in the formula
#' @param object count matrix with features in rows and samples in columns
#' @param probenum a vector of numbers of probes in each gene, default = rep(1, NROW(object))
#' @param features_all vector of all features to be run
#' @param sizefact size factor
#' @param sizefact_BG size factor for background
#' @param preci1 precision matrix for regression coefficients
#' @param threshold_mean average background level
#' @param preci2 precision for the background, default=10000
#' @param sizescalebythreshold whether to scale the size factor, default=TRUE
#' @param controlRandom list of random effect control parameters
#'
#' @importFrom lme4 lFormula
#'
#' @return a list with parameter estimation
#' #' \itemize{
#'   \item X, design matrix for fixed effect
#'   \item Z, design matrix for random effect
#'   \item rt, random effect terms
#'   \item para0, =NA
#'   \item para, estimated parameters, including regression coefficients, r and threshold in rows and features in columns
#'   \item sizefact, same as input sizefact
#'   \item sizefact0, NA
#'   \item preci1, input precision matrix for regression coefficients
#'   \item Im0, NA
#'   \item Im, Information matrix of parameters
#'   \item conv0, NA
#'   \item conv, vector of convergence, 0 converged, 1 not converged
#'   \item features_high, NA
#'   \item features_all, same as the input features_all
#'   \item theta, list of estimated random effect parameters(for relative covariance matrix)
#'   \item varcov, list of estimated variance covariance parameter estimation
#'   \item MAP random effect
#' }
#'
#' @rdname fitNBthmDE-methods
#' @aliases fitNBthmDE,matrix-method
setMethod(
    "fitNBthmDE", "matrix",
    function(form, annot, object, probenum = rep(1, NROW(object)),
    features_all, sizefact, sizefact_BG, preci1, threshold_mean = NULL,
    preci2=10000, sizescalebythreshold = TRUE, controlRandom = list()) {
        if (is.null(names(probenum))) names(probenum) <- rownames(object)

        n_feature <- length(features_all)

        cRandom <- list(
            nmh_s = 40, nmh_e = 200,
            thin_s = 8, thin_e = 20, useprior = TRUE, thetapri = 1,
            nu = NA, lower = NA, upper = NA, iterations = 300,
            lower0 = 0.01, upper0 = 20, preciu0 = 1
        )


        # Options for random effect
        cRandomNames <- names(cRandom)
        cRandom[(controlN <- names(controlRandom))] <- controlRandom
        if (length(unkwn <- controlN[!controlN %in% cRandomNames])) {
            warning("Unknown names in control: ", paste(unkwn, collapse = ", "))
        }


        annot$fake <- 1
        if (length(form) == 2) {
            form[[3]] <- form[[2]]
        }
        form[[2]] <- as.name("fake")
        resu <- lme4::lFormula(form, data = annot)

        X <- resu$X

        rt <- resu$reTrms
        Z <- t(as.matrix(rt$Zt))
        Lambdati <- rt$Lambdat
        mapping <- function(theta) theta[rt$Lind]

        if (missing(preci1)) {
            preci <- diag(1, ncol(X))
        }

        rl <- nrow(Lambdati)
        if (is.na(cRandom$preciu0)) {
            cRandom$preciu0 <- 1
        }



        cRandom$preciu <- diag(cRandom$preciu0, rl)


        lower <- rt$lower
        lower[lower == 0] <- cRandom$lower0
        cRandom$lower <- lower

        upper <- -rt$lower
        upper[upper == 0] <- cRandom$upper0
        cRandom$upper <- upper

        message(rl)
        message("----------")
        if (is.na(cRandom$nu)) {
            nu <- 4 + rl
            cRandom$nu <- nu
        }

        cluster_size <- sum(rt$Lind == 1)
        temp_size <- nrow(Lambdati) / cluster_size

        para <- matrix(0, nrow = (ncol(X) + 2), ncol = n_feature)
        theta <- matrix(0, nrow = max(rt$Lind), ncol = n_feature)
        varcov <- theta
        colnames(theta) <- features_all
        colnames(varcov) <- features_all
        colnames(para) <- features_all

        Im <- list()
        nmh_sq <- floor(seq(cRandom$nmh_s, cRandom$nmh_e, length.out = cRandom$iterations))
        thin_sq <- floor(seq(cRandom$thin_s, cRandom$thin_e, length.out = cRandom$iterations))
        cRandom$nmh_sq <- nmh_sq
        cRandom$thin_sq <- thin_sq


        Uvec <- matrix(0, ncol(Z), n_feature)
        colnames(Uvec) <- features_all
        conv <- numeric(n_feature)
        names(conv) <- features_all
        for (feature in features_all) {
            NBthmmodfeat <- fitNBthmDEfeat(
                t(object[feature, ]), probenum[feature], X, Z, sizefact, sizefact_BG, preci1, threshold_mean, preci2, Lambdati,
                mapping, cluster_size, temp_size, rl, rt, sizescalebythreshold, cRandom
            )
            para[, feature] <- NBthmmodfeat$para_fix
            theta[, feature] <- NBthmmodfeat$theta
            varcov[, feature] <- NBthmmodfeat$varcov
            Im[[feature]] <- NBthmmodfeat$Im
            Uvec[, feature] <- NBthmmodfeat$Uvec
            conv[feature] <- NBthmmodfeat$conv
            message(feature)
            message("----------")
        }
        names(Im) <- features_all
        names(theta) <- features_all

        paraname <- c(colnames(X), c("r", "threshold"))
        rownames(para) <- paraname

        return(list(
            X = X,
            Z = Z,
            rt = rt,
            para0 = NA,
            para = para,
            sizefact = sizefact,
            sizefact0 = NA,
            preci1 = preci1,
            conv0 = NA,
            conv = conv,
            Im0 = NA,
            Im = Im,
            features_high = NA,
            features_all = features_all,
            theta = theta,
            varcov = varcov,
            Uvec = Uvec
        ))
    }
)
