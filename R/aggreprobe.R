#' Generate aggregated counts of probes for the same target
#'
#' Generate Generate aggregated counts of probes for the same target, based on their score test results or correlation
#'
#' @param object a valid GeoMx S4 object
#' @param split indicator variable on whether it is for multiple slides (Yes, TRUE; No, FALSE)
#' @param use the method to determine outliers including score, cor, and both
#' @param corcutoff the cutoff value for correlation, default value: 0.85
#' @param ... additional argument list that might be used
#'
#' @return
#' \itemize{
#'   \item remain, the list of remaining probes of targets
#'   \item probenum, numerical vector of probe numbers of targets
#'   \item featuremat, the matrix of features
#' }
#'
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase experimentData
#' @importFrom Biobase phenoData
#' @importFrom Biobase protocolData
#' @importFrom Biobase annotation
#' @importFrom GeomxTools featureType
#' @importClassesFrom  GeomxTools NanoStringGeoMxSet
#'
#' @examples
#' data("demoData")
#' demoData <- aggreprobe(demoData, use = "cor")
#' @export
#' @docType methods
#' @rdname aggreprobe-methods
#'

setGeneric("aggreprobe",
    signature = c("object"),
    function(object, ...) standardGeneric("aggreprobe")
)

#' @rdname aggreprobe-methods
#' @aliases aggreprobe,NanoStringGeoMxSet-method
setMethod(
    "aggreprobe", "NanoStringGeoMxSet",
    function(object, split, use = c("score", "cor", "both"), corcutoff = 0.85, ...) {
        if(GeomxTools::featureType(object) == "Target") {
            stop("GeoMxSet object feature type is already target-level. ",
                 "No further aggregation can be performed.")
        }
        object_neg <- object[which(Biobase::fData(object)$CodeClass == "Negative"), ]
        countmat_neg <- Biobase::exprs(object_neg)
        object_nonneg <- object[which(Biobase::fData(object)$CodeClass != "Negative"), ]
        countmat_nonneg <- Biobase::exprs(object_nonneg)
        pDat <- Biobase::pData(object_neg)
        fDat <- Biobase::fData(object_neg)

        use <- match.arg(use)
        if (use == "score" | use == "both") {
            if (isFALSE(split)) {
                if (!any(c("sizefact" %in% colnames(pDat), "featfact" %in% colnames(fDat)))) {
                    stop("No sizefact and featfact is found. Please run `fitPoisBG` first. If you run `fitPoisBG` with `groupvar` before, please specify `split = TRUE` for `aggreprobe`.")
                }
                sizefact <- setNames(pDat[["sizefact"]], rownames(pDat))
                featfact <- setNames(fDat[["featfact"]], rownames(fDat))

                BGmod <- list(
                    sizefact = sizefact,
                    featfact = featfact,
                    countmat = countmat_nonneg
                )
            } else {
                if (!any(c("sizefact_sp" %in% colnames(pDat), "featfact_" %in% colnames(fDat)))) {
                    stop("No sizefact and featfact for multiple slides is found. Please run `fitPoisBG` first with `groupvar`.")
                }
                sizefact <- setNames(pDat[["sizefact_sp"]], rownames(pDat))
                featfact <- fDat[, grep("featfact_", colnames(fDat))]
                colnames(featfact) <- gsub("featfact_", "", colnames(featfact))
                idvar <- Biobase::notes(object)[["fitPoisBG_sp_var"]]
                id <- pDat[[idvar]]
                message(sprintf("The results are based on stored `groupvar`, %s", idvar))

                BGmod <- list(
                    sizefact = sizefact,
                    featfact = as.matrix(featfact),
                    countmat = countmat_nonneg,
                    id = id
                )
            }
        } else {
            BGmod <- NULL
        }

        result <- aggreprobe(
            object = countmat_nonneg,
            probenames = rownames(countmat_nonneg),
            featurenames = Biobase::fData(object_nonneg)$TargetName,
            negmod = BGmod,
            use = use,
            corcutoff = corcutoff, 
            ...
        )

        # assemble assay data for negative and non-negative
        # assay data for negative data
        targetCounts_neg <- do.call(rbind, NanoStringNCTools::esBy(object_neg,
            GROUP = "TargetName",
            FUN = function(x) {
                Biobase::esApply(x, 2, identity)
            }, simplify = FALSE
        ))
        # assay data for non-negative data
        targetCounts_nonneg <- result$featuremat

        # separate feature data by negative and non-negative
        # feature data for non-negative
        targetFeats <- Biobase::fData(object_nonneg)
        targetFeats <-
            targetFeats[!duplicated(targetFeats[["TargetName"]]), ]
        rownames(targetFeats) <- targetFeats[, "TargetName"]
        probeColumns <- c("RTS_ID", "QCFlags", "ProbeID")
        targetFeats_nonneg <-
            targetFeats[, !colnames(targetFeats) %in% probeColumns]
        targetFeats_nonneg$probenum <- result$probenum[match(rownames(targetFeats), names(result$probenum))]
        targetFeats_nonneg$proberemained <- cbind(result$remain)[match(rownames(targetFeats), names(result$remain)), ]

        # feature data for negative probes
        targetFeats <- Biobase::fData(object_neg)
        rownames(targetFeats) <- targetFeats[, "RTS_ID"]
        targetFeats[, "TargetName"] <- targetFeats[, "RTS_ID"]
        probeColumns <- c("RTS_ID", "QCFlags", "ProbeID")
        targetFeats_neg <-
            targetFeats[, !colnames(targetFeats) %in% probeColumns]
        targetFeats_neg$probenum <- 1
        targetFeats_neg$proberemained <- targetFeats_neg$TargetName

        # combine feature data
        targetFeats <-
            Biobase::AnnotatedDataFrame(rbind(
                targetFeats_nonneg[rownames(targetCounts_nonneg), ],
                targetFeats_neg[rownames(targetCounts_neg), ]
            ),
            dimLabels = c("featureNames", "featureColumns")
            )

        targetObject <- GeomxTools::NanoStringGeoMxSet(
            assayData = rbind(targetCounts_nonneg, targetCounts_neg),
            phenoData = Biobase::phenoData(object),
            featureData = targetFeats,
            experimentData = Biobase::experimentData(object),
            annotation = Biobase::annotation(object),
            protocolData = Biobase::protocolData(object),
            featureType = "Target",
            check = FALSE
        )
    }
)

#' Generate aggregated counts of probes for the same target
#'
#' Generate Generate aggregated counts of probes for the same target, based on their score test results or correlation
#'
#' @param object matrix of probes
#' @param probenames vector of names of probe
#' @param featurenames vector of names of features each probe corresponding to
#' @param negmod Poisson Background model object for negative probes
#' @param use the method to determine outliers including score, cor, and both
#' @param corcutoff the cutoff value for correlation 
#' @param ... additional argument list that might be used
#'
#' @return
#' \itemize{
#'   \item remain, the list of remaining probes of targets
#'   \item probenum, numerical vector of probe numbers of targets
#'   \item featuremat, the matrix of features
#' }
#'
#' @rdname aggreprobe-methods
#' @aliases aggreprobe,matrix-method
#'
setMethod(
    "aggreprobe", "matrix",
    function(object, probenames, featurenames, negmod, use = c("score", "cor", "both"), corcutoff=0.85, ...) {
        use <- match.arg(use)

        # Stop when there are probes with all 0 counts
        if(any(rowSums(object[probenames,])==0))
            stop("There are all 0 probes in the count matrix, remove them and rerun aggreprobe.")


        # select rows with probenames
        object <- object[probenames,]

        probemat_ls <- lapply(split(probenames, featurenames), function(x) object[x, ])

        remain <- split(probenames, featurenames)
        feature_useall <- names(which(sapply(remain, length) <= 2))
        feature_toagre <- names(which(sapply(remain, length) > 2))
        ## outlier by score test
        if (use == "score" | use == "both") {
            scores <- BGScoreTest(object, negmod, ...)$scores
            if (NCOL(scores) > 1) {
                scores <- apply(scores, 1, median)
            }
            scores_ls <- split(scores, featurenames)
            scores_ls <- scores_ls[feature_toagre]
            remainscores <- lapply(scores_ls, function(x) names(which(abs(x - median(x)) <= 1.96 * 1.25 * mad(x))))
            remainscores <- c(remain[feature_useall], remainscores)
        }
        ## outlier by correlation
        if (use == "cor" | use == "both") {
            cor_ls <- lapply(probemat_ls, function(x) cor(t(x)))[feature_toagre]
            remaincor <- remain[feature_useall]
            remaincortemp <- remain[feature_toagre]


            mean_cor_ls <- lapply(names(cor_ls), function(x) {
                apply(cor_ls[[x]][remaincortemp[[x]], remaincortemp[[x]]], 2, function(z) (sum(z) - 1) / (length(z) - 1))
            })
            names(mean_cor_ls) <- names(cor_ls)
            remaincortemp <- lapply(mean_cor_ls, function(x) {
                if (all(x > corcutoff)) y <- names(x) else y <- setdiff(names(x), names(which.min(x)))
                y
            })
            remaincor <- c(remaincor, remaincortemp[names(which(sapply(remaincortemp, length) <= 2))])
            remaincortemp <- remaincortemp[names(which(sapply(remaincortemp, length) > 2))]


            cor_ls <- cor_ls[names(remaincortemp)]


            remaincor <- c(remaincor, remaincortemp)
        }
        if (use == "both") {
            remainscores <- remainscores[unique(featurenames)]
            remaincor <- remaincor[unique(featurenames)]
            remain <- mapply(function(x, y) intersect(x, y), remaincor, remainscores)
        } else if (use == "score") {
            remain <- remainscores[unique(featurenames)]
        } else {
            remain <- remaincor[unique(featurenames)]
        }

        probenum <- sapply(remain, length)
        featuremat <- t(sapply(remain, function(x) colSums(object[x, , drop = FALSE])))
        rownames(featuremat) <- names(remain)
        return(list(
            remain = remain,
            probenum = probenum,
            featuremat = featuremat
        ))
    }
)
