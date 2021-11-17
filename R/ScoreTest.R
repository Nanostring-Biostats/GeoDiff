#library(GeoDiff)
library(magrittr)
library(GeomxTools)
library(Rcpp)
library(Rfast)
library(Matrix)
#' Testing for features above the background
#'
#' Testing for features above the background using Poisson background model as reference
#'
#' @param object a valid GeoMx S4 object with featfact and sizefact
#' @param split indicator variable on whether it is for multiple slides (Yes, TRUE; No, FALSE)
#' @param adj adjustment factor for the number of probes in each gene, default =1 i.e.
#'            each target only consists of one probe
#' @param removeoutlier whether to remove outlier
#' @param useprior whether to use the prior that the expression level of background follows a Beta distribution, l
#'                 eading to a more conservative test
#' @param ... additional argument list that might be used
#'
#' @return a valid GeoMx S4 object including the following items
#' \itemize{
#'   \item pvalues - Background score test pvalues, in featureData
#'   \item scores - Background score test statistics, in featureData
#' }
#'
#' if split is TRUE, a valid GeoMx S4 object including the following items
#' \itemize{
#'   \item pvalues_XX - Background score test pvalues vector, column name (denoted as XX) the same as slide names, in featureData
#'   \item scores_XX - Background score test statistics vector, column name (denoted as XX) the same as slide names, in featureData
#' }
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase varLabels
#' @importFrom Biobase fvarLabels
#' @importFrom Biobase notes
#' @importFrom Biobase sampleNames
#' @importFrom Biobase featureNames
#' @importFrom Biobase annotation
#'
#'
#'
#'
#' @examples
#'
#' data(demoData)
#' demoData <- fitPoisBG(demoData, size_scale = "sum")
#' demoData <- aggreprobe(demoData, use = "cor")
#' demoData <- BGScoreTest(demoData, adj = 1, useprior = FALSE)
#' demoData <- fitPoisBG(demoData, size_scale = "sum", groupvar = "slide name")
#' demoData <- BGScoreTest(demoData, adj = 1, useprior = TRUE, split = TRUE)
#'
#' @export
#' @docType methods
#' @rdname BGScoreTest-methods

setGeneric("BGScoreTest",
           signature = c("object"),
           function(object, ...) standardGeneric("BGScoreTest")
)

#' @rdname BGScoreTest-methods
#' @aliases BGScoreTest,NanoStringGeoMxSet-method
setMethod(
  "BGScoreTest", "NanoStringGeoMxSet",
  function(object, split = FALSE, adj = 1, removeoutlier = FALSE, useprior = FALSE) {
    posdat <- object[-which(Biobase::fData(object)$CodeClass == "Negative"), ]
    countmat <- Biobase::exprs(posdat)
    
    pDat <- Biobase::pData(object)
    fDat <- Biobase::fData(object)
    
    # calculate probenum for the dataset
    if ("probenum" %in% fvarLabels(posdat)) {
      probenum <- fData(posdat)[["probenum"]]
    } else {
      warning("No `probenum` is found. For targets with >1 probe, this ",
              "is allowed in order to run `aggreprobe` with `use=\"score\"`")
      probenum <- rep(1, nrow(posdat))
    }
    names(probenum) <- rownames(fData(posdat))
    
    if (isFALSE(split)) {
      if (!any(c("sizefact" %in% colnames(pDat), "featfact" %in% colnames(fDat)))) {
        stop("Please run `fitPoisBG` first. If you run `fitPoisBG` before, please specify `split = TRUE`.")
      }
      
      sizefact <- setNames(object[["sizefact"]], Biobase::sampleNames(object))
      featfact <- setNames(fDat[["featfact"]], Biobase::featureNames(object))
      featfact <- featfact[-which(is.na(featfact))]
      
      BGmod <- list(
        sizefact = sizefact,
        featfact = featfact,
        countmat = countmat
      )
      
      result <- BGScoreTest(
        object = countmat,
        BGmod = BGmod,
        probenum = probenum,
        adj = adj,
        removeoutlier = removeoutlier,
        useprior = useprior
      )
      
      if (any(c("pvalues", "scores") %in% Biobase::varLabels(object))) {
        warning("`pvalues` and `scores` exist in the phenodata. Those values are replaced.")
      }
      
      Biobase::fData(object)[["pvalues"]] <- NA
      Biobase::fData(object)[["pvalues"]][match(names(result$pvalues), Biobase::featureNames(object), nomatch = 0)] <- result$pvalues
      Biobase::fData(object)[["scores"]] <- NA
      Biobase::fData(object)[["scores"]][match(names(result$scores), Biobase::featureNames(object), nomatch = 0)] <- result$scores
    } else {
      if (!any(c("sizefact_sp" %in% colnames(pDat), "featfact_" %in% colnames(fDat)))) {
        stop("Please run `fitPoisBG` first with `groupvar`.")
      }
      
      idvar <- Biobase::notes(object)[["fitPoisBG_sp_var"]]
      id <- pDat[[idvar]]
      message(sprintf("The results are based on stored `groupvar`, %s", idvar))
      
      sizefact <- setNames(pDat[["sizefact_sp"]], Biobase::sampleNames(object))
      featfact <- fDat[, paste0("featfact_", unique(id))]
      featfact <- featfact[which(fDat$CodeClass == "Negative"), ]
      colnames(featfact) <- gsub("featfact_", "", colnames(featfact))
      
      BGmod <- list(
        sizefact = sizefact,
        featfact = as.matrix(featfact),
        countmat = countmat,
        id = id
      )
      
      result <- BGScoreTest_sp(
        object = countmat,
        BGmod = BGmod,
        probenum = probenum,
        adj = adj,
        removeoutlier = removeoutlier,
        useprior = useprior
      )
      
      if (length(c(grep("pvalues_", Biobase::fvarLabels(object)), grep("scores_", Biobase::fvarLabels(object)))) > 0) {
        warning("`pvalues_sp` and `scores_sp` exist in the phenodata. Those values are replaced.")
      }
      # append results to the object
      for (index in unique(id)) {
        Biobase::fData(object)[[paste0("pvalues_", index)]] <- NA
        Biobase::fData(object)[[paste0("pvalues_", index)]][match(rownames(result$pvalues), Biobase::featureNames(object), nomatch = 0)] <- result$pvalues[, index]
        Biobase::fData(object)[[paste0("scores_", index)]] <- NA
        Biobase::fData(object)[[paste0("scores_", index)]][match(rownames(result$scores_sp), Biobase::featureNames(object), nomatch = 0)] <- result$scores_sp[, index]
      }
    }
    return(object)
  }
)

#' Testing for features above the background
#'
#' Testing for features above the background using Poisson background model as reference
#'
#' @param object count matrix with features in rows and samples in columns
#' @param BGmod a list of sizefact, sizefact, and countmat
#' @param adj adjustment factor for the number of feature in each gene, default =1 i.e.
#'            each target only consists of one probe
#' @param probenum a vector of numbers of probes in each gene
#' @param removeoutlier whether to remove outlier
#' @param useprior whether to use the prior that the expression level of background follows a Beta distribution,
#'            leading to a more conservative test
#'
#' @importFrom graphics boxplot
#'
#' @return a list of following items
#' \itemize{
#'   \item pvalues - Background score test pvalues
#'   \item scores - Background score test statistics
#' }
#' @rdname BGScoreTest-methods
#' @aliases BGScoreTest,matrix-method

BGScoreTest_function = function(object, BGmod, adj = 1, probenum, removeoutlier = FALSE, useprior = FALSE) {
  if (removeoutlier == TRUE) {
    boxobj <- graphics::boxplot(BGmod$featfact, plot = FALSE)
    
    if (length(boxobj$out) > 0) {
      featfact <- BGmod$featfact[-which(BGmod$featfact %in% boxobj$out)]
    } else {
      featfact <- BGmod$featfact
    }
    message(sprintf("%s negative probes are removed prior to the score test.", length(boxobj$out)))
  } else {
    featfact <- BGmod$featfact
  }
  
  sizefact <- BGmod$sizefact
  
  if (useprior == FALSE) {
    if (missing(probenum)) {
      prodfact <- sizefact * mean(adj * featfact)
      
      scores <- apply(object, 1, function(x) sum(x - prodfact) / sqrt(sum(prodfact)))
    } else {
      if (is.null(names(probenum))) names(probenum) <- rownames(object)
      scores <- sapply(
        names(probenum),
        function(feat) {
          prodfact <- sizefact * mean(probenum[feat] * featfact)
          sum(object[feat, ] - prodfact) / sqrt(sum(prodfact))
        }
      )
    }
  } else {
    if (missing(probenum)) {
      featfact0 <- mean(adj * featfact)
      sigma <- var(adj * featfact) / (mean(adj * featfact))^2
      deno <- (sizefact * sigma * featfact0 + 1) * featfact0
      ### focus here!!! - issue with losing names attached to data... 
      #scores <- apply(object, 1, function(x) sum((x - sizefact * featfact0) / deno) / sqrt(sum(sizefact / deno)))
      new_object_numerator = ((object - sizefact * featfact0))
      new_object_denominator = 1/(sqrt(sum(sizefact / deno))*deno)
      new_object_denominator = as(new_object_denominator, "sparseMatrix") 
   
      # print(NROW(new_object_numerator))
      # print(NCOL(new_object_numerator))
      # print(NROW(new_object_denominator))
      # print(NCOL(new_object_denominator))
      quotient = new_object_numerator%*%new_object_denominator
      # print(NROW(quotient))
      # print(NCOL(quotient))
      scores_ned = quotient[,1]
      names(scores_ned) = rownames(object)
      scores = scores_ned
      #print(scores[[993]])
      # print(typeof(scores_ned))
      
    } else {
      if (is.null(names(probenum))) names(probenum) <- rownames(object)
      
      scores <- sapply(
        names(probenum),
        function(feat) {
          featfact0 <- mean(probenum[feat] * featfact)
          sigma <- var(probenum[feat] * featfact) / (mean(probenum[feat] * featfact))^2
          deno <- (sizefact * sigma * featfact0 + 1) * featfact0
          sum((object[feat, ] - sizefact * featfact0) / deno) / sqrt(sum(sizefact / deno))
        }
      )
    }
  }
  
  pvalues <- pnorm(scores, lower.tail = FALSE)
  
  return(list(
    pvalues = pvalues,
    scores = scores
  ))
}
setMethod(
  "BGScoreTest", "dgCMatrix",BGScoreTest_function
)
setMethod(
  "BGScoreTest", "matrix",BGScoreTest_function
)

#' Testing for features above the background, multiple slides case
#'
#' Testing for features above the background using Poisson background model as reference, multiple slides case
#'
#' @param object count matrix with features in rows and samples in columns
#' @param BGmod fitted background model, multiple slides case
#' @param probenum a vector of numbers of probes in each gene
#' @param adj adjustment factor for the number of probes in each feature, default =1 i.e.
#'            each target only consists of one probe
#' @param removeoutlier whether to remove outlier
#' @param useprior whether to use the prior that the expression level of background follows the Beta distribution,
#'                 leading to a more conservative test
#' @param ... additional argument list that might be used
#'
#' @importFrom graphics boxplot
#'
#' @return a list of following items
#' \itemize{
#'   \item pvalues - Background score test pvalues matrix, columns the same as slide names
#'   \item scores_sp - Background score test statistics matrix, columns the same as slide names
#' }
#'
#' @docType methods
#' @rdname BGScoreTest_sp-methods
#'
setGeneric("BGScoreTest_sp",
           signature = c("object"),
           function(object, ...) standardGeneric("BGScoreTest_sp")
)

#' @rdname BGScoreTest_sp-methods
#' @aliases BGScoreTest_sp,matrix-method
setMethod(
  "BGScoreTest_sp", "matrix",
  function(object, BGmod, adj = 1, probenum, removeoutlier = FALSE, useprior = FALSE) {
    id <- BGmod$id
    uniid <- unique(as.character(id))
    
    if (removeoutlier == TRUE) {
      
      #   boxobj <- apply(BGmod$featfact, 2, function(x) boxplot(x, plot = FALSE))
      
      featfact <- apply(BGmod$featfact, 2, function(x) {
        boxobj <- graphics::boxplot(x, plot = FALSE)
        message(sprintf("%s negative probes are removed prior to the score test.", length(boxobj$out)))
        x[which(x %in% boxobj$out)] <- NA
        x
      })
    } else {
      featfact <- BGmod$featfact
    }
    
    sizefact <- BGmod$sizefact
    
    
    if (useprior == FALSE) {
      if (missing(probenum)) {
        prodfact <- lapply(uniid, function(x) sizefact[x == id] * mean(adj * featfact[, x], na.rm = TRUE))
        names(prodfact) <- uniid
        #   scores <- apply(countmat, 2, function(x) sum(x - ab)/sqrt(sum(ab)))
        scores_sp <- sapply(uniid, function(x) apply(object[, x == id, drop = FALSE], 1, function(y) sum(y - prodfact[[x]]) / sqrt(sum(prodfact[[x]]))))
      } else {
        if (is.null(names(probenum))) names(probenum) <- rownames(object)
        scores_sp <- sapply(names(probenum), function(feat) {
          prodfact <- lapply(uniid, function(x) sizefact[x == id, drop = FALSE] * mean(probenum[feat] * featfact[, x], na.rm = TRUE))
          names(prodfact) <- uniid
          #   scores <- apply(countmat, 2, function(x) sum(x - ab)/sqrt(sum(ab)))
          sapply(uniid, function(x) sum(object[feat, x == id, drop = FALSE] - prodfact[[x]]) / sqrt(sum(prodfact[[x]])))
        })
        scores_sp <- t(scores_sp)
      }
    } else {
      if (missing(probenum)) {
        featfact0 <- Rfast::colmeans(adj * featfact)#, na.rm = TRUE)
        sigma <- apply(adj * featfact, 2, var, na.rm = TRUE) / featfact0^2
        deno <- lapply(uniid, function(x) (sizefact[x == id] * sigma[x] * featfact0[x] + 1) * featfact0[x])
        names(deno) <- uniid
        
        scores_sp <- sapply(uniid, function(x) apply(object[, x == id, drop = FALSE], 1, function(y) sum((y - sizefact[x == id] * featfact0[x]) / deno[[x]]) / sqrt(sum(sizefact[x == id] / deno[[x]]))))
        #  scores <- apply(scores2, 1, mean)
      } else {
        if (is.null(names(probenum))) names(probenum) <- rownames(object)
        scores_sp <- sapply(names(probenum), function(feat) {
          featfact0 <- Rfast::colmeans(probenum[feat] * featfact)#, na.rm = TRUE)
          sigma <- apply(probenum[feat] * featfact, 2, var, na.rm = TRUE) / featfact0^2
          deno <- lapply(uniid, function(x) (sizefact[x == id] * sigma[x] * featfact0[x] + 1) * featfact0[x])
          names(deno) <- uniid
          
          sapply(uniid, function(x) sum((object[feat, x == id, drop = FALSE] - sizefact[x == id] * featfact0[x]) / deno[[x]]) / sqrt(sum(sizefact[x == id] / deno[[x]])))
        })
        
        scores_sp <- t(scores_sp)
      }
    }
    
    pvalues <- pnorm(scores_sp, lower.tail = FALSE)
    
    return(list(
      pvalues = pvalues,
      scores_sp = scores_sp
    ))
  }
)