#library(R)
#' Estimate Poisson background model for either single slide or multiple slides
#'
#' Estimate Poisson background model for either single slide or multiple slides
#'
#' @param object a valid GeoMx S4 object
#' @param groupvar the group variable name for slide
#' @param iterations maximum iterations to be run, default=10
#' @param tol tolerance to determine convergence, default = 1e-3
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param ... additional argument list that might be used
#'
#' @importFrom Biobase exprs
#' @importFrom Biobase fData
#' @importFrom Biobase sampleNames
#' @importFrom Biobase featureNames
#'
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @return a valid GeoMx S4 object if split is FALSE
#' \itemize{
#'   \item sizefact - estimated size factor in phenoData
#'   \item featfact - estimated feature factor in featureData
#' }
#'
#' a valid GeoMx S4 object if split is TRUE,
#' \itemize{
#'   \item sizefact - estimated size factor in phenoData
#'   \item featfact_XX - estimated feature factor vector, column name (denoted as XX) the same as the slide id, in featureData for each unique slide
#'   \item fitPoisBG_sp_var - the column name for slide, in experimentData
#' }
#'
#' @examples
#'
#' data(demoData)
#' demoData <- fitPoisBG(demoData, size_scale = "sum")
#' data(demoData)
#' demoData <- fitPoisBG(demoData, groupvar = "slide name", size_scale = "sum")
#' @export
#' @docType methods
#' @rdname fitPoisBG-methods
setGeneric("fitPoisBG",
           signature = c("object"),
           function(object, ...) standardGeneric("fitPoisBG")
)

#' @rdname fitPoisBG-methods
#' @aliases fitPoisBG,NanoStringGeoMxSet-method
setMethod(
  "fitPoisBG", "NanoStringGeoMxSet",
  function(object, groupvar = NULL, iterations = 10, tol = 1e-3, size_scale = c("sum", "first"), ...) {
    # check on tol
    tol <- as.double(tol)
    stopifnot(length(tol) == 1)
    stopifnot(tol >= 0)
    
    # extract the negative probes matrix
    negdat <- object[which(Biobase::fData(object)$CodeClass == "Negative"), ]
    countmat <- Biobase::exprs(negdat)
    
    if (is.null(groupvar)) {
      # fit the model as a single slide
      split <- FALSE
    } else if (!(groupvar %in% varLabels(object))) {
      stop(sprintf("%s is not found in the S4 object.", groupvar))
    } else {
      # extract the groupvar
      id <- object[[groupvar]]
      
      if (length(unique(id)) == 1) {
        # this is equivalent as fitting the model as a single slide
        split <- FALSE
        warning(sprintf("`%s` has only one value, %s", groupvar, unique(id)))
      } else {
        # fit the model for multiple slides
        split <- TRUE
      }
    }
    
    if (isFALSE(split)) {
      # calling the fitPoisBG function when split is FALSE
      result <- fitPoisBG(
        object = countmat,
        iterations = iterations,
        tol = tol,
        size_scale = size_scale
      )
      
      object[["sizefact"]] <- result$sizefact[Biobase::sampleNames(object)]
      Biobase::fData(object)[["featfact"]] <- NA
      Biobase::fData(object)[["featfact"]][match(names(result$featfact), Biobase::featureNames(object), nomatch = 0)] <- result$featfact
    } else {
      # calling the fitPoisBG_sp function when split is TRUE
      result <- fitPoisBG_sp(
        object = countmat,
        id = id,
        iterations = iterations,
        tol = tol,
        size_scale = size_scale
      )
      
      # append results to the object
      object[["sizefact_sp"]] <- result$sizefact[Biobase::sampleNames(object)]
      for (index in unique(result$id)) {
        Biobase::fData(object)[[paste0("featfact_", index)]] <- NA
        Biobase::fData(object)[[paste0("featfact_", index)]][match(rownames(result$featfact), Biobase::featureNames(object), nomatch = 0)] <- result$featfact[, index]
      }
      
      Biobase::notes(object)$fitPoisBG_sp_var <- groupvar
    }
    
    return(object)
  }
)

#' Estimate Poisson background model
#'
#' Estimate Poisson background model:
#'
#' @param object count matrix with features in rows and samples in columns
#' @param iterations maximum iterations to be run, default=10
#' @param tol tolerance to determine convergence, default = 1e-3
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#'
#' @return a list of following items
#' \itemize{
#'   \item sizefact - estimated size factor
#'   \item featfact - estimated feature factor
#'   \item countmat - the input count matrix
#' }
#'
#' @rdname fitPoisBG-methods
#' @aliases fitPoisBG,matrix-method
#' 
#' 
fitPoisBG_function = function(object, iterations = 10, tol = 1e-3, size_scale = c("sum", "first")) {
  size_scale <- match.arg(size_scale)
  n_feature <- NROW(object)
  n_sample <- NCOL(object)
  ind_na <- which(is.na(object), arr.ind = TRUE)
  featfact0 <- numeric(n_feature)
  sizefact <- colMeans(object)
  
  if (size_scale == "first") {
    scale_fac <- sizefact[1]
  } else if (size_scale == "sum") {
    scale_fac <- sum(sizefact)
  }
  
  sizefact <- sizefact / scale_fac
  sizefact0 <- sizefact
  sizefact_mat <- matrix(rep(sizefact, n_feature), n_feature, n_sample, byrow = TRUE)
  sizefact_mat[ind_na] <- NA
  object_rowsum = rowSums(object)
  object_colsum = colSums(object)
  for (iter in seq_len(iterations)) {
    featfact <- object_rowsum / Rfast::rowsums(sizefact_mat)
    featfact_mat <- matrix(rep(featfact, n_sample), n_feature, n_sample)
    featfact_mat[ind_na] <- NA
    
    sizefact <- object_colsum / Rfast::colsums(featfact_mat)
    names(sizefact) = colnames(object)
    if (size_scale == "first") {
      scale_fac <- sizefact[1]
    } else if (size_scale == "sum") {
      scale_fac <- sum(sizefact)
    }
    
    sizefact <- sizefact / scale_fac
    
    sizefact_mat <- matrix(rep(sizefact, n_feature), n_feature, n_sample, byrow = TRUE)
    sizefact_mat[ind_na] <- NA
    
    message(sprintf(
      "Iteration = %s, squared error = %e",
      iter,
      sum((sizefact - sizefact0)^2) + sum((featfact - featfact0)^2)
    ))
    if (sum((sizefact - sizefact0)^2) + sum((featfact - featfact0)^2) < tol) {
      break
    }
    
    sizefact0 <- sizefact
    featfact0 <- featfact
  }
  message("Model converged.")
  
  return(list(
    sizefact = sizefact,
    featfact = featfact,
    countmat = object
  ))
}
setMethod(
  "fitPoisBG", "dgCMatrix",fitPoisBG_function
)
setMethod(
  "fitPoisBG", "matrix",fitPoisBG_function
)



#' Estimate Poisson background model for multiple slides
#'
#' Estimate Poisson background model for multiple slides:
#'
#' @param object count matrix with features in rows and samples in columns
#' @param id character vector same size as sample size representing slide names of each sample
#' @param iterations maximum iterations to be run, default=10
#' @param tol tolerance to determine convergence, default = 1e-3
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param ... additional argument list that might be used
#'
#' @return a list of following items
#' \itemize{
#'   \item sizefact - estimated size factor
#'   \item featfact - estimated feature factor matrix, column names the same as the slide id
#'   \item countmat - the input count matrix
#'   \item id - the input id
#' }
#'
#'
#' @docType methods
#' @rdname fitPoisBG_sp-methods
#'
setGeneric("fitPoisBG_sp",
           signature = c("object"),
           function(object, ...) standardGeneric("fitPoisBG_sp")
)

#' @rdname fitPoisBG_sp-methods
#' @aliases fitPoisBG_sp,matrix-method
setMethod(
  "fitPoisBG_sp", "matrix",
  function(object, id, iterations = 10, tol = 1e-3, size_scale = c("sum", "first")) {
    size_scale <- match.arg(size_scale)
    
    uniid <- unique(as.character(id))
    n_feature <- NROW(object)
    n_sample <- NCOL(object)
    ind_na <- which(is.na(object), arr.ind = TRUE)
    
    sizefact <- Rfast::colmeans(object, 2)
    
    if (size_scale == "first") {
      scale_fac <- sizefact[1]
    } else if (size_scale == "sum") {
      scale_fac <- sum(sizefact)
    }
    
    
    sizefact <- sizefact / scale_fac
    
    sizefact0 <- sizefact
    sizefact_mat <- matrix(rep(sizefact, n_feature), n_feature, n_sample, byrow = TRUE)
    sizefact_mat[ind_na] <- NA
    featfact0 <- matrix(0, n_feature, length(uniid))
    for (iter in seq_len(iterations)) {
      #featfact <- sapply(uniid, function(x) {
      # apply(object[, x == id, drop = FALSE], 1, sum, na.rm = TRUE) /
      #   apply(sizefact_mat[, x == id, drop = FALSE], 1, sum, na.rm = TRUE)
      featfact <- sapply(uniid, function(x) {
        Rfast::rowsums(object[, x == id, drop = FALSE]) /
          Rfast::rowsums(sizefact_mat[, x == id, drop = FALSE])
      })
      
      featfact_mat <- featfact[, id]
      featfact_mat[ind_na] <- NA
      
      #sizefact <- apply(object, 2, sum, na.rm = TRUE) / apply(featfact_mat, 2, sum, na.rm = TRUE)
      sizefact <- Rfast::colsums(object) / Rfast::colsums(featfact_mat)
      if (size_scale == "first") {
        scale_fac <- sizefact[1]
      } else if (size_scale == "sum") {
        scale_fac <- sum(sizefact)
      }
      
      sizefact <- sizefact / scale_fac
      
      
      sizefact_mat <- matrix(rep(sizefact, n_feature), n_feature, n_sample, byrow = TRUE)
      sizefact_mat[ind_na] <- NA
      
      message(sprintf(
        "Iteration = %s, squared error = %e",
        iter,
        sum((sizefact - sizefact0)^2) + sum((featfact - featfact0)^2)
      ))
      
      if (sum((sizefact - sizefact0)^2) + sum((featfact - featfact0)^2) < tol) {
        break
      }
      
      sizefact0 <- sizefact
      featfact0 <- featfact
    }
    message("Model converged.")
    
    return(list(
      sizefact = sizefact,
      featfact = featfact,
      countmat = object,
      id = id
    ))
  }
)


#' Perform diagnosis on Poisson background model
#'
#' Perform diagnosis on Poisson background model
#'
#'
#' @param object a valid GeoMx S4 object
#' @param split indicator variable on whether it is for multiple slides (Yes, TRUE; No, FALSE)
#' @param padj whether to adjust p value for outlier detection, default =TRUE
#' @param padj_method p value adjustment method, default ="BH"
#' @param cutoff p value(or adjusted p value) cutoff to determine outliers
#' @param generate_ppplot whether to generate ppplot, default =TRUE
#' @param ... additional argument list that might be used
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase notes
#' @importFrom Biobase assayDataElement
#'
#'
#' @return a valid S4 object
#' \itemize{
#'   \item lowtail - A matrix of lower tail probabilty, in assay slot
#'   \item uptail - A matrix of upper tail probability, in assay slot
#'   \item disper (or disper_sp if non single-valued groupvar is provided) - dispersion parameter in experimenetData
#'   \item low_outlier - A matrix to indicate lower outliers (0:False, 1:True) in assay slot
#'   \item upper_outlier - A matrix to indicate upper outliers (0:False, 1:True) in assay slot
#' }
#'
#' @examples
#' data(demoData)
#' demoData <- fitPoisBG(demoData, size_scale = "sum")
#' demoData <- diagPoisBG(demoData)
#' Biobase::notes(demoData)$disper
#' demoData <- fitPoisBG(demoData, groupvar = "slide name")
#' demoData <- diagPoisBG(demoData, split = TRUE)
#' Biobase::notes(demoData)$disper_sp
#' @export
#' @docType methods
#' @rdname diagPoisBG-methods


setGeneric("diagPoisBG",
           signature = c("object"),
           function(object, ...) standardGeneric("diagPoisBG")
)

#' @rdname diagPoisBG-methods
#' @aliases diagPoisBG,NanoStringGeoMxSet-method
setMethod(
  "diagPoisBG", "NanoStringGeoMxSet",
  function(object, split = FALSE, padj = FALSE, padj_method = "BH", cutoff = 1e-6, generate_ppplot = TRUE) {
    negdat <- object[which(Biobase::fData(object)$CodeClass == "Negative"), ]
    countmat <- Biobase::exprs(negdat)
    pDat <- Biobase::pData(negdat)
    fDat <- Biobase::fData(negdat)
    
    if (isFALSE(split)) {
      if (!any(c("sizefact" %in% colnames(pDat), "featfact" %in% colnames(fDat)))) {
        stop("Please run `fitPoisBG` first. If you run `fitPoisBG` before, please specify `split = TRUE`.")
      }
      sizefact <- setNames(pDat[["sizefact"]], rownames(pDat))
      featfact <- setNames(fDat[["featfact"]], rownames(fDat))
      
      BGmod <- list(
        sizefact = sizefact,
        featfact = featfact,
        countmat = countmat
      )
    } else {
      if (!any(c("sizefact_sp" %in% colnames(pDat), "featfact_" %in% colnames(fDat)))) {
        stop("Please run `fitPoisBG` first with `groupvar`.")
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
        countmat = countmat,
        id = id
      )
    }
    
    
    result <- diagPoisBG(BGmod,
                         padj = padj,
                         padj_method = padj_method,
                         cutoff = cutoff,
                         generate_ppplot = generate_ppplot
    )
    
    # add lowtail matrix
    lowtail_mat <- matrix(NA,
                          nrow = nrow(object), ncol = ncol(object),
                          dimnames = dimnames(object)
    )
    lowtail_mat[rownames(result$uptail_prob), colnames(result$uptail_prob)] <- result$uptail_prob
    Biobase::assayDataElement(object, "lowtail_prob") <- lowtail_mat
    
    # add uptail matrix
    uptail_mat <- matrix(NA,
                         nrow = nrow(object), ncol = ncol(object),
                         dimnames = dimnames(object)
    )
    uptail_mat[rownames(result$uptail_prob), colnames(result$uptail_prob)] <- result$uptail_prob
    Biobase::assayDataElement(object, "uptail_prob") <- uptail_mat
    
    # add disper parameter
    if (isTRUE(split)) {
      Biobase::notes(object)$disper_sp <- result$disper
    } else {
      Biobase::notes(object)$disper <- result$disper
    }
    
    # add upper tail outlier
    low_outlier_mat <- matrix(0,
                              nrow = nrow(object), ncol = ncol(object),
                              dimnames = dimnames(object)
    )
    for (i in seq_len(nrow(result$outlier$low_outlier))) {
      low_outlier_mat[rownames(result$outlier$low_outlier)[i], result$outlier$low_outlier[i, "col"]] <- 1
    }
    Biobase::assayDataElement(object, "low_outlier") <- low_outlier_mat
    
    # add lower tail outlier
    up_outlier_mat <- matrix(0,
                             nrow = nrow(object), ncol = ncol(object),
                             dimnames = dimnames(object)
    )
    for (i in seq_len(nrow(result$outlier$up_outlier))) {
      up_outlier_mat[rownames(result$outlier$up_outlier)[i], result$outlier$up_outlier[i, "col"]] <- 1
    }
    Biobase::assayDataElement(object, "up_outlier") <- up_outlier_mat
    
    return(object)
  }
)

#' Perform diagnosis on Poisson background model
#'
#' Perform diagnosis on Poisson background model
#'
#' @param object a list of sizefact, featfact, countmat, or id (if it is for mutliple slides)
#' @param padj whether to adjust p value for outlier detection, default =TRUE
#' @param padj_method p value adjustment method, default ="BH"
#' @param cutoff p value(or adjusted p value) cutoff to determine outliers
#' @param generate_ppplot whether to generate ppplot, default =TRUE
#'
#'
#' @importFrom graphics abline
#'
#' @return a list of following items
#' \itemize{
#'   \item lowtail - A matrix of lower tail probabilty
#'   \item uptail - A matrix of upper tail probability
#'   \item disper - dispersion parameter
#'   \item outlier - A list of coodinates of lower and upper outliers
#' }
#' @rdname diagPoisBG-methods
#' @aliases diagPoisBG,list-method

setMethod(
  "diagPoisBG", "list",
  function(object, padj = FALSE, padj_method = "BH", cutoff = 1e-6, generate_ppplot = TRUE) {
    countmat <- object$countmat
    if (NCOL(object$featfact) == 1) {
      countmat_expected <- (object$featfact %*% t(object$sizefact))
    } else {
      countmat_expected <- sweep(object$featfact[, object$id], 2, object$sizefact, FUN = "*")
    }
    
    countmat <- as.matrix(countmat)
    #conmat = rbind(as.vector(countmat),as.vector(countmat_expected))
    #unique_pairs = unique(conmat, MARGIN=2)
    #print(length(unique_pairs))
    lowtail_prob1 <- ppois(q = countmat, lambda = countmat_expected)
    lowtail_prob2 <- ppois(q = countmat - 1, lambda = countmat_expected)
    
    #print(length(unique()) )   
    lowtail_prob <- (lowtail_prob1 + lowtail_prob2) / 2
    
    uptail_prob <- 1 - lowtail_prob
    
    # simualte data (do it only once)
    #WHY?
    countmat_simu <- t(rpois(rep(1, ncol(countmat)), lambda = countmat_expected))
    lowtail_prob1_simu <- ppois(q = countmat_simu, lambda = countmat_expected)
    lowtail_prob2_simu <- ppois(q = countmat_simu - 1, lambda = countmat_expected)
    lowtail_prob_simu <- (lowtail_prob1_simu + lowtail_prob2_simu) / 2
    
    if (generate_ppplot) {
      y <- sort(lowtail_prob, na.last = TRUE)
      y_simu <- sort(lowtail_prob_simu, na.last = TRUE)
      plot(y_simu, y, ylim = c(0, 1), xlab = "Empircal CDF from simulated data", ylab = "Empircal CDF", main = "Poisson model")
      graphics::abline(a = 0, b = 1)
    }
    
    tmp_mat = (countmat - countmat_expected)^2 / countmat_expected
    tmp_mat[is.na(tmp_mat)] <- 0
    disper = sum(Rfast::colsums(tmp_mat))/length(unname(countmat))
    #disper <- mean((countmat - countmat_expected)^2 / countmat_expected, na.rm = TRUE)
    
    
    if (padj) {
      lowtail_prob <- matrix(p.adjust(lowtail_prob, method = padj_method), nrow = nrow(lowtail_prob), ncol = ncol(lowtail_prob))
      uptail_prob <- matrix(p.adjust(uptail_prob, method = padj_method), nrow = nrow(uptail_prob), ncol = ncol(uptail_prob))
    }
    
    
    
    low_outlier <- which(lowtail_prob < cutoff, arr.ind = TRUE)
    up_outlier <- which(uptail_prob < cutoff, arr.ind = TRUE)
    
    return(list(
      lowtail_prob = lowtail_prob,
      uptail_prob = uptail_prob,
      disper = disper,
      outlier = list(
        low_outlier = low_outlier,
        up_outlier = up_outlier
      )
    ))
  }
)

#' Compute Quantile Range
#'
#' Compute Quantile Range, a metric representing signal strength for QC purpose
#'
#' @param object a valid GeoMx S4 object
#' @param split indicator variable on whether it is for multiple slides
#' @param probs numeric vector of probabilities with values in [0,1] passed to quantile
#' @param removeoutlier indicator on whether to remove outliers, default: FALSE
#' @param ... additional argument list that might be used
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase notes
#' @importFrom Biobase annotation
#'
#'
#' @return a valid S4 object with probabilities in phenoData
#'
#' @examples
#' data(demoData)
#' demoData <- fitPoisBG(demoData, size_scale = "sum")
#' demoData <- diagPoisBG(demoData)
#' demoData <- aggreprobe(demoData, use = "cor")
#' Biobase::notes(demoData)$disper
#' demoData <- QuanRange(demoData, split = FALSE, probs = c(0.75, 0.8, 0.9, 0.95))
#'
#' data(demoData)
#' demoData <- fitPoisBG(demoData, groupvar = "slide name")
#' demoData <- diagPoisBG(demoData, split = TRUE)
#' demoData <- aggreprobe(demoData, use = "cor")
#' Biobase::notes(demoData)$disper_sp
#' demoData <- QuanRange(demoData, split = TRUE, probs = c(0.75, 0.8, 0.9, 0.95))
#' @export
#' @docType methods
#' @rdname QuanRange-methods

setGeneric("QuanRange",
           signature = c("object"),
           function(object, ...) standardGeneric("QuanRange")
)

#' @rdname QuanRange-methods
#' @aliases QuanRange,NanoStringGeoMxSet-method
setMethod(
  "QuanRange", "NanoStringGeoMxSet",
  function(object, split = FALSE, probs, removeoutlier = FALSE, ...) {
    fDat <- Biobase::fData(object)
    pDat <- Biobase::pData(object)
    # get the positive probes
    
    posdat <- object[-which(fDat$CodeClass == "Negative"), ]
    countmat <- Biobase::exprs(posdat)
    
    # get the negative probes
    negdat <- object[which(fDat$CodeClass == "Negative"), ]
    pDatNeg <- Biobase::pData(negdat)
    fDatNeg <- Biobase::fData(negdat)
    
    if (isFALSE(split)) {
      if (!any(c("sizefact" %in% colnames(pDat), "featfact" %in% colnames(fDat)))) {
        stop("Please run `fitPoisBG` first. If you run `fitPoisBG` before, please specify `split = TRUE`.")
      }
      sizefact <- setNames(pDatNeg[["sizefact"]], rownames(pDatNeg))
      featfact <- setNames(fDatNeg[["featfact"]], rownames(fDatNeg))
      
      BGmod <- list(
        sizefact = sizefact,
        featfact = featfact,
        countmat = countmat
      )
    } else {
      if (!any(c("sizefact_sp" %in% colnames(pDatNeg), "featfact_" %in% colnames(fDatNeg)))) {
        stop("Please run `fitPoisBG` first with `groupvar`.")
      }
      sizefact <- setNames(pDatNeg[["sizefact_sp"]], rownames(pDatNeg))
      featfact <- fDatNeg[, grep("featfact_", colnames(fDatNeg))]
      colnames(featfact) <- gsub("featfact_", "", colnames(featfact))
      idvar <- Biobase::notes(object)[["fitPoisBG_sp_var"]]
      id <- pDat[[idvar]]
      message(sprintf("The results are based on stored `groupvar`, %s", idvar))
      
      BGmod <- list(
        sizefact = sizefact,
        featfact = as.matrix(featfact),
        countmat = countmat,
        id = id
      )
    }
    
    # calculate probenum for the dataset
    if ("probenum" %in% fvarLabels(posdat)) {
      probenum <- fData(posdat)[["probenum"]]
    } else {
      stop("No `probenum` is found. Run `aggreprobe` first.")
    }
    names(probenum) <- rownames(fData(posdat))
    
    result <- QuanRange(
      object = countmat,
      probenum = probenum,
      BGmod = BGmod,
      probs = probs,
      removeoutlier = removeoutlier
    )
    
    for (prob in as.character(probs)) {
      object[[prob]] <- result[, prob]
    }
    return(object)
  }
)


#' Compute Quantile Range
#'
#' Compute Quantile Range, a metric representing signal strength for QC purpose
#'
#' @param object count matrix with features in rows and samples in columns
#' @param BGmod a list of sizefact, sizefact, countmat, and id (if it is for multiple slides)
#' @param probenum a vector of numbers of probes in each gene
#' @param probs numeric vector of probabilities with values in [0,1] passed to quantile
#' @param removeoutlier indicator on whether to remove outliers, default: FALSE
#'
#'
#' @importFrom graphics boxplot
#'
#' @return a matrix of quantile range in rows and probs in columns
#'
#' @rdname QuanRange-methods
#' @aliases QuanRange,matrix-method

setMethod(
  "QuanRange", "matrix",
  function(object, probenum, BGmod, probs, removeoutlier = FALSE) {
    if (NCOL(BGmod$featfact) == 1) {
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
      
      back <- mean(featfact) * BGmod$sizefact
    } else {
      if (removeoutlier == TRUE) {
        featfact <- apply(BGmod$featfact, 2, function(x) {
          boxobj <- graphics::boxplot(x, plot = FALSE)
          message(sprintf("%s negative probes are removed prior to the score test.", length(boxobj$out)))
          x[which(x %in% boxobj$out)] <- NA
          x
        })
      } else {
        featfact <- BGmod$featfact
      }
      
      back <- Rfast::colmeans(featfact[, BGmod$id], na.rm = TRUE) * BGmod$sizefact
    }
    
    
    
    
    object <- sweep(object, 1, probenum, FUN = "/")
    
    quan <- sapply(probs, function(y) apply(object, 2, function(x) quantile(x, probs = y)))
    
    colnames(quan) <- probs
    quanrange <- sweep(quan, 1, back, FUN = "-")
    
    return(quanrange)
  }
)

