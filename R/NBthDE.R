#' Negative Binomial threshold model for differential expression analysis
#'
#' Negative Binomial threshold model for differential expression analysis
#'
#' @param object a valid GeoMx S4 object
#' @param form model formula
#' @param split indicator variable on whether it is for multiple slides (Yes, TRUE; No, FALSE)
#' @param ROIs_high ROIs with high expressions defined based on featfact and featfact
#' @param features_high subset of features which are well above the background
#' @param features_all full list of features
#' @param sizefact_start initial value for size factors
#' @param sizefact_BG size factor for background
#' @param threshold_mean average threshold level
#' @param preci2 precision for the background, default=10000
#' @param lower_threshold lower limit for the threshold, default=0.01
#' @param prior_type empirical bayes prior type, choose from c("equal","byvar", "contrast")
#' @param sizefactrec whether to recalculate sizefact, default=TRUE
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param sizescalebythreshold XXXX, default = FALSE
#' @param iterations how many iterations need to run to get final results, default=2,
#'                   the first iteration apply the model only on features_high and construct the prior then refit the model using this prior for all genes.
#' @param covrob whether to use robust covariance in calculating covariance. default=FALSE
#' @param preci1con The user input constant term in specifying precision matrix 1, default=1/25
#' @param cutoff term in calculating precision matrix 1, default=10
#' @param confac The user input factor for contrast in precision matrix 1, default=1
#' @param ... additional argument list that might be used
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase annotation
#'
#'
#' @return a list of
#' \itemize{
#'   \item X, design matrix
#'   \item para0, estimated parameters for the first iteration, including regression coefficients, r and threshold in rows and features in columns
#'   \item para, estimated parameters, including regression coefficients, r and threshold in rows and features in columns
#'   \item sizefact, estimated sizefact
#'   \item sizefact0, estimated sizefact in iter=1
#'   \item preci1, precision matrix for regression coefficients estimated in iter=1
#'   \item Im0, Information matrix of parameters in iter=1
#'   \item Im, Information matrix of parameters in iter=2
#'   \item conv0, vector of convergence for iter=1, 0 converged, 1 not converged
#'   \item conv, vector of convergence for iter=2, 0 converged, 1 not converged
#'   \item features_high, same as the input features_high
#'   \item features_all, same as the input features_all
#' }
#'
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
#' NBthDEmod1 <- fitNBthDE(
#'     form = ~group,
#'     split = FALSE,
#'     object = demoData,
#'     ROIs_high = ROIs_high,
#'     features_high = features_high,
#'     features_all = features_all,
#'     sizefact_start = demoData[, ROIs_high][["sizefact_fitNBth"]],
#'     sizefact_BG = demoData[, ROIs_high][["sizefact"]],
#'     preci2 = 10000,
#'     prior_type = "contrast",
#'     covrob = FALSE,
#'     preci1con = 1/25,
#'     sizescalebythreshold = TRUE
#' )
#'
#'
#'
#'
#

#' @export
#' @docType methods
#' @rdname fitNBthDE-methods
#'
setGeneric("fitNBthDE",
    signature = c("object"),
    function(object, ...) standardGeneric("fitNBthDE")
)

#' @rdname fitNBthDE-methods
#' @aliases fitNBthDE,NanoStringGeoMxSet-method
setMethod(
    "fitNBthDE", "NanoStringGeoMxSet",
    function(object, form, split, ROIs_high = NULL,
    features_high = NULL, features_all = NULL,
    sizefact_start = NULL, sizefact_BG = NULL,
    threshold_mean = NULL, preci2 = 10000, lower_threshold = 0.01,
    prior_type = c("contrast", "equal"), sizefactrec = TRUE,
    size_scale = c("sum", "first"), sizescalebythreshold = FALSE,
    iterations = 2, covrob = FALSE, preci1con=1/25, cutoff = 10, confac = 1) {
        fDat <- Biobase::fData(object)
        pDat <- Biobase::pData(object)

        posdat <- object[-which(fDat$CodeClass == "Negative"), ]
        countmat <- Biobase::exprs(posdat)

        fDatNeg <- fDat[which(fDat$CodeClass == "Negative"), ]

        # only calculate thmean if any of the three params are missing
        if (any(c(is.null(sizefact_BG), is.null(sizefact_start), is.null(ROIs_high), is.null(threshold_mean)))) {
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

        # setting default value for sizefact_start
        if (is.null(sizefact_start)) {
            if (!("sizefact_fitNBth" %in% colnames(pDat))) {
                stop("Please run `fitNBth` first.")
            } else {
                sizefact_start <- pDat[ROIs_high, ][["sizefact_fitNBth"]]
            }
            names(sizefact_start) <- rownames(pDat[ROIs_high, ])
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

        # replace threshold_mean for WTA
        if (is.null(threshold_mean)) {
            if (any(grepl("WTA", toupper(Biobase::annotation(object))))) {
                threshold_mean <- Biobase::notes(object)[["threshold"]]
            } else if (any(grepl("CTA", toupper(Biobase::annotation(object))))) {
                threshold_mean <- thmean
            } else {
                stop("No information is found to determine the data type (CTA or WTA).")
            }
        }

        result <- fitNBthDE(
            form = form,
            annot = annot[ROIs_high, ],
            object = countmat[, ROIs_high],
            probenum = probenum,
            features_high = features_high,
            features_all = features_all,
            sizefact_start = sizefact_start,
            sizefact_BG = sizefact_BG,
            threshold_mean = threshold_mean,
            preci2 = preci2,
            prior_type = prior_type,
            sizefactrec = sizefactrec,
            sizescalebythreshold = sizescalebythreshold,
            iterations = iterations,
            covrob = covrob,
            preci1con = preci1con
        )

        return(result)
    }
)


#' Negative Binomial threshold model for differential expression analysis
#'
#' Negative Binomial threshold model for differential expression analysis
#'
#' @param form model formula
#' @param annot annotations files with variables in the formula
#' @param object count matrix with features in rows and samples in columns
#' @param probenum a vector of numbers of probes in each gene, default = rep(1, NROW(object))
#' @param features_high subset of features which are well above the background
#' @param features_all full list of features
#' @param sizefact_start initial value for size factors
#' @param sizefact_BG size factor for background
#' @param threshold_mean average threshold level
#' @param preci2 precision for the background, default=10000
#' @param lower_threshold lower limit for the threshold, default=0.01
#' @param prior_type empirical bayes prior type, choose from c("contrast", "equal")
#' @param sizefactrec whether to recalculate sizefact, default=TRUE
#' @param size_scale method to scale the sizefact, sum(sizefact)=1 when size_scale="sum", sizefact[1]=1 when size_scale="first"
#' @param sizescalebythreshold XXXX, default = FALSE
#' @param iterations how many iterations need to run to get final results, default=2,
#'                   the first iteration apply the model only on features_high and construct the prior then refit the model using this prior for all genes.
#' @param covrob whether to use robust covariance in calculating covariance. default=FALSE
#' @param preci1con The user input constant term in specifying precision matrix 1, default=1/25
#' @param cutoff term in calculating precision matrix 1, default=10
#' @param confac The user input factor for contrast in precision matrix 1, default=1
#' @param ... additional argument list that might be used
#'
#' @importFrom Matrix bdiag
#' @importFrom robust covRob
#'
#' @return a list of
#' \itemize{
#'   \item X, design matrix
#'   \item para0, estimated parameters for the first iteration, including regression coefficients, r and threshold in rows and features in columns
#'   \item para, estimated parameters, including regression coefficients, r and threshold in rows and features in columns
#'   \item sizefact, estimated sizefact
#'   \item sizefact0, estimated sizefact in iter=1
#'   \item preci1, precision matrix for regression coefficients estimated in iter=1
#'   \item Im0, Information matrix of parameters in iter=1
#'   \item Im, Information matrix of parameters in iter=2
#'   \item conv0, vector of convergence for iter=1, 0 converged, 1 not converged
#'   \item conv, vector of convergence for iter=2, 0 converged, 1 not converged
#'   \item features_high, same as the input features_high
#'   \item features_all, same as the input features_all
#' }
#'
#' @rdname fitNBthDE-methods
#' @aliases fitNBthDE,matrix-method
#' 
fitNBthDE_funct =     function(form, annot, object, probenum,
                               features_high, features_all, sizefact_start, sizefact_BG,
                               threshold_mean, preci2=10000, lower_threshold = 0.01,
                               prior_type = c("contrast", "equal"), sizefactrec = TRUE,
                               size_scale = c("sum", "first"), sizescalebythreshold = FALSE,
                               iterations = 2, covrob = FALSE, preci1con=1/25, cutoff = 10, confac = 1) {
  if (iterations == 1) {
    if (!setequal(features_high, features_all)) {
      warning("features_high and features_all need to be identical when iterations=1,
            assign features_high <- features_all")
    }
    features_high <- features_all
  }
  
  X <- model.matrix(form, data = annot)
  
  
  sizefact0 <- sizefact <- sizefact_start
  n_sample <- nrow(X)
  n_para <- ncol(X)
  n_var <- max(attributes(X)$assign)
  n_levels <- as.list(table(attributes(X)$assign))
  var_ind <- as.numeric(names(n_levels))
  
  prior_type <- match.arg(prior_type)
  if (prior_type == "equal") {
    preci1 <- preci1con * t(X) %*% diag(1, n_sample) %*% X
  } else if (prior_type == "contrast") {
    contrvec <- t(rep(1 / n_sample, n_sample)) %*% X
    
    preci1 <- (preci1con) * t(contrvec) %*% contrvec
  }
  
  
  
  
  if (is.null(names(probenum))) names(probenum) <- rownames(object)
  
  if (sizescalebythreshold) {
    startpara <- c(rep(0, ncol(X)), 1, 1.0)
  } else {
    startpara <- c(rep(0, ncol(X)), 1, threshold_mean)
  }
  #might lead to memory blowing up... need to figure out t(sparse_mat[subset,])
  object_mat = as.matrix(object)
  for (iter in seq_len(iterations)) {
    if (iter == 1) {
      result <- NBthDE_paraOptall(
        t(object_mat[features_high, ]), X, sizefact_BG, sizefact,
        preci1, threshold_mean * probenum[features_high], preci2,
        startpara, sizescalebythreshold, (iter == iterations)
      )
      para <- result$par
      colnames(para) <- features_high
      conv <- result$conv
      names(conv) <- features_high
      Im <- result$hes
      names(Im) <- features_high
      Im0 <- Im
      para0 <- para
      conv0 <- conv
    } else {
      result <- NBthDE_paraOptall(
        t(object_mat[features_all, ]), X, sizefact_BG, sizefact,
        preci1, threshold_mean * probenum[features_all], preci2,
        startpara, sizescalebythreshold, (iter == iterations)
      )
      para <- result$par
      colnames(para) <- features_all
      conv <- result$conv
      names(conv) <- features_all
      Im <- result$hes
      names(Im) <- features_all
    }
    
    features_remain <- NA
    if ((iterations > 1) & (iter == 1)) {
      features_remain <- names(which(colMeans(abs(para[2:n_para, , drop = FALSE])) < cutoff))
      if (prior_type == "equal") {
        if (covrob) {
          cov_mat <- robust::covRob(t(para[seq_len(n_para), features_remain]), na.action = na.omit)$cov
        } else {
          cov_mat <- cov(t(para[seq_len(n_para), features_remain, drop = FALSE]), use = "pairwise.complete.obs")
        }
        preci1 <- solve(cov_mat)
      } else if (prior_type == "contrast") {
        if (covrob) {
          if (n_para == 2) {
            cov_mat0 <- (robust::covRob(t(para[seq_len(n_para), features_remain]), na.action = na.omit)$cov)[2:n_para, 2:n_para, drop = FALSE]
          } else {
            cov_mat0 <- robust::covRob(t(para[2:n_para, features_remain]), na.action = na.omit)$cov
          }
        } else {
          cov_mat0 <- cov(t(para[2:n_para, features_remain, drop = FALSE]), use = "pairwise.complete.obs")
        }
        
        contrvec <- t(rep(1 / n_sample, n_sample)) %*% X
        contrmat <- rbind(contrvec, cbind(rep(0, (n_para - 1)), diag(1, (n_para - 1))))
        preci_list <- list(`0` = diag(confac * preci1con, nrow = 1), preci_mat = solve(cov_mat0))
        preci10 <- Matrix::bdiag(preci_list)
        preci1 <- as.matrix(t(contrmat) %*% preci10 %*% contrmat)
      }
    }
    
    
    if (sizefactrec) {
      # if((EB)&(iter==1)){
      #   genes_NA <- features_high[unique(which(is.na(para), arr.ind = TRUE)[,2])]
      #   features_remain <- setdiff(features_high, genes_NA)
      # } else {
      #   genes_NA <- features_all[unique(which(is.na(para), arr.ind = TRUE)[,2])]
      #   features_remain <- setdiff(features_all, genes_NA)
      # }
      size_scale <- match.arg(size_scale)
      features_remain <- names(which(colMeans(abs(para[2:n_para, , drop = FALSE])) < cutoff))
      
      
      for (i in seq_len(length(sizefact))) {
        fun <- NBthDE_scalenll(X[i, ], object[features_remain, i], probenum[features_remain], para[seq_len(n_para), features_remain], t(para[n_para + 1, features_remain]), sizefact_BG[i], para[n_para + 2, features_remain], sizescalebythreshold, threshold_mean)
        sizefact[i] <- optim(c(sizefact[i]), fun, lower = c(0), method = "L-BFGS-B")$par
      }
      
      if (size_scale == "first") {
        scale_fac <- sizefact[1]
      } else if (size_scale == "sum") {
        scale_fac <- sum(sizefact)
      }
      
      sizefact <- sizefact / scale_fac
      
      message(sprintf(
        "Iteration = %s, squared error = %e",
        iter,
        sum((sizefact - sizefact0)^2)
      ))
      
      if (iter == 1) {
        sizefact0 <- sizefact
      }
    }
  }
  
  
  paraname <- c(colnames(X), c("r", "threshold"))
  rownames(para0) <- rownames(para) <- paraname
  
  return(list(
    X = X,
    para0 = para0,
    para = para,
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

setMethod(
    "fitNBthDE", "dgCMatrix", fitNBthDE_funct
)
setMethod(
  "fitNBthDE", "matrix", fitNBthDE_funct
)


#' Generate list of Wald test inference results on model coefficients
#'
#' Generate list of Wald test inference results including parameter estimation and p value
#'
#' @param object DE model, output by fitNBthDE or fitNBthmDE
#' @param fullpara whether to generate results on all parameters
#' @param ... additional argument list that might be used
#'
#' @return
#' \itemize{
#'   \item estimate, coefficients estimate
#'   \item wald_stat, Wald test statistics
#'   \item p_value, p value of Wald test
#'   \item se, standard error
#' }
#'
#' @examples
#' data(NBthmDEmod2)
#' coeff <- coefNBth(NBthmDEmod2)
#' @export
#' @docType methods
#' @rdname coefNBth-methods

setGeneric("coefNBth",
    signature = c("object"),
    function(object, ...) standardGeneric("coefNBth")
)

#' @rdname coefNBth-methods
#' @aliases coefNBth,list-method
setMethod(
    "coefNBth", "list",
    function(object, fullpara = FALSE) {
        X <- object$X
        if (fullpara) {
            n_para <- ncol(X) + 2
            paraname <- c(colnames(X), c("r", "threshold"))
        } else {
            n_para <- ncol(X)
            paraname <- colnames(X)
        }
        estimate <- object$para[seq_len(n_para), ]
        se <- lapply(object$Im, function(x) tryCatch(
            {
                x_tmp <- diag(solve(x))
                x_tmp[which(x_tmp <0 )] <- NaN
                res <- sqrt(x_tmp)[seq_len(n_para)]
                return(res)
            },
            error = function(err) NA))

        wald_stat <- mapply(function(x, y) (x / y), as.list(as.data.frame(estimate)), se, SIMPLIFY = FALSE)


        p_value <- sapply(wald_stat, function(x) 2 * pnorm(-abs(x)))
        wald_stat <- sapply(wald_stat, function(x) x)
        se <- sapply(se, function(x) x)

        rownames(estimate) <- rownames(wald_stat) <- rownames(p_value) <- rownames(se) <- paraname

        return(list(
            estimate = estimate,
            wald_stat = wald_stat,
            p_value = p_value,
            se = se
        ))
    }
)


#' Generate list of Wald test inference results on user specified contrasts
#'
#' Generate list of Wald test inference results including contrast estimation and p value
#'
#' @param object DE model, output by fitNBthDE or fitNBthmDE
#' @param test type of statistical test, choose from c("two-sided", ">", "<")
#' @param method contrasts methods, only matrix of contrast vector is allowed for now, default=diag(1,ncol(object$X)), i.e. testing the regression coefficients
#' @param baseline testing baseline, default=0.
#' @param ... additional argument list that might be used
#'
#' @return
#' \itemize{
#'   \item estimate, contrasts estimate
#'   \item wald_stat, Wald test statistics
#'   \item p_value, p value of Wald test
#'   \item se, standard error
#' }
#'
#' @examples
#' data(NBthmDEmod2)
#' coeff <- contrastNBth(NBthmDEmod2)
#' @export
#'
#' @docType methods
#' @rdname contrastNBth-methods
#'
setGeneric("contrastNBth",
    signature = c("object"),
    function(object, ...) standardGeneric("contrastNBth")
)

#' @rdname contrastNBth-methods
#' @aliases contrastNBth,list-method
setMethod(
    "contrastNBth", "list",
    function(object, test = c("two-sided", ">", "<"),
    method = diag(1, ncol(object$X)),
    baseline = rep(0, ncol(method))) {
        test <- match.arg(test)
        X <- object$X
        n_para <- ncol(X)
        if (is.matrix(method)) {
            contrast <- method
        } else {
            stop("matrix is the only contrast method is allowed for now")
        }
        estimate <- t(contrast) %*% object$para[seq_len(n_para), ]
        estimate <- estimate - matrix(rep(baseline, ncol(estimate)), ncol = ncol(estimate))
        se <- lapply(object$Im, function(x) tryCatch(sqrt(diag(t(contrast) %*% solve(x)[seq_len(n_para), seq_len(n_para)] %*% contrast)), error = function(err) NA))
        wald_stat <- mapply(function(x, y) (x / y), as.list(as.data.frame(estimate)), se, SIMPLIFY = FALSE)
        if (test == "two-sided") {
            p_value <- sapply(wald_stat, function(x) 2 * pnorm(-abs(x)))
        } else if (test == ">") {
            p_value <- sapply(wald_stat, function(x) pnorm(x, lower.tail = FALSE))
        } else {
            p_value <- sapply(wald_stat, function(x) pnorm(x, lower.tail = TRUE))
        }

        wald_stat <- sapply(wald_stat, function(x) x)
        se <- sapply(se, function(x) x)
        rownames(estimate) <- rownames(wald_stat) <- rownames(p_value) <- rownames(se) <- colnames(contrast)

        return(list(
            estimate = estimate,
            wald_stat = wald_stat,
            p_value = p_value,
            se = se
        ))
    }
)



#' Generate DE table using the inference list generated by coefNBth or contrastNBth
#'
#' Generate DE table using the inference list generated by coefNBth or contrastNBth
#'
#' @param object inference list from coefNBth or contrastNBth
#' @param variable needed to construct
#' @param NAto1 whether to replace NA in pvalue by 1
#' @param padj whether to adjust p value
#' @param padj_method p value adjustment method, default="BH"
#' @param ... additional argument list that might be used
#'
#' @return DEtab, DE table
#'
#' @examples
#' data(NBthmDEmod2)
#' coeff <- coefNBth(NBthmDEmod2)
#' DEtab <- DENBth(coeff, variable = "regiontubule")
#' @export
#'
#' @docType methods
#' @rdname DENBth-methods
#'
setGeneric("DENBth",
    signature = c("object"),
    function(object, ...) standardGeneric("DENBth")
)

#' @rdname DENBth-methods
#' @aliases DENBth,list-method
setMethod(
    "DENBth", "list",
    function(object, variable, NAto1 = TRUE, padj = TRUE, padj_method = "BH") {
        DEtab <- as.data.frame(cbind(log2FC = object$estimate[variable, ], pvalue = object$p_value[variable, ]))
        if (NAto1) {
            DEtab$pvalue[is.na(DEtab$pvalue)] <- 1
        }

        if (padj) {
            DEtab$adjp <- p.adjust(DEtab$pvalue, method = padj_method)
        }


        return(DEtab)
    }
)


