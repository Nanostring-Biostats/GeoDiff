
NBth_paranll <- function(Y, probenum, sizefact, sizefact0, threshold) {
    loglik_ind <- function(x) {
        featfact <- x[1]
        r <- x[2]
        # tmp1 = 1/(1 + (sizefact0*min(featfact, threshold)+sizefact*max(featfact-threshold,0))/r)
        tmp1 <- sizefact0 * min(featfact, probenum * threshold) + sizefact * max(featfact - probenum * threshold, 0)

        # tmp3 = dnbinom(x = Y, size = r, prob = tmp1, log = TRUE)

        tmp3 <- dnbinom(x = Y, size = r, mu = tmp1, log = TRUE)

        # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
        -(sum(tmp3))
    }
}

NBth_scalenll <- function(Y, probenum, featfactvec, rvec, sizefact0, threshold) {
    loglik_ind <- function(sizefact) {


        # tmp1 = 1/(1 + (sizefact0*pmin(featfactvec, threshold)+sizefact*pmax(featfactvec-threshold,0))/rvec)
        tmp1 <- sizefact0 * pmin(featfactvec, probenum * threshold) + sizefact * pmax(featfactvec - probenum * threshold, 0)

        # tmp3 <- dnbinom(x = Y, size=rvec, prob = tmp1, log = TRUE)

        tmp3 <- dnbinom(x = Y, size = rvec, mu = tmp1, log = TRUE)

        # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
        -(sum(tmp3))
    }
}



NBth_thnll <- function(Y, probenum, sizefact, sizefact0, featfactvec, rvec) {
    loglik_ind <- function(threshold) {


        #  tmp1 = 1/(1 + sweep((sizefact0%*%(pmin(featfactvec, threshold))+sizefact%*%(pmax(featfactvec-threshold,0))), 2, rvec, FUN="/"))

        tmp1 <- sizefact0 %*% (pmin(featfactvec, probenum * threshold)) + sizefact %*% (pmax(featfactvec - probenum * threshold, 0))



        # tmp3 <- dnbinom(x = as.matrix(Y), size = matrix(rep(rvec, each=length(sizefact)), ncol=length(rvec)),
        #                 prob = tmp1, log = TRUE)

        tmp3 <- dnbinom(
            x = t(as.matrix(Y)), size = matrix(rep(rvec, each = length(sizefact)), ncol = length(rvec)),
            mu = tmp1, log = TRUE
        )



        # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
        -(sum(tmp3))
    }
}


NBth_paraopt <- function(countmat, probenum, sizefact, sizefact0, threshold, start = c(0.5, 0.5), lower = c(0.01, 0.01)) {
    if (is.null(names(probenum))) probenum <- rownames(countmat)
    para <- sapply(rownames(countmat), function(x) {
        fun1 <- NBth_paranll(countmat[x, ], probenum[x], sizefact, sizefact0, threshold)
        result <- tryCatch(optim(start, fun1, lower = lower, method = "L-BFGS-B")$par, error = function(err) c(NA, NA))
    })
    return(para)
}
