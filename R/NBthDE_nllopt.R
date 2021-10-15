






NBthDE_scalenll <- function(X, Y, probenum, regcoefmat, rvec, sizefact0, threshold, sizescale, threshold_mean) {
    tmp0 <- 2^(X %*% regcoefmat)
    if (sizescale) {
        loglik_ind <- function(sizefact) {
            tmp1 <- probenum * threshold_mean * sizefact0 * threshold + probenum * threshold_mean * sizefact * tmp0


            tmp3 <- dnbinom(x = Y, size = rvec, mu = tmp1, log = TRUE)

            # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
            -(sum(tmp3))
        }
    } else {
        loglik_ind <- function(sizefact) {
            tmp1 <- sizefact0 * probenum * threshold + sizefact * tmp0


            tmp3 <- dnbinom(x = Y, size = rvec, mu = tmp1, log = TRUE)

            # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
            -(sum(tmp3))
        }
    }
    return(loglik_ind)
}
