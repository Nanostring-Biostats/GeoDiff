PoisthNorm_scalenll <- function(X, Y, probenum, betamat, sizefact0, threshold, sizescale, threshold_mean) {
    tmp0 <- 2^(X %*% betamat)

    if (sizescale) {
        loglik_ind <- function(sizefact) {
            tmp1 <- probenum * sizefact0 * threshold_mean * threshold + probenum * sizefact * threshold_mean * tmp0


            tmp3 <- dpois(x = Y, lambda = tmp1, log = TRUE)

            # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
            -(sum(tmp3))
        }
    } else {
        loglik_ind <- function(sizefact) {
            tmp1 <- probenum * sizefact0 * threshold + sizefact * tmp0


            tmp3 <- dpois(x = Y, lambda = tmp1, log = TRUE)

            # - ((log(phi) - m0)^2)/(2 * (sigma^2)) - log(sigma)
            -(sum(tmp3))
        }
    }

    return(loglik_ind)
}
