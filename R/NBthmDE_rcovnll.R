NBthm_rcovnll <- function(Umat, Lambdati, mapping) {
    N <- nrow(Umat)
    loglik_fun <- function(theta) {
        Lambdati@x[] <- mapping(theta)
        Sigma <- as.matrix(Lambdati) %*% t(as.matrix(Lambdati))
        value <- -(N / 2) * determinant(Sigma, logarithm = TRUE)$modulus - 1 / 2 * (sum(apply(t(solve(Sigma, t(Umat))) * Umat, 1, sum)))
        return(-value)
    }
}



NBthm_rcovnlliWp <- function(Umat, Lambdati, mapping, nu, Lambda0) {
    N <- nrow(Umat)
    d <- ncol(Umat)
    loglik_fun <- function(theta) {
        Lambdati@x[] <- mapping(theta)
        Sigma <- as.matrix(Lambdati) %*% t(as.matrix(Lambdati))
        value <- -((N + nu + d + 1) / 2) * determinant(Sigma, logarithm = TRUE)$modulus -
            (1 / 2) * (sum(apply(t(solve(Sigma, t(Umat))) * Umat, 1, sum))) -
            (1 / 2) * sum(diag(solve(Sigma, Lambda0)))
        return(-value)
    }
}
