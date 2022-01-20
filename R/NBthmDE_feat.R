fitNBthmDEfeat <- function(Y, probenum, X, Z, sizefact, sizefact0, preci1, threshold_mean, preci2, Lambdati, mapping, cluster_size,
    temp_size, rl, rt, sizescalebythreshold, cRandom) {
    nmh_sq <- cRandom$nmh_sq
    thin_sq <- cRandom$thin_sq

    thetapri <- cRandom$thetapri
    nu <- cRandom$nu
    lower <- cRandom$lower
    upper <- cRandom$upper
    iterations <- cRandom$iterations
    useprior <- cRandom$useprior
    preciu <- cRandom$preciu


    para_fixmat <- matrix(0, iterations, ncol(X) + 2)
    theta_mat <- matrix(0, iterations, max(rt$Lind))
    if (missing(preci1)) {
        preci1 <- diag(1, ncol(X))
    }

    Im <- NA

    conv <- FALSE
    dive <- FALSE

    if (sizescalebythreshold) {
        sizefact <- sizefact * threshold_mean * probenum
        sizefact0 <- sizefact0 * threshold_mean * probenum
        threshold_mean <- 1.0
    } else {
        threshold_mean <- threshold_mean * probenum
    }

    para_fix0 <- para_fix <- c(numeric(ncol(X)), 1, threshold_mean)


    Sigma <- as.matrix(Lambdati) %*% t(as.matrix(Lambdati))

    for (iter in seq_len(iterations)) {
        Tem <- Sigma[seq_len(temp_size), seq_len(temp_size), drop = FALSE]


        U <- NBthmDE_uOpt(numeric(rl), X, Z, Y, sizefact0, sizefact, para_fix, preciu, 0)$par
        # print(U)
        ## Metropolis Hasting

        for (iter_conv in seq_len(5)) {
            if (!conv) {
                Umat <- NBthmDE_mh(
                    Tem, U, X, as.matrix(Z), Y, sizefact0,
                    sizefact, para_fix, nmh_sq[iter]
                )
            } else {
                Umat <- NBthmDE_mh(
                    Tem, U, X, as.matrix(Z), Y, sizefact0,
                    sizefact, para_fix, nmh_sq[iterations]
                )
            }
            #    print(colMeans(Umat))
            #
            ## MLE for covariance

            if (useprior) {
                covnll <- NBthm_rcovnlliWp(Umat, Lambdati, mapping, nu = nu, diag(thetapri, rl))
            } else {
                covnll <- NBthm_rcovnll(Umat, Lambdati, mapping)
            }

            result <- tryCatch(
                expr = {
                    optim(rt$theta, covnll, lower = lower, upper = upper, method = "L-BFGS-B")
                },
                error = function(e) {
                    message(sprintf("`optim` encountered singularity problems. Re-assigning the initial values (%s/5).", iter_conv))
                    return(NULL)
                }
            )

            if (!is.null(result)) {
                break
            }
        }

        if (is.null(result)) {
            return(list(
                para_fix = NA,
                theta = NA,
                Im = NA,
                conv = NA,
                Uvec = NA
            ))
        }


        theta <- result$par
        # print(theta)
        theta_mat[iter, ] <- theta
        # print(theta)
        Lambdati@x[] <- mapping(theta)

        Sigma <- as.matrix(Lambdati) %*% t(as.matrix(Lambdati))
        ## MLE for parameters for fixed effects

        if (!conv) {
            result <- NBthmDE_fparaOptfeat(X, Z, Y, sizefact0, sizefact, preci1, threshold_mean, preci2, Umat[((seq_len(thin_sq[iter])) * floor(nmh_sq[iter] / thin_sq[iter])), ], para_fix, iter == iterations)
        } else {
            result <- NBthmDE_fparaOptfeat(X, Z, Y, sizefact0, sizefact, preci1, threshold_mean, preci2, Umat[((seq_len(thin_sq[iterations])) * floor(nmh_sq[iterations] / thin_sq[iterations])), ], para_fix, 1)
        }


        para_fix <- result$par

        para_fixmat[iter, ] <- para_fix
        if (conv | (iter == iterations)) {
            hes <- result$hes
            gr0 <- NBthmDE_gradM(
                Y, X, Z, para_fix, Umat[((seq_len(thin_sq[iterations])) * floor(nmh_sq[iterations] / thin_sq[iterations])), ],
                sizefact0, sizefact, preci1, preci2, threshold_mean
            )
            grt <- plyr::alply(gr0, 2, function(x) x %*% t(x))
            a <- array(unlist(grt), c(ncol(X) + 2, ncol(X) + 2, length(grt)))
            grm <- apply(a, seq_len(2), mean)

            Im <- hes - grm
        }

        if (conv) break


        if (sum((para_fix - para_fix0)^2) < 5e-3) {
            conv <- TRUE
        }
        if (mean(diag(as.matrix(Lambdati))) > 19.98) {
            dive <- TRUE
        }

        # print(iter)
        if (dive) {
            Im <- NA
            break
        }

        para_fix0 <- para_fix
    }
    Tem <- Sigma[seq_len(temp_size), seq_len(temp_size), drop = FALSE]
    Temupper <- Tem
    Temupper[lower.tri(Tem)] <- 0
    varcov <- methods::as(Temupper, "sparseMatrix")@x[]
    # print(mconv)
    message(sprintf("conv = %s", conv))
    message(sprintf("theta = %s", paste0(round(theta, 4), collapse = " ")))
    message(sprintf("varcov = %s", paste0(round(varcov, 4), collapse = " ")))
    message(sprintf("Iteration = %s", iter))

    #  plot(theta_mat[200:300,1])

    return(list(
        para_fix = para_fix,
        theta = theta,
        Im = Im,
        conv = conv,
        varcov = varcov,
        Uvec = colMeans(Umat)
    ))
}
