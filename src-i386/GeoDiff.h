#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

arma::vec dnbinom_mu_vec(arma::vec x, double sz, arma::vec mu, int lg);
