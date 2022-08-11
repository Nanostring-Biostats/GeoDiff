#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec dnbinom_mu_vec(arma::vec x, double sz, arma::vec mu, int lg){
  int N = x.n_elem;
  arma::vec prob(N);
  //Rcpp::dnbinom_mu(x, sz, mu, lg)
  for(int i=0; i<N; i++)
    prob(i) = R::dnbinom_mu(x(i), sz, mu(i), lg);

  return(prob);
}
