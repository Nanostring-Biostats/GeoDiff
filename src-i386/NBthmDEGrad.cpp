#include <cmath>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace Rcpp;
using namespace roptim;

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec NBthmDE_grad(arma::vec& y, arma::mat& X, arma::mat& Z, arma::vec& x, arma::vec& u, arma::vec& alpha0,
               arma::vec& alpha, arma::mat& preci1, double preci2, double threshold0){
  int n = X.n_cols;
  arma::vec gr = arma::zeros<arma::vec>(n+2);

  arma::vec beta = x(arma::span(0,n-1));
  double r = x(n);
  double threshold = x(n+1);

  arma::vec tmp0 = arma::exp2(X*beta+Z*u);
  arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
  arma::vec tmp2 = (y/tmp1-1.0)/(1.0+tmp1/r);

  gr(arma::span(0,n-1)) = (-log(2.0)*(tmp2%alpha%tmp0).t()*X+beta.t()*preci1).t();

  arma::vec pLr = -arma::log(1.0+tmp1/r);
  for(int i = 0; i < y.n_elem; i++){
    for(int j = 0; j < y(i); j++){
      pLr(i) += 1.0/(j+r);
    }
  }
  pLr += -(y-tmp1)/(r+tmp1);
  gr(n) = -arma::sum(pLr);

  arma::mat tmp3 = tmp2.t()*alpha0;
  gr(n+1) = -tmp3(0,0)+(threshold-threshold0)*preci2;
  return gr;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat NBthmDE_gradM(arma::vec& y, arma::mat& X, arma::mat& Z, arma::vec& x, arma::mat& u_mat, arma::vec& alpha0,
               arma::vec& alpha, arma::mat& preci1, double preci2, double threshold0){

    int n = X.n_cols;
    int nmh = u_mat.n_rows;
    int m = y.n_elem;



    arma::vec beta = x(arma::span(0,n-1));
    double r = x(n);
    double threshold = x(n+1);

    arma::mat tmp0 = Z*u_mat.t();

    arma::vec tmp00 = X*beta;

    tmp0.each_col() += tmp00;

    tmp0 = arma::exp2(tmp0);

    arma::mat tmp1=tmp0;

    tmp1.each_col() %= alpha;

    arma::vec tmp11 = alpha0*threshold;

    tmp1.each_col() += tmp11;

    arma::mat tmp2(m, nmh, arma::fill::zeros);

    arma::mat grm(n+2, nmh, arma::fill::zeros);

    arma::vec pLr(m);
    arma::mat tmp3(1,1);
    for(int i=0; i < nmh; i++){
      tmp2.col(i) = (y/tmp1.col(i)-1.0)/(1.0+tmp1.col(i)/r);
      grm(arma::span(0,n-1), i) = (-log(2.0)*(tmp2.col(i)%alpha%tmp0.col(i)).t()*X+beta.t()*preci1).t();
      pLr = -arma::log(1.0+tmp1.col(i)/r);
      for(int k = 0; k < y.n_elem; k++){
        for(int j = 0; j < y(k); j++){
          pLr(k) += 1.0/(j+r);
        }
      }
      pLr += -(y-tmp1.col(i))/(r+tmp1.col(i));
      grm(n, i) = -arma::sum(pLr);
      tmp3 = tmp2.col(i).t()*alpha0;
      grm(n+1, i) = -tmp3(0,0)+(threshold-threshold0)*preci2;
    }

    return grm;

}


/*** R

*/
