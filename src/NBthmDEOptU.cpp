#include <cmath>
#include "GeoDiff.h"
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



class NBthDE_unll : public Functor {
public:
  arma::mat X;
  arma::mat Z;
  arma::vec y;
  arma::vec alpha0;
  arma::vec alpha;
  arma::vec x;
  arma::mat preciu;

  double operator()(const arma::vec &u) override {
    int n = X.n_cols;
    arma::vec beta = x(arma::span(0,n-1));
    double r = x(n);
    double threshold = x(n+1);
    arma::vec tmp = arma::exp2(X*beta+Z*u);
    arma::vec mu = alpha0*threshold+alpha%tmp;
    arma::mat pen10 = u.t()*preciu*u;
    double pen1 = pen10(0,0)/2.0;
    return(-arma::sum(dnbinom_mu_vec(y, r, mu, 1))+pen1);
  }
};

// // [[Rcpp::export]]
// double haha(arma::mat X, arma::mat Z, arma::vec y, arma::vec alpha0,
//             arma::vec alpha, arma::vec x, arma::mat preciu, arma::vec u){
//   pos_u f;
//   f.X=X;
//   f.Z=Z;
//   f.y=y;
//   f.alpha0 = alpha0;
//   f.alpha=alpha;
//   f.x=x;
//   f.preciu=preciu;
//   return(f(u));
// }






// [[Rcpp::export]]
List NBthmDE_uOpt(arma::vec& u0, arma::mat& X, arma::mat& Z, arma::vec& y, arma::vec& alpha0,
          arma::vec& alpha, arma::vec& x, arma::mat& preciu, bool calhes) {
  int len_u = u0.n_elem;
  NBthDE_unll f;
  f.X=X;
  f.Z=Z;
  f.y=y;
  f.alpha0 = alpha0;
  f.alpha = alpha;
  f.x=x;
  f.preciu = preciu;


  arma::vec lower = arma::ones<arma::vec>(len_u) * (-100);
  arma::vec upper = arma::ones<arma::vec>(len_u) * 100;


  Roptim<NBthDE_unll> opt("L-BFGS-B");
  opt.set_lower(lower);
  opt.set_upper(upper);
  //opt.control.maxit = maxit;

  opt.set_hessian(calhes);

 // opt.set_hessian(true);
  opt.control.pgtol=1e-3;
  arma::vec u = u0;
   // arma::zeros<arma::vec>(n+2);
  // x(arma::span(0,n-1))=arma::solve(X, arma::log2(y/alpha + 0.001));
  // x(n) = 1;
  // x(n+1)=threshold0;

  opt.minimize(f, u);

  // arma::mat hes = opt.hessian();
  // double hes_det = arma::log_det(hes);
  // double hes_det1 = arma::log_det(hes(arma::span(1,n-1), arma::span(1,n-1)));

  return List::create(Named("par") = opt.par(),
                      Named("conv") = opt.convergence(),
                      Named("hes") = opt.hessian());


}



/*** R
#save(mat, Z, alpha0, alpha, para_fix, Umat, threshold0, preci2, Y, Lambdati, file="testData.Rdata")
#load("testData.Rdata")
#
# preciu <- as.matrix(solve(Lambdati))
#
# haha(mat, Z, Y, alpha0,
#      alpha, para_fix, preciu, u)

#lik_fun2 <- lik_probe2(mat, Z, Y, alpha0=alpha0, alpha=alpha,  Umat[((1:10)*20),], threshold0, preci2)

#lik_fun2(para_fix0)
#
# microbenchmark::microbenchmark(test(mat, Z, Y, alpha0=alpha0, alpha=alpha, x=para_fix, threshold0, preci2, Umat[((1:20)*20),]), lik_fun2(para_fix))
#
#
#
# system.time(result1<- optim(para_fix0, lik_fun2, lower=c(rep(-Inf,ncol(mat)), 0.01,0.01),
#                 method="L-BFGS-B"))
#

#system.time(result2 <- mleprobe2(mat, Z, Y, alpha0, alpha, preci1, threshold0, preci2, Umat[((1:20)*10),], para_fix0, FALSE))

*/
