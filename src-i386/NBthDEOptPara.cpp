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






class NBthDE_paranll : public Functor {
public:
  arma::vec y;
  arma::mat X;
  arma::vec alpha0;
  arma::vec alpha;
  arma::mat preci1;
  double preci2;
  double threshold0;


  double operator()(const arma::vec &x) override {
    int n = X.n_cols;
    int m = y.n_elem;
    arma::vec beta = x(arma::span(0,n-1));
    double r = x(n);
    double threshold = x(n+1);


    arma::vec tmp0 = arma::exp2(X*beta);
    arma::vec tmp1 = alpha0*threshold+alpha%tmp0;


    arma::vec llk = dnbinom_mu_vec(y, r, tmp1, 1);

 //   Rcout << "mean(llk.row(0))" << arma::mean(llk.row(0)) << "\n";
    //arma::mat pen
    arma::mat pen10 = beta.t()*preci1*beta;
    double pen1 = pen10(0,0)/2.0;
    //+nmh*(1.0/2.0)*pow((threshold-threshold0),2)*preci2
    return(-arma::sum(llk)+pen1+(1.0/2.0)*pow((threshold-threshold0),2)*preci2);
  }



  void Gradient(const arma::vec &x, arma::vec &gr) override {
    int n = X.n_cols;
    int m = y.n_elem;


    gr = arma::zeros<arma::vec>(n+2);

    arma::vec beta = x(arma::span(0,n-1));

    double r = x(n);
    double threshold = x(n+1);


    arma::vec tmp0 = arma::exp2(X*beta);

 //   Rcout << "tmp0:" << tmp0(arma::span(0,4)) << "\n";

    arma::vec tmp1 = alpha0*threshold+alpha%tmp0;

 //   Rcout << "tmp1:" << tmp1(arma::span(0,4)) << "\n";

    arma::vec tmp2 = (y/tmp1-1.0)/(1.0+tmp1/r);

 //   Rcout << "tmp2:" << tmp2(arma::span(0,4)) << "\n";

    gr(arma::span(0,n-1)) = (-log(2.0)*(tmp2%alpha%tmp0).t()*X+beta.t()*preci1).t();

    arma::vec pLr = -arma::log(1.0+tmp1/r);

    for(int k = 0; k < y.n_elem; k++){
      for(int j = 0; j < y(k); j++){
        pLr(k) += 1.0/(j+r);
      }
    }

    pLr += -(y-tmp1)/(r+tmp1);

    gr(n) = -arma::sum(pLr);

    arma::mat tmp3 = tmp2.t()*alpha0;
    gr(n+1) = -tmp3(0,0)+(threshold-threshold0)*preci2;

//    Rcout << "gr:" << gr << "\n";
    }

};

// // [[Rcpp::export]]
// double test(arma::mat X, arma::vec y, arma::vec& alpha0,
//             arma::vec alpha, arma::vec x, arma::mat preci1,
//             double threshold0, double preci2){
//   NBthDE_paranll f;
//   f.X=X;
//   f.y=y;
//   f.alpha0 = alpha0;
//   f.alpha = alpha;
//   f.preci1=preci1;
//   f.threshold0=threshold0;
//   f.preci2=preci2;
//
//   return(f(x));
// }
//
// // [[Rcpp::export]]
// arma::vec test2(arma::mat X, arma::vec y, arma::vec alpha0,
//             arma::vec alpha, arma::vec x, arma::mat preci1,
//             double threshold0, double preci2){
//   int n = X.n_cols;
//   int m = y.n_elem;
//
//
//   arma::vec gr = arma::zeros<arma::vec>(n+2);
//
//   arma::vec beta = x(arma::span(0,n-1));
//
//   double r = x(n);
//   double threshold = x(n+1);
//
//
//   arma::vec tmp0 = arma::exp2(X*beta);
//
//   //   Rcout << "tmp0:" << tmp0(arma::span(0,4)) << "\n";
//
//   arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
//
//   //   Rcout << "tmp1:" << tmp1(arma::span(0,4)) << "\n";
//
//   arma::vec tmp2 = (y/tmp1-1.0)/(1.0+tmp1/r);
//
//   //   Rcout << "tmp2:" << tmp2(arma::span(0,4)) << "\n";
//
//   gr(arma::span(0,n-1)) = (-log(2.0)*(tmp2%alpha%tmp0).t()*X+beta.t()*preci1).t();
//
//   arma::vec pLr = -arma::log(1.0+tmp1/r);
//
//   for(int k = 0; k < y.n_elem; k++){
//     for(int j = 0; j < y(k); j++){
//       pLr(k) += 1.0/(j+r);
//     }
//   }
//
//   pLr += -(y-tmp1)/(r+tmp1);
//
//   gr(n) = -arma::sum(pLr);
//
//   arma::mat tmp3 = tmp2.t()*alpha0;
//   gr(n+1) = -tmp3(0,0)+(threshold-threshold0)*preci2;
//
//   Rcout << "gr:" << gr << "\n";
//
//   return(gr);
// }



// [[Rcpp::export]]
List NBthDE_paraOptfeat(arma::mat& X,
               arma::vec y,
               arma::vec alpha0,
               arma::vec alpha,
               arma::mat& preci1,
               double threshold0,
               double preci2,
               arma::vec& x0,
               bool calhes) {
  NBthDE_paranll f;
  f.X=X;
  f.y=y;
  f.alpha0 = alpha0;
  f.alpha = alpha;
  f.preci1=preci1;
  f.threshold0=threshold0;
  f.preci2=preci2;

  int n = X.n_cols;


  arma::vec lower = arma::ones<arma::vec>(n+2) * (-100);
  lower(n) = 0.01;
  lower(n+1) = 0.01;
  arma::vec upper = arma::ones<arma::vec>(n+2) * 100;
  upper(n) = 1000;
  upper(n+1) = 1000000;

  Roptim<NBthDE_paranll> opt("L-BFGS-B");
  opt.set_lower(lower);
  opt.set_upper(upper);
  //opt.control.maxit = maxit;

  opt.set_hessian(calhes);

 // opt.set_hessian(true);
  opt.control.pgtol=1e-3;
  arma::vec x = x0;
   // arma::zeros<arma::vec>(n+2);
  // x(arma::span(0,n-1))=arma::solve(X, arma::log2(y/alpha + 0.001));
  // x(n) = 1;
  // x(n+1)=threshold0;

  opt.minimize(f, x);

  // arma::mat hes = opt.hessian();
  // double hes_det = arma::log_det(hes);
  // double hes_det1 = arma::log_det(hes(arma::span(1,n-1), arma::span(1,n-1)));

  return List::create(Named("par") = opt.par(),
                      Named("conv") = opt.convergence(),
                      Named("hes") = opt.hessian());


}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List NBthDE_paraOptall(arma::mat& Y,
                   arma::mat& X,
                   arma::vec& alpha0,
                   arma::vec& alpha,
                   arma::mat& preci1,
                   arma::vec& threshold0,
                   double preci2,
                   arma::vec& x0,
                   bool sizescale,
                   bool calhes){

  int n = X.n_cols;
  int m = Y.n_cols;
  arma::mat par(n+2,m);
  List hes(m);


  arma::vec conv(m);

  if(sizescale){
    for(int i=0; i < m; i++){
      List result = NBthDE_paraOptfeat(X, Y.col(i),
                                       threshold0(i)*alpha0, threshold0(i)*alpha,
                                       preci1, 1.0, preci2, x0, calhes);

      par.col(i) = (as<arma::vec>(result["par"]));
      hes[i] = result["hes"];
      conv(i) = result["conv"];
    }

  } else {
    for(int i=0; i < m; i++){
      List result = NBthDE_paraOptfeat(X, Y.col(i),
                             alpha0, alpha,
                             preci1, threshold0(i), preci2, x0, calhes);

      par.col(i) = (as<arma::vec>(result["par"]));
      hes[i] = result["hes"];
      conv(i) = result["conv"];
    }
  }
  return List::create(Named("par") = par,
                      Named("conv") = conv,
                      Named("hes") = hes);
}



/*** R
#save(mat, Z, alpha0, alpha, para_fix ,Umat, threshold0, preci2, file="testData.Rdata")
#load("testData.Rdata")
#test(mat, Z, Y, alpha0=alpha0,
#     alpha=alpha, x=para_fix0, NBmod$preci1, threshold0, preci2, Umat[((1:10)*20),])


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
