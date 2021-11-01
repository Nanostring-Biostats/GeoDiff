#include <cmath>  // std::pow

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace Rcpp;
using namespace roptim;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double poisprobe0(NumericVector y, arma::mat X, arma::vec x,arma::vec alpha0,
//                   arma::vec alpha, arma::mat preci1, double preci2, double threshold0)
// {  int N=y.length();
//   arma::vec tmp(N);
//   NumericVector z(1);
//
//
//   int n = X.n_cols;
//   arma::vec beta = x(arma::span(0,n-1));
//   double threshold = x(n);
//   arma::vec tmp0 = arma::exp2(X*beta);
//   arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
//
//   for(int i=0; i < N; i++){
//     z[0] = y[i];
//
//     tmp(i) = dpois(z, tmp1(i), true)[0];
//   }
//
//   arma::mat pen10 = beta.t()*preci1*beta;
//   double pen1 = pen10(0,0);
//   return -arma::sum(tmp)+pen1/2+pow((threshold-threshold0),2)*preci2/2;}
//
//
//
//
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double poisprobe00(arma::vec y, arma::mat X, arma::vec x,arma::vec alpha0,
//                    arma::vec alpha, arma::mat preci1, double preci2, double threshold0)
// {  int N=y.n_elem;
//   arma::vec tmp(N);
//
//
//
//   int n = X.n_cols;
//   arma::vec beta = x(arma::span(0,n-1));
//   double threshold = x(n);
//   arma::vec tmp0 = arma::exp2(X*beta);
//   arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
//
//   tmp = y%log(tmp1)-tmp1;
//
//   arma::mat pen10 = beta.t()*preci1*beta;
//   double pen1 = pen10(0,0);
//   return -arma::sum(tmp)+pen1/2+pow((threshold-threshold0),2)*preci2/2;}


// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::vec poisprobegr(arma::vec y, arma::mat X, arma::vec x,arma::vec alpha0,
//                       arma::vec alpha, arma::mat preci1, double preci2, double threshold0)
// { int n = X.n_cols;
//   arma::vec gr = arma::zeros<arma::vec>(n+1);
//
//   arma::vec beta = x(arma::span(0,n-1));
//   double threshold = x(n);
//
//   arma::vec tmp0 = arma::exp2(X*beta);
//   arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
//   arma::vec tmp2 = (y/tmp1-1);
//
//   gr(arma::span(0,n-1)) = (-log(2.0)*(tmp2%alpha%tmp0).t()*X+beta.t()*preci1).t();
//
//   arma::mat tmp3 = tmp2.t()*alpha0;
//   gr(n) = -tmp3(0,0)+(threshold-threshold0)*preci2;
//   return gr;}




class PoisthNorm_paranll : public Functor {
public:
  arma::vec y;
  arma::mat X;
  arma::vec alpha0;
  arma::vec alpha;
  arma::mat preci1;
  double preci2;
  double threshold0;

  double operator()(const arma::vec &x) override {
    int N=y.n_elem;
    arma::vec tmp(N);



    int n = X.n_cols;
    arma::vec beta = x(arma::span(0,n-1));
    double threshold = x(n);
    arma::vec tmp0 = arma::exp2(X*beta);
    arma::vec tmp1 = alpha0*threshold+alpha%tmp0;

    tmp = y%log(tmp1)-tmp1;

    arma::mat pen10 = beta.t()*preci1*beta;
    double pen1 = pen10(0,0);
    return -arma::sum(tmp)+pen1/2.0+pow((threshold-threshold0),2)*preci2/2.0;
  }

  void Gradient(const arma::vec &x, arma::vec &gr) override {
    int n = X.n_cols;
    gr = arma::zeros<arma::vec>(n+1);

    arma::vec beta = x(arma::span(0,n-1));
    double threshold = x(n);

    arma::vec tmp0 = arma::exp2(X*beta);
    arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
    arma::vec tmp2 = (y/tmp1-1);
    gr(arma::span(0,n-1)) = (-log(2.0)*(tmp2%alpha%tmp0).t()*X+beta.t()*preci1).t();

    arma::mat tmp3 = tmp2.t()*alpha0;
    gr(n) = -tmp3(0,0)+(threshold-threshold0)*preci2;

  }

};




// [[Rcpp::export]]
List PoisthNorm_paraOptfeat(arma::vec y,
                  arma::mat& X,
                  arma::vec alpha0,
                  arma::vec alpha,
                  arma::mat& preci1,
                  double preci2,
                  double threshold0,
                  bool calhes) {
  PoisthNorm_paranll f;
  f.X = X;
  f.y = y;
  f.alpha = alpha;
  f.alpha0 = alpha0;
  f.preci1 = preci1;
  f.preci2 = preci2;
  f.threshold0 = threshold0;

  int n = X.n_cols;


  arma::vec lower = arma::ones<arma::vec>(n+1) * (-50);
  lower(n) = 0.01;
  arma::vec upper = arma::ones<arma::vec>(n+1) * 50;
  upper(n) = 1000000;

  Roptim<PoisthNorm_paranll> opt("L-BFGS-B");
  opt.set_lower(lower);
  opt.set_upper(upper);
  opt.control.maxit = 500;
  opt.set_hessian(calhes);
  opt.control.pgtol=1e-8;
  arma::vec x = arma::zeros<arma::vec>(n+1);


  x(arma::span(0,n-1))=arma::solve(X, arma::log2(y/alpha + 0.001));

  x(n)=threshold0;

  opt.minimize(f, x);

  // double hes_det = arma::log_det(hes);
  // double hes_det1 = arma::log_det(hes(arma::span(1,n-1), arma::span(1,n-1)));

  return List::create(Named("par") = opt.par(),
                      Named("conv") = opt.convergence(),
                      Named("hes") = opt.hessian());
//                      Named("hes_det1") = arma::log_det(hes(arma::span(1,n-1), arma::span(1,n-1))));


}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List PoisthNorm_paraOptall(arma::mat& Y,
                     arma::mat& X,
                     arma::vec& alpha0,
                     arma::vec& alpha,
                     arma::mat& preci1,
                     arma::vec& threshold0,
                     double preci2,
                     bool sizescale,
                     bool calhes){

  int n = X.n_cols;
  int m = Y.n_cols;
  arma::mat par(n+1,m);
  // arma::vec hes_det(m);
  // arma::vec hes_det1(m);
  arma::vec conv(m);
  List hes(m);
  if(sizescale) {
    for(int i=0; i < m; i++){
      List result = PoisthNorm_paraOptfeat(Y.col(i),X,
                                           threshold0(i)*alpha0,
                                           threshold0(i)*alpha,
                                           preci1,
                                           preci2,
                                           1.0,
                                           calhes);
      par.col(i) = (as<arma::vec>(result["par"]));
      // hes_det(i) = result["hes_det"];
      // hes_det1(i) = result["hes_det1"];
      hes[i] = result["hes"];
      conv(i) = result["conv"];
    }
  } else {
    for(int i=0; i < m; i++){
      List result = PoisthNorm_paraOptfeat(Y.col(i),X,
                                           alpha0,
                                           alpha,
                                           preci1,
                                           preci2,
                                           threshold0(i),
                                           calhes);
      par.col(i) = (as<arma::vec>(result["par"]));
      // hes_det(i) = result["hes_det"];
      // hes_det1(i) = result["hes_det1"];
      hes[i] = result["hes"];
      conv(i) = result["conv"];
  }


  }
  return List::create(Named("par") = par,
                      Named("conv") = conv,
                      Named("hes") = hes);
                      // Named("hes_det") = hes_det,
                      // Named("hes_det1") = hes_det1);
}


// // NumericVector conv(arma::vec y){
// //   return as<NumericVector>(wrap(y));
// // }
//
//
