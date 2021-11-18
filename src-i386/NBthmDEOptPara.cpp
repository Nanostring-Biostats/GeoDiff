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






class NBthmDE_fparanll : public Functor {
public:
  arma::mat X;
  arma::mat Z;
  arma::vec y;
  arma::vec alpha0;
  arma::vec alpha;
  arma::mat u_mat;
  arma::mat preci1;
  double preci2;
  double threshold0;

  double operator()(const arma::vec &x) override {
    int n = X.n_cols;
    int nmh = u_mat.n_rows;
    int m = y.n_elem;
    arma::vec beta = x(arma::span(0,n-1));
    double r = x(n);
    double threshold = x(n+1);
    arma::mat tmp0 = Z*u_mat.t();

    arma::vec tmp1 = X*beta;

    tmp0.each_col() += tmp1;

    arma::mat tmp = arma::exp2(tmp0);

    tmp.each_col() %= alpha;

    arma::vec tmp2 = alpha0*threshold;

    tmp.each_col() += tmp2;


    arma::mat llk(m,nmh, arma::fill::zeros);

    for(int i=0; i < nmh; i++){
      llk.col(i) = dnbinom_mu_vec(y, r, tmp.col(i), 1);
    }
 //   Rcout << "mean(llk.row(0))" << arma::mean(llk.row(0)) << "\n";
    //arma::mat pen
    arma::mat pen10 = beta.t()*preci1*beta;
    double pen1 = pen10(0,0)/2.0;
    //+nmh*(1.0/2.0)*pow((threshold-threshold0),2)*preci2
    return(-arma::accu(llk)/static_cast<double>(nmh)+pen1+(1.0/2.0)*pow((threshold-threshold0),2)*preci2);
  }



  void Gradient(const arma::vec &x, arma::vec &gr) override {
    int n = X.n_cols;
    int nmh = u_mat.n_rows;
    int m = y.n_elem;

    gr = arma::zeros<arma::vec>(n+2);

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
      gr += grm.col(i)/static_cast<double>(nmh);
    }


  }
};

// // [[Rcpp::export]]
// double test(arma::mat X, arma::mat Z, arma::vec y, arma::vec alpha0,
//             arma::vec alpha, arma::vec x, arma::mat preci1,
//             double threshold0, double preci2, arma::mat u_mat){
//   lik_probe2_mh f;
//   f.X=X;
//   f.Z=Z;
//   f.y=y;
//   f.alpha0 = alpha0;
//   f.alpha = alpha;
//   f.u_mat=u_mat;
//   f.preci1=preci1;
//   f.threshold0=threshold0;
//   f.preci2=preci2;
//
//   return(f(x));
// }



// [[Rcpp::export]]
List NBthmDE_fparaOptfeat(arma::mat& X,
               arma::mat& Z,
               arma::vec& y,
               arma::vec& alpha0,
               arma::vec& alpha,
               arma::mat& preci1,
               double threshold0,
               double preci2,
               arma::mat& u_mat,
               arma::vec& x0,
               bool calhes) {
  NBthmDE_fparanll f;
  f.X=X;
  f.Z=Z;
  f.y=y;
  f.alpha0 = alpha0;
  f.alpha = alpha;
  f.u_mat=u_mat;
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

  Roptim<NBthmDE_fparanll> opt("L-BFGS-B");
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
