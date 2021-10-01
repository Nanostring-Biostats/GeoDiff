#include <cmath>  // std::pow
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





class NBthmDE_unll : public Functor {
public:
  arma::mat X;
  arma::mat Z;
  arma::vec y;
  arma::vec alpha0;
  arma::vec alpha;
  arma::vec x;

  double operator()(const arma::vec &u) override {
    int n = X.n_cols;
    arma::vec beta = x(arma::span(0,n-1));
    double r = x(n);
    double threshold = x(n+1);
    arma::vec tmp = arma::exp2(X*beta+Z*u);
    arma::vec mu = alpha0*threshold+alpha%tmp;

    return(arma::sum(dnbinom_mu_vec(y, r, mu, 1)));
  }
};

// // [[Rcpp::export]]
// double haha(arma::mat X, arma::mat Z, arma::vec y, arma::vec alpha0,
//             arma::vec alpha, arma::vec x, arma::vec u){
//   lik_probe_mh f;
//   f.X=X;
//   f.Z=Z;
//   f.y=y;
//   f.alpha0 = alpha0;
//   f.alpha=alpha;
//   f.x=x;
//
//   return(f(u));
// }

// [[Rcpp::export]]
List condi_u(arma::mat& Tem, int ind, arma::vec& u,
           int temp_size){
  float mu;
  float Sig;

  if(temp_size==1) {
    mu = 0;
    Sig = Tem(0,0);
  } else {
    int cluster_ind = (ind+temp_size-1)/temp_size;
    //Rcout << "cluster_ind:" << cluster_ind << "\n";
    int temp_ind = ind - (cluster_ind-1)*temp_size;
    //Rcout << "temp_ind:" << temp_ind << "\n";
    arma::vec u_temp = u(arma::span(((cluster_ind-1)*temp_size),(cluster_ind*temp_size-1)));
    //Rcout << "u_temp:" << u_temp << "\n";
    arma::vec ind_vec = arma::linspace(0,(temp_size-1), temp_size);
    arma::uvec lef = arma::find(ind_vec!=(temp_ind-1));
    arma::uvec temp_ind_v(1);
    temp_ind_v.fill(temp_ind-1);
    //Rcout << "lef:" << lef << "\n";
    //arma::mat mu = Tem.submat(temp_ind_v, lef)*arma::solve(Tem.submat(lef,lef), u_temp.elem(lef)));
    //arma::mat Sig = Tem.submat(temp_ind_v,temp_ind_v)+Tem.submat(temp_ind_v, lef)*arma::solve(Tem.submat(lef,lef), Tem.submat(lef, temp_ind_v));
    arma::mat mu0 = Tem.submat(temp_ind_v, lef)*arma::solve(Tem.submat(lef,lef), u_temp.elem(lef));
    arma::mat Sig0 = Tem.submat(temp_ind_v,temp_ind_v)+Tem.submat(temp_ind_v, lef)*arma::solve(Tem.submat(lef,lef), Tem.submat(lef, temp_ind_v));
    mu = mu0(0,0);
    Sig = Sig0(0,0);
  }

  //Rcout << "mu:" << mu << "\n";
  //Rcout << "Sig:" << Sig << "\n";
  // arma::rowvec mu2 = mu(lef);
  // arma::vec mu3 = mu2.elem(lef);
  // Rcout << "mu2:" << mu2 << "\n";
  // //
  return List::create(Named("mu") = mu,
                      Named("Sig") = Sig);

}



// [[Rcpp::export]]
arma::mat NBthmDE_mh(arma::mat& Tem, arma::vec& u, arma::mat& X, arma::mat& Z, arma::vec& y, arma::vec& alpha0,
                 arma::vec& alpha, arma::vec& x, int nmh){
  int len_u = u.n_elem;
  int temp_size = Tem.n_rows;
  float v=0;
  arma::vec u_tmp = u;
  NBthmDE_unll f;
  f.X=X;
  f.Z=Z;
  f.y=y;
  f.alpha0 = alpha0;
  f.alpha=alpha;
  f.x=x;

  arma::mat u_mat(nmh, len_u, arma::fill::zeros);

  for(int i=0; i < nmh; i++){
    for (int j=0; j < len_u; j++){
      u_tmp = u;
      List condidis = condi_u(Tem, (j+1) ,u, temp_size);
      float mu = condidis["mu"];
//      Rcout << "mu" << mu << "\n";
      float Sig = condidis["Sig"];
//      Rcout << "Sig" << Sig << "\n";
      u_tmp(j) = R::rnorm(mu, sqrt(Sig));
//      Rcout << "u_tmp(j):" << u_tmp(j) << "\n";

      float ratio = exp(f(u_tmp)-f(u));
      if(ratio>1) {
        u(j) = u_tmp(j);
      } else{
          v=R::runif(0,1);
          if(v<ratio)
            u(j) = u_tmp(j);
        }

    }
    u_mat.row(i) = u.t();
  }


    return(u_mat);

}



