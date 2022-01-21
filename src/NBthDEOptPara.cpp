
#include <cmath>
#include "GeoDiff.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <chrono>
#include <iostream>

using namespace std::chrono;

// some of the parameters are not clearly named (Ned doesn't know what biological significance they have)
// (see NBthDE_model_description.pdf in this repo).
// this code defines an optimisation function (NBthDE_paranll)
// defines a function, given some data, that returns the optimised parameters of that function (NBthDE_paraOptfeat)
// and defines a function to iterate through all the columns of a given 

// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace Rcpp;
using namespace roptim;

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

arma::vec approx_dnbinom_mu_vec(const arma::vec &x_vec, 
                             const double s, 
                             const arma::vec &mu_vec){
  //function based on reading the #$%^& R documents of dnbinom!!!
  //everything is already logarithmic - separates out what can be
  //precalculated
  // prob(x,s,mu) = ((s+x-1) choose (s-1))*p^s *(1-p)^x
  
  int N = x_vec.n_elem;
  arma::vec prob(N);
  double s_binom = lgamma(s);
  double p =0;
  for(int i=0; i<N; i++){
    if (x_vec(i)==0 && mu_vec(i)==0){
      prob(i) = 0;
    }
    else{
      p = s/(s+mu_vec(i));
      prob(i) = lgamma(s+x_vec(i)) - s_binom - lgamma(x_vec(i)+1) + s*log(p) + x_vec(i)*log(1-p);
    }
    if(std::isinf(prob(i))){ // in some instances prob can go to -inf, this catches that
      prob(i) = 0; //check this doesn't wreck everything
    }
    if(std::isnan(prob(i))){ // in some instances prob can go to nan, this catches that
      prob(i) = 0; //check this doesn't wreck everything
    }
  }
  return(prob);
}

class NBthDE_paranll : public Functor {
public:
  arma::vec y; //column of t(object[features_high, ]) in R code (fitNBthDE_funct in NBthDE.R)
  arma::mat X; // X = model.matrix(form, data = annot) in R code (fitNBthDE_funct in NBthDE.R)
  arma::vec alpha0; // sizefact_BG in R code
  arma::vec alpha; // sizefact in R code
  arma::mat preci1; // related to preci1con (see NBthDE.R lines 311-317)
  double preci2; // (preci2 in fitNBthDE_funct in NBthDE.R (defaults to 10000))
  double threshold0; // threshold_mean * probenum[features_high]
  arma::vec beta; //this vector is shared between the operator function and the gradient
  arma::vec tmp0; // this vector is shared between the operator function and the gradient 
  arma::vec tmp1; // vector is shared between the operator function and the gradient
  
  double operator()(const arma::vec &x) override {
    int n = X.n_cols;
    beta = x(arma::span(0,n-1)); //subset of vecotr x
    double r = x(n);
    
    if(std::isnan(r)){ // in some instances r can go to nan, this catches that
      throw 20;
    }
    
    double threshold = x(n+1);
    
    arma::vec tmpneg1 = X*beta; // calculate tmp0_i = 2^(X_ij * beta_j) in stages to speed it up 
    tmp0 = arma::zeros<arma::vec>(tmpneg1.n_elem);
    for(int l=0;l<tmpneg1.n_elem;l++){
      tmp0(l) = pow(2.0, tmpneg1(l));
    }
    tmp1 = alpha0*threshold+alpha%tmp0; // % here is element-wise multiplication
    
    arma::vec llk = approx_dnbinom_mu_vec(y, r, tmp1); //rate-limiting step

    if (std::isinf(sum(llk))){
      throw 20;
    }
    
    if (std::isnan(sum(llk))){
      throw 20;
    }
    arma::mat pen10 = beta.t()*preci1*beta;
    
    double pen1 = pen10(0,0)/2.0;
    
    double result = -arma::sum(llk)+pen1+(1.0/2.0)*pow((threshold-threshold0),2)*preci2;
    
    return(result);
  }
 
  void Gradient(const arma::vec &x, arma::vec &gr) override {
    int n = X.n_cols;
    int m = y.n_elem;
    
    gr = arma::zeros<arma::vec>(n+2);
    
    double r = x(n);
    double threshold = x(n+1);
    
    arma::vec tmp2 = (y/tmp1-1.0)/(1.0+tmp1/r);
    
    gr(arma::span(0,n-1)) = -log(2.0)*X.t()*(tmp2%alpha%tmp0)+preci1.t()*beta;  //rate-limiting step
    
    arma::vec tmp4 = 1 + tmp1/r;//rate-limiting step
    arma::vec pLr = arma::zeros<arma::vec>(tmp4.n_elem);
    for(int k = 0; k < m; k++){
      pLr(k) = -log(tmp4(k));
      for(int j = 0; j < y(k); j++){
        pLr(k) += 1.0/(j+r);
      }
    }

    
    pLr += -(y-tmp1)/(r+tmp1);
    gr(n) = -arma::sum(pLr);
    arma::mat tmp3 = tmp2.t()*alpha0;
    
    gr(n+1) = -tmp3(0,0)+(threshold-threshold0)*preci2;
    
  }
  
};
                         
// [[Rcpp::export]]
List NBthDE_paraOptfeat(arma::mat& X, //X = model.matrix(form, data = annot) in R code (fitNBthDE_funct in NBthDE.R)
                        arma::vec y, //data vector part of t(object[features_high, ]) in R code (fitNBthDE_funct in NBthDE.R)
                        arma::vec alpha0, //sizefact_BG in R code
                        arma::vec alpha,//sizefact in R code
                        arma::mat& preci1, //related to preci1con (see NBthDE.R lines 311-317)
                        double threshold0, //threshold_mean * probenum[features_high]
                        double preci2, //(preci2 in fitNBthDE_funct in NBthDE.R (defaults to 10000))
                        arma::vec& x0, //startpara <- c(rep(0, ncol(X)), 1, (1.0 or threshold_mean))
                        bool calhes) { //(iter == iterations) (see fitNBthDE_funct in NBthDE.R for details)
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
  
  
  opt.set_hessian(calhes);
  
  opt.control.pgtol=1e-3;
  arma::vec x = x0;
  
  opt.minimize(f, x);
  
  
  return List::create(Named("par") = opt.par(),
                      Named("conv") = opt.convergence(),
                      Named("hes") = opt.hessian());
  
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List NBthDE_paraOptall(arma::sp_mat &Y, //t(object[features_high, ]) in R code (fitNBthDE_funct in NBthDE.R)
                       arma::mat &X, // X = model.matrix(form, data = annot) in R code (fitNBthDE_funct in NBthDE.R)
                       arma::vec &alpha0, //sizefact_BG in R code
                       arma::vec &alpha, //sizefact in R code
                       arma::mat &preci1, //related to preci1con (see NBthDE.R lines 311-317)
                       arma::vec &threshold0, //threshold_mean * probenum[features_high]
                       double preci2, //(preci2 in fitNBthDE_funct in NBthDE.R (defaults to 10000))
                       arma::vec &x0, //startpara <- c(rep(0, ncol(X)), 1, (1.0 or threshold_mean))
                       bool sizescale, //sizescalebythreshold (defaults to FALSE)
                       bool calhes){ //(iter == iterations) (see fitNBthDE_funct in NBthDE.R for details)
  
  int n = X.n_cols;
  int m = Y.n_cols;
  int n_rows = Y.n_rows;
  arma::mat par(n+2,m);
  List hes(m);
  
  
  arma::vec conv(m);
  Rcout << "number of columns: "<< m << " \n";
  int failcount = 0;
  if(sizescale){
    for(int i=0; i < m; i++){
      try{
        arma::vec Ycol(n_rows);
        for(int k=0; k<n_rows; k++){
          if(Y(k,i)!=0){
            Ycol(k)=Y(k,i);
          }
        }
        
        List result = NBthDE_paraOptfeat(X, Ycol,
                                         threshold0(i)*alpha0, threshold0(i)*alpha,
                                         preci1, 1.0, preci2, x0, calhes);
        
        par.col(i) = (as<arma::vec>(result["par"]));
        hes[i] = result["hes"];
        conv(i) = result["conv"];
      }
      catch (...){
        failcount++;
      }
    }

  } else {
    for(int i=0; i < m; i++){
      try{
        arma::vec Ycol(n_rows);
        for(int k=0; k<n_rows; k++){
          Ycol(k)=Y(k,i);
        }
        List result = NBthDE_paraOptfeat(X, Ycol,
                                         alpha0, alpha,
                                         preci1, threshold0(i), preci2, x0, calhes);
        
        par.col(i) = (as<arma::vec>(result["par"]));
        hes[i] = result["hes"];
        conv(i) = result["conv"];
      }
      catch (...){
        failcount++;
      }
    }
  }
  Rcout << "failed columns: " << failcount << "\n";
  return List::create(Named("par") = par,
                      Named("conv") = conv,
                      Named("hes") = hes);
}
