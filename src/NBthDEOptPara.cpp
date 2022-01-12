#include <cmath>
#include "GeoDiff.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <chrono>
#include <iostream>


using namespace std::chrono;


// Use auto keyword to avoid typing long
// type definitions to get the timepoint
// at this instant use function now()

// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace Rcpp;
using namespace roptim;

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// trying to make the vectorized 
arma::vec ref_dnbinom_mu_vec(const arma::vec &y, 
                             const double r, 
                             const arma::vec &tmp1){
  int N = y.n_elem;
  arma::vec prob(N);

  for(int i=0; i<N; i++)
    prob(i) = R::dnbinom_mu(y(i), r, tmp1(i), 1);
  return(prob);
}

class NBthDE_paranll : public Functor {
public:
  arma::vec y; // wa column of another matrix?
  arma::mat X; // our data matrix?
  arma::vec alpha0; // what is alpha0?
  arma::vec alpha; // what is alpha0?
  arma::mat preci1; // what is alpha0?
  double preci2; // what is alpha0?
  double threshold0; // what this a threashold for?
  arma::vec beta;
  arma::vec tmp0;
  arma::vec tmp1;
  
  double operator()(const arma::vec &x) override {
    int n = X.n_cols;
    beta = x(arma::span(0,n-1)); //subset of vecotr x
    //Rcout << "optimiser iteration \n";
    double r = x(n);
    
    if(std::isnan(r)){ // in some instances r can go to nan, this catches that
      throw 20;
    }
    
    double threshold = x(n+1);
    

    Rcout << "duration llk: "<< duration.count() << " \n";
    //auto starttmp0 = high_resolution_clock::now();
    //tmp0 = arma::exp2(X*beta);
    //auto stoptmp0 = high_resolution_clock::now();
    //auto durationtmp0 = duration_cast<microseconds>(stoptmp0 - starttmp0);
    //Rcout << "duration tmp0: "<< durationtmp0.count() << " \n";
    
    //auto starttmp01 = high_resolution_clock::now();
    arma::vec tmpneg1 = X*beta;
    tmp0 = arma::zeros<arma::vec>(tmpneg1.n_elem);
    for(int l=0;l<tmpneg1.n_elem;l++){
      tmp0(l) = pow(2.0, tmpneg1(l));
    }
    //auto stoptmp01 = high_resolution_clock::now();
    //auto durationtmp01 = duration_cast<microseconds>(stoptmp01 - starttmp01);
    //Rcout << "duration tmp01: "<< durationtmp01.count() << " \n";
    
    
    //auto starttmp1 = high_resolution_clock::now();
    tmp1 = alpha0*threshold+alpha%tmp0; // % here is element-wise multiplication
    //auto stoptmp1 = high_resolution_clock::now();
    //auto durationtmp1 = duration_cast<microseconds>(stoptmp1 - starttmp1);
    //Rcout << "duration tmp1: "<< durationtmp1.count() << " \n";
    
    //auto start = high_resolution_clock::now();
    arma::vec llk = ref_dnbinom_mu_vec(y, r, tmp1); //rate-limiting step
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(stop - start);
    //Rcout << "duration llk: "<< duration.count() << " \n";
    
    if (std::isinf(sum(llk))){
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
    
    //beta = x(arma::span(0,n-1));
    
    double r = x(n);
    double threshold = x(n+1);
    
    //auto starttmp2 = high_resolution_clock::now();
    arma::vec tmp2 = (y/tmp1-1.0)/(1.0+tmp1/r);
    //auto stoptmp2 = high_resolution_clock::now();
    //auto durationtmp2 = duration_cast<microseconds>(stoptmp2 - starttmp2);
    //Rcout << "duration tmp2: "<< durationtmp2.count() << " \n";
    
    //auto start1 = high_resolution_clock::now();
    gr(arma::span(0,n-1)) = -log(2.0)*X.t()*(tmp2%alpha%tmp0)+preci1.t()*beta;  //rate-limiting step
    //auto stop1 = high_resolution_clock::now();
    //auto duration1 = duration_cast<microseconds>(stop1 - start1);
    //Rcout << "duration transpose: "<< duration1.count() << " \n";

    //gr(arma::span(0,n-1)) = -log(2.0)*X.t()*(tmp2%alpha%tmp0)+preci1.t()*beta;  //rate-limiting step
    
    
    // auto start2 = high_resolution_clock::now();
    // arma::vec tmp4 = 1 + tmp1/r;//rate-limiting step

    // arma::vec pLr = arma::zeros<arma::vec>(tmp4.n_elem);
    // for(int k = 0; k < m; k++){
    //   pLr(k) = -log(tmp4(k));
    // }
    // auto stop2 = high_resolution_clock::now();
    // auto duration2 = duration_cast<microseconds>(stop2 - start2);
    // Rcout << "duration log matrix: "<< duration2.count() << " \n";

    // for(int k = 0; k < m; k++){
    //   for(int j = 0; j < y(k); j++){
    //     pLr(k) += 1.0/(j+r);
    //   }
    // }

    arma::vec tmp5 = 1 + tmp1/r;//rate-limiting step
    arma::vec pLr = arma::zeros<arma::vec>(tmp5.n_elem);
    for(int k = 0; k < m; k++){
      pLr(k) = -log(tmp5(k));
      for(int j = 0; j < y(k); j++){
        pLr(k) += 1.0/(j+r);
      }
    }
    pLr += -(y-tmp1)/(r+tmp1);
    gr(n) = -arma::sum(pLr);
    //auto stop4 = high_resolution_clock::now();
    //auto duration4 = duration_cast<microseconds>(stop4 - start4);
    //Rcout << "duration plr fancy: "<< duration4.count() << " \n";
    
    //auto start5 = high_resolution_clock::now();
    arma::vec tmp52 = 1 + tmp1/r;//rate-limiting step
    double running_total = 0;
    for(int k = 0; k < m; k++){
      running_total += -log(tmp52(k));
      for(int j = 0; j < y(k); j++){
        running_total += 1.0/(j+r);
      }
      running_total += -(y(k)-tmp1(k))/(r+tmp1(k));
    }
    
    //auto stop5 = high_resolution_clock::now();
    //auto duration5 = duration_cast<microseconds>(stop5 - start5);
    //Rcout << "duration just sum: "<< duration5.count() << " \n";
    //Rcout << "gr(n): "<< gr(n) << " \n";
    //Rcout << "-running_total: "<< -running_total << " \n";
    
    //auto start4tmp3 = high_resolution_clock::now();
    arma::mat tmp3 = tmp2.t()*alpha0;
    
    gr(n+1) = -tmp3(0,0)+(threshold-threshold0)*preci2;
    //auto stop4tmp3 = high_resolution_clock::now();
    //auto duration4tmp3 = duration_cast<microseconds>(stop4tmp3 - start4tmp3);
    //Rcout << "duration plr fancy: "<< duration4tmp3.count() << " \n";

    
  }
  
};

// [[Rcpp::export]]
List NBthDE_paraOptfeat(arma::mat &X, //define these terms
                        const arma::vec &y,
                        arma::vec alpha0,
                        arma::vec alpha,
                        arma::mat &preci1,
                        double threshold0,
                        double preci2,
                        arma::vec &x0,
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
List NBthDE_paraOptall(arma::sp_mat &Y,
                       arma::mat &X,
                       arma::vec &alpha0,
                       arma::vec &alpha,
                       arma::mat &preci1,
                       arma::vec &threshold0,
                       double preci2,
                       arma::vec &x0,
                       bool sizescale,
                       bool calhes){
  
  int n = X.n_cols;
  int m = Y.n_cols;
  int n_rows = Y.n_rows;
  arma::mat par(n+2,m);
  List hes(m);
  
  int vector_time = 0;
  int optim_time = 0;
  
  arma::vec conv(m);
  Rcout << "number of columns: "<< m << " \n";
  int failcount = 0;
  if(sizescale){
    for(int i=0; i < m; i++){
      try{
        //Rcout << i << "\n";
        auto start = high_resolution_clock::now();
        arma::vec Ycol(n_rows);
        for(int k=0; k<n_rows; k++){
          if(Y(k,i)!=0){
            Ycol(k)=Y(k,i);
          }
        }
        //Rcout << "Ycol \n";
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        vector_time += duration.count();
        auto start1 = high_resolution_clock::now();
        List result = NBthDE_paraOptfeat(X, Ycol,
                                         threshold0(i)*alpha0, threshold0(i)*alpha,
                                         preci1, 1.0, preci2, x0, calhes);
        auto stop1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(stop1 - start1);
        optim_time += duration1.count();
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
