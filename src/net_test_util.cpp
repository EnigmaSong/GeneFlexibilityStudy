#include <RcppArmadillo.h>
using namespace Rcpp;

// Collection of Utility functions for net_test()
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.E_uutA_uppertri)]]
double E_uutA_uppertri(const arma::vec &u, const arma::mat &A){
  const int n = u.n_elem;
  const double n_double = (double) n;
  double res=0.0;
  
  for(int i=0; i< n-1; i++){
    for(int j = i+1; j<n; j++){
      res += u(i)*u(j)*A(i,j);
    }
  }
  res /= (n_double*(n_double-1.0)/2.0);
  
  return res;
}

// [[Rcpp::export(.E_a4a1)]]
double E_a4a1(const arma::vec &u, const arma::mat &A){
  const int n = u.n_elem;
  const double n_double = (double) n;
  double res=0.0;
  
  for(int i=0; i< n; i++){
    for(int j=0; j<n; j++){
      res += A(i,j)*u(j);
    }
    res -= A(i,i)*u(i);
  }
  res /= (n_double*(n_double-1.0));
  
  return res;
}
