#include <RcppArmadillo.h>
using namespace Rcpp;

// Rcpp implmentation of Cai and Liu (2015), Large-Scale Multiple Testing of Correlations
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int min(int x, int y){
  if(x<y) return x;
  else return y;
}

// Cross-correlation implementation in armadillo
//' @rdname corFamily
//' @export
// [[Rcpp::export]]
arma::mat cor_arma(const arma::mat& X,
                   const arma::mat& Y) {
  arma::mat out(X.n_cols, Y.n_cols);
  
  // Degenerate case first
  if (X.n_cols == 0 || Y.n_cols == 0) {
    return out;
  } else if (X.n_rows == 0 || X.n_rows == 1 || Y.n_rows == 0 || Y.n_rows == 1) {
    out.fill(Rcpp::NumericVector::get_na());
  } else {
    //sample covariance
    out = arma::cor(X, Y, 0);
  }
  
  return out;
}

// [[Rcpp::export]]
arma::mat CL_stat_twosample(arma::mat& X, arma::mat& Y, double& kappa_x, double& kappa_y, int b_size = 100){
  double n_X = (double)X.n_rows;
  double n_Y = (double)Y.n_rows;
  int p_X = X.n_cols;
  double temp_a, temp_b;

  b_size = min(b_size, p_X);

  double threshold_X = 2*sqrt(log(p_X)/n_X);
  double threshold_Y = 2*sqrt(log(p_X)/n_Y);
  
  //for debug
  // Rcout << n_X << "\n";
  // Rcout << n_Y << "\n";
  // Rcout << p_X << "\n";
  // Rcout << b_size << "\n";
  // Rcout << threshold_X << "\n";
  // Rcout << threshold_Y << "\n";
  
  arma::mat T_stat = arma::zeros(X.n_cols,X.n_cols);
  arma::mat rho_1, rho_2;
  arma::mat trho_1, trho_2, trho_sq;

  // centering the data
  // Compute colmun-wise means.
  // X.each_row() -= mean(X);
  // Y.each_row() -= mean(Y);

  arma::mat x1, x2, y1, y2;

  int loop_len = (p_X - 1)/ b_size+1;
  //To reduce memory burden, use for loop to compute block-wise covariance (e.g. 100 by 100)
  // instead of computing the full corelation matrix.
  // It sacrifices speed, but runs.
  // 
  for(int i = 0; i < loop_len; i++){
    x1 = X.cols(b_size*i, min(b_size*(i+1), p_X)-1);
    y1 = Y.cols(b_size*i, min(b_size*(i+1), p_X)-1);
    for(int j = i ; j < loop_len; j++){
      x2 = X.cols(b_size*j, min(b_size*(j+1), p_X)-1);
      y2 = Y.cols(b_size*j, min(b_size*(j+1), p_X)-1);
      //empirical correlation
      rho_1 = cor_arma(x1,x2);
      rho_2 = cor_arma(y1,y2);
      // Rcout<< rho_1(0,1) << "\n";
      // Rcout<< rho_1(3,4) << "\n";
    
      //thresholded correlation
      trho_1 = rho_1 %(abs(rho_1)/sqrt(kappa_x/n_X*(1-rho_1%rho_1)%(1-rho_1%rho_1))> threshold_X);
      trho_2 = rho_2 %(abs(rho_2)/sqrt(kappa_y/n_Y*(1-rho_2%rho_2)%(1-rho_2%rho_2))> threshold_Y);
      //save the maximum of square of trho_1 and trho_2
      trho_sq = trho_1;// First create a matrix with the same size of trho_1
      for(int k = 0; k< min(b_size, p_X-b_size*i); k++){
        for(int l = 0; l< min(b_size, p_X-b_size*j); l++){
          temp_a = trho_1(k,l)*trho_1(k,l);
          temp_b = trho_2(k,l)*trho_2(k,l);
          if(temp_a> temp_b){
            trho_sq(k,l) = temp_a;
          }else{
            trho_sq(k,l) = temp_b;
          }
        }
      }
      // Rcout << trho_sq(0,1) << "\n";
      // Rcout << trho_sq(1,4) << "\n";
      // Rcout << trho_sq(2,3) << "\n";
      // Rcout << trho_sq(3,4) << "\n";
      //computing test statistics
      // T_stat= (cor_x-cor_y)/sqrt((kappa_x/n_x + kappa_y/n_y)*(1-pmax(tcor_x^2,tcor_y^2))^2)
      // T_stat((b_size*i):(b_size*i + min(b_size, p_X-b_size*i)),(b_size*j):(b_size*j + min(b_size, p_X-b_size*j)))= (rho_1 - rho_2)/sqrt((kappa_x/n_X+kappa_y/n_Y)*(1-trho%trho)%(1-trho%trho));
      for(int k = 0; k< min(b_size, p_X-b_size*i); k++){
        for(int l = 0; l< min(b_size, p_X-b_size*j); l++){
          T_stat(b_size*i + k, b_size*j + l) = (rho_1(k,l)-rho_2(k,l))/sqrt((kappa_x/n_X+kappa_y/n_Y)*(1-trho_sq(k,l))*(1-trho_sq(k,l)));
          // T_stat(b_size*j + l, b_size*i + k) = T_stat(b_size*i+k, b_size*j + l);
        }
      }
    }
  }
  
  //T_stat is a symmetric matrix, diagonal is 0
  for(int k = 0; k < p_X; k++){
    T_stat(k,k) = 0;
    for(int l = k+1; l < p_X; l++){
      T_stat(l,k) = T_stat(k,l);
    }
  }
  
  return T_stat;
}

// [[Rcpp::export]]
arma::mat CL_stat_onesample(arma::mat& X, double& kappa_x, int b_size = 100){
  double n_X = (double)X.n_rows;
  int p_X = X.n_cols;
  
  b_size = min(b_size, p_X);
  
  double threshold_X = 2*sqrt(log(p_X)/n_X);
  
  //for debug
  // Rcout << n_X << "\n";
  // Rcout << n_Y << "\n";
  // Rcout << p_X << "\n";
  // Rcout << b_size << "\n";
  // Rcout << threshold_X << "\n";
  // Rcout << threshold_Y << "\n";
  
  arma::mat T_stat = arma::zeros(X.n_cols,X.n_cols);
  arma::mat rho_1;
  arma::mat trho_1, trho_sq;
  
  // centering the data
  // Compute colmun-wise means.
  // X.each_row() -= mean(X);
  // Y.each_row() -= mean(Y);
  
  arma::mat x1, x2;
  
  int loop_len = (p_X - 1)/ b_size+1;
  //To reduce memory burden, use for loop to compute block-wise covariance (e.g. 100 by 100)
  // instead of computing the full corelation matrix.
  // It sacrifices speed, but runs.
  // 
  for(int i = 0; i < loop_len; i++){
    x1 = X.cols(b_size*i, min(b_size*(i+1), p_X)-1);
    for(int j = i ; j < loop_len; j++){
      x2 = X.cols(b_size*j, min(b_size*(j+1), p_X)-1);
      //empirical correlation
      rho_1 = cor_arma(x1,x2);
      // Rcout<< rho_1(0,1) << "\n";
      // Rcout<< rho_1(3,4) << "\n";
      
      //thresholded correlation
      trho_1 = rho_1 %(abs(rho_1)/sqrt(kappa_x/n_X*(1-rho_1%rho_1)%(1-rho_1%rho_1))> threshold_X);
      //save the maximum of square of trho_1 and trho_2
      trho_sq = trho_1%trho_1;
      
      // Rcout << trho_sq(0,1) << "\n";
      // Rcout << trho_sq(1,4) << "\n";
      // Rcout << trho_sq(2,3) << "\n";
      // Rcout << trho_sq(3,4) << "\n";
      //computing test statistics
      // T_stat= (cor_x-cor_y)/sqrt((kappa_x/n_x + kappa_y/n_y)*(1-pmax(tcor_x^2,tcor_y^2))^2)
      // T_stat((b_size*i):(b_size*i + min(b_size, p_X-b_size*i)),(b_size*j):(b_size*j + min(b_size, p_X-b_size*j)))= (rho_1 - rho_2)/sqrt((kappa_x/n_X+kappa_y/n_Y)*(1-trho%trho)%(1-trho%trho));
      for(int k = 0; k< min(b_size, p_X-b_size*i); k++){
        for(int l = 0; l< min(b_size, p_X-b_size*j); l++){
          T_stat(b_size*i + k, b_size*j + l) = abs(rho_1(k,l))/(sqrt(kappa_x/n_X)*(1-trho_sq(k,l)));
          // T_stat(b_size*j + l, b_size*i + k) = T_stat(b_size*i+k, b_size*j + l);
        }
      }
    }
  }
  
  //T_stat is a symmetric matrix, diagonal is 0
  for(int k = 0; k < p_X; k++){
    T_stat(k,k) = 0;
    for(int l = k+1; l < p_X; l++){
      T_stat(l,k) = T_stat(k,l);
    }
  }
  
  return T_stat;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# set.seed(2020)
# X = matrix(rnorm(50000),100,500)
# Y = matrix(rnorm(50000),100,500)
# res_cpp = cor_CL(X,Y,verbose=TRUE)
# res_r = cor_CL_R(X,Y,verbose=TRUE)
# 
# round(res_r$stat[1:5,101:105],4)
# round(res_cpp$stat[1:5,101:105],4)
# 
# summary(res_r$stat[upper.tri(res_r$stat)]-res_cpp$stat[upper.tri(res_cpp$stat)])
# sum(res_r$stat[upper.tri(res_r$stat)]!=res_cpp$stat[upper.tri(res_cpp$stat)])
# 
# res_r$threshold
# res_cpp$threshold
# 
# sum(res_r$decision[upper.tri(res_r$decision)],na.rm=TRUE)
# sum(res_cpp$decision[upper.tri(res_cpp$decision)],na.rm=TRUE)
*/
