
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export(VAR_sim)]]
arma::mat VAR_sim(int n, arma::vec mu, arma::mat Sigma, arma::field<arma::mat> coeffs, std::string error_dist, arma::mat P1 , arma::mat Q1, int df = 1)
{
  int d = coeffs(0).n_rows;
  arma::mat E(n ,d); arma::mat X(n,d);
  if (error_dist == "normal"){
    E = mvnrnd( mu, Sigma, n ).t();
  }
  if (error_dist == "t"){
    E = mvnrnd( zeros(d), Sigma, n ).t();
    arma::vec u = chi2rnd(df, n) ;
    u /=  df;
    E.each_col() /= sqrt(u);
    E.each_row() += mu.t();
  }
  if (error_dist == "garch"){ //BEKK(1,1)
    arma::mat Z = mvnrnd( zeros(d), eye(d,d), n ).t();
    arma::mat Sigma_t = Sigma; arma::mat Sigma_t1 = Sigma;
    E.row(0) = Z.row(0) * chol(Sigma_t).t();
    for(int ii = 1; ii < n; ii++ ) {
      Sigma_t = Sigma + P1 * Z.row(ii-1).t() * Z.row(ii-1)  * P1.t() + Q1 * Sigma_t1  * Q1.t();
      E.row(ii) = Z.row(ii) * chol(Sigma_t).t();
      Sigma_t1 = Sigma_t;
    }
  }

  int p = size(coeffs)(0);

  X = E;
  for (int row = p; row < n; row ++) {
// X.row(row) =  E.row(row);
      for(int ii = 1; ii < p+1; ii++ ) {
        X.row(row) = X.row(row) +   X.row(row-ii) * coeffs(ii-1).t() ;
       }
     }

  return X;

}

