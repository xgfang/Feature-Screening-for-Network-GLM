#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat mvrnormArma(int n, vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

vec ArmaEigenValues(mat M) {
  return eig_sym(M);
}

arma::mat rmvnorm(
    int n,
    const arma::vec& mu,
    const arma::mat& Sigma) {
  unsigned int p = Sigma.n_cols;
  
  // First draw N x P values from a N(0,1)
  Rcpp::NumericVector draw = Rcpp::rnorm(n*p);
  
  // Instantiate an Armadillo matrix with the
  // drawn values using advanced constructor
  // to reuse allocated memory
  arma::mat Z = arma::mat(
    draw.begin(), n, p, false, true);
  // Simpler, less performant alternative
  // arma::mat Z = Rcpp::as<arma::mat>(draw);
  
  // Generate a sample from the Transformed
  // Multivariate Normal
  arma::mat Y = arma::repmat(mu, 1, n).t() +
    Z * arma::chol(Sigma);
  return Y;
}