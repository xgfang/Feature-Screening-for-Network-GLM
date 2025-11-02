#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
mat Network_Cpp(int n) {
  IntegerVector indices = sample(50, n, true);  // Generate random indices with replacement
  vec B = rbinom(n*n,1,0.6);
  
  mat A = arma::mat(B.begin(),n,n,false,true);
  A.diag().zeros();
  
  mat idx_mat(n,n);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      idx_mat(i,j) = indices(i)==indices(j);
    }
  }
  
  A = idx_mat % A;
  
  mat W = A;
  W.each_col() /= sum(A,1);
  
  W.elem(find_nan(W)).zeros();
  
  return W;
}

// [[Rcpp::export]]
double modeRcpp(NumericVector x) {
  std::map<double, int> counts;
  
  // Count the occurrences of each unique value
  for (int i = 0; i < x.size(); ++i) {
    counts[x[i]]++;
  }
  
  // Find the value with the highest count (mode)
  double mode = x[0];  // Initialize mode with the first value
  int maxCount = 0;
  
  for (const auto& pair : counts) {
    if (pair.second > maxCount) {
      maxCount = pair.second;
      mode = pair.first;
    }
  }
  
  return mode;
}

// [[Rcpp::export]]
vec mode(mat mat1) {
  int nCols = mat1.n_cols;
  vec modes(nCols);
  
  for (int i = 0; i < nCols; ++i) {
    NumericVector col_i = wrap(mat1.col(i));
    double mode = modeRcpp(col_i);
    modes(i) = mode;
  }
  return modes;
}

// [[Rcpp::export]]
mat mvrnormArma(int n, vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

// [[Rcpp::export]]
List dt_gen_rcpp(int n, int p, int m, int m0, mat W, double rho, mat beta, std::string family){
  
  NumericVector fill1 = rnorm(n*p,0,1);
  mat X1 = mat(fill1.begin(), n, p, false, true);
  
  mat Sigma2 = mat(p,p,fill::ones)*0.5;
  Sigma2.diag().ones();
  vec mu(p,fill::zeros);
  mat X2 = mvrnormArma(n,mu,Sigma2);
  
  mat sigma3 = mat(p,p,fill::ones)*0.3;
  for(int i = 0; i<4; i++){
    for(int j = 0; j<4; j++){
      sigma3(i,j) = 0.15;
    }
  }
  sigma3.diag().ones();
  mat X3 = mvrnormArma(n, mu,sigma3);
  
  vec Y_t1 = vec(n,fill::zeros), Y_t2 = vec(n,fill::zeros), Y_t3 = vec(n,fill::zeros);
  vec theta1(n), theta2(n), theta3(n);
  mat Y1 = mat(m,n,fill::zeros), Y2 = mat(m,n,fill::zeros), Y3 = mat(m,n,fill::zeros);
  
  vec YY1 = X1*beta.col(0);
  vec YY2 = X2*beta.col(1);
  vec YY3 = X3*beta.col(2);
  
  List res,X,Y;
  X["X1"] = X1;
  X["X2"] = X2;
  X["X3"] = X3;
  
  if(family == "gaussian")
  {
    for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
        theta1(j) = YY1(j)+rho*sum(Y_t1%W.row(j).t());
        Y_t1(j) = as<double>(rnorm(1,theta1(j),1));
        theta2(j) = YY2(j)+rho*sum(Y_t2%W.row(j).t());
        Y_t2(j) = as<double>(rnorm(1,theta2(j),1));
        theta3(j) = YY3(j)+rho*sum(Y_t3%W.row(j).t());
        Y_t3(j) = as<double>(rnorm(1,theta3(j),1));
        if(j==(n-1)){
          Y1.row(i) = Y_t1.t();
          Y2.row(i) = Y_t2.t();
          Y3.row(i) = Y_t3.t();
        }
      }
    }
    Y1 = sum(Y1.rows(m0,m-1),0)/(m-m0);
    Y2 = sum(Y2.rows(m0,m-1),0)/(m-m0);
    Y3 = sum(Y3.rows(m0,m-1),0)/(m-m0);
    
    Y["Y1"] = Y1.t();
    Y["Y2"] = Y2.t();
    Y["Y3"] = Y3.t();
  } else if(family == "binomial"){
    for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
        theta1(j) = YY1(j)+rho*sum(Y_t1%W.row(j).t());
        double invlog1 = (exp(theta1(j))/(1+exp(theta1(j))));
        Y_t1(j) = as<int>(rbinom(1,1,invlog1));
        theta2(j) = YY2(j)+rho*sum(Y_t2%W.row(j).t());
        double invlog2 = (exp(theta2(j))/(1+exp(theta2(j))));
        Y_t2(j) = as<int>(rbinom(1,1,invlog2));
        theta3(j) = YY3(j)+rho*sum(Y_t3%W.row(j).t());
        double invlog3 = (exp(theta3(j))/(1+exp(theta3(j))));
        Y_t3(j) = as<int>(rbinom(1,1,invlog3));
        if(j==(n-1)){
          Y1.row(i) = Y_t1.t();
          Y2.row(i) = Y_t2.t();
          Y3.row(i) = Y_t3.t();
        }
      }
    }
    Y1 = mode(Y1.rows(m0,m-1));
    Y2 = mode(Y2.rows(m0,m-1));
    Y3 = mode(Y3.rows(m0,m-1));
    Y["Y1"] = Y1;
    Y["Y2"] = Y2;
    Y["Y3"] = Y3;
  } else if(family == "poisson"){
    for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
        theta1(j) = YY1(j)+rho*sum(Y_t1%W.row(j).t());
        double exp1 = exp(theta1(j));
        Y_t1(j) = as<int>(rpois(1,exp1));
        theta2(j) = YY2(j)+rho*sum(Y_t2%W.row(j).t());
        double exp2 = exp(theta2(j));
        Y_t2(j) = as<int>(rpois(1,exp2));
        theta3(j) = YY3(j)+rho*sum(Y_t3%W.row(j).t());
        double exp3 = exp(theta3(j));
        Y_t3(j) = as<int>(rpois(1,exp3));
        if(j==(n-1)){
          Y1.row(i) = Y_t1.t();
          Y2.row(i) = Y_t2.t();
          Y3.row(i) = Y_t3.t();
        }
      }
    }
    Y1 = mode(Y1.rows(m0,m-1));
    Y2 = mode(Y2.rows(m0,m-1));
    Y3 = mode(Y3.rows(m0,m-1));
    Y["Y1"] = Y1;
    Y["Y2"] = Y2;
    Y["Y3"] = Y3;
  }
  res["X"] = X;
  res["Y"] = Y;
  
  return res;
  
}