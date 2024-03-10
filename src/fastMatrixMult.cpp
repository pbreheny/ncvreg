#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd fastMatrixMultiply(Eigen::MatrixXd A, Eigen::MatrixXd B) {
  Eigen::MatrixXd result = A * B;
  return result;
}
