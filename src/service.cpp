#include <RcppArmadillo.h>
#include "service.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//' signC
//'
//' Finds the sign of a double
//'
//' @param x Number to find the sign of
//' @export
// [[Rcpp::export]]
double signC(double x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

//' absC
//'
//' Finds the absolute value of a double
//'
//' @param x Number to take abs value of
//' @export
// [[Rcpp::export]]
double absC(double x) {
  return std::abs(x);
}

//' absC2
//'
//' Finds the absolute value of each element in a numeric vector
//'
//' @param x Vector to take abs value of
//' @export
// [[Rcpp::export]]
arma::vec absC2(arma::vec x) {
  arma::vec out = arma::vec(x.size());
  for (int i = 0; i < x.size(); i++) {
    out(i) = absC(x(i));
  }
  return out;
}

//' softC
//'
//' Thresholds x by lambda
//'
//' @param x Number to soft threshold
//' @param lam Parameter by which to soft threshold
//' @export
// [[Rcpp::export]]
double softC(double x, double lam) {
  double out = signC(x) * std::max(absC(x) - lam, 0.0);
  return out;
}

//' softC2
//'
//' Soft thresholds each element of a vector x by lambda
//'
//' @param x Vector whose elements should be soft thresholded
//' @param lam Parameter by which to soft threshold
//' @export
// [[Rcpp::export]]
arma::vec softC2(arma::vec x, double lam){
  arma::vec xal = abs(x) - lam;
  arma::uvec sparse = arma::uvec(xal > 0);
  return xal % sparse % sign(x);
}

//' sqrtC
//'
//' Finds the square root of a double
//'
//' @param x Number to take sqrt of
//' @export
// [[Rcpp::export]]
double sqrtC(double x) {
  return std::sqrt(x);
}

