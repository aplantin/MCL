#include "service.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' mclC
//'
//' Runs my multilevel compositional lasso
//'
//' @param Z OTU matrix (group level)
//' @param W OTU matrix (within-group level)
//' @param y Vector of outcomes
//' @param facZ Scaling factors for Z
//' @param facW Scaling factors for W
//' @param groups Vector of group membership
//' @param groupIdx List of feature indices in each group
//' @param mu Lagrange scaling parameter
//' @param lam1 Lambda for group-level lasso penalty
//' @param lam2 Lambda for within-group lasso penalty
//' @param thresh Parameter defining convergence
//' @param maxit Maximum number of iterations allowed
//' @export
// [[Rcpp::export]]
Rcpp::List mclC(arma::mat Z, arma::mat W, arma::vec y,
                arma::vec facZ, arma::vec facW,
                arma::vec groups, Rcpp::List groupIdx,
                double mu, double lam1, double lam2,
                double thresh, double maxit) {
  //int n = Z.n_rows;
  int q = Z.n_cols;
  int p = W.n_cols;
  int J = groupIdx.size();

  arma::vec beta = arma::zeros<arma::vec>(q);
  arma::vec gamma = arma::zeros<arma::vec>(p);
  arma::vec xi = arma::zeros<arma::vec>(J + 1);

  int count = 0;
  double diff = thresh + 1;
  double mu1 = 1 / (1 + mu);

  while (count < maxit && diff > thresh) {
    diff = 0.0;
    for (int i = 0; i < q; i++) {
      double old = beta(i);
      double t1 = as_scalar(Z.col(i).t() * (y - W * gamma - Z * beta + Z.col(i) * beta(i))
                              - mu * ((accu(beta) - beta(i)) + xi(J))); // scale
      beta(i) = mu1 * softC(t1, lam1);
      diff += facZ(i) * absC(beta(i) - old)/p; // scale
    }
    for (int j = 0; j < p; j++) {
      int thisgrp = groups(j);
      arma::uvec grp = groupIdx[thisgrp];
      double old = gamma(j);
      double t2 = as_scalar(W.col(j).t() * (y - Z * beta - W * gamma + W.col(j) * gamma(j))
                              - mu * ((accu(gamma(grp)) - gamma(j)) + xi(thisgrp)) ); // scale
      gamma(j) = mu1 * softC(t2, lam2);
      diff += facW(j) * absC(gamma(j) - old); // scale
    }
    arma::vec xi_update = arma::zeros<arma::vec>(J + 1);
    for (int j = 0; j < J; j++) {
      arma::uvec grp = groupIdx[j];
      xi_update(j) = accu(gamma(grp)); // scale
    }
    xi_update(J) = as_scalar(facZ.t()*beta); // scale
    xi += mu * xi_update;
    count++;
  }
  arma::vec beta_scaled = beta % facZ;
  arma::vec gamma_scaled = gamma % facW;
  return Rcpp::List::create(Rcpp::Named("beta") = beta_scaled, Rcpp::Named("gamma") = gamma_scaled);
}

//' findMaxLams
//'
//' Finds maximum lambda value (smallest lambda such
//' that all coefficients are zero) for both lam1 and lam2
//'
//' @param Z OTU matrix (group level)
//' @param W OTU matrix (within-group level)
//' @param y Vector of outcomes
//' @param facZ Scaling factors for Z
//' @param facW Scaling factors for W
//' @param groups Vector indicating group membership for each feature
//' @param groupIdx List indicating feature indices in each group
//' @param mu Lagrange scaling parameter
//' @param maxit Maximum number of iterations allowed
//' @param thresh Parameter defining convergence
//' @export
// [[Rcpp::export]]
arma::vec findMaxLams(arma::mat Z, arma::mat W, arma::vec y,
                      arma::vec facZ, arma::vec facW,
                      arma::vec groups, Rcpp::List groupIdx,
                      double mu, double maxit, double thresh) {

  // first, lambda for Z/beta
  double maxlam1 = max(absC2(Z.t() * y));
  double alt2 = max(absC2(W.t() * y));
  bool all_zero = true;
  while (all_zero == true) {
    maxlam1 *= 0.95;
    Rcpp::List res = mclC(Z, W, y, facZ, facW, groups, groupIdx, mu,
                          maxlam1, alt2, thresh, maxit);
    arma::vec bet = as<arma::vec>(res[0]);
    if (any(bet != 0)) {
      all_zero = false;
    }
  }
  maxlam1 = maxlam1 / 0.95;

  // now, lambda for W/gamma
  double alt1 = max(absC2(Z.t() * y));
  double maxlam2 = max(absC2(W.t() * y));
  all_zero = true;
  while (all_zero == true) {
    maxlam2 *= 0.95;
    Rcpp::List res = mclC(Z, W, y, facZ, facW, groups, groupIdx, mu,
                          alt1, maxlam2, thresh, maxit);
    arma::vec gam = as<arma::vec>(res[1]);
    if (any(gam != 0)) {
      all_zero = false;
    }
  }
  maxlam2 = maxlam2 / 0.95;
  arma::vec out(2);
  out(0) = maxlam1;
  out(1) = maxlam2;
  return out;
}


