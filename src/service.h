#ifndef SERVICE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define SERVICE_H

#include <RcppArmadillo.h>

double signC(double x);
double absC(double x);
arma::vec absC2(arma::vec x);
double softC(double x, double lam);
arma::vec softC2(arma::vec x, double lam);
double sqrtC(double x);

#endif
