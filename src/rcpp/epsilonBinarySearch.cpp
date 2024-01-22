#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Binary search for maximum epsilon such that Bayesian FDR can be controlled
//' under threshold delta with T positive results
//'
//' @param d matrix of non-negative values to evaluate against epsilon
//' @param eta numeric in [0, 1], maximum tolerable Bayesian FDR
//' @param n_positive integer, number of positive labels
//' @param lower numeric, initial lower bound
//' @param upper numeric, initial upper bound
//' @param toler numeric, error tolerance for convergence criteria
//'
//' @return maximum epsilon such that Bayesian FDR controlled under delta with T
//'  positive results
// [[Rcpp::export]]
double epsilonBinarySearch(arma::mat & diffs,
                           double eta,
                           int n_positive,
                           double lower,
                           double upper,
                           double toler = 0.00001) {
  int n_row = diffs.n_rows;
  int n_col = diffs.n_cols;
  double initial_bound = upper - lower;
  double eps;
  double bayesian_fdr;
  double n_iter = 0;
  Rcpp::NumericVector vij(n_col);
  while ((initial_bound / pow(2.0, n_iter)) > toler) {
    eps = (upper + lower) / 2;
    for (int j = 0; j < n_col; j++) {
      double counter = 0.0;
      for (int i = 0; i < n_row; i++) {
        if (diffs(i, j) > eps) {
          counter += 1.0;
        }
      }
      vij[j] = counter / n_row;
    }
    vij.sort(true);
    bayesian_fdr = 0.0;
    int current_positives = 0;
    for (int i = 0; i < n_col; i++) {
      if (vij[i] >= vij[n_positive - 1]) {
        bayesian_fdr += 1.0 - vij[i];
        current_positives++;
      } else {
        break;
      }
    }
    bayesian_fdr /= current_positives;
    if (bayesian_fdr < eta) {
      lower = eps;
    } else {
      upper = eps;
    }
    n_iter++;
  }
  return eps;
}
