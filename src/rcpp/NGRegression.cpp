#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

class NGRegression {
public:
  arma::mat X;
  arma::vec y;
  arma::mat Q;
  double rho;
  // priors
  arma::mat M_0inv;
  arma::vec m_0;
  double a_0;
  double b_0;
  double a_n;
  double b_n;
  double logLBayes;
  NGRegression(arma::mat X_,
                   arma::vec y_,
                   arma::mat Q_):
    X(X_), y(y_), Q(Q_) {
    N = X.n_rows;
    p = X.n_cols;
    m_nstar = zeros(p);
    M_nstar_inv = zeros(p, p);
    vec_p = zeros(p);
    vec_N = zeros(N);
  }

  void SetPriors(arma::mat M_0inv_,
                 arma::vec m_0_,
                 double a_0_,
                 double b_0_) {
    M_0inv = M_0inv_;
    m_0 = m_0_;
    a_0 = a_0_;
    b_0 = b_0_;
    m_nstar = m_0 + X.t() * Q * y;
    M_nstar_inv = M_0inv + X.t() * Q * X;
    M_nstar_inv = chol(M_nstar_inv, "upper");
    M_nstar_chol.eye(p, p);
    M_nstar_chol = solve(trimatu(M_nstar_inv), M_nstar_chol);
    beta_mean = trimatl(M_nstar_chol.t()) * m_nstar;
    b_n = dot(beta_mean, beta_mean);
    beta_mean = trimatu(M_nstar_chol) * beta_mean;
    a_n = a_0 + (double) N / 2.0;
    M_0inv_chol = chol(M_0inv, "upper");
    vec_p = solve(trimatu(M_0inv_chol), m_0);
    b_n = dot(vec_p, vec_p) - b_n;
    b_n += dot(y, Q * y);
    b_n = b_0 + b_n / 2.0;
    // pre-compute quantities for WAIC
    arma::vec mu = X * beta_mean;
    double sigma2_mean = b_n / (a_n - 1);
    vec_N = y - mu;
    // log p(y|theta_bayes)
    logLBayes = dot(vec_N, Q * vec_N) / (-2 * sigma2_mean);
    logLBayes += log_det_sympd(Q) / 2;
    logLBayes -= (double) N / 2 * log(2 * datum::pi);
  }

  List ExactSample(int n_samples) {
    double sigma2;
    arma::vec beta = zeros(p);
    arma::mat beta_sim = zeros(n_samples, p);
    arma::vec sigma2_sim = zeros(n_samples);
    arma::mat YFit_sim = zeros(n_samples, N);
    arma::vec mu = zeros(N);
    // log-likelihood for complete data
    arma::vec ll = zeros(n_samples);
    // pointwise-predicitive densities
    arma::mat log_ppd = zeros(n_samples, N);
    for (int i = 0; i < n_samples; i++) {
      sigma2 = 1.0 / randg(distr_param(a_n, 1.0 / b_n));
      beta = randn(p, distr_param(0, 1));
      beta = trimatu(M_nstar_chol) * beta;
      beta = beta * sqrt(sigma2);
      beta += beta_mean;
      beta_sim.row(i) = beta.t();
      sigma2_sim(i) = sigma2;
      mu = X * beta;
      for (int j = 0; j < N; j++) {
        YFit_sim(i, j) = randn(distr_param(mu[j], sqrt(sigma2)));
      }
      vec_N.fill(pow(sigma2, 0.5));
      vec_N = log_normpdf(y, mu, vec_N);
      ll(i) = sum(vec_N);
      log_ppd.row(i) = vec_N.t();
    }
    // compute DIC
    double p_DIC = 2 * (logLBayes - mean(ll));
    double DIC = -2 * (logLBayes - p_DIC);
    // compute WAIC
    // column-wise mean
    arma::mat temp = mean(exp(log_ppd));
    double lppd = accu(log(temp));
    temp = var(log_ppd);
    double p_WAIC2 = accu(temp);
    double WAIC = -2 *(lppd - p_WAIC2);
    // predictive loss from Gelfand and Ghosh (1998)
    for (int i = 0; i < N; i++) {
      vec_N(i) = y(i) - mean(YFit_sim.col(i));
    }
    double G = dot(vec_N, vec_N);
    double P = accu(var(YFit_sim));
    List out = List::create(_["beta"] = beta_sim,
                            _["sigma2"] = sigma2_sim,
                            _["YFit"] = YFit_sim,
                            _["DIC"] = DIC,
                            _["p_DIC"] = p_DIC,
                            _["WAIC"] = WAIC,
                            _["p_WAIC2"] = p_WAIC2,
                            _["pred_G"] = G,
                            _["pred_P"] = P,
                            _["pred_D"] = G + P);
    return out;
  }

private:
  int N;
  int p;
  arma::mat M_nstar_inv;
  arma::vec m_nstar;
  arma::mat M_nstar_chol;
  arma::mat beta_mean;
  arma::mat M_0inv_chol;
  arma::vec vec_p;
  arma::vec vec_N;
};

RCPP_MODULE(NGRegression_module) {

  class_<NGRegression>("NGRegression")
  .constructor<arma::mat, arma::vec, arma::mat>()

  .field("X", &NGRegression::X)
  .field("y", &NGRegression::y)
  .field("Q", &NGRegression::Q)
  .field("rho", &NGRegression::rho)
  .field("M_0inv", &NGRegression::M_0inv)
  .field("m_0", &NGRegression::m_0)
  .field("a_0", &NGRegression::a_0)
  .field("b_0", &NGRegression::b_0)
  .field("a_n", &NGRegression::a_n)
  .field("b_n", &NGRegression::b_n)

  .method("ExactSample", &NGRegression::ExactSample)
  .method("SetPriors", &NGRegression::SetPriors)
  ;
}
