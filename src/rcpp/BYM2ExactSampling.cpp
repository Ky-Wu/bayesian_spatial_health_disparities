#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

class BYM2ExactSampler {
public:
  arma::mat X;
  arma::vec y;
  arma::mat Q;
  double rho;
  // priors
  List priors;
  arma::mat M_0inv;
  arma::vec m_0;
  arma::mat M_nstar_inv;
  arma::vec m_nstar;
  double a_0;
  double b_0;
  double a_n;
  double b_n;
  arma::vec gammastar_mean;
  double llBayes;
  double logCollapsedBayes;
  BYM2ExactSampler(arma::mat X_,
                   arma::vec y_,
                   arma::mat Q_,
                   double rho_):
    X(X_), y(y_), Q(Q_), rho(rho_) {
    N = X.n_rows;
    p = X.n_cols;
    XtX = X.t() * X;
    Xty = X.t() * y;
    m_nstar = zeros(N + p);
    M_nstar_inv = zeros(N + p, N + p);
    gammastar_mean = zeros(N + p);
    gammastar = zeros(N + p);
    vec_p = zeros(p);
    vec_N = zeros(N);
    eig_sym(Lambda, P, Q);
  }

  void SetPriors(arma::mat M_0inv_,
                 arma::vec m_0_,
                 double a_0_,
                 double b_0_) {
    M_0inv = M_0inv_;
    m_0 = m_0_;
    a_0 = a_0_;
    b_0 = b_0_;
    m_nstar.subvec(0, p - 1) = Xty / (1.0 - rho) + m_0;
    m_nstar.subvec(p, N + p - 1) = y / (1.0 - rho);
    // M_nstar_inv = X_star' V_yinv X_star
    M_nstar_inv.submat(0, 0, p - 1, p - 1) = XtX / (1.0 - rho) + M_0inv;
    M_nstar_inv.submat(0, p, p - 1, N + p - 1) = X.t() / (1.0 - rho);
    M_nstar_inv.submat(p, 0, N + p - 1, p - 1) = X / (1.0 - rho);
    M_nstar_inv.submat(p, p, N + p - 1, N + p - 1) = Q / rho;
    for (int i = 0; i < N; i++) {
      M_nstar_inv(p + i, p + i) += 1.0 / (1.0 - rho);
    }
    M_nstar_inv = chol(M_nstar_inv, "upper");
    M_nstar_chol.eye(N + p, N + p);
    M_nstar_chol = solve(trimatu(M_nstar_inv), M_nstar_chol);
    gammastar_mean = trimatl(M_nstar_chol.t()) * m_nstar;
    a_n = a_0 + (double) N / 2.0;
    b_n = dot(gammastar_mean, gammastar_mean);
    gammastar_mean = trimatu(M_nstar_chol) * gammastar_mean;
    M_0inv_chol = chol(M_0inv, "upper");
    vec_p = solve(trimatu(M_0inv_chol), m_0);
    b_n = dot(vec_p, vec_p) - b_n;
    b_n += dot(y, y) / (1 - rho);
    b_n = b_0 + b_n / 2.0;
    // pre-compute quantities for WAIC
    arma::vec mu = X * gammastar_mean.subvec(0, p - 1) +
      gammastar_mean.subvec(p, N + p - 1);
    double sigma2_mean = b_n / (a_n - 1);
    vec_N.fill(pow(sigma2_mean * (1 - rho), 0.5));
    vec_N = log_normpdf(y, mu, vec_N);
    // log p(y|theta_bayes)
    llBayes = sum(vec_N);
    vec_p = gammastar_mean.subvec(0, p - 1);
    logCollapsedBayes = collapsedLikelihood(vec_p, sigma2_mean);
  }

  List ExactSample(int n_samples) {
    double sigma2;
    arma::vec gamma_star = zeros(N + p);
    arma::mat beta_sim = zeros(n_samples, p);
    arma::mat gamma_sim = zeros(n_samples, N);
    arma::vec sigma2_sim = zeros(n_samples);
    arma::mat YFit_sim = zeros(n_samples, N);
    arma::vec mu = zeros(N);
    // log-likelihood for complete data
    arma::vec ll = zeros(n_samples);
    // pointwise-predicitive densities
    arma::mat log_ppd = zeros(n_samples, N);
    // non-spatial likelihood
    arma::vec collapsedLL = zeros(n_samples);
    for (int i = 0; i < n_samples; i++) {
      sigma2 = 1.0 / randg(distr_param(a_n, 1.0 / b_n));
      gammastar = randn(N + p, distr_param(0, 1));
      gammastar = trimatu(M_nstar_chol) * gammastar;
      gammastar = gammastar * sqrt(sigma2);
      gammastar += gammastar_mean;
      beta_sim.row(i) = gammastar.subvec(0, p - 1).t();
      gamma_sim.row(i) = gammastar.subvec(p, N + p - 1).t();
      sigma2_sim(i) = sigma2;
      mu = X * gammastar.subvec(0, p - 1) + gammastar.subvec(p, N + p - 1);
      for (int j = 0; j < N; j++) {
        YFit_sim(i, j) = randn(distr_param(mu[j], sqrt(sigma2 * (1 - rho))));
      }
      vec_N.fill(pow(sigma2 * (1 - rho), 0.5));
      vec_N = log_normpdf(y, mu, vec_N);
      ll(i) = sum(vec_N);
      log_ppd.row(i) = vec_N.t();
      vec_p = gammastar.subvec(0, p - 1);
      collapsedLL(i) = collapsedLikelihood(vec_p, sigma2);
    }
    // compute DIC
    double p_DIC = 2 * (llBayes - mean(ll));
    double p_DIC2 = 2 * var(ll);
    double DIC = -2 * (llBayes - p_DIC);
    // compute WAIC
    // column-wise mean
    arma::mat temp = mean(exp(log_ppd));
    double lppd = accu(log(temp));
    double p_WAIC = 2 * accu(log(temp) - mean(log_ppd));
    temp = var(log_ppd);
    double p_WAIC2 = accu(temp);
    double WAIC = -2 *(lppd - p_WAIC2);
    // predictive loss from Gelfand and Ghosh (1998)
    for (int i = 0; i < N; i++) {
      vec_N(i) = y(i) - mean(YFit_sim.col(i));
    }
    double G = dot(vec_N, vec_N);
    double P = accu(var(YFit_sim));
    // collapsed likelihood DIC
    double p_DIC_collapsed = 2 * (logCollapsedBayes - mean(collapsedLL));
    double DIC_collapsed = -2 * (logCollapsedBayes - p_DIC_collapsed);
    List out = List::create(_["beta"] = beta_sim,
                            _["gamma"] = gamma_sim,
                            _["sigma2"] = sigma2_sim,
                            _["YFit"] = YFit_sim,
                            _["DIC"] = DIC,
                            _["p_DIC"] = p_DIC,
                            _["p_DIC2"] = p_DIC2,
                            _["WAIC"] = WAIC,
                            _["p_WAIC"] = p_WAIC,
                            _["p_WAIC2"] = p_WAIC2,
                            _["p_DIC_collapsed"] = p_DIC_collapsed,
                            _["DIC_collapsed"] = DIC_collapsed,
                            _["pred_G"] = G,
                            _["pred_P"] = P,
                            _["pred_D"] = G + P);
    return out;
  }

  double collapsedLikelihood(arma::vec & beta,
                              double sigma2) {
    double out = (double) N / -2.0 * log(2 * datum::pi * sigma2);
    out -= (sum(rho / Lambda) + N * (1 - rho)) / 2;
    vec_N = y - X * beta;
    vec_N = P.t() * vec_N;
    vec_N = pow(rho / Lambda + (1 - rho), -0.5) % vec_N;
    out -= dot(vec_N, vec_N) / (2 * sigma2);
    return out;
  }

private:
  int N;
  int p;
  arma::mat XtX;
  arma::vec Xty;
  arma::mat M_nstar_chol;
  arma::mat M_0inv_chol;
  arma::vec gammastar;
  arma::mat P;
  arma::vec Lambda;
  arma::vec vec_p;
  arma::vec vec_N;
};

RCPP_MODULE(BYM2ExactSampler_module) {

  class_<BYM2ExactSampler>("BYM2ExactSampler")
  .constructor<arma::mat, arma::vec, arma::mat, double>()

  .field("X", &BYM2ExactSampler::X)
  .field("y", &BYM2ExactSampler::y)
  .field("Q", &BYM2ExactSampler::Q)
  .field("rho", &BYM2ExactSampler::rho)
  .field("M_0inv", &BYM2ExactSampler::M_0inv)
  .field("m_0", &BYM2ExactSampler::m_0)
  .field("a_0", &BYM2ExactSampler::a_0)
  .field("b_0", &BYM2ExactSampler::b_0)
  .field("a_n", &BYM2ExactSampler::a_n)
  .field("b_n", &BYM2ExactSampler::b_n)
  .field("gammastar_mean", &BYM2ExactSampler::gammastar_mean)
  .field("M_nstar_inv", &BYM2ExactSampler::M_nstar_inv)
  .field("m_nstar", &BYM2ExactSampler::m_nstar)

  .method("ExactSample", &BYM2ExactSampler::ExactSample)
  .method("SetPriors", &BYM2ExactSampler::SetPriors)
  ;
}
