#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

class BYM2FlatBetaMCMC {
public:
  arma::mat X;
  arma::vec y;
  arma::mat Q;
  arma::vec beta;
  arma::vec gamma;
  double rho;
  double sigma2;
  // priors
  double a_sigma;
  double b_sigma;
  double astar_sigma;
  double bstar_sigma;
  double a_rho;
  double b_rho;
  double astar_rho;
  double bstar_rho;
  double lower_rho;
  double upper_rho;
  arma::vec vec_N;
  arma::vec Lambda;
  arma::mat P;
  arma::mat C;
  arma::mat U;
  arma::vec D;
  arma::mat matrix_pp;
  arma::vec vec_p;
  arma::vec vec2_p;

  BYM2FlatBetaMCMC(arma::mat X_,
                   arma::vec y_,
                   arma::mat Q_):
    X(X_), y(y_), Q(Q_) {
    N = X.n_rows;
    p = X.n_cols;
    XtX = X.t() * X;
    XtX_R = chol(XtX, "upper");
    beta_ols = X.t() * y;
    beta_ols = solve(trimatl(XtX_R.t()), beta_ols);
    beta_ols = solve(trimatu(XtX_R), beta_ols);
    // Q = P Lambda P'
    eig_sym(Lambda, P, Q);
    // store intermediate Q^{-1/2} temporarily in U
    U = P * diagmat(pow(Lambda, -0.5)) * P.t();
    // matrix_NN = I - H
    matrix_pN = solve(trimatl(XtX_R.t()), X.t());
    matrix_NN = -1 * matrix_pN.t() * matrix_pN;
    matrix_NN.diag() += 1;
    // simultaneous diagonalization of Q and I - H
    // Q^0.5 * (I - H) * Q^0.5
    matrix_NN = U * matrix_NN * U;
    eig_sym(D, C, matrix_NN);
    U = U * C;
    UtX = U.t() * X;
    // e := y - X * beta_ols
    // r := y - X * beta - gamma
    e = y - X * beta_ols;
    Ute = U.t() * e;
    vec_p = zeros(p);
    vec2_p = zeros(p);
    vec_N = zeros(N);
    matrix_pp = zeros(p, p);
    matrix_pN = zeros(p, N);
    matrix_NN = zeros(N, N);
  }

  // priors moved out of constructor since >6 args in constructor not supported
  // in Rcpp modules
  void setPriors(double a_sigma_,
                 double b_sigma_,
                 double a_rho_ = 0.0,
                 double b_rho_ = 0.0,
                 double lower_rho_ = 0.0,
                 double upper_rho_ = 1.0) {
    a_sigma = a_sigma_;
    b_sigma = b_sigma_;
    a_rho = a_rho_;
    b_rho = b_rho_;
    lower_rho = lower_rho_;
    upper_rho = upper_rho_;
    astar_sigma = a_sigma + (double) N / 2.0;
    astar_rho = a_sigma + (double) N / 2.0;
  }

  void initOLS() {
    // initialize beta to beta_ols
    beta = beta_ols;
    // initialize gamma to 0
    gamma = zeros(N);
    // rtr = ||y - X * beta - gamma||^2
    rtr = dot(e, e);
    // initialize rho
    rho = 0.5;
    // initialize sigma2
    sigma2 = rtr / (N - 1);
  }

  void initZeros() {
    beta = zeros(p);
    gamma = zeros(N);
    rtr = dot(y, y);
    sigma2 = rtr / (N - 1);
    updateRho(false);
  }

  void initRandom() {
    beta = randn(p, distr_param(0, 1));
    gamma = randn(N, distr_param(0.0, 1.0 / N));
    vec_N = y - X * beta - gamma;
    rtr = dot(vec_N, vec_N);
    sigma2 = rtr / (N - 1);
    updateRho(false);
  }

  List MCMCSample(int n_samples) {
    arma::mat beta_sim = zeros(n_samples, p);
    arma::mat gamma_sim = zeros(n_samples, N);
    arma::vec sigma2_sim = zeros(n_samples);
    arma::vec rho_sim = zeros(n_samples);
    arma::mat YFit_sim = zeros(n_samples, N);
    arma::vec log_postd = zeros(n_samples);
    arma::vec mu = zeros(N);
    // log-likelihood for complete data
    arma::vec ll = zeros(n_samples);
    // pointwise-predicitive densities
    arma::mat log_ppd = zeros(n_samples, N);
    for (int i = 0; i < n_samples; i++) {
      updateAllParams();
      beta_sim.row(i) = beta.t();
      sigma2_sim(i) = sigma2;
      gamma_sim.row(i) = gamma.t();
      rho_sim(i) = rho;
      mu = X * beta + gamma;
      for (int j = 0; j < N; j++) {
        YFit_sim(i, j) = randn(distr_param(mu[j], sqrt(sigma2 * (1 - rho))));
      }
      log_postd(i) = logPostd(false);
      vec_N.fill(pow(sigma2 * (1 - rho), 0.5));
      vec_N = log_normpdf(y, mu, vec_N);
      ll(i) = sum(vec_N);
      log_ppd.row(i) = vec_N.t();
    }
    // estimated posterior mean
    mu = (X * arma::mean(beta_sim, 0).t() + arma::mean(gamma_sim, 0).t()).as_col();
    vec_N = zeros(N);
    vec_N.fill(pow(mean(sigma2_sim) * (1 - mean(rho_sim)), 0.5));
    vec_N = log_normpdf(y, mu, vec_N);
    double logLBayes = sum(vec_N);
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
    List out = List::create(_["beta"] = beta_sim,
                            _["gamma"] = gamma_sim,
                            _["sigma2"] = sigma2_sim,
                            _["rho"] = rho_sim,
                            _["YFit"] = YFit_sim,
                            _["log_postd"] = log_postd,
                            _["DIC"] = DIC,
                            _["p_DIC"] = p_DIC,
                            _["WAIC"] = WAIC,
                            _["p_WAIC2"] = p_WAIC2);
    return out;
  }

  void burnMCMCSample(int n_iter) {
    for (int i = 0; i < n_iter; i++) {
      updateAllParams();
    }
  }

  void updateAllParams() {
    updateBeta();
    updateGamma();
    updateSigma2(true);
    updateRho(false);
  }

  double logPostd(bool update_rtr = true) {
    if (update_rtr) {
      vec_N = y - X * beta - gamma;
      rtr = dot(vec_N, vec_N);
    }
    // likelihood
    double out = (double) -N / 2.0 * log(2.0 * datum::pi * sigma2 * (1 - rho));
    out -= rtr / (2 * sigma2 * (1 - rho));
    // CAR prior on phi
    out -= (double) N / 2.0 * log(2 * datum::pi * sigma2 * rho);
    out -= dot(gamma, gamma) / (2 * sigma2 * rho);
    // IG prior on sigma2
    out -= (a_sigma + 1) * log(sigma2) + 1 / (b_sigma * sigma2);
    // truncated IG prior on 1 - rho
    out -= (a_rho + 1) * log(1 - rho) + 1 / (b_rho * (1 - rho));
    return out;
  }

  DataFrame rhoChainAnalysis(int n_chains,
                     int n_burnin,
                     int n_samples) {
    int total_samps = n_chains * n_samples;
    int indx;
    double init_rho;
    arma::vec sigma2_sim = zeros(total_samps);
    arma::vec rho_sim = zeros(total_samps);
    arma::vec log_postd = zeros(total_samps);
    arma::vec chain_indx = zeros(total_samps);
    arma::vec init_rhos = zeros(total_samps);
    arma::vec rtr_values = zeros(total_samps);
    for (int i = 0; i < n_chains; i++) {
      initOLS();
      init_rho = randu(distr_param(0, 1));
      rho = init_rho;
      burnMCMCSample(n_burnin);
      for (int j = 0; j < n_samples; j++) {
        indx = i * n_samples + j;
        updateBeta();
        updateGamma();
        updateSigma2(true);
        updateRho(false);
        sigma2_sim(indx) = sigma2;
        rho_sim(indx) = rho;
        log_postd(indx) = logPostd(false);
        chain_indx(indx) = i;
        init_rhos(indx) = init_rho;
        rtr_values(indx) = rtr;
      }
    }
    DataFrame out = DataFrame::create(_["chain"] = chain_indx,
                                      _["init_rho"] = init_rhos,
                                      _["sigma2"] = sigma2_sim,
                                      _["rho"] = rho_sim,
                                      _["log_postd"] = log_postd,
                                      _["rtr"] = rtr_values);
    return out;
  }

  void updateBeta() {
    // E(beta | y, phi, sigma2, rho) = beta_ols - (X'X)^(-1)X'gamma
    // = (X'X)^(-1)X'(y - gamma)
    vec_N = y - gamma;
    // use beta temporarily
    beta = X.t() * vec_N;
    vec_p = solve(trimatl(XtX_R.t()), beta);
    // var(beta | y, phi, sigma2, rho) = sigma2 * (1 - rho) (X'X)^{-1}
    vec_p += randn(p, distr_param(0.0, pow(sigma2 * (1.0 - rho), 0.5)));
    beta = solve(trimatu(XtX_R), vec_p);
  }

  void updateGamma() {
    // B = (1 - rho) / rho * Q + (I - H) = U[(1 - rho) / rho I_n + D]U'
    // E(gamma | y, beta, sigma2, rho) =
    // B^{-1}[y - XC^{-1}(X'Xbeta + X'B^(-1)y)]
    // C = X'X + X'B^(-1)X
    // vec_N = B^(-1)y
    vec_N = U.t() * y;
    vec_N = (1 / ((1 - rho) / rho + D)) % vec_N;
    vec_N = U * vec_N;
    vec_p = X.t() * vec_N + XtX * beta;
    // matrix_pp = chol(C) = chol(X'X + X'B^{-1}X)
    matrix_pN = UtX.t() * diagmat(pow((1 - rho) / rho + D, -0.5));
    matrix_pp =  matrix_pN * matrix_pN.t();
    matrix_pp += XtX;
    matrix_pp = chol(matrix_pp, "upper");
    vec2_p = solve(trimatl(matrix_pp.t()), vec_p);
    vec_p = solve(trimatu(matrix_pp), vec2_p);
    gamma = y - X * vec_p;
    gamma = U.t() * gamma;
    gamma = (1.0 / ((1.0 - rho) / rho + D)) % gamma;
    gamma = U * gamma;
    // var(gamma | y, beta, sigma2, rho) = sigma2 * [1 / (1 - rho) * I + Q / rho]
    // sigma2 * P[1 / (1 - rho) * I + 1 / rho * Lambda]^{-1}P^T
    vec_N = Lambda / rho;
    vec_N += 1.0 / (1.0 - rho);
    vec_N = pow(vec_N, -0.5);
    vec_N = vec_N % randn(N, distr_param(0, 1));
    vec_N = P * vec_N;
    vec_N *= pow(sigma2, 0.5);
    gamma = gamma + vec_N;
  }

  void updateSigma2(bool update_rtr = true) {
    if (update_rtr) {
      vec_N = y - X * beta - gamma;
      rtr = dot(vec_N, vec_N);
    }
    bstar_sigma = b_sigma + rtr / (2.0 * (1.0 - rho));
    sigma2 = 1.0 / randg(distr_param(astar_sigma, 1.0 / bstar_sigma));
  }

  void updateRho(bool update_rtr = true) {
    if (update_rtr) {
      vec_N = y - X * beta - gamma;
      rtr = dot(vec_N, vec_N);
    }
    // p(rho) = TruncIG(1 - rho | a_rho, b_rho, [0, 1])
    // p(rho | y) = TruncIG(1 - rho | astar_rho, bstar_rho, [0, 1])
    bstar_rho = b_rho + rtr / (2.0 * sigma2);
    double rho_star = 2.0;
    while (rho_star > upper_rho || rho_star < lower_rho) {
      rho_star = 1.0 - 1.0 / randg(distr_param(astar_rho, 1.0 / bstar_rho));
    }
    rho = rho_star;
  }

private:
  int N;
  int p;
  double rtr;
  arma::vec e;
  arma::vec Ute;
  arma::mat XtX;
  arma::mat XtX_R;
  arma::vec beta_ols;
  arma::mat UtX;
  arma::mat matrix_pN;
  arma::mat matrix_NN;
};

RCPP_MODULE(BYM2FlatBetaMCMC_module) {

  class_<BYM2FlatBetaMCMC>("BYM2FlatBetaMCMC")
  .constructor<arma::mat, arma::vec, arma::mat>()

  .field("X", &BYM2FlatBetaMCMC::X)
  .field("y", &BYM2FlatBetaMCMC::y)
  .field("Q", &BYM2FlatBetaMCMC::Q)
  .field("beta", &BYM2FlatBetaMCMC::beta)
  .field("gamma", &BYM2FlatBetaMCMC::gamma)
  .field("rho", &BYM2FlatBetaMCMC::rho)
  .field("sigma2", &BYM2FlatBetaMCMC::sigma2)
  .field("astar_sigma", &BYM2FlatBetaMCMC::astar_sigma)
  .field("bstar_sigma", &BYM2FlatBetaMCMC::bstar_sigma)
  .field("astar_rho", &BYM2FlatBetaMCMC::astar_rho)
  .field("bstar_rho", &BYM2FlatBetaMCMC::bstar_rho)
  .field("vec_N", &BYM2FlatBetaMCMC::vec_N)
  .field("Lambda", &BYM2FlatBetaMCMC::Lambda)
  .field("P", &BYM2FlatBetaMCMC::P)
  .field("matrix_pp", &BYM2FlatBetaMCMC::matrix_pp)
  .field("vec_p", &BYM2FlatBetaMCMC::vec_p)
  .field("vec2_p", &BYM2FlatBetaMCMC::vec2_p)
  .field("U", &BYM2FlatBetaMCMC::U)
  .field("D", &BYM2FlatBetaMCMC::D)

  .method("setPriors", &BYM2FlatBetaMCMC::setPriors)
  .method("initOLS", &BYM2FlatBetaMCMC::initOLS)
  .method("initZeros", &BYM2FlatBetaMCMC::initZeros)
  .method("initRandom", &BYM2FlatBetaMCMC::initRandom)
  .method("MCMCSample", &BYM2FlatBetaMCMC::MCMCSample)
  .method("burnMCMCSample", &BYM2FlatBetaMCMC::burnMCMCSample)
  .method("updateBeta", &BYM2FlatBetaMCMC::updateBeta)
  .method("updateGamma", &BYM2FlatBetaMCMC::updateGamma)
  .method("updateRho", &BYM2FlatBetaMCMC::updateRho)
  .method("updateSigma2", &BYM2FlatBetaMCMC::updateSigma2)
  .method("logPostd", &BYM2FlatBetaMCMC::logPostd)
  .method("rhoChainAnalysis", &BYM2FlatBetaMCMC::rhoChainAnalysis)
  ;
}
