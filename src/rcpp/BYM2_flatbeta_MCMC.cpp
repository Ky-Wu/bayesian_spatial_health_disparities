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
  double rho, sigma2;
  // priors
  double a_sigma, b_sigma;
  double astar_sigma, bstar_sigma;
  std::string rho_prior_type;
  // Truncated inverse gamma prior on 1 - rho
  double a_rho, b_rho;
  double astar_rho, bstar_rho;
  // PC prior on rho
  double lambda_rho;
  double lower_rho, upper_rho;
  arma::vec Lambda;
  arma::mat P;
  arma::mat C;
  arma::mat U;
  arma::vec D;
  arma::mat matrix_pp;
  arma::vec vec_p;
  arma::vec vec2_p;
  arma::vec vec_N;

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
    // Q^{-1/2} * (I - H) * Q^{-1/2}
    matrix_NN = U * matrix_NN * U;
    eig_sym(D, C, matrix_NN);
    U = U * C;
    Uty = U.t() * y;
    UtX = U.t() * X;
    Xty = X.t() * y;
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
                 std::string rho_prior_type_,
                 double a_rho_,
                 double b_rho_,
                 double lower_rho_,
                 double upper_rho_,
                 double lambda_rho_) {
    a_sigma = a_sigma_;
    b_sigma = b_sigma_;
    if (rho_prior_type_ == "pc") {
      lambda_rho = lambda_rho_;
      a_rho = 0.0;
      b_rho = 0.0;
    } else if (rho_prior_type_ == "truncig") {
      a_rho = a_rho_;
      b_rho = b_rho_;
      lambda_rho = -1.0;
    } else {
      throw std::runtime_error("allowed rho prior types: 'pc', 'truncig'");
    }
    rho_prior_type = rho_prior_type_;
    a_rho = a_rho_;
    b_rho = b_rho_;
    lambda_rho = lambda_rho_;
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
    double p_DIC2 = 2 * var(ll);
    double DIC = -2 * (logLBayes - p_DIC);
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
    List out = List::create(_["beta"] = beta_sim,
                            _["gamma"] = gamma_sim,
                            _["sigma2"] = sigma2_sim,
                            _["rho"] = rho_sim,
                            _["YFit"] = YFit_sim,
                            _["log_postd"] = log_postd,
                            _["DIC"] = DIC,
                            _["p_DIC"] = p_DIC,
                            _["p_DIC2"] = p_DIC2,
                            _["WAIC"] = WAIC,
                            _["p_WAIC"] = p_WAIC,
                            _["p_WAIC2"] = p_WAIC2,
                            _["pred_G"] = G,
                            _["pred_P"] = P,
                            _["pred_D"] = G + P);
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

  void updateBeta() {
    // E(beta | y, phi, sigma2, rho) = beta_ols - (X'X)^(-1)X'gamma
    // = (X'X)^(-1)X'(y - gamma)
    // use beta temporarily
    beta = Xty - X.t() * gamma;
    vec_p = solve(trimatl(XtX_R.t()), beta);
    // var(beta | y, phi, sigma2, rho) = sigma2 * (1 - rho) (X'X)^{-1}
    vec_p += randn(p, distr_param(0.0, pow(sigma2 * (1.0 - rho), 0.5)));
    beta = solve(trimatu(XtX_R), vec_p);
  }

  void updateGamma() {
    // B = (1 - rho) / rho * Q + (I - H) = U[(1 - rho) / rho I_n + D]U'
    // E(gamma | y, beta, sigma2, rho) =
    // ((1 - rho / rho) Q + I_n)^{-1}(y - X %*% \beta)
    vec_N = y - X * beta;
    vec_N = P.t() * vec_N;
    vec_N = (1.0 / ((1.0 - rho) / rho * Lambda + 1.0)) % vec_N;
    gamma = vec_N;
    // var(gamma | y, beta, sigma2, rho) = sigma2 * [1 / (1 - rho) * I + Q / rho]^{-1}
    // sigma2 * P[1 / (1 - rho) * I + 1 / rho * Lambda]^{-1}P^T
    vec_N = Lambda / rho;
    vec_N += 1.0 / (1.0 - rho);
    vec_N = pow(vec_N, -0.5);
    vec_N = vec_N % randn(N, distr_param(0, 1));
    vec_N *= pow(sigma2, 0.5);
    gamma = P * (gamma + vec_N);
  }

  arma::vec gammaMean() {
    // B = (1 - rho) / rho * Q + (I - H) = U[(1 - rho) / rho I_n + D]U'
    // E(gamma | y, beta, sigma2, rho) =
    // B^{-1}[y - XC^{-1}(X'Xbeta + X'B^(-1)y)]
    // C = X'X + X'B^(-1)X
    // vec_N = B^(-1)y
    vec_N = Uty;
    vec_N = (1 / ((1 - rho) / rho + D)) % vec_N;
    vec_p = UtX.t() * vec_N + XtX * beta;
    // matrix_pp = chol(C) = chol(X'X + X'B^{-1}X)
    matrix_pN = UtX.t() * diagmat(pow((1 - rho) / rho + D, -0.5));
    matrix_pp =  matrix_pN * matrix_pN.t();
    matrix_pp += XtX;
    matrix_pp = chol(matrix_pp, "upper");
    vec2_p = solve(trimatl(matrix_pp.t()), vec_p);
    vec_p = solve(trimatu(matrix_pp), vec2_p);
    arma:;vec gamma_mean = Uty - UtX * vec_p;
    gamma_mean = (1.0 / ((1.0 - rho) / rho + D)) % gamma_mean;
    gamma_mean = U * gamma_mean;
    return gamma_mean;
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
    bstar_rho = b_rho + rtr / (2.0 * sigma2);
    double rho_star = 2.0;
    double logMH;
    double logu;
    while (rho_star > upper_rho || rho_star < lower_rho) {
      rho_star = 1.0 - 1.0 / randg(distr_param(astar_rho, 1.0 / bstar_rho));
    }
    if (rho_prior_type == "pc") {
      // p(rho)
      logMH = pow(2.0 * KLD(rho), 0.5) - pow(2.0 * KLD(rho_star), 0.5);
      logMH = logMH * lambda_rho;
      // p(gamma | sigma2, rho)
      logMH -= (double) N / 2.0 * log(rho_star / rho);
      vec_N = Q * gamma;
      logMH -= dot(gamma, vec_N) / (2.0 * sigma2) * (1.0 / rho_star - 1.0 / rho);
      logu = log(randu(distr_param(0, 1)));
      rho_star = logu < logMH ? rho_star : rho;
    }
    rho = rho_star;
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
    // prior on rho
    if (rho_prior_type == "pc") {
      out -= rhoPCPriorlogd(rho);
    } else {
      out -= (a_rho + 1) * log(1 - rho) + 1 / (b_rho * (1 - rho));
    }
    return out;
  }

  double KLD(double rho_star) {
    double out = (double) -N * rho_star / 2.0;
    out -= sum(log(rho_star / Lambda + (1 - rho_star))) / 2;
    out += sum(rho_star / Lambda) / 2;
    return out;
  }

  double rhoPCPriorlogd(double rho_star) {
    double d = sqrt(2.0 * KLD(rho_star));
    double logd = log(lambda_rho) - lambda_rho * d;
    return logd;
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
  arma::vec Uty;
  arma::vec Xty;
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
  .field("lower_rho", &BYM2FlatBetaMCMC::lower_rho)
  .field("upper_rho", &BYM2FlatBetaMCMC::upper_rho)
  .field("lambda_rho", &BYM2FlatBetaMCMC::lambda_rho)
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
  .method("gammaMean", &BYM2FlatBetaMCMC::gammaMean)
  ;
}
