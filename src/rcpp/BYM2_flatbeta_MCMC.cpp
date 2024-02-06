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
  double lambda_rho;
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
                   arma::mat Q_,
                   double a_sigma_,
                   double b_sigma_,
                   double lambda_rho_):
    X(X_), y(y_), Q(Q_),
    a_sigma(a_sigma_), b_sigma(b_sigma_), lambda_rho(lambda_rho_) {
    N = X.n_rows;
    astar_sigma = a_sigma + (double) N / 2.0;
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

  List MCMCSample(int n_samples) {
    arma::mat beta_sim = zeros(n_samples, p);
    arma::mat gamma_sim = zeros(n_samples, N);
    arma::vec sigma2_sim = zeros(n_samples);
    arma::vec rho_sim = zeros(n_samples);
    arma::mat YFit_sim = zeros(n_samples, N);
    arma::vec mu = zeros(N);
    for (int i = 0; i < n_samples; i++) {
      updateBeta();
      beta_sim.row(i) = beta.t();
      updateGamma();
      gamma_sim.row(i) = gamma.t();
      updateSigma2(true);
      sigma2_sim(i) = sigma2;
      updateRho(false);
      rho_sim(i) = rho;
      mu = X * beta + gamma;
      for (int j = 0; j < N; j++) {
        YFit_sim(i, j) = randn(distr_param(mu[j], sqrt(sigma2 * (1 - rho))));
      }
    }
    List out = List::create(_["beta"] = beta_sim,
                            _["gamma"] = gamma_sim,
                            _["sigma2"] = sigma2_sim,
                            _["rho"] = rho_sim,
                            _["YFit"] = YFit_sim);
    return out;
  }

  void burnMCMCSample(int n_iter) {
    for (int i = 0; i < n_iter; i++) {
      updateBeta();
      updateGamma();
      updateSigma2(true);
      updateRho(false);
    }
  }

  double logPostd(bool update_rtr = true) {
    if (update_rtr) {
      vec_N = y - X * beta - gamma;
      rtr = dot(vec_N, vec_N);
    }
    double out = (double) -N / 2.0 * log(2.0 * datum::pi * sigma2 * (1 - rho));
    out -= rtr / (2 * sigma2 * (1 - rho));
    // IG prior on sigma2
    out -= (a_sigma + 1) * log(sigma2) + 1 / (b_sigma * sigma2);
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
    double theta = randg(distr_param(1.0, 1.0 / lambda_rho));
    double rho_star = exp(-1.0 * theta);
    double MH = (double) N / 2.0 * log((1.0 - rho) / (1.0 - rho_star));
    MH += rtr / (2.0 * sigma2) * (1.0 / (1.0 - rho) - 1.0 / (1.0 - rho_star));
    MH += (lambda_rho - 1) * log(rho / rho_star);
    if (log(randu(distr_param(0, 1))) <= MH) {
      rho = rho_star;
    }
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
  .constructor<arma::mat, arma::vec, arma::mat, double, double, double>()

  .field("X", &BYM2FlatBetaMCMC::X)
  .field("y", &BYM2FlatBetaMCMC::y)
  .field("Q", &BYM2FlatBetaMCMC::Q)
  .field("beta", &BYM2FlatBetaMCMC::beta)
  .field("gamma", &BYM2FlatBetaMCMC::gamma)
  .field("rho", &BYM2FlatBetaMCMC::rho)
  .field("sigma2", &BYM2FlatBetaMCMC::sigma2)
  .field("lambda_rho", &BYM2FlatBetaMCMC::lambda_rho)
  .field("astar_sigma", &BYM2FlatBetaMCMC::astar_sigma)
  .field("bstar_sigma", &BYM2FlatBetaMCMC::bstar_sigma)
  .field("vec_N", &BYM2FlatBetaMCMC::vec_N)
  .field("Lambda", &BYM2FlatBetaMCMC::Lambda)
  .field("P", &BYM2FlatBetaMCMC::P)
  .field("matrix_pp", &BYM2FlatBetaMCMC::matrix_pp)
  .field("vec_p", &BYM2FlatBetaMCMC::vec_p)
  .field("vec2_p", &BYM2FlatBetaMCMC::vec2_p)
  .field("U", &BYM2FlatBetaMCMC::U)
  .field("D", &BYM2FlatBetaMCMC::D)

  .method("initOLS", &BYM2FlatBetaMCMC::initOLS)
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
