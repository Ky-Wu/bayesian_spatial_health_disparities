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
    double delta_rho;
    BYM2FlatBetaMCMC(arma::mat X_,
                  arma::vec y_,
                  arma::mat Q_,
                  double a_sigma_,
                  double b_sigma_,
                  double delta_rho_ = 0.1):
      X(X_), y(y_), Q(Q_),
      a_sigma(a_sigma_), b_sigma(b_sigma_), delta_rho(delta_rho_) {
      N = X.n_rows;
      p = X.n_cols;
      XtX = X.t() * X;
      XtX_R = chol(XtX, "upper");
      beta_ols = X.t() * y;
      beta_ols = solve(trimatl(XtX_R.t()), beta_ols);
      beta_ols = solve(trimatu(XtX_R), beta_ols);
      // simultaneous diagonalization of Q and I - H
      eig_sym(Lambda, P, Q);
      // store intermediate Q^{-1/2} temporarily in U
      U = P * diagmat(pow(Lambda, -0.5)) * P.t();
      // matrix_NN = I - H
      matrix_pN = solve(trimatl(XtX_R.t()), X.t());
      matrix_NN = -1 * matrix_pN.t() * matrix_pN;
      matrix_NN.diag() += 1;
      // Q^{-1/2} * (I - H) * Q^{-1/2}
      matrix_NN = U * matrix_NN * U;
      eig_sym(D, P, matrix_NN);
      U = U * P;
      UtX = U.t() * X;
      // e := y - X * beta
      // r := y - X * beta - gamma
      e = y - X * beta_ols;
      Ute = U.t() * e;
      vec_p = zeros(p);
      vec_N = zeros(N);
      matrix_pp = zeros(p, p);
      matrix_pN = zeros(p, N);
      matrix_NN = zeros(N, N);
      rho_lower = std::max({0.0, rho - delta_rho});
      rho_upper = std::min({1.0, rho + delta_rho});
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

    arma::mat MCMCSample(int n_iter) {
      arma::mat out = zeros(n_iter, N + p + 2);
      for (int i = 0; i < n_iter; i++) {
        updateBeta();
        out.submat(i, 0, i, p - 1) = beta;
        updateGamma();
        out.submat(i, p, i, N + p - 1) = gamma;
        updateSigma2();
        out.submat(i, N + p + 1, i, N + p + 1) = sigma2;
        updateRho();
        out.submat(i, N + p, i, N + p) = rho;

      }
      return out;
    }

    void BurnMCMCSample(int n_iter) {
      for (int i = 0; i < n_iter; i++) {
        updateBeta();
        updateGamma();
        updateRho();
        updateSigma2();
      }
    }
    void updateBeta() {
      // vec_N = 1 / (1 - rho) B^(-1) e
      vec_N = diagmat(1.0 / ((1.0 - rho) / rho + D)) * Ute;
      vec_N = U * vec_N;
      // vec_N = beta_ols - (X'X)^{-1}X'[gamma - (1 - rho)^{-1} B^(-1) e]
      vec_N = gamma - vec_N;
      vec_N = y - vec_N;
      beta = X.t() * vec_N;
      beta = solve(trimatl(XtX_R.t()), beta);
      // var(beta | y, phi, sigma2, rho) = sigma2 * (1 - rho) (X'X)^{-1}
      beta += randn(p, distr_param(0.0, pow(sigma2 * (1.0 - rho), 0.5)));
      beta = solve(trimatu(XtX_R), beta);
    }

    void updateGamma() {
      // var(gamma | y, beta, sigma2, rho) =
      // sigma2 * P[1 / (1 - rho) * I + 1 / rho * Lambda]^{-1}P^T
      vec_N = Lambda / rho;
      vec_N += 1.0 / (1.0 - rho);
      vec_N = pow(vec_N, -0.5);
      vec_N = vec_N % randn(N, distr_param(0, 1));
      vec_N = P * vec_N;
      vec_N *= pow(sigma2, 0.5);
      // E(gamma | y, beta, sigma2, rho) =
      // (1 - rho)^{-1}B^{-1}[e - XC^{-1}X'X(beta - beta_ols)]
      vec_p = beta - beta_ols;
      vec_p = XtX * vec_p;
      // matrix_pp = chol(C) = chol(X'X + (1 - rho)^{-1}X'B^{-1}X)
      matrix_pN = UtX.t() * diagmat(1 / ((1 - rho) / rho + D));
      matrix_pp =  matrix_pN * matrix_pN.t();
      matrix_pp += XtX;
      matrix_pp = chol(matrix_pp, "upper");
      vec_p = solve(trimatl(matrix_pp.t()), vec_p);
      vec_p = solve(trimatu(matrix_pp), vec_p);
      gamma = X * vec_p;
      gamma = e - gamma;
      gamma = U.t() * gamma;
      gamma = (1.0 / ((1.0 - rho) / rho + D)) % gamma;
      gamma = U * gamma + vec_N;
    }

    void updateSigma2(bool update_rtr = true) {
      double bstar_sigma;
      if (update_rtr) {
        vec_N = y - X * beta - gamma;
        rtr = dot(vec_N, vec_N);
      }
      bstar_sigma = b_sigma + rtr / (2 * (1 - rho));
      sigma2 = 1.0 / randg(distr_param(a_sigma + (double) N / 2, bstar_sigma));
    }

    void updateRho(bool update_rtr = true) {
      double rho_star = randu(distr_param(rho_lower, rho_upper));
      double new_lower = std::max({0.0, rho_star - delta_rho});
      double new_upper = std::min({1.0, rho_star + delta_rho});
      double MH = (double) N / 2 * log((1 - rho) / (1 - rho_star));
      MH += rtr / (2 * sigma2) * (1 / (1 - rho) - 1 / (1 - rho_star));
      MH += log((new_upper - new_lower) / (rho_upper - rho_lower));
      if (log(randu(distr_param(0, 1))) <= MH) {
        rho = rho_star;
        rho_lower = new_lower;
        rho_upper = new_upper;
      }
    }

  private:
    int N;
    int p;
    double rtr;
    double rho_lower;
    double rho_upper;
    arma::vec e;
    arma::vec Ute;
    arma::mat XtX;
    arma::mat XtX_R;
    arma::vec beta_ols;
    arma::vec Lambda;
    arma::mat P;
    arma::mat U;
    arma::mat UtX;
    arma::vec D;
    arma::vec vec_p;
    arma::vec vec_N;
    arma::mat matrix_pp;
    arma::mat matrix_pN;
    arma::mat matrix_NN;
};

RCPP_MODULE(BYM2FlatBetaMCMC_module) {

  class_<BYM2FlatBetaMCMC>("BYM2FlatBetaMCMC")
  .constructor<arma::mat, arma::vec, arma::mat, double, double>()

  .field("X", &BYM2FlatBetaMCMC::X)
  .field("y", &BYM2FlatBetaMCMC::y)
  .field("Q", &BYM2FlatBetaMCMC::Q)
  .field("beta", &BYM2FlatBetaMCMC::beta)
  .field("gamma", &BYM2FlatBetaMCMC::gamma)
  .field("rho", &BYM2FlatBetaMCMC::rho)
  .field("sigma2", &BYM2FlatBetaMCMC::sigma2)
  .field("delta_rho", &BYM2FlatBetaMCMC::delta_rho)

  .method("initOLS", &BYM2FlatBetaMCMC::initOLS)
  .method("MCMCSample", &BYM2FlatBetaMCMC::MCMCSample)
  .method("BurnMCMCSample", &BYM2FlatBetaMCMC::BurnMCMCSample)
  .method("updateBeta", &BYM2FlatBetaMCMC::updateBeta)
  .method("updateGamma", &BYM2FlatBetaMCMC::updateGamma)
  ;
}
