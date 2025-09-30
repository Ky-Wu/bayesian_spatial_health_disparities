functions{
  real rhoPCPrior_log(real x, int N, vector L, real l) {
    real kld = -N / 2;
    kld = kld + sum(x / L + (1 - x)) / 2;
    kld = kld - sum(log(x / L + (1 - x))) / 2;
    real d = sqrt(2 * kld);
    d =  log(l) - (l * d);
    return d;
  }
}
data {
  int<lower=0> N;
  // matrix[N, N] Sigma_scaled; // var-cov matrix of phi
  matrix[N, N] Sigma_chol;
  //vector[N] mu_phi;
  int<lower=0> Y[N]; // response
  vector<lower=0>[N] E;
  int<lower=1> p; // num covariates
  matrix[N, p] X; // design matrix
  real<lower=0> a0_sigma;
  real<lower=0> b0_sigma;
  real<lower=0> lambda_rho;
  vector<lower=0>[N] Lambda;
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  vector[p] beta; // covariates
  real<lower=0> sigma2; // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  vector[N] eta; // spatial effects
  vector[N] eps; // unstructered error
}
transformed parameters {
  vector[N] phi;
  phi = Sigma_chol * eta;
  vector[N] alpha;
  alpha = (sqrt(1 - rho) * eps + sqrt(rho) * phi) * sqrt(sigma2) + log_E;
}
model {
  Y ~ poisson_log_glm(X, alpha, beta);
  //beta ~ normal(0.0, 100.0);
  eta ~ normal(0.0, 1.0);
  eps ~ normal(0.0, 1.0);
  //phi ~ multi_normal(mu_phi, Sigma_scaled);
  sigma2 ~ inv_gamma(a0_sigma, b0_sigma);
  //rho ~ beta(1, 1);
  rho ~ rhoPCPrior(N, Lambda, lambda_rho);
}
generated quantities {
  real logit_rho = log(rho / (1.0 - rho));
  vector[N] log_mu = alpha + X * beta;
}
