data {
  int<lower=0> N;
  // matrix[N, N] Sigma_scaled; // var-cov matrix of phi
  matrix[N, N] Sigma_chol;
  //vector[N] mu_phi;
  real Y[N]; // response
  int<lower=1> p; // num covariates
  matrix[N, p] X; // design matrix
  real<lower=0> a0_sigma;
  real<lower=0> b0_sigma;
}
parameters {
  vector[p] beta; // covariates
  real<lower=0> sigma2; // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  vector[N] eta; // for spatial error
}
transformed parameters {
  vector[N] phi;
  phi = Sigma_chol * eta;
  vector[N] gamma;
  gamma = sqrt(rho * sigma2) * phi;
  vector [N] mu;
  mu = X * beta + gamma;
}
model {
  Y ~ normal(mu, sqrt(sigma2 * (1 - rho)));
  //Y ~ normal_id_glm(X, gamma, beta, sqrt(sigma2 * (1 - rho)));
  //beta ~ normal(0.0, 100.0);
  eta ~ normal(0.0, 1.0);
  //phi ~ multi_normal(mu_phi, Sigma_scaled);
  sigma2 ~ inv_gamma(a0_sigma, b0_sigma);
  rho ~ beta(1, 1);
}
generated quantities {
  real logit_rho = log(rho / (1.0 - rho));
  real YFit[N] = normal_rng(mu, sqrt(sigma2 * (1 - rho)));
}
