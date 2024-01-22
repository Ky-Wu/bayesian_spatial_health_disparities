data {
  int<lower=0> N;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2; // and node1[i] < node2[i]
  matrix[N, N] Sigma_scaled; // var-cov matrix of phi
  vector[N] mu_phi;
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

  vector[N] phi; // spatial effects
}
transformed parameters {
  real sigmaSp;
  real sigmaNSp;
  sigmaSp = sqrt(sigma2 * rho);
  sigmaNSp = sqrt((1 - rho) * sigma2);
}
model {
  Y ~ normal(X * beta + phi * sigmaSp, sigmaNSp);
  phi ~ multi_normal(mu_phi, Sigma_scaled);
  beta ~ normal(0.0, 1.0);
  sigma2 ~ inv_gamma(a0_sigma, b0_sigma);
  rho ~ beta(0.5, 0.5);
}
//generated quantities {
  //real logit_rho = log(rho / (1.0 - rho));
  //real YFit[N] = normal_rng(X * beta + phi * sigmaSp, sigmaNSp);
//}
