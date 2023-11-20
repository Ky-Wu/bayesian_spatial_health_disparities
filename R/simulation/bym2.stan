data {
  int<lower=0> N;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2; // and node1[i] < node2[i]

  array[N] real Y; // continuous outcomes
  int<lower=1> p; // num covariates
  matrix[N, p] X; // design matrix

  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}
parameters {
  real beta0; // intercept
  vector[p] beta; // covariates

  real<lower=0> sigma; // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N] phi; // spatial effects
}
model {
  Y ~ normal(beta0 + X * beta + sigma * sqrt(rho / scaling_factor) * phi, sigma * sqrt(1 - rho));
  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  // soft sum-to-zero constraint on phi
  sum(phi) ~ normal(0, 0.001 * N); // equivalent to mean(phi) ~ normal(0,0.001)

  beta0 ~ normal(0.0, 1.0);
  beta ~ normal(0.0, 1.0);
  sigma ~ normal(0, 1.0);
  rho ~ beta(0.5, 0.5);
}
generated quantities {
  real logit_rho = log(rho / (1.0 - rho));
}
