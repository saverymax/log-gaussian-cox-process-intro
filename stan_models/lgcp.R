# https://github.com/mbjoseph/gpp-speed-test/blob/master/stan/full_pois.stan
log_gaussian_cox_process <- '
data {
  int<lower = 1> N;
  vector[N] X;
  vector[N] Z;
  array[N] int y;
  matrix[N, N] D;
}

parameters {
  vector[N] r;
  real<lower=0> eta;
  real<lower=0> phi;
  real<lower=0> sigma;
  real alpha;
  real beta;
  real gamma;
  real delta;
}

transformed parameters {
  cov_matrix[N] Sigma;
  vector[N] w;
  real<lower = 0> eta_sq;
  real<lower = 0> sig_sq;
  
  eta_sq = pow(eta, 2);
  sig_sq = pow(sigma, 2);
  
  for (i in 1:(N-1)) {
    for (j in (i + 1):N) {
      Sigma[i, j] = eta_sq * exp(-D[i, j] * phi);
      Sigma[j, i] = Sigma[i, j];
    }
  }
  // For efficiency use CD: https://mc-stan.org/docs/2_18/stan-users-guide/simulating-from-a-gaussian-process.html
  for (k in 1:N) Sigma[k, k] = eta_sq + sig_sq;
  w = cholesky_decompose(Sigma) * r;
}

model {
  eta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  phi ~ normal(0, 5);
  r ~ normal(0, 1);
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  delta ~ normal(0, 1);
  y ~ poisson_log(alpha + beta * X + gamma + delta * Z + w);
}
generated quantities{
  vector[N] lambda_rep;
  vector[N] b_rep;
  vector[N] lambda_bias;
  array[N] int y_rep;
  lambda_rep = exp(alpha + beta * X);
  b_rep = exp(gamma + delta * Z);
  // Thin the process
  lambda_bias = exp(alpha + beta * X + gamma + delta * Z);
  y_rep = poisson_log_rng(alpha + beta * X + gamma + delta * Z + w);
}
'