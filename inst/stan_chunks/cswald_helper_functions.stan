// log(exp(a) - exp(b)), returning -inf when b >= a
real swald_log_diff_exp(real a, real b) {
  if (b >= a) return negative_infinity();
  return log_diff_exp(a, b);
}

// log of the standard normal CDF, stable in the far lower tail
real swald_log_Phi(real z) {
  real p = Phi(z);
  return p >= 1e-300 ? log(p) : std_normal_lcdf(z | );
}

// log survivor function of the shifted Wald distribution from its standardised arguments
real swald_log_surv(real z1, real z2, real log_c) {
  real surv = Phi(-z1) - exp(log_c + swald_log_Phi(z2));
  if (surv < 1e-300) {
    return swald_log_diff_exp(std_normal_lcdf(-z1 | ),
                              log_c + std_normal_lcdf(z2 | ));
  }
  return log(surv);
}

// log density of the shifted Wald distribution
real swald_lpdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0 || bound <= 0) return negative_infinity();

  real sigma_sq = square(sigma);
  real log_t = log(t_shifted);

  real log_norm = log(bound) - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log_t);

  real residual = bound - drift * t_shifted;
  real log_kernel = -0.5 * square(residual) / (sigma_sq * t_shifted);

  return log_norm + log_kernel;
}

// log survivor function (1 - CDF) of the shifted Wald distribution
real swald_lccdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return 0;
  if (bound <= 0) return negative_infinity();

  real sigma_sqrt_t = sigma * sqrt(t_shifted);
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;
  real log_c = 2 * bound * drift / square(sigma);

  return swald_log_surv(z1, z2, log_c);
}

// shifted Wald log density for a vector of decision times t = rt - ndt
vector swald_log_dens_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  if (min(t) <= 0 || min(bound) <= 0) return rep_vector(negative_infinity(), n);
  vector[n] residual = bound - drift .* t;
  return log(bound)
         - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log(t))
         - 0.5 * square(residual) ./ (square(sigma) .* t);
}

// shifted Wald log survivor function for a vector of decision times t = rt - ndt
vector swald_log_surv_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  if (min(bound) <= 0) return rep_vector(negative_infinity(), n);
  vector[n] denom = sigma .* sqrt(t);
  vector[n] dxt = drift .* t;
  vector[n] z1 = (dxt - bound) ./ denom;
  vector[n] z2 = -(dxt + bound) ./ denom;
  vector[n] log_c = 2 * (bound .* drift) ./ square(sigma);
  vector[n] p2 = Phi(z2);
  vector[n] surv = Phi(-z1) - exp(fmin(log_c + log(fmax(p2, 1e-300)), 0));
  vector[n] out = log(fmax(surv, 1e-300));
  for (k in 1:n) {
    if (p2[k] < 1e-300 || surv[k] < 1e-300) {
      out[k] = swald_log_surv(z1[k], z2[k], log_c[k]);
    }
  }
  return out;
}
