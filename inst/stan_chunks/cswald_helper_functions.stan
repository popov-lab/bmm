// log-PDF of the shifted Wald distribution
// Optimized to compute entirely in log-space for numerical stability
real swald_lpdf(real rt, real drift, real bound, real ndt, real sigma) {
  // compute shifted response time
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return negative_infinity();

  // pre-compute common terms
  real sigma_sq = square(sigma);
  real log_t = log(t_shifted);

  // log-normalization: log(bound / sqrt(2*pi*sigma^2*t^3))
  // = log(bound) - 0.5*(log(2*pi) + 2*log(sigma) + 3*log(t))
  real log_norm = log(bound) - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log_t);

  // log-kernel: -((bound - drift*t)^2) / (2*sigma^2*t)
  real residual = bound - drift * t_shifted;
  real log_kernel = -0.5 * square(residual) / (sigma_sq * t_shifted);

  return log_norm + log_kernel;
}

// log shifted Wald survivor function (complementary CDF)
// Optimized for numerical stability using Stan's log-space CDF functions
// and log_diff_exp to avoid overflow/underflow issues
real swald_lccdf(real rt, real drift, real bound, real ndt, real sigma) {
  // compute shifted response time
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return 0;  // process hasn't started, survival = 1

  // pre-compute common terms
  real sqrt_t = sqrt(t_shifted);
  real sigma_sqrt_t = sigma * sqrt_t;
  real sigma_sq = square(sigma);

  // standardized arguments for the normal CDF
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;

  // log of the scaling constant (kept in log-space to avoid overflow)
  real log_c = 2 * bound * drift / sigma_sq;

  // Survival function: S(t) = 1 - Phi(z1) - exp(c)*Phi(z2)
  //                         = Phi(-z1) - exp(c)*Phi(z2)
  //
  // Compute in log-space for numerical stability:
  // log(S(t)) = log(exp(log_Phi(-z1)) - exp(log_c + log_Phi(z2)))
  //           = log_diff_exp(std_normal_lccdf(z1), log_c + std_normal_lcdf(z2))
  real log_term1 = std_normal_lccdf(z1 | );  // log(1 - Phi(z1)) = log(Phi(-z1))
  real log_term2 = log_c + std_normal_lcdf(z2 | );  // log(exp(c) * Phi(z2))

  // log_diff_exp(a, b) = log(exp(a) - exp(b)), stable when a > b
  return log_diff_exp(log_term1, log_term2);
}
