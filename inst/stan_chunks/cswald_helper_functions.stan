// log(exp(a) - exp(b)) that returns -inf instead of NaN when rounding makes
// b >= a, so one underflowing survivor cannot poison the whole gradient
real swald_log_diff_exp(real a, real b) {
  if (b >= a) return negative_infinity();
  return log_diff_exp(a, b);
}

// log(Phi(z)) with an exact derivative: Phi's derivative is the normal density,
// whereas std_normal_lcdf's gradient is an approximation (relative error up to
// ~1e-4). The lcdf is used only where Phi rounds to 0 (z < ~-37.5)
real swald_log_Phi(real z) {
  real p = Phi(z);
  return p > 0 ? log(p) : std_normal_lcdf(z | );
}

// log shifted Wald survivor S = Phi(-z1) - exp(log_c) * Phi(z2), with
// z1 = (drift*t - bound) / (sigma*sqrt(t)), z2 = -(drift*t + bound) / (sigma*sqrt(t)),
// log_c = 2*bound*drift/sigma^2. Computed in probability space for the exact
// derivative (see swald_log_Phi); exp(log_c)*Phi(z2) <= Phi(-z1) <= 1, so the
// exp cannot overflow. The log-space form takes over when the difference
// underflows or rounds to <= 0
real swald_log_surv(real z1, real z2, real log_c) {
  real surv = Phi(-z1) - exp(log_c + swald_log_Phi(z2));
  if (surv <= 0) {
    return swald_log_diff_exp(std_normal_lcdf(-z1 | ),
                              log_c + std_normal_lcdf(z2 | ));
  }
  return log(surv);
}

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
real swald_lccdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return 0;  // process hasn't started, survival = 1

  real sigma_sqrt_t = sigma * sqrt(t_shifted);
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;
  real log_c = 2 * bound * drift / square(sigma);

  return swald_log_surv(z1, z2, log_c);
}

// Vectorized forms of swald_lpdf / swald_lccdf for the loop = FALSE families.
// Both take the already-shifted time t = rt - ndt. Keep them in step with the
// scalar forms above and with .dwald() / .pwald() in R/distributions.R

// t <= 0 anywhere makes the summed target -inf, exactly as the scalar sum would
vector swald_log_dens_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  if (min(t) <= 0) return rep_vector(negative_infinity(), n);
  vector[n] residual = bound - drift .* t;
  return log(bound)
         - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log(t))
         - 0.5 * square(residual) ./ (square(sigma) .* t);
}

// requires t > 0 elementwise. Elements whose probability-space survivor
// underflows are recomputed by the scalar swald_log_surv below. The fmax/fmin
// clamps are required: without them a log(0) or exp(inf) in a discarded
// element injects NaN adjoints into the shared parameters even though its
// value is overwritten. Every clamped element is one the loop recomputes
vector swald_log_surv_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
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
