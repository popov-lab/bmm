// Single-accumulator log-PDF for LBA with gamma drift
real lba_gamma_single_lpdf(real t, real v, real b, real A, real s) {
  if (A < 1e-6) {
    return log(b) - 2 * log(t) + gamma_lpdf(b / t | v, s);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  real v1 = v + 1;
  real pdf_val = (v / (s * A)) * (
    gamma_cdf(hi | v1, s) - gamma_cdf(lo | v1, s));
  return log(fmax(pdf_val, 1e-10));
}

// Single-accumulator log-survival for LBA with gamma drift
// Survival is clamped to [1e-10, 1] using fmin/fmax for gradient-safe bounds
real lba_gamma_single_lccdf(real t, real v, real b, real A, real s) {
  if (A < 1e-6) {
    return gamma_lccdf(b / t | v, s);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  real v1 = v + 1;
  real M = (v / s) * (gamma_cdf(hi | v1, s) - gamma_cdf(lo | v1, s));
  real surv = (b * gamma_cdf(hi | v, s) - (b - A) * gamma_cdf(lo | v, s)
               - t * M) / A;
  return fmin(log(fmax(surv, 1e-10)), 0.0);
}
