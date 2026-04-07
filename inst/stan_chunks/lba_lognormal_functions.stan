// Single-accumulator log-PDF for LBA with lognormal drift
// Uses Phi_approx (inv_logit approximation) for speed
real lba_lognormal_single_lpdf(real t, real v, real b, real A, real s) {
  if (A < 1e-6) {
    return log(b) - 2 * log(t) + lognormal_lpdf(b / t | v, s);
  }
  real c = 1.702;
  real hi = b / t;
  real lo = (b - A) / t;
  real s2 = square(s);
  real log_hi = log(hi);
  real log_lo = log(lo);
  real mu_s2 = exp(v + s2 / 2);
  real pdf_val = (mu_s2 / A) * (
    inv_logit(c * (log_hi - v - s2) / s) -
    inv_logit(c * (log_lo - v - s2) / s));
  return log(fmax(pdf_val, 1e-10));
}

// Single-accumulator log-survival for LBA with lognormal drift
// Survival is clamped to [1e-10, 1] using fmin/fmax for gradient-safe bounds
real lba_lognormal_single_lccdf(real t, real v, real b, real A, real s) {
  real c = 1.702;
  if (A < 1e-6) {
    return log1m_inv_logit(c * (log(b / t) - v) / s);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  real s2 = square(s);
  real log_hi = log(hi);
  real log_lo = log(lo);
  real mu_s2 = exp(v + s2 / 2);
  real M = mu_s2 * (
    inv_logit(c * (log_hi - v - s2) / s) -
    inv_logit(c * (log_lo - v - s2) / s));
  real F_hi = inv_logit(c * (log_hi - v) / s);
  real F_lo = inv_logit(c * (log_lo - v) / s);
  real surv = (b * F_hi - (b - A) * F_lo - t * M) / A;
  return fmin(log(fmax(surv, 1e-10)), 0.0);
}
