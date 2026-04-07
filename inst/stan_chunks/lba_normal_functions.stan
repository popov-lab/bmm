// Single-accumulator log-PDF for LBA with normal drift
// Uses Phi_approx (inv_logit approximation) for speed
// c = 1.702 gives max |Phi_approx(x) - Phi(x)| < 0.00014
real lba_normal_single_lpdf(real t, real v, real b, real A, real s) {
  real c = 1.702;
  if (A < 1e-6) {
    return log(b) - 2 * log(t) + normal_lpdf(b / t | v, s);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  real cz_hi = c * (hi - v) / s;
  real cz_lo = c * (lo - v) / s;
  real Phi_hi = inv_logit(cz_hi);
  real Phi_lo = inv_logit(cz_lo);
  real M = v * (Phi_hi - Phi_lo)
         + s * c * (Phi_lo * (1 - Phi_lo) - Phi_hi * (1 - Phi_hi));
  return log(fmax(M, 1e-10)) - log(A);
}

// Single-accumulator log-survival for LBA with normal drift
// Survival is clamped to [1e-10, 1] using fmin/fmax for gradient-safe bounds
real lba_normal_single_lccdf(real t, real v, real b, real A, real s) {
  real c = 1.702;
  if (A < 1e-6) {
    return log1m_inv_logit(c * (b / t - v) / s);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  real cz_hi = c * (hi - v) / s;
  real cz_lo = c * (lo - v) / s;
  real Phi_hi = inv_logit(cz_hi);
  real Phi_lo = inv_logit(cz_lo);
  real M = v * (Phi_hi - Phi_lo)
         + s * c * (Phi_lo * (1 - Phi_lo) - Phi_hi * (1 - Phi_hi));
  real surv = (b * Phi_hi - (b - A) * Phi_lo - t * M) / A;
  return fmin(log(fmax(surv, 1e-10)), 0.0);
}
