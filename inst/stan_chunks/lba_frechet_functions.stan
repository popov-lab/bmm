// Single-accumulator log-PDF for LBA with Frechet drift
// Uses 20-point trapezoidal quadrature for M integral when A >= 1e-6
real lba_frechet_single_lpdf(real t, real v, real b, real A, real s) {
  if (A < 1e-6) {
    real z = b / (t * s);
    return log(b) - 2 * log(t) + log(v) - log(s)
           - (1 + v) * log(z) - pow(z, -v);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  int n_quad = 20;
  real h = (hi - lo) / n_quad;
  real M = 0;
  for (q in 0:n_quad) {
    real u = lo + q * h;
    real log_z = log(u / s);
    real f_u = (v / s) * exp(-(1 + v) * log_z - exp(-v * log_z));
    real w = (q == 0 || q == n_quad) ? 0.5 : 1.0;
    M += w * u * f_u;
  }
  M *= h;
  return log(fmax(M, 1e-10)) - log(A);
}

// Single-accumulator log-survival for LBA with Frechet drift
// Survival is clamped to [1e-10, 1] using fmin/fmax for gradient-safe bounds
real lba_frechet_single_lccdf(real t, real v, real b, real A, real s) {
  if (A < 1e-6) {
    real log_F = -pow(b / (t * s), -v);
    return log1m_exp(log_F);
  }
  real hi = b / t;
  real lo = (b - A) / t;
  int n_quad = 20;
  real h = (hi - lo) / n_quad;
  real M = 0;
  for (q in 0:n_quad) {
    real u = lo + q * h;
    real log_z = log(u / s);
    real f_u = (v / s) * exp(-(1 + v) * log_z - exp(-v * log_z));
    real w = (q == 0 || q == n_quad) ? 0.5 : 1.0;
    M += w * u * f_u;
  }
  M *= h;
  real F_hi = exp(-pow(hi / s, -v));
  real F_lo = exp(-pow(lo / s, -v));
  real surv = (b * F_hi - (b - A) * F_lo - t * M) / A;
  return fmin(log(fmax(surv, 1e-10)), 0.0);
}
