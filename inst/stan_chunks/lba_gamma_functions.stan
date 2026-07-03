// Single-accumulator log-PDF for LBA with gamma drift
// Exact formulas follow the R reference implementation in distributions.R
real lba_gamma_single_lpdf(real t, real v, real b, real A, real s) {
  real hi = b / t;
  real lo = (b - A) / t;
  real log_M = log(v) - log(s) +
    log_diff_exp(gamma_lcdf(hi | v + 1, s), gamma_lcdf(lo | v + 1, s));
  return log_M - log(A);
}

// Single-accumulator log-survival for LBA with gamma drift
real lba_gamma_single_lccdf(real t, real v, real b, real A, real s) {
  real hi = b / t;
  real lo = (b - A) / t;
  real log_M = log(v) - log(s) +
    log_diff_exp(gamma_lcdf(hi | v + 1, s), gamma_lcdf(lo | v + 1, s));
  real log_u = log_diff_exp(
    log(b) + gamma_lcdf(hi | v, s),
    log(b - A) + gamma_lcdf(lo | v, s)
  );
  // survival numerator (u - t*M) in log space: log_diff_exp is exact where
  // log_u > log(t) + log_M; the cancellation tail (survival numerically 0) floors
  // to the same penalty as lba_log_clip (log(1e-300)).
  real log_tM = log(t) + log_M;
  return (log_u > log_tM ? log_diff_exp(log_u, log_tM) : log(1e-300)) - log(A);
}
