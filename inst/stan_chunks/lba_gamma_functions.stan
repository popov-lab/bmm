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
  real surv_num = exp(log_u) - t * exp(log_M);
  return lba_log_positive(surv_num) - log(A);
}
