// Single-accumulator log-PDF for LBA with lognormal drift
// Exact formulas follow the R reference implementation in distributions.R
real lba_lognormal_single_lpdf(real t, real v, real b, real A, real s) {
  real hi = b / t;
  real lo = (b - A) / t;
  real s2 = square(s);
  real log_M = v + s2 / 2 +
    log_diff_exp(
      std_normal_lcdf((log(hi) - v - s2) / s | ),
      std_normal_lcdf((log(lo) - v - s2) / s | )
    );
  return log_M - log(A);
}

// Single-accumulator log-survival for LBA with lognormal drift
real lba_lognormal_single_lccdf(real t, real v, real b, real A, real s) {
  real hi = b / t;
  real lo = (b - A) / t;
  real s2 = square(s);
  real log_M = v + s2 / 2 +
    log_diff_exp(
      std_normal_lcdf((log(hi) - v - s2) / s | ),
      std_normal_lcdf((log(lo) - v - s2) / s | )
    );
  real log_u = log_diff_exp(
    log(b) + std_normal_lcdf((log(hi) - v) / s | ),
    log(b - A) + std_normal_lcdf((log(lo) - v) / s | )
  );
  real surv_num = exp(log_u) - t * exp(log_M);
  return lba_log_positive(surv_num) - log(A);
}
