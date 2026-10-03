// Gamma drift (shape v, rate s); why the tails are split: .lba_gamma_log_dF() in R/distributions.R

// log(F(hi) - F(lo)) of Gamma(alpha, beta), hi > lo: midpoint, lccdf pair above the mean, lcdf pair below
real lba_gamma_log_dF(real lo, real hi, real alpha, real beta) {
  real du = hi - lo;
  real u_m = 0.5 * (lo + hi);
  // the midpoint rule's expansion parameter, du |d log f / du| at the midpoint
  if (du * (abs(alpha - 1) / u_m + beta) < 1e-3) {
    real r = (alpha - 1) / u_m - beta;
    return log(du) + gamma_lpdf(u_m | alpha, beta)
      + log1p(square(du) / 24 * (square(r) - (alpha - 1) / square(u_m)));
  }
  if (lo * beta > alpha) {
    return log_diff_exp(gamma_lccdf(lo | alpha, beta), gamma_lccdf(hi | alpha, beta));
  }
  return log_diff_exp(gamma_lcdf(hi | alpha, beta), gamma_lcdf(lo | alpha, beta));
}

// Single-accumulator log-PDF for LBA with gamma drift
real lba_gamma_single_lpdf(real t, real v, real b, real A, real s) {
  return log(v) - log(s) + lba_gamma_log_dF((b - A) / t, b / t, v + 1, s) - log(A);
}

// Single-accumulator log-survival for LBA with gamma drift
real lba_gamma_single_lccdf(real t, real v, real b, real A, real s) {
  real lo = (b - A) / t;
  real hi = b / t;
  real log_tM = log(t) + log(v) - log(s) + lba_gamma_log_dF(lo, hi, v + 1, s);
  real log_u = log_sum_exp(log(A) + gamma_lcdf(lo | v, s), log(b) + lba_gamma_log_dF(lo, hi, v, s));
  // u_num - t M > 0; guard before log_diff_exp, never clamp after it
  return (log_u > log_tM ? log_diff_exp(log_u, log_tM) : lba_log_floor()) - log(A);
}
