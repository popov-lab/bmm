// Gamma drift d ~ Gamma(shape v, rate s). M = (v / s) [F(hi; v + 1) - F(lo; v + 1)].
//
// Two facts about Stan Math's gamma CDFs bind here (stan-dev/math #3408,
// measured on 5.3.0 and 5.4.0): gamma_lcdf rounds to 0 once log Q < -37, so
// a difference of two lcdf values is -Inf for every fast response with
// b / t >= 10, and its shape partial is off by 2.5e-3 at shape 5, by 0.1 to
// 50 % for shape >= 10 in the lower tail, NaN once P saturates (shape >= 12,
// s b / t above ~1e3) and an exception above 5e4. gamma_lccdf is finite far
// beyond exp underflow. So the difference is taken from the upper tail
// (lccdf pair) beyond the mean, from the lower tail (lcdf pair) below it,
// and where the interval is narrow the midpoint rule through gamma_lpdf
// (exact partials) replaces both. The shape is v + 1 with v an autodiff
// parameter, and s b / t grows without bound as ndt approaches the fastest
// response: the gamma drift is the least robust of the four there.

// log(F(hi) - F(lo)) of Gamma(alpha, beta), hi > lo
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
  // u_num - t M = t int_{lo}^{hi} F(u) du > 0; the difference loses digits
  // only where the survivor itself is negligible
  return (log_u > log_tM ? log_diff_exp(log_u, log_tM) : lba_log_floor()) - log(A);
}
