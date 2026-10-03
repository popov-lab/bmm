// noise-distribution dispatchers shared by the SDT models
// dist_type follows .sdt_dists in R/distributions.R: 1 = normal, 2 = gumbel_min, 3 = gumbel_max, 4 = logistic

// CDF dispatch: F(eta)
real sdt_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return Phi(eta);
  if (dist_type == 2) return 1 - exp(-exp(eta));          // gumbel_min
  if (dist_type == 3) return exp(-exp(-eta));             // gumbel_max
  if (dist_type == 4) return inv_logit(eta);              // logistic
  reject("sdt_cumprob: unknown dist_type: ", dist_type);
}

// Log-CDF dispatch: log(F(eta))
real sdt_log_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(eta);
  if (dist_type == 2) return log1m_exp(-exp(eta));       // gumbel_min
  if (dist_type == 3) return -exp(-eta);                  // gumbel_max
  if (dist_type == 4) return -log1p_exp(-eta);            // logistic
  reject("sdt_log_cumprob: unknown dist_type: ", dist_type);
}

// converts d_a to noise-SD units (1 at sdratio = 1); takes sdratio on the natural scale
real sdt_rms_scale(real sdratio) {
  return sqrt((1 + square(sdratio)) / 2);
}

// Log complementary CDF: log(1 - F(eta))
real sdt_log_one_minus_cumprob(real eta, int dist_type) {
  // std_normal_lccdf(eta) underflows to -inf from eta ~ 8.3, std_normal_lcdf(-eta) does not
  if (dist_type == 1) return std_normal_lcdf(-eta);
  if (dist_type == 2) return -exp(eta);                   // gumbel_min
  if (dist_type == 3) return log1m_exp(-exp(-eta));       // gumbel_max
  if (dist_type == 4) return -eta - log1p_exp(-eta);      // logistic
  reject("sdt_log_one_minus_cumprob: unknown dist_type: ", dist_type);
}

// Quantile (inverse-CDF) dispatch: F^{-1}(u), u in (0, 1)
real sdt_quantile(real u, int dist_type) {
  if (dist_type == 1) return inv_Phi(u);
  if (dist_type == 2) return log(-log1m(u));             // gumbel_min
  if (dist_type == 3) return -log(-log(u));              // gumbel_max
  if (dist_type == 4) return logit(u);                   // logistic
  reject("sdt_quantile: unknown dist_type: ", dist_type);
}
