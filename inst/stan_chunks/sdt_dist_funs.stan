// Shared noise-distribution dispatchers for SDT models
// dist_type codes follow the position in the .sdt_dists registry
// (R/distributions.R): 1=Gaussian, 2=Gumbel_min, 3=Gumbel_max, 4=Logistic

// Log-CDF dispatch: log(F(eta))
// Uses numerically stable Stan primitives for each distribution
real sdt_log_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(eta);
  if (dist_type == 2) return -exp(-eta);                 // gumbel_min
  if (dist_type == 3) return log1m_exp(-exp(eta));       // gumbel_max
  return -log1p_exp(-eta);                                // logistic
}

// Log complementary CDF: log(1 - F(eta))
real sdt_log_one_minus_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return std_normal_lccdf(eta);
  if (dist_type == 2) return log1m_exp(-exp(-eta));      // gumbel_min
  if (dist_type == 3) return -exp(eta);                   // gumbel_max
  return -eta - log1p_exp(-eta);                          // logistic
}
