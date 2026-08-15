// Shared noise-distribution dispatchers for SDT models
// dist_type codes follow the position in the .sdt_dists registry
// (R/distributions.R): 1=Gaussian, 2=Gumbel_min, 3=Gumbel_max, 4=Logistic
// Gumbel_min is the smallest-extreme-value distribution (cdf 1 - exp(-exp(x))),
// Gumbel_max the largest (cdf exp(-exp(-x))).

// Log-CDF dispatch: log(F(eta))
// Uses numerically stable Stan primitives for each distribution
real sdt_log_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(eta);
  if (dist_type == 2) return log1m_exp(-exp(eta));       // gumbel_min
  if (dist_type == 3) return -exp(-eta);                  // gumbel_max
  if (dist_type == 4) return -log1p_exp(-eta);            // logistic
  reject("sdt_log_cumprob: unknown dist_type: ", dist_type);
}

// Log complementary CDF: log(1 - F(eta))
real sdt_log_one_minus_cumprob(real eta, int dist_type) {
  // std_normal_lcdf(-eta), not std_normal_lccdf(eta): the two are equal by
  // symmetry, but lccdf loses accuracy from eta ~ 8 and underflows to -inf
  // from eta ~ 10, which makes the likelihood log(0) for any cell with
  // y < trials. lcdf stays exact past eta = 29.
  if (dist_type == 1) return std_normal_lcdf(-eta);
  if (dist_type == 2) return -exp(eta);                   // gumbel_min
  if (dist_type == 3) return log1m_exp(-exp(-eta));       // gumbel_max
  if (dist_type == 4) return -eta - log1p_exp(-eta);      // logistic
  reject("sdt_log_one_minus_cumprob: unknown dist_type: ", dist_type);
}
