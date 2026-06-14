// Shared SDT noise-distribution dispatch
// dist_type: 1=Gaussian, 2=Gumbel_min, 3=Gumbel_max, 4=Logistic
// Used by every SDT model that needs the noise CDF, log-CDF, log-CCDF, or
// quantile on the latent evidence scale.

// CDF dispatch: F(eta)
real sdt_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return Phi(eta);
  if (dist_type == 2) return exp(-exp(-eta));
  if (dist_type == 3) return 1 - exp(-exp(eta));
  return inv_logit(eta);
}

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

// Quantile (inverse-CDF) dispatch: F^{-1}(u), u in (0, 1)
real sdt_quantile(real u, int dist_type) {
  if (dist_type == 1) return inv_Phi(u);
  if (dist_type == 2) return -log(-log(u));              // gumbel_min
  if (dist_type == 3) return log(-log1m(u));             // gumbel_max
  return logit(u);                                        // logistic
}
