// CDF dispatch for noise distribution selection
// dist_type: 1=Gaussian, 2=Gumbel_min, 3=Gumbel_max, 4=Logistic
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

// Binary SDT likelihood with binomial on aggregated data
// Parameters:
//   y:         count of "old"/"signal" responses
//   mu:        internal parameter (fixed to 0, required by brms)
//   dprime:    sensitivity (on link scale, transformed by brms)
//   criterion: response bias (identity link)
//   stimulus:  0 = noise/new, 1 = signal/old
//   dist_type: noise distribution (1-4, see sdt_cumprob)
//   trials:    total number of trials in this cell
real sdt_binary_lpmf(int y, real mu, real dprime, real criterion,
                     int stimulus, int dist_type, int trials) {
  real eta = dprime / 2.0 * (2 * stimulus - 1) - criterion;
  real p = sdt_cumprob(eta, dist_type);
  return binomial_lpmf(y | trials, p);
}
