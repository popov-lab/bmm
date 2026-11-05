// log probability of one category of an ordered-threshold model; needs sdt_dist_functions.stan
// R twin: .ordinal_category_log_probs() in R/distributions.R

// thresholds are on the noise-standardized axis; the latent variable has mean `location` and SD `scale`
real ordinal_log_prob_cat(int cat, vector thresholds, real location, real scale,
                          int dist_type) {
  int K_full = num_elements(thresholds) + 1;
  if (cat == 1) {
    return sdt_log_cumprob((thresholds[1] - location) / scale, dist_type);
  }
  if (cat == K_full) {
    return sdt_log_one_minus_cumprob(
      (thresholds[K_full - 1] - location) / scale, dist_type
    );
  }

  // an interval above 0 is taken from the upper tails, where both log cdf values would round to 0
  real eta_lo = (thresholds[cat - 1] - location) / scale;
  real eta_hi = (thresholds[cat] - location) / scale;
  if (eta_lo > 0) {
    return log_diff_exp(sdt_log_one_minus_cumprob(eta_lo, dist_type),
                        sdt_log_one_minus_cumprob(eta_hi, dist_type));
  }
  return log_diff_exp(sdt_log_cumprob(eta_hi, dist_type),
                      sdt_log_cumprob(eta_lo, dist_type));
}
