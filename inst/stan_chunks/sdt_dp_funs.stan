// DPSDT log-probability for a single rating observation.
// Requires sdt_rating_logprob_category() from sdt_rating_funs.stan.

real sdt_dp_row_lpmf(int y, int category, int stimulus,
                     int n_ratings, int dist_type, int thresh_type,
                     real dprime, real criterion, real spacing,
                     real Ro, real Rn, real sdratio,
                     array[] real deltas) {
  if (y == 0) return 0;

  vector[n_ratings - 1] thresholds = sdt_make_thresholds_rating(
    criterion, spacing, deltas, n_ratings, thresh_type
  );
  int K_full = num_elements(thresholds) + 1;
  real log_base = sdt_rating_logprob_category(
    category, stimulus, thresholds, dprime, criterion, sdratio, dist_type
  );

  if (stimulus == 1) {
    if (category == K_full) {
      return y * log_sum_exp(
        log1m_inv_logit(Ro) + log_base,
        log_inv_logit(Ro)
      );
    }
    return y * (log1m_inv_logit(Ro) + log_base);
  }

  if (category == 1) {
    return y * log_sum_exp(
      log1m_inv_logit(Rn) + log_base,
      log_inv_logit(Rn)
    );
  }
  return y * (log1m_inv_logit(Rn) + log_base);
}
