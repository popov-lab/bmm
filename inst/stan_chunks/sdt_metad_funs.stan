// Meta-d' log-probability for a single rating observation.
// Requires sdt_rating_logprob_category() from sdt_rating_funs.stan.

real sdt_metad_row_lpmf(int y, int category, int stimulus,
                        int n_ratings, int dist_type, int thresh_type,
                        real dprime, real criterion, real spacing,
                        real metad, real sdratio,
                        array[] real deltas) {
  if (y == 0) return 0;

  int mid = (n_ratings - 1) / 2 + 1;
  vector[n_ratings - 1] thresholds = sdt_make_thresholds_rating(
    criterion, spacing, deltas, n_ratings, thresh_type
  );
  real scale = stimulus == 1 ? exp(sdratio) : 1.0;
  real d_shift = dprime / 2.0 * (2 * stimulus - 1);
  real metad_shift = metad / 2.0 * (2 * stimulus - 1);
  real crit = thresholds[mid];
  real log_norm;

  if (category <= mid) {
    log_norm = sdt_log_cumprob((crit - d_shift) / scale, dist_type) -
      sdt_log_cumprob((crit - metad_shift) / scale, dist_type);
  } else {
    log_norm = sdt_log_one_minus_cumprob((crit - d_shift) / scale, dist_type) -
      sdt_log_one_minus_cumprob((crit - metad_shift) / scale, dist_type);
  }

  return y * (
    sdt_rating_logprob_category(
      category, stimulus, thresholds, metad, criterion, sdratio, dist_type
    ) + log_norm
  );
}
