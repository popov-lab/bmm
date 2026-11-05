// log probability of one Likert response under the cumulative model; needs sdt_dist_functions.stan,
// threshold_functions.stan and ordinal_functions.stan. R twin: dlikert_cumulative()

// the K - 1 thresholds are built around 0 and the latent variable has mean `activation` and SD 1;
// threshold types without a spacing get 0, which they ignore
real likert_cumulative_log_prob(int y, real activation, real spacing, array[] real deltas,
                                int K, int dist_type, int thresh_type) {
  vector[K - 1] thr = make_thresholds(0, spacing, deltas, K, thresh_type);
  return ordinal_log_prob_cat(y, thr, activation, 1, dist_type);
}
